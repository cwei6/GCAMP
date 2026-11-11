"""Inter-session pre-alignment: bring every run onto the first run's frame
before suite2p, so the same afferents land on the same pixels across runs.

suite2p's registration corrects frame-to-frame motion (translation + small
nonrigid deformation, limited to maxregshift ~10 % of the frame and NO
rotation). Between sessions (pre vs post an hour later, or different days) the
field is re-found by hand and can be offset by tens of pixels and rotated by a
few degrees — too much for suite2p alone. Here we:

  1. build a STRUCTURAL image per run (red vascular channel if every run has
     it, else the 10th-percentile of the green movie = mCitrine + resting GCaMP)
  2. search rotation (±max_angle_deg) x translation (FFT phase correlation) of
     each run onto run 0, keeping the transform with the best correlation
  3. report it; if it exceeds min_shift_px / min_angle_deg, the convert step
     applies it (rotation about the frame centre, then shift) to every frame of
     that run (both channels) — suite2p then only has residual motion to fix

Outputs prealign.json (transform + achieved correlation per run) and
qc/prealign.png (structural images before/after).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, rotate as nd_rotate, affine_transform


def structural_image(sbx_path: Path, use_red: bool, n_frames: int = 600) -> np.ndarray:
    from sbxreader import sbx_memmap          # only the suite2p route reads .sbx here
    mm = sbx_memmap(str(sbx_path))
    T = mm.shape[0]
    step = max(1, T // n_frames)
    if use_red and mm.shape[2] >= 2:
        img = np.asarray(mm[::step, 0, 1]).astype(np.float32).mean(axis=0)
    else:
        sub = np.asarray(mm[::step, 0, 0]).astype(np.float32)
        img = np.percentile(sub, 10, axis=0) if sub.shape[0] > 20 else sub.mean(axis=0)
    del mm
    return img


def _prep(img: np.ndarray) -> np.ndarray:
    hp = img - gaussian_filter(img, 15)          # remove vignetting / low-frequency
    hp = gaussian_filter(hp, 1.0)
    hp = hp - hp.mean()
    return hp / (hp.std() + 1e-9)


def _phase_corr(ref: np.ndarray, img: np.ndarray, max_shift: int):
    """(dy, dx, ncc) so that shifting img by (dy, dx) best matches ref."""
    F = np.fft.rfft2(ref) * np.conj(np.fft.rfft2(img))
    cc = np.fft.fftshift(np.fft.irfft2(F, s=ref.shape))
    cy, cx = ref.shape[0] // 2, ref.shape[1] // 2
    my, mx = min(cy, max_shift), min(cx, max_shift)
    win = cc[cy - my:cy + my + 1, cx - mx:cx + mx + 1]
    iy, ix = np.unravel_index(int(np.argmax(win)), win.shape)
    dy, dx = iy - my, ix - mx
    ncc = float(win[iy, ix] / ref.size)
    return int(dy), int(dx), ncc


def estimate_transform(ref: np.ndarray, img: np.ndarray, max_angle_deg: float = 6.0,
                       max_shift_frac: float = 0.3):
    """Search rotation x translation of img onto ref. Returns dict."""
    r = _prep(ref)
    max_shift = int(max_shift_frac * min(ref.shape))
    best = None
    coarse = np.arange(-max_angle_deg, max_angle_deg + 1e-9, 1.0)
    for stage, angles in (("coarse", coarse), ("fine", None)):
        if angles is None:
            a0 = best["angle_deg"]
            angles = np.arange(a0 - 0.75, a0 + 0.76, 0.25)
        for ang in angles:
            im = nd_rotate(img, ang, reshape=False, order=1, mode="nearest") if ang else img
            dy, dx, ncc = _phase_corr(r, _prep(im), max_shift)
            if best is None or ncc > best["ncc"]:
                best = {"angle_deg": float(ang), "dy": dy, "dx": dx, "ncc": ncc}
    # reference correlation with no transform, for the report
    _, _, ncc0 = _phase_corr(r, _prep(img), 0)
    best["ncc_before"] = ncc0
    return best


def apply_transform(frames: np.ndarray, angle_deg: float, dy: int, dx: int,
                    cval: float | None = None) -> np.ndarray:
    """Rotate (about centre) then shift a (T, H, W) block. Linear interpolation."""
    if not angle_deg and not dy and not dx:
        return frames
    T, H, W = frames.shape
    th = np.deg2rad(angle_deg)
    c, s = np.cos(th), np.sin(th)
    # output(y, x) = input(M @ (y - cy, x - cx) + centre - shift)
    M = np.array([[c, -s], [s, c]])
    centre = np.array([H / 2.0, W / 2.0])
    offset = centre - M @ (centre + np.array([dy, dx]))
    fill = float(np.median(frames[0])) if cval is None else cval
    out = np.empty_like(frames)
    for i in range(T):
        out[i] = affine_transform(frames[i].astype(np.float32), M, offset=offset,
                                  order=1, mode="constant", cval=fill).astype(frames.dtype)
    return out


def prealign_runs(specs, use_red: bool, out_dir: Path, max_angle_deg: float = 6.0,
                  min_shift_px: float = 2.0, min_angle_deg: float = 0.5,
                  enabled: bool = True) -> list[dict]:
    """Estimate a transform per run (relative to run 0). Returns a list of dicts
    with keys angle_deg, dy, dx, ncc, ncc_before, apply (bool)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    imgs = [structural_image(s.sbx_path, use_red) for s in specs]
    results = [{"angle_deg": 0.0, "dy": 0, "dx": 0, "ncc": 1.0, "ncc_before": 1.0, "apply": False}]
    print("Pre-alignment of runs onto the first run "
          f"({'red vascular channel' if use_red else 'green structural image'}):")
    print(f"  {specs[0].name}: reference")
    for s, img in zip(specs[1:], imgs[1:]):
        tr = estimate_transform(imgs[0], img, max_angle_deg=max_angle_deg)
        tr["apply"] = enabled and (abs(tr["dy"]) >= min_shift_px or abs(tr["dx"]) >= min_shift_px
                                   or abs(tr["angle_deg"]) >= min_angle_deg)
        results.append(tr)
        print(f"  {s.name}: shift (dy={tr['dy']}, dx={tr['dx']}) px, rotation {tr['angle_deg']:+.2f} deg, "
              f"structural correlation {tr['ncc_before']:.2f} -> {tr['ncc']:.2f}"
              + ("  [applied before suite2p]" if tr["apply"] else "  [small — left to suite2p]"))
        if tr["ncc"] < 0.3:
            print("    WARNING: low correlation even after alignment — this run may not be the same field of view")
    # figure: structural images before / after
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(specs)
        fig, axes = plt.subplots(2, n, figsize=(3.6 * n, 7), squeeze=False)
        ref = _prep(imgs[0])
        for i, (s, img, tr) in enumerate(zip(specs, imgs, results)):
            aligned = apply_transform(img[None].astype(np.float32), tr["angle_deg"], tr["dy"], tr["dx"])[0] \
                if i else img
            for row, im, title in ((0, img, "as recorded"), (1, aligned, f"aligned (r={tr['ncc']:.2f})")):
                a = _prep(im); rgb = np.dstack([np.clip(ref / 4 + .5, 0, 1), np.clip(a / 4 + .5, 0, 1), np.clip(ref / 4 + .5, 0, 1)])
                axes[row][i].imshow(rgb if i else np.clip(ref / 4 + .5, 0, 1), cmap="gray")
                axes[row][i].set_title(f"{s.condition}: {s.name}\n{title if i else 'reference'}", fontsize=8)
                axes[row][i].axis("off")
        fig.suptitle("Pre-alignment: magenta = reference run, green = this run (white = overlap)")
        fig.tight_layout(); fig.savefig(out_dir / "prealign.png", dpi=130); plt.close(fig)
    except Exception as e:  # pragma: no cover
        print(f"  prealign: figure skipped ({e})")
    (out_dir.parent / "prealign.json").write_text(json.dumps(
        [{"run": s.name, **r} for s, r in zip(specs, results)], indent=2), encoding="utf-8")
    return results
