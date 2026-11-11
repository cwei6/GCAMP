"""fullreg — rigid registration of the RAW frames (15.6 Hz) before binning,
in Python, for the registration comparison (reg_compare.py) and as an
optional first step of the unpack. It does in one translation-only pass what
the lab's legacy MATLAB chain does in two (ProcessScanboxData_concatFirst.m,
dataSet 'Afferents', method 'rigid': CorrectExpt removes the sub-pixel jitter
of every raw frame against the mean of its own 15-frame chunk by DFT
cross-correlation, then AlignPlanes / RegisterSBX aligns every frame to the
mean of the longest still epoch with MATLAB's intensity-based 'rigid'
imregtform on edge-cropped, 2x-binned frames, then GenerateExptProjections
bins). The legacy chain itself, run by MATLAB, is legacy_register_run.m
(variant D of reg_compare.py); this module is the Python stand-in and the
pipeline's own option.

What it does, per recording
  1. reference: raw frames from the longest still stretch of the run (the
     binned behaviour csv next to the .sbx; without one, the middle of the
     run), 300 of them spread over the stretch, averaged, then each aligned to
     that average and averaged again (two rounds) — a sharp still-period
     template, as DetermineReference / suite2p's refImg do;
  2. every raw frame: rigid shift against the reference by phase correlation
     (frames lightly smoothed, slow background removed, 2x downsampled for
     speed, parabolic sub-pixel peak; afferent2p.motion._Matcher — the same
     measure the motion check uses), applied as a whole-pixel translation (no
     interpolation blur); the border a shift uncovers is filled with the
     reference so the edge stays static;
  3. the registered frames are averaged into bins of N frames (the same bins
     as unpack_bin.py), written as uint16 ImageJ tiff, cropped by the 99.5th
     percentile of |shift| on every side so the filled border is gone;
  4. sidecars like the unpack: <stem>_<tag>.json (with the registration
     facts: reference frames, shift statistics, crop), _meanproj.tif, the
     behaviour csv (copied — same bins), plus <stem>_<tag>_fullreg_shifts.csv
     (dy, dx, corr per raw frame) and _fullreg.png (shift traces, speed).

Differences from the legacy chain, all small for 2-s bins: the shifts are
whole-pixel (legacy: sub-pixel with bilinear interpolation — after averaging
31 frames the difference is a blur below 0.3 px RMS); no rotation (the legacy
imregtform 'rigid' allows it; within one 16-min run it is nil, and prealign
handles it between runs); the estimate uses the whole frame under a Hann
window instead of cropping [80 80 40 40] px off the edges; the reference is
the still stretch of this run (legacy: of the concatenated day), trimmed
like DetermineReference. Not done by either: non-rigid warps, z-motion.
"""
from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import numpy as np

from .motion import _Matcher, whole_shift

REF_N = 300            # raw frames in the reference
REF_ROUNDS = 2         # align-and-average rounds
CROP_PCT = 99.5        # |shift| percentile cropped off every side


def _still_stretch(beh, N: int, T: int, min_bins: int = 6, trim_bins: int = 2):
    """(first raw frame, last raw frame) of the longest still stretch in the
    binned behaviour csv (still_fraction >= 0.9), or None. A stretch long
    enough loses `trim_bins` bins at each end (the legacy DetermineReference
    trims 5 s: the animal settling / starting is not reference material)."""
    if beh is None or "still_fraction" not in beh.columns:
        return None
    sf = beh["still_fraction"].to_numpy(float)
    still = sf >= 0.9
    best, cur_start = None, None
    for i, s in enumerate(np.r_[still, False]):
        if s and cur_start is None:
            cur_start = i
        elif not s and cur_start is not None:
            if best is None or i - cur_start > best[1] - best[0]:
                best = (cur_start, i)
            cur_start = None
    if best is None or best[1] - best[0] < min_bins:
        return None
    b0, b1 = best
    if b1 - b0 >= min_bins + 2 * trim_bins:
        b0 += trim_bins; b1 -= trim_bins
    a, b = b0 * N, min(b1 * N, T)
    return int(a), int(b)


def build_reference(frames: list, max_shift: float = 40.0, down: int = 2, sigma: float = 1.0, hp: float = 30.0,
                    rounds: int = REF_ROUNDS) -> tuple[np.ndarray, dict]:
    """Mean of the frames, then each frame aligned to it (whole-pixel) and
    averaged again, `rounds` times. Returns (reference float32, info)."""
    stack = np.stack([np.asarray(f, np.float32) for f in frames])
    ref = stack.mean(axis=0)
    info = {"n_frames": int(stack.shape[0]), "rounds": []}
    for r in range(rounds):
        M = _Matcher(ref, max_shift=max_shift, sigma=sigma, hp=hp, down=down)
        acc = np.zeros_like(ref); n = 0; corrs = []
        for f in stack:
            dy, dx, c = M.shift(f)
            corrs.append(c)
            acc += shift_frame(f, dy, dx, fill=ref); n += 1
        ref = acc / max(n, 1)
        info["rounds"].append({"corr_median": float(np.median(corrs)), "corr_min": float(np.min(corrs))})
    return ref.astype(np.float32), info


def shift_frame(frame: np.ndarray, dy: float, dx: float, fill=None) -> np.ndarray:
    """Undo a measured shift (the content sits dy rows / dx columns from the
    reference) by a whole-pixel translation; the uncovered border takes
    `fill` (an image of the same shape) or 0."""
    dy, dx = int(round(dy)), int(round(dx))
    H, W = frame.shape
    out = np.array(fill, np.float32, copy=True) if fill is not None else np.zeros((H, W), np.float32)
    ys, ye = max(0, dy), min(H, H + dy)          # rows of the source that land inside
    xs, xe = max(0, dx), min(W, W + dx)
    if ye > ys and xe > xs:
        out[ys - dy:ye - dy, xs - dx:xe - dx] = frame[ys:ye, xs:xe]
    return out


def register_bin_sbx(sbx_path: Path, N: int, out_dir: Path, chan: int = 1, tag: str | None = None,
                     max_shift: float = 40.0, down: int = 2, sigma: float = 1.0, hp: float = 30.0,
                     smooth_time: float = 1.0, smooth_shifts: int = 0, crop_pct: float = CROP_PCT,
                     behavior_csv: Path | None = None, global_fix: bool = True, overwrite: bool = False, log=print) -> Path:
    """Register every raw frame of channel `chan` of the .sbx and write the
    N-frame binned movie into out_dir (named like unpack_bin's output).
    smooth_time: Gaussian sigma (frames) applied along time to the frames the
    shifts are ESTIMATED on (dim raw frames; the suite2p route used ~2 for
    this data) — the frames that are shifted and binned stay unsmoothed.
    Returns the tiff path."""
    import tifffile
    from sbxreader import sbx_memmap, sbx_get_metadata
    sbx_path = Path(sbx_path); out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    meta = sbx_get_metadata(str(sbx_path))
    fs = float(meta["frame_rate"]); nchan = int(meta["num_channels"])
    if chan > nchan:
        raise ValueError(f"{sbx_path.name}: channel {chan} but the file has {nchan}")
    sec = N / fs
    tag = tag or f"bin{N}f_{sec:.2f}s".replace(".", "p")
    stem = f"{sbx_path.stem}_{tag}"
    out_tif = out_dir / f"{stem}.tif"
    if out_tif.exists() and not overwrite:
        log(f"  exists, skipping: {out_tif.name}")
        return out_tif
    mm = sbx_memmap(str(sbx_path))
    T, H, W = int(mm.shape[0]), int(mm.shape[3]), int(mm.shape[4])
    n_out = T // N
    if n_out == 0:
        raise ValueError(f"{sbx_path.name}: only {T} frames, fewer than the bin of {N}")
    t_start = time.time()
    # ---- behaviour (for the still stretch and the figure)
    beh = None
    if behavior_csv is None:
        cands = sorted(sbx_path.parent.glob(f"{sbx_path.stem}_{tag}_behavior.csv")) or sorted(sbx_path.parent.glob(f"{sbx_path.stem}_bin*_behavior.csv"))
        behavior_csv = cands[0] if cands else None
    if behavior_csv is not None and Path(behavior_csv).exists():
        try:
            import pandas as pd
            beh = pd.read_csv(behavior_csv)
        except Exception:
            beh = None
    # ---- reference from the stillest stretch
    stretch = _still_stretch(beh, N, T)
    if stretch is None:
        mid = T // 2; stretch = (max(0, mid - 20 * N), min(T, mid + 20 * N))
        why = "middle of the run (no still stretch found in the behaviour csv)"
    else:
        why = f"longest still stretch of the behaviour csv ({(stretch[1] - stretch[0]) / fs:.0f} s)"
    idx = np.unique(np.linspace(stretch[0], stretch[1] - 1, min(REF_N, stretch[1] - stretch[0])).astype(int))
    log(f"  {sbx_path.name}: {T} frames @ {fs:.2f} Hz, {H}x{W}; reference = {len(idx)} frames from the {why}")
    frames = [np.asarray(mm[int(i), 0, chan - 1], np.float32) for i in idx]
    ref, ref_info = build_reference(frames, max_shift=max_shift, down=down, sigma=sigma, hp=hp)
    del frames
    M = _Matcher(ref, max_shift=max_shift, sigma=sigma, hp=hp, down=down)
    # ---- one pass: shifts, whole-pixel correction, binning
    gb = n_out * H * W * 2 / 1e9
    use_imagej = gb < 3.9
    kw = dict(imagej=True, metadata={"axes": "TYX", "finterval": sec}) if use_imagej else dict(bigtiff=True)
    part = out_tif.with_name(out_tif.name + ".partial")
    dst = tifffile.memmap(str(part), shape=(n_out, H, W), dtype=np.uint16, **kw)
    dy_all = np.zeros(T, np.float32); dx_all = np.zeros(T, np.float32); cc_all = np.zeros(T, np.float32)
    acc_proj = np.zeros((H, W), np.float64)
    bins_per_chunk = max(1, 248 // N)          # ~250 raw frames per chunk keeps the memory under 1 GB
    t_next = time.time()
    for b0 in range(0, n_out, bins_per_chunk):
        b1 = min(b0 + bins_per_chunk, n_out)
        raw = np.asarray(mm[b0 * N:b1 * N, 0, chan - 1]).astype(np.float32)
        est = raw
        if smooth_time and smooth_time > 0 and raw.shape[0] > 2:
            from scipy.ndimage import gaussian_filter1d
            est = gaussian_filter1d(raw, float(smooth_time), axis=0, mode="nearest")
        sh = np.zeros((raw.shape[0], 3), np.float32)
        for k in range(raw.shape[0]):
            fp = M.prep(est[k])
            sh[k] = whole_shift(M, fp)
        del est
        if smooth_shifts and smooth_shifts > 1:
            from scipy.ndimage import median_filter
            sh[:, 0] = median_filter(sh[:, 0], smooth_shifts, mode="nearest")
            sh[:, 1] = median_filter(sh[:, 1], smooth_shifts, mode="nearest")
        dy_all[b0 * N:b1 * N] = sh[:, 0]; dx_all[b0 * N:b1 * N] = sh[:, 1]; cc_all[b0 * N:b1 * N] = sh[:, 2]
        reg = np.stack([shift_frame(raw[k], sh[k, 0], sh[k, 1], fill=ref) for k in range(raw.shape[0])])
        binned = reg.reshape(b1 - b0, N, H, W).mean(axis=1)
        dst[b0:b1] = np.clip(np.round(binned), 0, 65535).astype(np.uint16)
        acc_proj += binned.sum(axis=0)
        if time.time() - t_next > 20:
            log(f"    {b1}/{n_out} bins ({time.time() - t_start:.0f} s)")
            t_next = time.time()
    dst.flush(); del dst; del mm
    # ---- crop the filled border
    amp = np.hypot(dy_all, dx_all)
    crop = int(np.ceil(np.percentile(np.abs(np.r_[dy_all, dx_all]), crop_pct)))
    crop = max(crop, 0)
    arr = tifffile.memmap(str(part), mode="r")
    if crop > 0:
        final = out_tif.with_name(out_tif.name + ".crop")
        d2 = tifffile.memmap(str(final), shape=(n_out, H - 2 * crop, W - 2 * crop), dtype=np.uint16, **kw)
        for b0 in range(0, n_out, 64):
            b1 = min(b0 + 64, n_out)
            d2[b0:b1] = arr[b0:b1, crop:H - crop, crop:W - crop]
        d2.flush(); del d2; del arr
        part.unlink(); final.replace(out_tif)
    else:
        del arr
        part.replace(out_tif)
    Hc, Wc = H - 2 * crop, W - 2 * crop
    proj = np.clip(np.round(acc_proj / n_out), 0, 65535).astype(np.uint16)[crop:H - crop, crop:W - crop]
    tifffile.imwrite(str(out_dir / f"{stem}_meanproj.tif"), proj, imagej=True)
    # ---- sidecars
    still_raw = None
    if beh is not None and "still_fraction" in beh.columns:
        sf = beh["still_fraction"].to_numpy(float)[:n_out]
        still_raw = np.repeat(sf >= 0.5, N)[:T]
    stats = {"shift_median_px": float(np.median(amp)), "shift_p95_px": float(np.percentile(amp, 95)), "shift_max_px": float(amp.max()),
             "dy_range": [float(dy_all.min()), float(dy_all.max())], "dx_range": [float(dx_all.min()), float(dx_all.max())],
             "corr_median": float(np.median(cc_all)), "corr_p05": float(np.percentile(cc_all, 5)),
             "jump_p95_px": float(np.percentile(np.hypot(np.diff(dy_all), np.diff(dx_all)), 95)),
             "within_bin_range_median_px": float(np.median(np.maximum(np.ptp(dy_all[:n_out * N].reshape(n_out, N), axis=1),
                                                                      np.ptp(dx_all[:n_out * N].reshape(n_out, N), axis=1))))}
    if still_raw is not None and still_raw.size == T and (~still_raw).sum() >= 10 and still_raw.sum() >= 10:
        stats["shift_sd_still_px"] = float(np.hypot(dy_all[still_raw].std(), dx_all[still_raw].std()))
        stats["shift_sd_running_px"] = float(np.hypot(dy_all[~still_raw].std(), dx_all[~still_raw].std()))
    info = {"source_sbx": str(sbx_path), "channel": chan, "channels_in_file": nchan, "frames_in": int(T), "fs_in_hz": fs,
            "bin_frames": N, "seconds_per_frame": sec, "fs_out_hz": fs / N, "frames_out": int(n_out),
            "crop_LRTB_px": [crop, crop, crop, crop], "height": Hc, "width": Wc, "dtype": "uint16",
            "values": "mean of (65535 - raw) over each bin, every raw frame shifted to the reference first",
            "registration": {"method": "rigid whole-pixel translation per raw frame (phase correlation, 2x downsampled search, parabolic sub-pixel peak, rounded)",
                             "reference": f"{len(idx)} raw frames {int(idx[0])}-{int(idx[-1])} from the {why}, {REF_ROUNDS} align-and-average rounds",
                             "reference_rounds": ref_info["rounds"], "smooth_time_sigma_frames": float(smooth_time), "smooth_shifts_frames": int(smooth_shifts),
                             "crop_in_raw_px": crop, "filled_border": "reference image", "shifts_csv": f"{stem}_fullreg_shifts.csv",
                             "stats": stats, "seconds": round(time.time() - t_start, 1)}}
    (out_dir / f"{stem}.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    import pandas as pd
    pd.DataFrame({"frame": np.arange(T), "time_s": np.arange(T) / fs, "dy": dy_all, "dx": dx_all, "corr": cc_all,
                  "running": (~still_raw).astype(int) if still_raw is not None and still_raw.size == T else -1}
                 ).to_csv(out_dir / f"{stem}_fullreg_shifts.csv", index=False, float_format="%.3f")
    if behavior_csv is not None and Path(behavior_csv).exists():
        shutil.copy2(behavior_csv, out_dir / f"{stem}_behavior.csv")
    try:
        _shift_figure(dy_all, dx_all, cc_all, fs, beh, N, out_dir / f"{stem}_fullreg.png", stem, stats)
    except Exception as e:
        log(f"    figure failed ({e})")
    log(f"  wrote {out_tif.name} ({n_out} bins, {Hc}x{Wc} after cropping {crop} px; shift median {stats['shift_median_px']:.1f} px, "
        f"max {stats['shift_max_px']:.1f}, within-bin range median {stats['within_bin_range_median_px']:.1f} px; {time.time() - t_start:.0f} s)")
    if global_fix:
        try:
            from .globalfix import apply as global_apply
            global_apply(out_tif, spf=sec, log=log)
        except Exception as e:
            log(f"  global correction FAILED ({e}) — the tiff is the plain registered, binned movie")
    return out_tif


def _shift_figure(dy, dx, cc, fs, beh, N, out_png, title, stats):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.arange(dy.size) / fs / 60.0
    fig, axes = plt.subplots(3, 1, figsize=(11, 5.2), sharex=True, gridspec_kw=dict(hspace=0.12, left=0.07, right=0.99, top=0.9, bottom=0.1))
    axes[0].plot(t, dy, lw=0.5, color="#1f77b4", label="dy (rows)"); axes[0].plot(t, dx, lw=0.5, color="#ff7f0e", label="dx (columns)")
    axes[0].set_ylabel("shift of the raw\nframe (px)", fontsize=8); axes[0].legend(fontsize=7, loc="upper right", frameon=False, ncol=2)
    axes[1].plot(t, cc, lw=0.5, color="#555"); axes[1].set_ylabel("correlation to\nthe reference", fontsize=8)
    if beh is not None and "velocity_mean" in beh.columns:
        tb = (np.arange(len(beh)) + 0.5) * N / fs / 60.0
        axes[2].fill_between(tb, 0, np.abs(beh["velocity_mean"].to_numpy(float)), color="#1f77b4", alpha=0.5, lw=0, step="mid")
    axes[2].set_ylabel("running\nspeed", fontsize=8); axes[2].set_xlabel("minutes", fontsize=8)
    fig.suptitle(f"{title} — full-rate registration: shift median {stats['shift_median_px']:.1f} px, max {stats['shift_max_px']:.1f}, "
                 f"within-bin range median {stats['within_bin_range_median_px']:.1f} px", fontsize=9, x=0.01, ha="left")
    fig.savefig(out_png, dpi=100); plt.close(fig)
