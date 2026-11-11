"""Compare the Python-registered movie against a MATLAB-pipeline output.

    python run_pipeline.py <config> --steps compare --ref <matlab file>

<matlab file> can be a 2-D projection (e.g. the mean image UnpackSBX or
WriteSbxProjection wrote) or a 3-D registered movie tif (its mean is used).
Produces qc/compare_matlab.png and qc/compare_matlab.txt with:
  * Pearson correlation of the two mean images after removing any global
    x/y offset (so a different crop or a constant shift does not count)
  * a sharpness score for each (variance of the Laplacian / mean^2) —
    better registration = sharper mean image
Correlation >~0.9 and similar-or-better sharpness means the two pipelines
did the same job.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import laplace, gaussian_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_ref(path: Path) -> np.ndarray:
    if path.suffix.lower() in (".png", ".jpg", ".jpeg"):
        img = plt.imread(str(path))
        if img.ndim == 3:
            img = img[..., :3].mean(axis=-1)
        return img.astype(np.float32)
    with tifffile.TiffFile(str(path)) as tf:
        n = len(tf.pages)
        if n == 1:
            arr = tf.pages[0].asarray()
            if arr.ndim == 3:            # RGB projection
                arr = arr.mean(axis=-1)
            return arr.astype(np.float32)
        step = max(1, n // 2000)
        acc, cnt = None, 0
        for i in range(0, n, step):
            fr = tf.pages[i].asarray().astype(np.float32)
            if fr.ndim == 3:
                fr = fr.mean(axis=-1)
            acc = fr if acc is None else acc + fr
            cnt += 1
        return acc / cnt


def _norm(img):
    lo, hi = np.percentile(img, [0.5, 99.7])
    return np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1)


def _zscore(img):
    img = img - img.mean()
    return img / (img.std() + 1e-9)


def _pearson_at(a: np.ndarray, b: np.ndarray, dy: int, dx: int, margin: int = 20) -> float:
    bs = np.roll(np.roll(b, dy, axis=0), dx, axis=1)
    aa, bb = a[margin:-margin, margin:-margin], bs[margin:-margin, margin:-margin]
    return float(np.corrcoef(_zscore(aa).ravel(), _zscore(bb).ravel())[0, 1])


def _best_offset(a: np.ndarray, b: np.ndarray, max_shift: int = 80, k: int = 8):
    """Integer (dy, dx) that best aligns b onto a. Candidates come from FFT
    cross-correlation of high-passed images; the winner is the candidate with
    the highest Pearson r (robust to periodic structure / multiple peaks)."""
    ha = _zscore(a - gaussian_filter(a, 20))
    hb = _zscore(b - gaussian_filter(b, 20))
    cc = np.fft.irfft2(np.fft.rfft2(ha) * np.conj(np.fft.rfft2(hb)), s=a.shape)
    cc = np.fft.fftshift(cc)
    cy, cx = a.shape[0] // 2, a.shape[1] // 2
    my, mx = min(cy, max_shift), min(cx, max_shift)
    win = cc[cy - my:cy + my + 1, cx - mx:cx + mx + 1]
    flat = np.argsort(win.ravel())[::-1][:k]
    cands = {(0, 0)}
    for f in flat:
        iy, ix = np.unravel_index(f, win.shape)
        cands.add((int(iy - my), int(ix - mx)))
    best = max(cands, key=lambda c: _pearson_at(a, b, *c))
    return best


def _center_crop(img, Ly, Lx):
    y0 = (img.shape[0] - Ly) // 2
    x0 = (img.shape[1] - Lx) // 2
    return img[y0:y0 + Ly, x0:x0 + Lx]


def _sharpness(img):
    s = gaussian_filter(img, 1.0)
    return float(laplace(s).var() / (s.mean() ** 2 + 1e-9))


def compare_to_reference(plane0: Path, out_dir: Path, total_frames: int,
                         Ly: int, Lx: int, ref_path: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    mov = np.memmap(plane0 / "data.bin", dtype=np.int16, mode="r",
                    shape=(total_frames, Ly, Lx))
    step = max(1, total_frames // 3000)
    py_mean = np.asarray(mov[::step]).astype(np.float32).mean(axis=0)
    del mov

    ref = _load_ref(Path(ref_path))
    print(f"compare: python mean {py_mean.shape}, reference {ref.shape} from {ref_path}")

    # bring to a common size: assume any size difference is a symmetric crop
    # (MATLAB projections are cropped by projParam.edge) or a 2x downsample
    if ref.shape != py_mean.shape:
        ry = py_mean.shape[0] / ref.shape[0]
        if 1.8 < ry < 2.2 and 1.8 < py_mean.shape[1] / ref.shape[1] < 2.2:
            py_mean = py_mean[:py_mean.shape[0] // 2 * 2, :py_mean.shape[1] // 2 * 2]
            py_mean = py_mean.reshape(py_mean.shape[0] // 2, 2, py_mean.shape[1] // 2, 2).mean((1, 3))
            print("  reference looks 2x downsampled — matched")
        Ly2, Lx2 = min(py_mean.shape[0], ref.shape[0]), min(py_mean.shape[1], ref.shape[1])
        py_c, ref_c = _center_crop(py_mean, Ly2, Lx2), _center_crop(ref, Ly2, Lx2)
    else:
        py_c, ref_c = py_mean, ref

    dy, dx = _best_offset(py_c, ref_c)
    ref_s = np.roll(np.roll(ref_c, dy, axis=0), dx, axis=1)
    m = 20
    a, b = py_c[m:-m, m:-m], ref_s[m:-m, m:-m]
    r = float(np.corrcoef(_zscore(a).ravel(), _zscore(b).ravel())[0, 1])
    sh_py, sh_ref = _sharpness(a), _sharpness(b)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    axes[0].imshow(_norm(a), cmap="gray"); axes[0].set_title("afferent2p (suite2p) registered mean")
    axes[1].imshow(_norm(b), cmap="gray"); axes[1].set_title(f"reference (MATLAB): {Path(ref_path).name}")
    rgb = np.dstack([_norm(a), _norm(b), _norm(a)])
    axes[2].imshow(rgb); axes[2].set_title("overlay (magenta=python, green=MATLAB)\nwhite = agreement")
    for ax in axes:
        ax.axis("off")
    fig.suptitle(f"r = {r:.3f} after offset (dy={dy}, dx={dx})   |   sharpness python {sh_py:.2e} vs MATLAB {sh_ref:.2e}")
    fig.tight_layout()
    fig.savefig(out_dir / "compare_matlab.png", dpi=150)
    plt.close(fig)

    verdict = ("SAME FIELD, comparable registration" if r >= 0.9 else
               "similar but check overlay" if r >= 0.7 else
               "LOW agreement — inspect: different FOV/crop, or one registration failed")
    lines = [
        "afferent2p vs MATLAB comparison", "=" * 40,
        f"reference: {ref_path}",
        f"global offset removed: dy={dy}, dx={dx} px",
        f"Pearson r of mean images: {r:.3f}   -> {verdict}",
        f"sharpness (Laplacian var / mean^2): python {sh_py:.3e}, MATLAB {sh_ref:.3e}",
        "  (higher = crisper mean = better motion correction; a RAW, unregistered",
        "   MATLAB projection should be LESS sharp than a registered one)",
    ]
    (out_dir / "compare_matlab.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
