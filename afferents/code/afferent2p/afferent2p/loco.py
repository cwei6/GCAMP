"""Locomotion (wheel quadrature) loading and still/moving state detection.

NOTE: the exact .quad file format on your rig should be verified against one
real file — this loader handles the common cases (a MATLAB .mat containing one
numeric vector, e.g. 'quad_data', or a raw binary vector) and will tell you
what it found. If it mis-reads your files, send one .quad file and the loader
can be pinned to the exact format.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np


def _load_quad_vector(path: Path) -> np.ndarray:
    """Best-effort read of a quadrature counts vector."""
    # 1) MATLAB file?
    try:
        from scipy.io import loadmat
        d = loadmat(str(path), squeeze_me=True)
        arrays = {
            k: np.asarray(v).ravel()
            for k, v in d.items()
            if not k.startswith("__") and np.asarray(v).size > 10
            and np.issubdtype(np.asarray(v).dtype, np.number)
        }
        if arrays:
            # prefer a variable with 'quad' in the name
            for k in arrays:
                if "quad" in k.lower():
                    print(f"  loco: read '{k}' ({arrays[k].size} samples) from {path.name}")
                    return arrays[k].astype(np.float64)
            k = max(arrays, key=lambda k: arrays[k].size)
            print(f"  loco: read '{k}' ({arrays[k].size} samples) from {path.name}")
            return arrays[k].astype(np.float64)
    except Exception:
        pass
    # 2) raw binary fallbacks
    for dtype in (np.int32, np.float64, np.int16, np.uint32):
        try:
            v = np.fromfile(str(path), dtype=dtype)
            if v.size > 10 and np.all(np.isfinite(v.astype(np.float64))):
                print(f"  loco: read raw {np.dtype(dtype).name} vector "
                      f"({v.size} samples) from {path.name} — VERIFY this is correct")
                return v.astype(np.float64)
        except Exception:
            continue
    raise ValueError(f"Could not parse quadrature file: {path}")


def find_quad_file(sbx_path: Path, pattern: str = "*quad*") -> Path | None:
    cands = [p for p in sbx_path.parent.glob(pattern) if p.suffix != ".sbx"]
    # prefer files sharing the run's stem
    stem_matches = [p for p in cands if sbx_path.stem in p.name]
    cands = stem_matches or cands
    return cands[0] if cands else None


def loco_for_run(
    sbx_path: Path,
    nframes: int,
    fs: float,
    quad_glob: str = "*quad*",
    counts_per_cm: float | None = None,
    still_thresh: float = 1.0,
    min_still_s: float = 5.0,
    smooth_s: float = 0.5,
    quad_kind: str = "cumulative",
) -> dict:
    """Return velocity per imaging frame and a boolean 'still' mask.

    If no quad file is found, everything is treated as still=True with a warning
    (so the pipeline still runs, but state-gating is a no-op for that run).
    """
    quad_path = find_quad_file(sbx_path, quad_glob)
    if quad_path is None:
        print(f"  loco: WARNING no quadrature file for {sbx_path.name}; "
              "treating whole run as still")
        return {
            "velocity": np.zeros(nframes),
            "still": np.ones(nframes, bool),
            "quad_file": None,
            "calibrated": False,
        }

    counts = _load_quad_vector(quad_path)
    # Scanbox quadrature files (<run>_quadrature.mat, variable 'quad_data',
    # int32 cumulative encoder position) may hold a different number of samples
    # than imaging frames; resample the cumulative position onto frames
    if counts.size != nframes:
        ratio = counts.size / nframes
        print(f"  loco: {counts.size} quad samples for {nframes} frames "
              f"({ratio:.2f} per frame) — resampling onto frames")
        x_old = np.linspace(0, 1, counts.size)
        x_new = np.linspace(0, 1, nframes)
        counts = np.interp(x_new, x_old, counts)

    # velocity from the counts vector.
    # quad_kind: 'cumulative' (default; encoder position sampled per frame ->
    # differentiate) or 'rate' (file already stores per-frame speed).
    if quad_kind == "rate":
        vel = counts.astype(float)
    else:
        vel = np.gradient(counts) * fs
    if counts_per_cm:
        vel = vel / float(counts_per_cm)

    # smooth
    w = max(1, int(round(smooth_s * fs)))
    kernel = np.ones(w) / w
    vel_s = np.convolve(vel, kernel, mode="same")

    still = np.abs(vel_s) < still_thresh
    # drop still bouts shorter than min_still_s
    min_len = int(round(min_still_s * fs))
    still_clean = still.copy()
    edges = np.flatnonzero(np.diff(np.r_[0, still.astype(int), 0]))
    for s, e in zip(edges[::2], edges[1::2]):
        if (e - s) < min_len:
            still_clean[s:e] = False

    frac = 100 * still_clean.mean()
    print(f"  loco: {sbx_path.name}: {frac:.0f}% of frames in still bouts "
          f">= {min_still_s}s (thresh {still_thresh})")
    return {
        "velocity": vel_s,
        "still": still_clean,
        "quad_file": str(quad_path),
        "calibrated": counts_per_cm is not None,
    }
