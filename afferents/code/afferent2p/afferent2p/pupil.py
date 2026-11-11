"""Pupil from DeepLabCut output (<run>_eye_*DLC_*.csv), aligned to imaging frames.

Format (pinned to the lab's files): DLC CSV with 3 header rows (scorer,
bodyparts, coords); bodyparts pupil_1..pupil_8 = points around the pupil edge
(45 deg apart), pupil_9 = centre; one row per eye-camera frame, and the eye
camera is triggered by the 2P frame clock, so row i = imaging frame i (same as
the quadrature file; verified on a paired recording: pupil vs running r = 0.65,
pupil lagging running by ~0.8 s).

The model was designed so that the centre point carries the geometry: per frame
  diam_px   : 2 x median distance from the centre point to the edge points with
              likelihood >= lik_min (needs the centre and >= min_edge_points
              edge points; occluded/mistracked edge points simply drop out)
  area_px   : pi * (diam/2)^2
  ell_major/ell_minor_px : ellipse fit through the valid edge points when >= 5
              are available (gaze-angle-robust diameter; supplementary)
  cx, cy    : centre point (gaze proxy)
  valid     : usable frame
Frames failing the point requirements are blinks/occlusions -> NaN. Short gaps
(<= max_gap_s) are linearly interpolated; long gaps stay NaN.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np


# ----------------------------------------------------------------------------- ellipse fit

def _fit_ellipse(x: np.ndarray, y: np.ndarray):
    """Direct least-squares ellipse fit (Halir & Flusser). Returns
    (cx, cy, a, b) with a >= b semi-axes, or None if degenerate."""
    mx, my = x.mean(), y.mean()
    xs, ys = x - mx, y - my
    D1 = np.column_stack([xs * xs, xs * ys, ys * ys])
    D2 = np.column_stack([xs, ys, np.ones_like(xs)])
    S1, S2, S3 = D1.T @ D1, D1.T @ D2, D2.T @ D2
    try:
        T = -np.linalg.solve(S3, S2.T)
        M = S1 + S2 @ T
        M = np.array([M[2] / 2, -M[1], M[0] / 2])
        w, v = np.linalg.eig(M)
        cond = 4 * v[0] * v[2] - v[1] ** 2
        idx = np.flatnonzero(np.real(cond) > 0)
        if idx.size == 0:
            return None
        a1 = np.real(v[:, idx[0]])
        A, B, C = a1
        D, E, F = T @ a1
    except np.linalg.LinAlgError:
        return None
    den = B * B - 4 * A * C
    if den >= 0:
        return None
    cx = (2 * C * D - B * E) / den
    cy = (2 * A * E - B * D) / den
    num = 2 * (A * E * E + C * D * D - B * D * E + den * F)
    root = np.sqrt((A - C) ** 2 + B * B)
    try:
        a = -np.sqrt(num * (A + C + root)) / den
        b = -np.sqrt(num * (A + C - root)) / den
    except FloatingPointError:
        return None
    if not (np.isfinite(a) and np.isfinite(b)) or a <= 0 or b <= 0:
        return None
    a, b = max(a, b), min(a, b)
    return cx + mx, cy + my, a, b


# ----------------------------------------------------------------------------- loading

def find_pupil_file(sbx_path: Path, pattern: str = "*eye*DLC*.csv") -> Path | None:
    cands = sorted(sbx_path.parent.glob(pattern))
    stem = [p for p in cands if sbx_path.stem in p.name]
    cands = stem or cands
    filt = [p for p in cands if "filtered" in p.name]
    cands = filt or cands
    return cands[0] if cands else None


def load_dlc_csv(path: Path):
    """Returns (X, Y, L) arrays of shape (frames, points) and the bodypart names."""
    import csv
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) < 4 or rows[1][0] != "bodyparts":
        raise ValueError(f"{path.name} does not look like a DeepLabCut CSV")
    bps_row, coords_row = rows[1][1:], rows[2][1:]
    data = np.array([[float(v) if v else np.nan for v in r[1:]] for r in rows[3:]], np.float64)
    bps = list(dict.fromkeys(bps_row))
    X = np.column_stack([data[:, i] for i, (b, c) in enumerate(zip(bps_row, coords_row)) if c == "x"])
    Y = np.column_stack([data[:, i] for i, (b, c) in enumerate(zip(bps_row, coords_row)) if c == "y"])
    L = np.column_stack([data[:, i] for i, (b, c) in enumerate(zip(bps_row, coords_row)) if c == "likelihood"])
    return X, Y, L, bps


def _guess_center(X, Y, L, bps):
    """The centre point: named like one, else the point nearest the median
    position of all points (pupil_9 in the lab's 9-point model)."""
    for i, b in enumerate(bps):
        if any(k in b.lower() for k in ("center", "centre", "middle")):
            return i
    mx, my = np.nanmedian(X), np.nanmedian(Y)
    d = np.sqrt((np.nanmedian(X, 0) - mx) ** 2 + (np.nanmedian(Y, 0) - my) ** 2)
    return int(np.argmin(d)) if len(bps) >= 5 else None


def pupil_from_dlc(X, Y, L, bps, lik_min: float = 0.6, edge_points=None,
                   center_point=None, min_edge_points: int = 3):
    """Per-frame pupil metrics from DLC landmarks (centre + edge points)."""
    n_pts = X.shape[1]
    if center_point is None:
        center_point = _guess_center(X, Y, L, bps)
    if edge_points is None:
        edge_points = [i for i in range(n_pts) if i != center_point]
    T = X.shape[0]
    out = {k: np.full(T, np.nan) for k in
           ("diam_px", "area_px", "ell_major_px", "ell_minor_px", "cx", "cy")}
    n_valid = np.zeros(T, int)
    edge = np.array(edge_points)
    okL = L[:, edge] >= lik_min
    if center_point is not None:
        okC = (L[:, center_point] >= lik_min) & np.isfinite(X[:, center_point])
        cx, cy = X[:, center_point], Y[:, center_point]
    else:   # no centre tracked: use the centroid of valid edge points
        okC = np.ones(T, bool)
        cx = np.nanmean(np.where(okL, X[:, edge], np.nan), 1)
        cy = np.nanmean(np.where(okL, Y[:, edge], np.nan), 1)
    r = np.sqrt((X[:, edge] - cx[:, None]) ** 2 + (Y[:, edge] - cy[:, None]) ** 2)
    r = np.where(okL, r, np.nan)
    n_valid = okL.sum(1)
    good = okC & (n_valid >= min_edge_points)
    import warnings
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rad = np.nanmedian(r, 1)
    out["diam_px"][good] = 2 * rad[good]
    out["area_px"][good] = np.pi * rad[good] ** 2
    out["cx"][good] = cx[good]; out["cy"][good] = cy[good]
    # supplementary ellipse fit where enough edge points are valid
    with np.errstate(all="ignore"):
        for t in np.flatnonzero(good & (n_valid >= 5)):
            ok = edge[okL[t]]
            fit = _fit_ellipse(X[t, ok], Y[t, ok])
            if fit is not None:
                out["ell_major_px"][t] = 2 * fit[2]; out["ell_minor_px"][t] = 2 * fit[3]
    return out, n_valid, edge_points, center_point


def _clean(diam: np.ndarray, fs: float, max_gap_s: float, jump_mad: float = 5.0):
    """Reject implausible jumps, interpolate short gaps."""
    d = diam.copy()
    # rolling-median outlier rejection (blinks that DLC still 'fit')
    from scipy.ndimage import median_filter
    win = max(3, int(round(5 * fs)) | 1)
    filled = np.where(np.isfinite(d), d, np.nanmedian(d))
    med = median_filter(filled, size=win, mode="nearest")
    resid = d - med
    mad = np.nanmedian(np.abs(resid)) * 1.4826 + 1e-9
    d[np.abs(resid) > jump_mad * mad] = np.nan
    valid = np.isfinite(d)
    # interpolate short gaps only
    idx = np.arange(d.size)
    if valid.sum() >= 2:
        interp = np.interp(idx, idx[valid], d[valid])
        gap_len = np.zeros(d.size, int)
        edges = np.flatnonzero(np.diff(np.r_[1, valid.astype(int), 1]))
        for s, e in zip(edges[::2], edges[1::2]):
            gap_len[s:e] = e - s
        short = (~valid) & (gap_len <= max(1, int(round(max_gap_s * fs))))
        d[short] = interp[short]
    return d, valid


def pupil_for_run(sbx_path: Path, nframes: int, fs: float, pattern: str = "*eye*DLC*.csv",
                  lik_min: float = 0.6, max_gap_s: float = 0.5, smooth_s: float = 0.25,
                  min_edge_points: int = 3, center_point: int | None = None) -> dict | None:
    """Pupil metrics aligned to the imaging frames of one run, or None if no file."""
    f = find_pupil_file(sbx_path, pattern)
    if f is None:
        print(f"  pupil: no DLC eye file for {sbx_path.name} (pattern {pattern})")
        return None
    X, Y, L, bps = load_dlc_csv(f)
    raw, n_valid, edge, center = pupil_from_dlc(X, Y, L, bps, lik_min=lik_min,
                                                 center_point=center_point,
                                                 min_edge_points=min_edge_points)
    n = X.shape[0]
    if n != nframes:
        print(f"  pupil: {n} eye frames for {nframes} imaging frames ({n / nframes:.2f} per frame) — resampling")
        xo, xn = np.linspace(0, 1, n), np.linspace(0, 1, nframes)
        for k in raw:
            v = raw[k]; ok = np.isfinite(v)
            raw[k] = np.interp(xn, xo[ok], v[ok]) if ok.sum() > 1 else np.full(nframes, np.nan)
            # keep NaNs where the source was NaN
            raw[k][np.interp(xn, xo, (~ok).astype(float)) > 0.5] = np.nan
    diam, valid = _clean(raw["diam_px"], fs, max_gap_s)
    if smooth_s > 0:
        from scipy.ndimage import gaussian_filter1d
        filled = np.where(np.isfinite(diam), diam, np.nanmedian(diam))
        diam_s = gaussian_filter1d(filled, smooth_s * fs)
        diam_s[~np.isfinite(diam)] = np.nan
    else:
        diam_s = diam
    frac = 100 * np.isfinite(diam_s).mean()
    print(f"  pupil: {f.name}: {len(bps)} points (edge {[bps[i] for i in edge]}, "
          f"centre {bps[center] if center is not None else 'none'}); diameter median "
          f"{np.nanmedian(diam_s):.1f} px ({np.nanpercentile(diam_s, 5):.0f}-{np.nanpercentile(diam_s, 95):.0f}), "
          f"{frac:.1f}% frames valid")
    return {
        "file": str(f), "diam_px": diam_s, "diam_raw_px": raw["diam_px"],
        "area_px": raw["area_px"], "ell_major_px": raw["ell_major_px"],
        "ell_minor_px": raw["ell_minor_px"],
        "cx": raw["cx"], "cy": raw["cy"], "valid": np.isfinite(diam_s),
        "n_valid_points": n_valid if n == nframes else None,
        "center_point": bps[center] if center is not None else None,
    }
