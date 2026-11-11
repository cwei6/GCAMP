"""Condition comparison (pre/post, or across days) at ROI level AND afferent level.

Pipeline (suite2p or manual ROIs):
  1. neuropil-corrected F  ->  moving-percentile dF/F  ->  z-score
  2. auto-curation of ROIs (suite2p's classifier is soma-trained and rejects
     axons; we keep ROIs with skewed traces / detected events instead)
  3. ROIs -> AFFERENTS: hierarchical clustering of ROIs by activity
     correlation (Blaeser et al., eLife 2024: r > 0.7 during quiet
     wakefulness); fragments of one axon share a trace, so they merge
  4. per-afferent trace = npix-weighted average of member ROIs (better SNR)
  5. per condition (pre / post), optionally still-gated:
       events/min, mean peak dF/F, mean peak dF (absolute), fraction of time
       active, deconvolved spike rate (suite2p spks), baseline F0
  6. summary with paired within-FOV statistics + group-level guidance

Event definition: z > event_z for >= event_min_dur_s with peak dF/F >= event_min_dff.
"""
from __future__ import annotations
import csv
from pathlib import Path

import numpy as np
from scipy.ndimage import percentile_filter
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import skew as sp_skew, wilcoxon

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------- basics

def moving_percentile_baseline(F: np.ndarray, fs: float, win_s: float, q: float,
                               stride: int | None = None) -> np.ndarray:
    """Moving percentile baseline. Computed on a temporally subsampled trace
    (~2 Hz) and linearly interpolated back: a 10th percentile over 32 s does not
    change on a sub-second scale, and this is ~50x faster (4000 ROIs x 60k
    frames: 20 s instead of 17 min)."""
    F = np.asarray(F, np.float32)
    T = F.shape[1]
    if stride is None:
        stride = max(1, int(fs / 2.0))
    if stride > 1 and T > 4 * stride:
        sub = F[:, ::stride]
        win = max(3, int(round(win_s * fs / stride)))
        if win % 2 == 0:
            win += 1
        base_sub = percentile_filter(sub, q, size=(1, win), mode="nearest")
        t_sub = np.arange(0, T, stride, dtype=np.float32)
        t = np.arange(T, dtype=np.float32)
        out = np.empty_like(F)
        for r in range(F.shape[0]):
            out[r] = np.interp(t, t_sub, base_sub[r])
        return out
    win = max(3, int(round(win_s * fs)))
    if win % 2 == 0:
        win += 1
    return percentile_filter(F, q, size=(1, win), mode="nearest")


def dff_and_z(Fc: np.ndarray, fs: float, base_win: float, base_q: float):
    Fc = np.asarray(Fc, np.float32)                     # float32: thousands of ROIs x 100k frames
    Fc = Fc - min(0.0, float(Fc.min())) + 1.0          # guard non-positive baselines
    F0 = moving_percentile_baseline(Fc, fs, base_win, base_q)
    dff = (Fc - F0) / F0
    med = np.median(dff, axis=1, keepdims=True)
    # robust SD (MAD) so events do not inflate their own threshold
    mad = np.median(np.abs(dff - med), axis=1, keepdims=True) * 1.4826 + 1e-9
    z = (dff - med) / mad
    return Fc, F0, dff, z


def detect_events(z: np.ndarray, dff: np.ndarray, fs: float,
                  z_thresh: float, min_dur_s: float, min_dff: float,
                  valid: np.ndarray | None = None):
    """Return list of (start, stop) index pairs for one trace."""
    above = z > z_thresh
    if valid is not None:
        above &= valid
    edges = np.flatnonzero(np.diff(np.r_[0, above.astype(int), 0]))
    events = []
    min_len = int(round(min_dur_s * fs))
    for s, e in zip(edges[::2], edges[1::2]):
        if (e - s) >= min_len and np.nanmax(dff[s:e]) >= min_dff:
            events.append((int(s), int(e)))
    return events


# ----------------------------------------------------------------------------- loading

def load_rois(plane0: Path):
    F = np.load(plane0 / "F.npy")
    Fneu = np.load(plane0 / "Fneu.npy")
    spks = np.load(plane0 / "spks.npy") if (plane0 / "spks.npy").exists() else None
    stat = (np.load(plane0 / "stat.npy", allow_pickle=True)
            if (plane0 / "stat.npy").exists() else None)
    iscell = np.load(plane0 / "iscell.npy") if (plane0 / "iscell.npy").exists() \
        else np.ones((F.shape[0], 2), np.float32)
    names = None
    if (plane0 / "roi_names.txt").exists():
        names = (plane0 / "roi_names.txt").read_text(encoding="utf-8").splitlines()
    return F, Fneu, spks, stat, iscell, names


def red_mean_image(plane0: Path):
    """suite2p's mean image of channel 2 (vascular tracer), if it was recorded."""
    try:
        ops = np.load(plane0 / "ops.npy", allow_pickle=True).item()
        img = ops.get("meanImg_chan2")
        return np.asarray(img, np.float32) if img is not None else None
    except Exception:
        return None


def auto_curate(dff: np.ndarray, z: np.ndarray, fs: float, a: dict, **kw) -> np.ndarray:
    """Keep ROIs whose traces look like real calcium transients (see curate.py)."""
    from .curate import curate
    return curate(dff, fs, a, **kw)[0]


def curate_only(plane0: Path, out_dir: Path, fs: float, cfg: dict) -> np.ndarray:
    """Re-run ROI curation alone (fast) and overwrite iscell.npy — for tuning
    analysis.curate_* on real data without repeating the whole analysis."""
    from .curate import curate
    out_dir.mkdir(parents=True, exist_ok=True)
    a = cfg.get("analysis", {})
    F, Fneu, spks, stat, iscell, names = load_rois(plane0)
    print(f"curate: dF/F for {F.shape[0]} ROIs x {F.shape[1]} frames ...", flush=True)
    _, _, dff, _ = dff_and_z(F - float(a.get("neuropil_coeff", 0.7)) * Fneu, fs,
                             float(a.get("baseline_win_s", 32.0)), float(a.get("baseline_percentile", 10)))
    redcell = np.load(plane0 / "redcell.npy") if (plane0 / "redcell.npy").exists() else None
    cls_path = plane0 / "iscell_classifier.npy"
    if not cls_path.exists() and (plane0 / "iscell.npy").exists():
        np.save(cls_path, iscell)
    keep, score, _ = curate(dff, fs, a, stat=stat, F=F, Fneu=Fneu, redcell=redcell,
                            report_dir=out_dir, names=names, red_img=red_mean_image(plane0))
    new = np.column_stack([keep.astype(np.float32), score.astype(np.float32)])
    np.save(plane0 / "iscell.npy", new); np.save(plane0 / "iscell_auto.npy", new)
    print(f"curate: kept {int(keep.sum())}/{len(keep)} ROIs -> {plane0 / 'iscell.npy'} "
          f"(report: {out_dir / 'roi_curation.png'})")
    return keep


# ----------------------------------------------------------------------------- clustering

def _corr_matrix(X: np.ndarray) -> np.ndarray:
    X = X - X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True) + 1e-9
    R = (X / sd) @ (X / sd).T / X.shape[1]
    np.fill_diagonal(R, 1.0)
    return np.clip(R, -1, 1)


def cluster_into_afferents(dff: np.ndarray, valid: np.ndarray, fs: float,
                           r_thresh="auto", r_min: float = 0.3, smooth_s: float = 0.5,
                           min_valid_s: float = 30.0, n_null: int = 200,
                           centroids: np.ndarray | None = None, max_dist_px: float = 0.0,
                           seed: int = 0):
    """ROIs -> afferents by activity correlation.

    Fragments of one axon carry the same calcium signal but independent noise,
    so their correlation is limited by per-ROI SNR. To make this robust:
      * traces are Gaussian-smoothed (smooth_s) before correlating — noise
        averages out, 0.5-2 s transients do not
      * threshold 'auto' = max(r_min, 99.9th percentile of a circular-shift
        null of all pairs), so it adapts to the recording's SNR; or a number
      * optional spatial gate: ROIs further apart than max_dist_px never merge
    Hierarchical clustering, average linkage on 1 - r.
    Returns labels (1..K), R, threshold used, null 99.9th percentile.
    """
    n, T = dff.shape
    if n == 1:
        return np.array([1]), np.ones((1, 1)), 1.0, np.nan, np.array([])
    use = valid if valid.sum() >= min_valid_s * fs else np.ones(T, bool)
    if use is not valid:
        print(f"  cluster: only {valid.sum() / fs:.0f} s of still data — using all frames")
    X = dff[:, use]
    if smooth_s > 0:
        from scipy.ndimage import gaussian_filter1d
        X = gaussian_filter1d(X, smooth_s * fs, axis=1)
    R = _corr_matrix(X)

    # circular-shift null: what correlation do unrelated traces reach by chance?
    rng = np.random.default_rng(seed)
    null = []
    Xn = X - X.mean(axis=1, keepdims=True); Xn /= (Xn.std(axis=1, keepdims=True) + 1e-9)
    for _ in range(n_null):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        s = rng.integers(int(0.1 * X.shape[1]), int(0.9 * X.shape[1]))
        null.append(float(Xn[i] @ np.roll(Xn[j], s) / X.shape[1]))
    null99 = float(np.percentile(null, 99.9)) if null else np.nan
    if r_thresh in (None, "auto"):
        thr = max(float(r_min), null99 if np.isfinite(null99) else r_min)
    else:
        thr = float(r_thresh)

    D = 1.0 - R
    if centroids is not None and max_dist_px and max_dist_px > 0:
        d = np.sqrt(((centroids[:, None, :] - centroids[None, :, :]) ** 2).sum(-1))
        D[d > max_dist_px] = 2.0
    D = (D + D.T) / 2
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    labels = fcluster(Z, t=1.0 - thr, criterion="distance")
    return labels, R, thr, null99, np.asarray(null)


# ----------------------------------------------------------------------------- metrics

def condition_metrics(dff, z, Fc, F0, spks, fs, cond_mask, valid, a,
                      vel=None, pupil=None):
    """Per-trace metrics for one condition (dicts of arrays over traces)."""
    z_thresh = float(a.get("event_z", 1.0)); min_dur = float(a.get("event_min_dur_s", 1.0))
    min_dff = float(a.get("event_min_dff", 0.05))
    m = cond_mask & valid
    minutes = m.sum() / fs / 60.0
    n = dff.shape[0]
    out = {k: np.full(n, np.nan) for k in
           ("events_per_min", "mean_peak_dff", "mean_peak_dF", "active_frac",
            "spks_per_min", "F0", "r_pupil", "r_loco", "pupil_median")}
    out["minutes"] = np.full(n, round(minutes, 2))
    if minutes == 0:
        return out, [[] for _ in range(n)]
    # coupling to behaviour: correlation of 0.5 s-smoothed dF/F with pupil size
    # and running speed over the WHOLE condition (not still-gated — the
    # coupling needs the behavioural variance), lag 0
    from scipy.ndimage import gaussian_filter1d
    cm = cond_mask
    sm = gaussian_filter1d(dff, 0.5 * fs, axis=1)
    for name, beh in (("r_pupil", pupil), ("r_loco", vel)):
        if beh is None:
            continue
        ok = cm & np.isfinite(beh)
        if ok.sum() < 30 * fs or np.nanstd(beh[ok]) == 0:
            continue
        b = beh[ok] - beh[ok].mean(); b /= (b.std() + 1e-9)
        X = sm[:, ok]; X = X - X.mean(1, keepdims=True); X /= (X.std(1, keepdims=True) + 1e-9)
        out[name] = X @ b / ok.sum()
    if pupil is not None and np.isfinite(pupil[m]).any():
        out["pupil_median"][:] = np.nanmedian(pupil[m])
    all_events = []
    for r in range(n):
        evts = detect_events(z[r], dff[r], fs, z_thresh, min_dur, min_dff, valid=m)
        all_events.append(evts)
        out["events_per_min"][r] = len(evts) / minutes
        if evts:
            out["mean_peak_dff"][r] = np.mean([np.nanmax(dff[r, s:e]) for s, e in evts])
            out["mean_peak_dF"][r] = np.mean([np.nanmax(Fc[r, s:e] - F0[r, s:e]) for s, e in evts])
            out["active_frac"][r] = sum(e - s for s, e in evts) / m.sum()
        else:
            out["mean_peak_dff"][r] = 0.0; out["mean_peak_dF"][r] = 0.0; out["active_frac"][r] = 0.0
        out["F0"][r] = np.median(F0[r, cond_mask])
        if spks is not None:
            out["spks_per_min"][r] = spks[r, m].sum() / minutes
    return out, all_events


def _write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def _fmt(v, nd=4):
    return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else round(float(v), nd)


# ----------------------------------------------------------------------------- main

def analyze_prepost(
    plane0: Path,
    out_dir: Path,
    run_conditions: list[str],
    run_nframes: list[int],
    fs: float,
    loco: list[dict] | None,
    cfg: dict,
    mean_image: np.ndarray | None = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    a = cfg.get("analysis", {})
    neu_coeff = float(a.get("neuropil_coeff", 0.7))
    base_q = float(a.get("baseline_percentile", 10))
    base_win = float(a.get("baseline_win_s", 32.0))
    gate = a.get("loco_gate", "still")
    r_cfg = a.get("cluster_r", "auto")
    r_thresh = "auto" if str(r_cfg).lower() == "auto" else float(r_cfg)
    r_min = float(a.get("cluster_r_min", 0.3))
    smooth_s = float(a.get("cluster_smooth_s", 0.5))
    max_dist = float(a.get("cluster_max_dist_px", 0) or 0)
    min_valid_s = float(a.get("cluster_min_still_s", 30.0))

    F, Fneu, spks, stat, iscell, names = load_rois(plane0)
    n_all = F.shape[0]
    print(f"analysis: dF/F for {n_all} ROIs x {F.shape[1]} frames ...", flush=True)
    Fc_all, F0_all, dff_all, z_all = dff_and_z(F - neu_coeff * Fneu, fs, base_win, base_q)

    # ---- curation
    # suite2p's built-in classifier is trained on somata and rejects axonal
    # ROIs, so its iscell.npy is not usable here. Rule:
    #   * first run: curate by trace statistics (skew / detected events), write
    #     iscell.npy and keep a copy as iscell_auto.npy
    #   * later runs: if iscell.npy still equals iscell_auto.npy nobody has
    #     curated in the GUI -> recompute (settings may have changed);
    #     if it differs, the user curated in the GUI -> respect that
    accepted = iscell[:, 0].astype(bool)
    auto_path, cls_path = plane0 / "iscell_auto.npy", plane0 / "iscell_classifier.npy"
    user_curated = auto_path.exists() and not np.array_equal(
        np.load(auto_path)[:, 0], iscell[:, 0])
    if bool(a.get("auto_curate", True)) and not user_curated:
        if not cls_path.exists() and (plane0 / "iscell.npy").exists():
            np.save(cls_path, iscell)                      # keep the classifier's verdict
        from .curate import curate
        redcell = (np.load(plane0 / "redcell.npy") if (plane0 / "redcell.npy").exists() else None)
        accepted, score, _ = curate(dff_all, fs, a, stat=stat, F=F, Fneu=Fneu, redcell=redcell,
                                    report_dir=out_dir, names=names, red_img=red_mean_image(plane0))
        # score -> iscell[:, 1] so the suite2p GUI probability slider (0.5) reproduces this
        new = np.column_stack([accepted.astype(np.float32), score.astype(np.float32)])
        np.save(plane0 / "iscell.npy", new); np.save(auto_path, new)
        print(f"analysis: suite2p classifier accepted {int(iscell[:, 0].sum())}/{n_all} ROIs "
              f"(soma-trained, not meaningful for axons); response-profile curation kept "
              f"{int(accepted.sum())}/{n_all} — see {out_dir / 'roi_curation.png'} / .csv / .txt; "
              "review in the suite2p GUI if you like (GUI edits are detected and take "
              "precedence on re-runs)")
    elif user_curated:
        print(f"analysis: using GUI-curated ROIs ({int(accepted.sum())}/{n_all})")
    if accepted.sum() == 0:
        accepted[:] = True
        print("analysis: WARNING no ROI passed curation — using all")
    keep = np.flatnonzero(accepted)
    Fc, F0, dff, z = Fc_all[keep], F0_all[keep], dff_all[keep], z_all[keep]
    spks_k = spks[keep] if spks is not None else None
    npix = np.array([int(stat[i]["npix"]) for i in keep]) if stat is not None else np.ones(len(keep))

    # ---- frame masks
    bounds = np.r_[0, np.cumsum(run_nframes)]
    total = bounds[-1]
    # conditions = distinct labels in run order (pre/post, or dates d260911, ...);
    # the baseline is analysis.baseline_condition, else 'pre' if present, else
    # the first condition; every other condition is compared against it
    present = list(dict.fromkeys(run_conditions))
    base = a.get("baseline_condition") or ("pre" if "pre" in present else present[0])
    if base not in present:
        print(f"analysis: baseline_condition {base!r} not among {present} — using {present[0]!r}")
        base = present[0]
    present = [base] + [c for c in present if c != base]
    others = present[1:]
    cond_mask = {c: np.zeros(total, bool) for c in present}
    still_all = np.ones(total, bool)
    for i, cond in enumerate(run_conditions):
        cond_mask[cond][bounds[i]:bounds[i + 1]] = True
        if loco is not None:
            still_all[bounds[i]:bounds[i + 1]] = loco[i]["still"][: run_nframes[i]]
    valid = still_all if (gate == "still" and loco is not None) else np.ones(total, bool)

    # behaviour vectors on the concatenated timeline (for coupling / gating)
    vel_all = np.full(total, np.nan); pupil_all = np.full(total, np.nan)
    if loco is not None:
        for i in range(len(run_conditions)):
            sl = slice(bounds[i], bounds[i + 1])
            vel_all[sl] = np.abs(loco[i]["velocity"][: run_nframes[i]])
            if loco[i].get("pupil") is not None:
                pupil_all[sl] = loco[i]["pupil"][: run_nframes[i]]
    have_pupil = np.isfinite(pupil_all).sum() > 100
    pupil_gate = str(a.get("pupil_gate", "none")).lower()
    pupil_note = ""
    if have_pupil and pupil_gate == "match" and len(present) >= 2:
        # compare conditions at matched arousal: keep frames whose pupil lies in
        # the overlap of all conditions' central (10th-90th pct) ranges
        lo = max(np.nanpercentile(pupil_all[cond_mask[c] & valid], 10) for c in present)
        hi = min(np.nanpercentile(pupil_all[cond_mask[c] & valid], 90) for c in present)
        if hi > lo:
            in_range = (pupil_all >= lo) & (pupil_all <= hi)
            valid = valid & in_range
            pupil_note = f"pupil-matched gating: frames with pupil in [{lo:.1f}, {hi:.1f}] px"
        else:
            pupil_note = "pupil-matched gating requested but the conditions' pupil ranges do not overlap — NOT applied"
    elif pupil_gate == "match" and not have_pupil:
        pupil_note = "pupil_gate: match requested but no pupil data found"

    # ---- ROIs -> afferents
    centroids = (np.array([[float(np.mean(stat[i]["ypix"])), float(np.mean(stat[i]["xpix"]))]
                           for i in keep]) if stat is not None else None)
    labels, R, thr_used, null99, null_r = cluster_into_afferents(
        dff, valid, fs, r_thresh=r_thresh, r_min=r_min, smooth_s=smooth_s,
        min_valid_s=min_valid_s, centroids=centroids, max_dist_px=max_dist)
    K = int(labels.max())
    aff_F = np.zeros((K, total)); aff_spks = np.zeros((K, total)) if spks_k is not None else None
    members = []
    for k in range(1, K + 1):
        idx = np.flatnonzero(labels == k)
        members.append(idx)
        w = npix[idx] / npix[idx].sum()
        aff_F[k - 1] = (Fc[idx] * w[:, None]).sum(axis=0)
        if aff_spks is not None:
            aff_spks[k - 1] = (spks_k[idx] * w[:, None]).sum(axis=0)
    aff_Fc, aff_F0, aff_dff, aff_z = dff_and_z(aff_F, fs, base_win, base_q)
    print(f"analysis: {len(keep)} ROIs -> {K} afferents (r > {thr_used:.2f}"
          f"{' auto' if r_thresh == 'auto' else ''}; shuffle-null 99.9% = {null99:.2f}; "
          f"{sum(len(m) > 1 for m in members)} multi-ROI clusters)")

    # ---- metrics
    beh_kw = dict(vel=vel_all if loco is not None else None,
                  pupil=pupil_all if have_pupil else None)
    roi_metrics = {c: condition_metrics(dff, z, Fc, F0, spks_k, fs, cond_mask[c], valid, a, **beh_kw)[0]
                   for c in present}
    aff_metrics, aff_events = {}, {}
    for c in present:
        aff_metrics[c], aff_events[c] = condition_metrics(
            aff_dff, aff_z, aff_Fc, aff_F0, aff_spks, fs, cond_mask[c], valid, a, **beh_kw)

    keys = ("minutes", "events_per_min", "mean_peak_dff", "mean_peak_dF",
            "active_frac", "spks_per_min", "F0", "r_pupil", "r_loco", "pupil_median")
    roi_rows = []
    for r in range(len(keep)):
        row = {"roi": int(keep[r]), "afferent": int(labels[r]), "npix": int(npix[r])}
        if names:
            row["name"] = names[keep[r]] if keep[r] < len(names) else ""
        for c in present:
            for k in keys:
                row[f"{c}_{k}"] = _fmt(roi_metrics[c][k][r])
        roi_rows.append(row)
    _write_csv(out_dir / "prepost_per_roi.csv", roi_rows)

    aff_rows = []
    for k in range(K):
        row = {"afferent": k + 1, "n_rois": len(members[k]),
               "rois": ";".join(str(int(keep[i])) for i in members[k]),
               "total_npix": int(npix[members[k]].sum())}
        for c in present:
            for kk in keys:
                row[f"{c}_{kk}"] = _fmt(aff_metrics[c][kk][k])
        for c in others:
            pre_r, post_r = aff_metrics[base]["events_per_min"][k], aff_metrics[c]["events_per_min"][k]
            row[f"{c}_vs_{base}_rate_change_per_min"] = _fmt(post_r - pre_r)
            row[f"{c}_vs_{base}_rate_ratio"] = _fmt((post_r + 0.01) / (pre_r + 0.01), 3)
            f0p, f0q = aff_metrics[base]["F0"][k], aff_metrics[c]["F0"][k]
            row[f"{c}_vs_{base}_F0_ratio"] = _fmt(f0q / f0p if f0p else np.nan, 3)
        aff_rows.append(row)
    _write_csv(out_dir / "prepost_per_afferent.csv", aff_rows)
    frame_condition = np.concatenate([np.full(f, c) for c, f in zip(run_conditions, run_nframes)])
    np.savez(out_dir / "afferent_traces.npz", dff=aff_dff, F0=aff_F0, labels=labels,
             keep=keep, fs=fs, bounds=bounds, valid=valid,
             run_conditions=np.array(run_conditions), run_nframes=np.array(run_nframes),
             frame_condition=frame_condition)

    # ---- summary
    lines = ["afferent2p per-afferent summary", "=" * 44,
             f"ROIs analyzed: {len(keep)} / {n_all}   ->  afferents: {K}",
             f"Gating: {gate}" + ("" if loco is not None else " (no loco data — ungated)"),
             f"Cluster rule: pairwise r > {thr_used:.2f} ({'auto: max(r_min, shuffle null)' if r_thresh == 'auto' else 'fixed'}; "
             f"null 99.9% = {null99:.2f}) on {smooth_s:g}s-smoothed dF/F during "
             + ("still frames" if valid.sum() < total else "all frames"),
             f"  {sum(len(m) > 1 for m in members)} afferents merge >1 ROI; see roi_correlation.png "
             "— if observed and null r overlap, fragments cannot be merged at this SNR"]
    if pupil_note:
        lines.append(pupil_note)
    if have_pupil:
        for c in present:
            mm = cond_mask[c] & valid
            lines.append(f"[{c}] pupil median {np.nanmedian(pupil_all[mm]):.1f} px in analysed frames "
                         f"(whole run {np.nanmedian(pupil_all[cond_mask[c]]):.1f} px); "
                         f"afferent coupling median r(pupil) = {np.nanmedian(aff_metrics[c]['r_pupil']):.2f}, "
                         f"r(|vel|) = {np.nanmedian(aff_metrics[c]['r_loco']):.2f}")
        pp = np.nanmedian(pupil_all[cond_mask[base] & valid])
        for c in others:
            pq = np.nanmedian(pupil_all[cond_mask[c] & valid])
            if pp > 0 and abs(pq - pp) / pp > 0.15:
                lines.append(f"  NOTE: pupil differs {c} vs {base} by {100 * (pq - pp) / pp:+.0f}% — arousal "
                             "state differs; consider analysis.pupil_gate: match")
    for c in present:
        m = aff_metrics[c]
        lines.append(f"[{c}] {m['minutes'][0]:.1f} min analysed | median events/min "
                     f"{np.nanmedian(m['events_per_min']):.3f} | median peak dF/F "
                     f"{np.nanmedian(m['mean_peak_dff']):.3f} | median spks/min "
                     f"{np.nanmedian(m['spks_per_min']):.2f} | median F0 {np.nanmedian(m['F0']):.0f}")
    if others:
        # afferents with events in only ONE condition: real silencing/appearance,
        # or a session-specific ROI from imperfect cross-session alignment
        ev = np.array([[aff_metrics[c]["events_per_min"][k] for c in present] for k in range(K)], float)
        one_only = int(np.sum((np.nan_to_num(ev) > 0).sum(1) == 1))
        lines.append(f"Baseline condition: {base}")
        lines.append(f"Afferents with events in only one condition: {one_only} / {K}"
                     + ("   <- if many, check qc/prealign.png and qc_report.txt: session-specific ROIs "
                        "from imperfect alignment inflate this" if one_only > K / 3 else ""))
        for c in others:
            pre_r = aff_metrics[base]["events_per_min"]; post_r = aff_metrics[c]["events_per_min"]
            ok = np.isfinite(pre_r) & np.isfinite(post_r)
            up = int(np.sum(post_r[ok] > pre_r[ok])); down = int(np.sum(post_r[ok] < pre_r[ok]))
            line = f"[{c} vs {base}] afferents with rate increase: {up} / {int(ok.sum())} (decrease: {down})"
            if ok.sum() >= 5 and np.any(post_r[ok] != pre_r[ok]):
                try:
                    pv = wilcoxon(post_r[ok], pre_r[ok]).pvalue
                    line += f"; Wilcoxon signed-rank (within-FOV, exploratory) p = {pv:.3g}"
                except ValueError:
                    pass
            f0r = np.array([r.get(f"{c}_vs_{base}_F0_ratio") or np.nan for r in aff_rows], float)
            line += f"; median F0 ratio {np.nanmedian(f0r):.3f}"
            lines.append(line)
        lines += ["  (F0 ratio < 0.9 = bleaching; then prefer events/min, spks/min and peak dF "
                  "over peak dF/F for amplitude)",
                  "", "Group statistics: one value per animal/FOV (e.g. median afferent "
                  "rate change vs baseline), treatment vs control mice."]
    else:
        lines += ["", "Single-condition run — no comparison made."]
    (out_dir / "prepost_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))

    # ---- figures
    _figures(out_dir, present, aff_metrics, aff_dff, aff_events, labels, members, keep,
             stat, mean_image, R, bounds, fs, valid, thr_used, null_r)
    return out_dir / "prepost_per_afferent.csv"


def _figures(out_dir, present, aff_metrics, aff_dff, aff_events, labels, members, keep,
             stat, mean_image, R, bounds, fs, valid, thr_used=0.7, null_r=None):
    K = len(members)
    cmap = plt.get_cmap("tab20")

    # 1. cluster map on mean image
    if stat is not None:
        fig, ax = plt.subplots(figsize=(10, 7))
        if mean_image is not None:
            lo, hi = np.percentile(mean_image, [0.5, 99.8])
            ax.imshow(np.clip((mean_image - lo) / max(hi - lo, 1e-9), 0, 1), cmap="gray")
        else:
            Ly = max(int(np.max(s["ypix"])) for s in stat) + 1
            Lx = max(int(np.max(s["xpix"])) for s in stat) + 1
            ax.imshow(np.zeros((Ly, Lx)), cmap="gray")
        for k, idx in enumerate(members):
            col = cmap(k % 20)
            for i in idx:
                s = stat[keep[i]]
                ax.scatter(s["xpix"], s["ypix"], s=0.4, color=col, alpha=0.6)
            s0 = stat[keep[idx[0]]]
            ax.text(np.mean(s0["xpix"]), np.mean(s0["ypix"]), str(k + 1), color=col,
                    fontsize=8, weight="bold")
        ax.set_title(f"{len(keep)} ROIs grouped into {K} afferents (colour = afferent)")
        ax.axis("off"); fig.tight_layout()
        fig.savefig(out_dir / "afferent_map.png", dpi=150); plt.close(fig)

    # 2. correlation matrix sorted by afferent + observed-vs-null histogram
    order = np.argsort(labels, kind="stable")
    fig, (ax, ah) = plt.subplots(1, 2, figsize=(12, 5))
    im = ax.imshow(R[np.ix_(order, order)], vmin=-0.2, vmax=1, cmap="viridis")
    edges = np.flatnonzero(np.diff(labels[order])) + 0.5
    for e in edges:
        ax.axhline(e, color="w", lw=0.5); ax.axvline(e, color="w", lw=0.5)
    ax.set_title("ROI-ROI correlation (smoothed dF/F), sorted by afferent")
    fig.colorbar(im, ax=ax, label="r")
    iu = np.triu_indices(R.shape[0], k=1)
    obs = R[iu]
    if obs.size:
        bins = np.linspace(-0.3, 1, 53)
        ah.hist(obs, bins=bins, alpha=0.7, label=f"observed pairs (n={obs.size})", density=True)
        if null_r is not None and len(null_r):
            ah.hist(null_r, bins=bins, alpha=0.5, label="circular-shift null", density=True)
        ah.axvline(thr_used, color="k", ls="--", label=f"merge threshold {thr_used:.2f}")
        ah.set_xlabel("pairwise r"); ah.set_ylabel("density"); ah.legend(fontsize=8)
        ah.set_title("Do fragments separate from chance? (pairs right of the line merge)")
    fig.tight_layout()
    fig.savefig(out_dir / "roi_correlation.png", dpi=150); plt.close(fig)

    # 3. per-afferent paired plot across conditions + traces
    nc = len(present)
    ncol = 3 if nc >= 2 else 1
    fig, axes = plt.subplots(1, ncol, figsize=(5.5 * ncol, 4.8))
    axes = np.atleast_1d(axes)
    if nc >= 2:
        xs = np.arange(nc)
        for ax, key, lab in zip(axes[:2], ("events_per_min", "mean_peak_dff"),
                                ("events / min", "mean peak dF/F")):
            for k in range(K):
                ys = [aff_metrics[c][key][k] for c in present]
                ax.plot(xs, ys, "-o", color=cmap(k % 20), ms=4, lw=1, alpha=0.8)
            ax.set_xticks(xs); ax.set_xticklabels(present, rotation=20 if nc > 2 else 0)
            ax.set_xlim(-0.4, nc - 0.6)
            ax.set_ylabel(lab); ax.set_title(f"per afferent (n={K})")
        tax = axes[2]
    else:
        tax = axes[0]
    t = np.arange(aff_dff.shape[1]) / fs
    off = 0.0
    for k in range(min(K, 12)):
        tr = aff_dff[k]
        tax.plot(t, tr + off, lw=0.4, color=cmap(k % 20))
        for c in present:
            for s, e in aff_events[c][k]:
                tax.plot(t[s:e], tr[s:e] + off, lw=1.2, color=cmap(k % 20))
        tax.text(t[0], off + 0.05, f"aff {k + 1}", fontsize=7, color=cmap(k % 20))
        off += max(0.6, float(np.nanpercentile(tr, 99.8)) * 1.2)
    for b in bounds[1:-1]:
        tax.axvline(b / fs, color="k", ls="--", lw=0.8)
    if valid.sum() < len(valid):
        tax.fill_between(t, -0.3, off, where=~valid, color="orange", alpha=0.08, lw=0)
    tax.set_xlabel("time (s)"); tax.set_ylabel("dF/F (offset)")
    tax.set_title("afferent traces (thick = detected events; shaded = moving)")
    fig.tight_layout()
    fig.savefig(out_dir / "prepost_overview.png", dpi=150); plt.close(fig)
