"""Afferent ROI curation: keep suite2p ROIs whose traces look like real calcium
transients, reject background speckle / noise ROIs.

suite2p's ROI detection on axonal FOVs returns thousands of ROIs, most of them
small patches of background whose "signal" is photon noise. Its built-in
classifier is soma-trained and useless here, so we classify from the RESPONSE
PROFILE of each ROI (plus a few spatial sanity checks), in the spirit of how
AQuA decides that an event is real (noise estimated from frame-to-frame
differences, peak significance relative to that noise, minimum duration):

  n_sig_events number of significant transients: 0.5-s smoothed z > z_on for
               >= min duration, with peak smoothed z >= z_peak (AQuA-like
               significance: noise-scaled peak + minimum duration).  Rejects
               smooth-but-eventless traces; lenient by default (>= 1).
  asymmetry    area after / area before the peak (1 s each) of the raw event-
               triggered average — a
               calcium transient rises fast and decays slowly (jGCaMP8s), a noise
               or motion bump is symmetric.  Rejects motion/vessel artefacts and
               slow background excursions that pass the event test.
  skew         skewness of dF/F (transients are positive-going); weak on its own
  snr          signal sd of the smoothed dF/F in excess of the white-noise
               expectation (noise sd from first differences / sqrt(w)); reported,
               off by default
  autocorr     lag-1 autocorrelation of the high-passed F-0.7Fneu (suite2p's
               stat 'snr'): the fraction of variance that is temporally smooth;
               photon noise ~0, real afferents 0.1-0.7 on real data.  THE main
               criterion; threshold 'auto' sits 4 robust SDs above the noise bulk.
  red_ratio    ROI intensity in the red (tracer) mean image / local annulus:
               > 2 = the ROI sits on a vessel (vessel walls make fake transients)
  npix         minimum ROI size (single-pixel speckle)
  red          optional: reject ROIs suite2p flags as red (vascular tracer) —
               vessel walls produce motion-driven "transients"

Each criterion c gives a soft score s_c = clip(value / threshold, 0, 2) / 2 (0.5
at threshold); the ROI score is min(s_c).  keep <=> score >= 0.5.  The score is
written into iscell.npy[:, 1], so the suite2p GUI's probability slider at 0.5
reproduces the selection, and raising it makes the selection stricter.
"""
from __future__ import annotations
import csv
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.stats import skew as sp_skew

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULTS = dict(
    curate_min_autocorr="auto",  # lag-1 autocorrelation of high-passed F-0.7Fneu (= suite2p 'snr');
                                 # 'auto' = max(0.05, median + 4 robust-SD of all ROIs) -> the noise bulk
    curate_min_sig_events=1,     # significant transients per whole recording (all runs together)
    curate_sig_z_on=3.0,         # smoothed z to start an event
    curate_sig_z_peak=4.0,       # smoothed z the peak must reach
    curate_sig_min_dur_s=0.5,
    curate_smooth_s=0.5,
    curate_min_npix=6,
    curate_max_red_ratio=2.0,    # 0 = off; ROI red (chan2) intensity / local red background -> vessel
    curate_reject_red=False,     # suite2p's redcell.npy flag (unreliable: flagged 100% on real data)
    curate_min_asymmetry=0.0,    # 0 = off; after/before-peak area ratio of the event-triggered average
    curate_asym_min_events=1,
    curate_snr_min=0.0,          # 0 = off; excess signal sd / noise sd after 0.5-s smoothing
    curate_skew_min=0.0,         # 0 = off; skewness of dF/F (Poisson noise is skewed too -> useless on dim FOVs)
    curate_min_aspect=0.0,       # 0 = off  (suite2p stat aspect_ratio)
    curate_max_neuropil_r=0.0,   # 0 = off  (corr of raw F with Fneu)
)


def _robust_sd(x: np.ndarray) -> np.ndarray:
    med = np.median(x, axis=-1, keepdims=True)
    return np.median(np.abs(x - med), axis=-1) * 1.4826 + 1e-12


def _noise_sd(x: np.ndarray) -> np.ndarray:
    """Photon/white noise level from first differences (robust)."""
    return np.median(np.abs(np.diff(x, axis=-1)), axis=-1) * 1.4826 / np.sqrt(2) + 1e-12


def _events_1d(zs: np.ndarray, z_on: float, z_peak: float, min_len: int):
    above = zs > z_on
    edges = np.flatnonzero(np.diff(np.r_[0, above.astype(int), 0]))
    out = []
    for s, e in zip(edges[::2], edges[1::2]):
        if (e - s) >= min_len and zs[s:e].max() >= z_peak:
            out.append((int(s), int(e)))
    return out


def _asymmetry(x: np.ndarray, align: np.ndarray, events, fs: float,
               pre_s: float = 2.0, post_s: float = 4.0, win_s: float = 1.0):
    """Shape of the EVENT-TRIGGERED AVERAGE (raw dF/F, peaks located on the
    lightly smoothed trace `align`).

    Returns (area_asymmetry, rise_half_s, decay_half_s):
      area_asymmetry = (area in the 1 s AFTER the peak) / (area in the 1 s BEFORE)
      A calcium transient rises within a frame or two and decays over hundreds
      of ms to seconds -> ratio >> 1; a noise or motion bump is symmetric -> ~1.
    Averaging the events first gives a far cleaner shape than judging single
    noisy transients. nan if not measurable.
    """
    pre, post = int(round(pre_s * fs)), int(round(post_s * fs))
    n = x.size
    wins = []
    for s, e in events:
        p = s + int(np.argmax(align[s:e]))
        lo, hi = p - pre, p + post
        w = np.full(pre + post, np.nan)
        w[max(0, -lo):pre + post - max(0, hi - n)] = x[max(lo, 0):min(hi, n)]
        wins.append(w)
    if not wins:
        return np.nan, np.nan, np.nan
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)     # all-NaN columns at recording edges
        eta = np.nanmean(np.array(wins), axis=0)
        eta = eta - np.nanmedian(eta[: pre // 2])
    pk = eta[pre]
    if not np.isfinite(pk) or pk <= 0:
        return np.nan, np.nan, np.nan
    w = int(round(win_s * fs))
    after = np.nansum(np.clip(eta[pre + 1:pre + 1 + w], 0, None))
    before = np.nansum(np.clip(eta[pre - w:pre], 0, None))
    ratio = (after + 0.5 * pk) / (before + 0.5 * pk)
    half = 0.5 * pk
    i = pre
    while i > 0 and eta[i] > half:
        i -= 1
    j = pre
    while j < eta.size - 1 and eta[j] > half:
        j += 1
    return float(ratio), float(max(pre - i, 0.5) / fs), float(max(j - pre, 0.5) / fs)


def _highpass_autocorr(Fc: np.ndarray, fs: float, win_s: float = 300.0) -> np.ndarray:
    """Lag-1 autocorrelation of F - 0.7 Fneu after removing a slow (5 min) running
    mean — identical in spirit to suite2p's stat['snr'] but insensitive to
    bleaching and to brightness steps between runs.  Computed on the raw
    corrected fluorescence, NOT on dF/F: the moving-percentile baseline used for
    dF/F adds its own slow variability and inflates the autocorrelation of pure
    noise (0.05 instead of 0.015 on real data)."""
    Fc = np.asarray(Fc, np.float32)
    n, T = Fc.shape
    stride = max(1, int(fs))                      # 1 Hz is plenty for a 5-min trend
    sub = Fc[:, ::stride]
    w = max(3, int(round(win_s * fs / stride)))
    trend_sub = uniform_filter1d(sub, size=w, axis=1, mode="nearest")
    t_sub = np.arange(0, T, stride, dtype=np.float32); t = np.arange(T, dtype=np.float32)
    out = np.zeros(n)
    for r in range(n):
        x = Fc[r] - np.interp(t, t_sub, trend_sub[r])
        v = x.var()
        out[r] = 1.0 - 0.5 * np.diff(x).var() / v if v > 0 else 0.0
    return out


def red_ratio(stat, red_img: np.ndarray, inner_px: int = 3, outer_px: int = 12) -> np.ndarray:
    """ROI mean intensity in the red (chan2, vascular tracer) mean image divided
    by the mean in a surrounding annulus — > ~2 means the ROI sits on a vessel."""
    from scipy.ndimage import binary_dilation
    Ly, Lx = red_img.shape
    out = np.full(len(stat), np.nan)
    st_in = np.ones((2 * inner_px + 1, 2 * inner_px + 1), bool)
    st_out = np.ones((2 * outer_px + 1, 2 * outer_px + 1), bool)
    for i, s in enumerate(stat):
        yp, xp = np.asarray(s["ypix"]), np.asarray(s["xpix"])
        ok = (yp >= 0) & (yp < Ly) & (xp >= 0) & (xp < Lx)
        if ok.sum() == 0:
            continue
        yp, xp = yp[ok], xp[ok]
        y0, y1 = max(0, yp.min() - outer_px), min(Ly, yp.max() + outer_px + 1)
        x0, x1 = max(0, xp.min() - outer_px), min(Lx, xp.max() + outer_px + 1)
        m = np.zeros((y1 - y0, x1 - x0), bool); m[yp - y0, xp - x0] = True
        ring = binary_dilation(m, st_out) & ~binary_dilation(m, st_in)
        crop = red_img[y0:y1, x0:x1]
        bg = crop[ring].mean() if ring.any() else np.nan
        out[i] = crop[m].mean() / bg if bg and bg > 0 else np.nan
    return out


def roi_features(dff: np.ndarray, fs: float, stat=None, F=None, Fneu=None,
                 redcell=None, smooth_s: float = 0.5, z_on: float = 3.0,
                 z_peak: float = 4.0, min_dur_s: float = 0.5, neu_coeff: float = 0.7,
                 red_img: np.ndarray | None = None) -> dict:
    n, T = dff.shape
    w = max(2, int(round(smooth_s * fs)))
    min_len = max(2, int(round(min_dur_s * fs)))
    lag = max(1, int(round(0.2 * fs)))
    snr = np.zeros(n); noise = np.zeros(n); skw = np.zeros(n); peak = np.zeros(n); ac = np.zeros(n); ac1 = np.zeros(n)
    n_ev = np.zeros(n, int); asym = np.full(n, np.nan); ev_dur = np.full(n, np.nan)
    rise = np.full(n, np.nan); decay = np.full(n, np.nan)
    events = []
    # process ROIs in chunks: with thousands of ROIs x 100k frames the smoothed
    # copies would otherwise need many GB
    step = max(1, int(2e8 // max(T, 1)))
    for c0 in range(0, n, step):
        c1 = min(n, c0 + step)
        if n > 500:
            print(f"  curation features: ROI {c0}/{n}", flush=True)
        d = dff[c0:c1].astype(np.float32, copy=False)
        sm = uniform_filter1d(d, size=w, axis=1, mode="nearest")
        sm3 = uniform_filter1d(d, size=3, axis=1, mode="nearest")
        nz = _noise_sd(d)                        # white-noise sd of raw dF/F
        sm_noise = nz / np.sqrt(w)               # what that noise becomes after smoothing
        # signal sd in excess of what white noise alone would leave after smoothing
        snr[c0:c1] = np.sqrt(np.maximum(sm.var(axis=1) - sm_noise ** 2, 0)) / sm_noise
        noise[c0:c1] = nz
        peak[c0:c1] = np.nanmax(sm, axis=1)
        skw[c0:c1] = sp_skew(d, axis=1, nan_policy="omit")
        zs = (sm - np.median(sm, axis=1, keepdims=True)) / _robust_sd(sm)[:, None]
        dm = d - d.mean(axis=1, keepdims=True)
        ac[c0:c1] = (dm[:, :-lag] * dm[:, lag:]).sum(1) / (np.sqrt((dm[:, :-lag] ** 2).sum(1) * (dm[:, lag:] ** 2).sum(1)) + 1e-12)
        # lag-1 autocorrelation == suite2p's stat['snr'] (1 - 0.5 var(diff)/var):
        # fraction of the variance that is temporally smooth; white noise -> ~0
        ac1[c0:c1] = 1.0 - 0.5 * np.diff(d, axis=1).var(axis=1) / (d.var(axis=1) + 1e-12)
        for k in range(c1 - c0):
            r = c0 + k
            ev = _events_1d(zs[k], z_on, z_peak, min_len)
            events.append(ev)
            n_ev[r] = len(ev)
            if ev:
                asym[r], rise[r], decay[r] = _asymmetry(d[k], sm3[k], ev, fs)
                ev_dur[r] = float(np.median([(e - s) / fs for s, e in ev]))
    feats = {
        "snr": snr, "noise_sd": noise, "skew": skw,
        "n_sig_events": n_ev, "events_per_min": n_ev / (T / fs / 60.0),
        "asymmetry": asym, "rise_half_s": rise, "decay_half_s": decay,
        "median_event_dur_s": ev_dur, "autocorr_dff_lag1": ac1, "autocorr_0.2s": ac,
        "peak_dff": peak,
    }
    if F is not None and Fneu is not None:
        feats["autocorr_lag1"] = _highpass_autocorr(F - neu_coeff * Fneu, fs)
    else:
        feats["autocorr_lag1"] = ac1
    if stat is not None:
        for key in ("npix", "aspect_ratio", "compact", "skew", "snr", "radius"):
            if key in stat[0]:
                feats[f"s2p_{key}"] = np.array([float(s.get(key, np.nan)) for s in stat])
    if F is not None and Fneu is not None:
        Fc = F - F.mean(axis=1, keepdims=True); Nc = Fneu - Fneu.mean(axis=1, keepdims=True)
        feats["neuropil_r"] = (Fc * Nc).sum(1) / (np.sqrt((Fc ** 2).sum(1) * (Nc ** 2).sum(1)) + 1e-12)
    if redcell is not None:
        feats["red"] = redcell[:, 0].astype(float)
        feats["red_prob"] = redcell[:, 1].astype(float)
    if red_img is not None and stat is not None:
        feats["red_ratio"] = red_ratio(stat, np.asarray(red_img, np.float32))
    return feats, events


def curate_rois(feats: dict, a: dict):
    """Return keep mask, score (0..1), per-criterion soft scores, thresholds."""
    p = {k: a.get(k, v) for k, v in DEFAULTS.items()}
    n = feats["snr"].size

    def soft_min(v, thr):                 # value must be >= thr
        v = np.nan_to_num(np.asarray(v, float), nan=0.0)
        if thr <= 0:
            return np.ones(n)
        return np.clip(v / thr, 0, 2) / 2

    crit = {
        "snr": soft_min(feats["snr"], float(p["curate_snr_min"])),
        "skew": soft_min(feats["skew"], float(p["curate_skew_min"])),
        "sig_events": soft_min(feats["n_sig_events"], float(p["curate_min_sig_events"])),
    }
    ac = np.nan_to_num(feats["autocorr_lag1"], nan=0.0)
    thr_ac = p["curate_min_autocorr"]
    if isinstance(thr_ac, str) and thr_ac.lower() == "auto":
        med = float(np.median(ac)); rsd = 1.4826 * float(np.median(np.abs(ac - med)))
        thr_ac = max(0.05, med + 4 * rsd)
        p["curate_min_autocorr"] = f"auto -> {thr_ac:.3f} (median {med:.3f} + 4 x robust SD {rsd:.3f})"
    else:
        thr_ac = float(thr_ac)
    p["_autocorr_thr"] = thr_ac
    if thr_ac > 0:
        crit["autocorr"] = soft_min(ac, thr_ac)
    if "s2p_npix" in feats:
        crit["npix"] = soft_min(feats["s2p_npix"], float(p["curate_min_npix"]))
    if float(p["curate_min_aspect"]) > 0 and "s2p_aspect_ratio" in feats:
        crit["aspect"] = soft_min(feats["s2p_aspect_ratio"], float(p["curate_min_aspect"]))
    # asymmetry: only judged when enough events; otherwise neutral (0.5 -> pass)
    asym_thr = float(p["curate_min_asymmetry"])
    if asym_thr > 0:
        s = soft_min(feats["asymmetry"], asym_thr)
        enough = feats["n_sig_events"] >= int(p["curate_asym_min_events"])
        crit["asymmetry"] = np.where(enough, s, np.maximum(s, 0.5))
    max_r = float(p["curate_max_neuropil_r"])
    if max_r > 0 and "neuropil_r" in feats:      # value must be <= thr
        crit["neuropil_r"] = np.where(feats["neuropil_r"] <= max_r, 1.0, 0.0)
    max_rr = float(p["curate_max_red_ratio"])
    if max_rr > 0 and "red_ratio" in feats:
        rr = np.nan_to_num(feats["red_ratio"], nan=1.0)
        crit["not_vessel"] = np.clip(max_rr / np.maximum(rr, 1e-6), 0, 2) / 2   # 0.5 at threshold
    if bool(p["curate_reject_red"]) and "red" in feats:
        frac = float(np.mean(feats["red"] > 0))
        if frac > 0.5:
            print(f"  curate: suite2p flagged {frac:.0%} of ROIs as red — ignoring redcell.npy "
                  "(degenerate); use curate_max_red_ratio instead")
        else:
            crit["not_red"] = np.where(feats["red"] > 0, 0.0, 1.0)
    score = np.min(np.column_stack(list(crit.values())), axis=1)
    keep = score >= 0.5
    return keep, score, crit, p


def write_curation_report(out_dir: Path, feats: dict, crit: dict, score: np.ndarray,
                          keep: np.ndarray, dff: np.ndarray, fs: float, p: dict,
                          events=None, names=None):
    out_dir.mkdir(parents=True, exist_ok=True)
    n = keep.size
    # ---- table
    cols = {"roi": np.arange(n), "keep": keep.astype(int), "score": np.round(score, 3)}
    if names is not None:
        cols["name"] = np.asarray(names)
    for k, v in feats.items():
        cols[k] = np.round(v, 4) if np.issubdtype(np.asarray(v).dtype, np.floating) else v
    for k, v in crit.items():
        cols[f"pass_{k}"] = (v >= 0.5).astype(int)
    with open(out_dir / "roi_curation.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(cols.keys())
        w.writerows(zip(*[np.asarray(v) for v in cols.values()]))
    # ---- text summary
    lines = [f"ROI curation: kept {int(keep.sum())} / {n}", "thresholds used:"]
    lines += [f"  {k}: {v}" for k, v in p.items() if not k.startswith("_")]
    lines.append("ROIs failing each criterion (an ROI can fail several):")
    for k, v in crit.items():
        lines.append(f"  {k:12s} {int((v < 0.5).sum())}")
    if keep.any():
        lines.append("feature medians, kept vs rejected:")
        for k in ("autocorr_lag1", "n_sig_events", "red_ratio", "asymmetry", "snr", "skew", "s2p_npix",
                  "s2p_aspect_ratio", "s2p_compact", "neuropil_r"):
            if k in feats:
                a_, b_ = feats[k][keep], feats[k][~keep]
                lines.append(f"  {k:16s} kept {np.nanmedian(a_):8.3f}   rejected "
                             f"{(np.nanmedian(b_) if b_.size else np.nan):8.3f}")
    lines.append("To loosen/tighten, edit analysis.curate_* in the config and re-run "
                 "--steps analyze; or drag the probability slider in the suite2p GUI "
                 "(0.5 = this selection).")
    (out_dir / "roi_curation.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---- figure: feature scatter + example traces
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(feats["skew"][~keep], feats["snr"][~keep], s=4, c="0.6", label=f"rejected ({int((~keep).sum())})")
    ax.scatter(feats["skew"][keep], feats["snr"][keep], s=8, c="tab:red", label=f"kept ({int(keep.sum())})")
    ax.axhline(p["curate_snr_min"], ls="--", c="k", lw=0.8); ax.axvline(p["curate_skew_min"], ls="--", c="k", lw=0.8)
    ax.set_xlabel("skewness of dF/F"); ax.set_ylabel("snr (smoothed signal / smoothed noise)")
    ax.set_yscale("symlog", linthresh=1); ax.legend(fontsize=8); ax.set_title("trace features")
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(feats["n_sig_events"][~keep] + 0.1, np.nan_to_num(feats["asymmetry"][~keep], nan=0.1), s=4, c="0.6")
    ax.scatter(feats["n_sig_events"][keep] + 0.1, np.nan_to_num(feats["asymmetry"][keep], nan=0.1), s=8, c="tab:red")
    ax.axhline(p["curate_min_asymmetry"], ls="--", c="k", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("significant transients (+0.1)"); ax.set_ylabel("asymmetry (area after/before peak)")
    ax.set_title("transient shape")
    ax = fig.add_subplot(gs[0, 2])
    ac1 = np.nan_to_num(feats["autocorr_lag1"], nan=0.0)
    bins = np.linspace(min(-0.05, ac1.min()), max(0.3, min(1.0, np.percentile(ac1, 99.9) + 0.05)), 60)
    ax.hist(ac1[~keep], bins=bins, color="0.6", label="rejected")
    ax.hist(ac1[keep], bins=bins, color="tab:red", label="kept")
    if p.get("_autocorr_thr", 0) > 0:
        ax.axvline(p["_autocorr_thr"], ls="--", c="k", lw=0.8)
    ax.set_yscale("log"); ax.set_xlabel("lag-1 autocorrelation of F-0.7Fneu (= suite2p 'snr')")
    ax.set_ylabel("ROIs"); ax.set_title("temporal smoothness (white noise ~ 0)"); ax.legend(fontsize=8)

    # example traces: best kept, borderline rejected, worst rejected
    order = np.argsort(-score)
    kept_sorted = [i for i in order if keep[i]]
    kept_idx = kept_sorted[:3] + (kept_sorted[-3:] if len(kept_sorted) > 3 else [])
    rej_sorted = [i for i in order if not keep[i]]
    border_idx = rej_sorted[:6]
    worst_idx = rej_sorted[-6:] if len(rej_sorted) > 6 else []
    t = np.arange(dff.shape[1]) / fs / 60.0
    seg = slice(0, min(dff.shape[1], int(10 * 60 * fs)))    # first 10 min
    for col, (title, idx) in enumerate([("kept: 3 highest + 3 lowest score", kept_idx),
                                        ("rejected, closest to threshold", border_idx),
                                        ("rejected, lowest score", worst_idx)]):
        ax = fig.add_subplot(gs[1:, col])
        off = 0.0
        for i in idx:
            x = dff[i, seg]; sc = np.nanpercentile(x, 99.5) - np.nanpercentile(x, 0.5) + 1e-6
            ax.plot(t[seg], (x - np.nanpercentile(x, 0.5)) / sc + off, lw=0.5, c="tab:red" if keep[i] else "0.3")
            ax.text(t[seg][0], off + 0.62, f"roi {i}  score {score[i]:.2f}  ev {feats['n_sig_events'][i]}  "
                    f"asym {feats['asymmetry'][i]:.1f}  snr {feats['snr'][i]:.1f}  skew {feats['skew'][i]:.1f}",
                    fontsize=7, va="bottom", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
            off += 1.2
        ax.set_yticks([]); ax.set_xlabel("time (min, first 10 min)"); ax.set_title(title, fontsize=10)
    fig.suptitle(f"ROI curation: kept {int(keep.sum())} / {n}  (thresholds dashed; see roi_curation.csv)")
    fig.tight_layout()
    fig.savefig(out_dir / "roi_curation.png", dpi=110)
    plt.close(fig)


def curate(dff: np.ndarray, fs: float, a: dict, stat=None, F=None, Fneu=None,
           redcell=None, report_dir: Path | None = None, names=None,
           red_img: np.ndarray | None = None):
    """Full curation: features -> keep/score (+ report). Returns keep, score, feats."""
    p = {k: a.get(k, v) for k, v in DEFAULTS.items()}
    feats, events = roi_features(dff, fs, stat=stat, F=F, Fneu=Fneu, redcell=redcell,
                                 smooth_s=float(p["curate_smooth_s"]),
                                 z_on=float(p["curate_sig_z_on"]),
                                 z_peak=float(p["curate_sig_z_peak"]),
                                 min_dur_s=float(p["curate_sig_min_dur_s"]),
                                 neu_coeff=float(a.get("neuropil_coeff", 0.7)), red_img=red_img)
    keep, score, crit, p = curate_rois(feats, a)
    if report_dir is not None:
        write_curation_report(report_dir, feats, crit, score, keep, dff, fs, p,
                              events=events, names=names)
    return keep, score, feats
