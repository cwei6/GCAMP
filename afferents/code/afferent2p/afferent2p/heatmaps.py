"""Heatmap activity comparison (slides 21 / 27 of LabMeeting_261012):

  H1  CFU x time heatmap per experiment — the runs concatenated in clock order
      (RearChr: pre r1 | pre r2 ‖ post r1 | post r2; Treg: day 1 | day 2 | ...),
      rows = the experiment's CFU identities (aqua_analyze.match_cfus: the same
      unit across runs), each row's trace read from EVERY run's binned movie
      under the identity's footprint (mapped with the match transform and that
      run's AQuA2 crop offset), dF/F with F0 = 20th percentile per run, then
      min-max normalised 0-1 per CFU over all segments; rows grouped by k-means
      (k = 3) on the segment means and ordered by their post - pre change; the
      event raster and the running speed underneath. Without movies the rows
      are AQuA2's own CFU dF/F where the CFU was detected and blank elsewhere.
  H2  up / down / unchanged per CFU: 20-s bin means of the normalised trace,
      condition vs baseline, Mann-Whitney p < 0.05 + direction; pie with
      "Total = N CFUs" and cfu_updown_<exp>.csv; per-CFU event rate scatter.
  H3  the 25-um tile maps already made by afferent2p.spatial (<exp>_maps.png).
  H4  the peri-locomotion-onset heatmap is made by afferent2p.locotrials.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .locotrials import _style

BIN_FRAMES_20S = 10


# ----------------------------------------------------------------------------- traces

def _footprint_in_run(mask_home: np.ndarray, tr_home: dict, tr_run: dict) -> np.ndarray:
    """A footprint found in its home run, carried into another run's crop frame:
    home -> run 1 (the home run's match transform) -> run (inverse of that
    run's transform). Small angles, so the inverse is applied as -angle,
    -shift; the footprint is dilated by 1 px to tolerate the residual."""
    from .prealign import apply_transform
    from scipy.ndimage import binary_dilation
    m = mask_home.astype(np.float32)[None]
    if tr_home.get("applied"):
        m = apply_transform(m, tr_home["angle_deg"], tr_home["dy"], tr_home["dx"], cval=0.0)
    if tr_run.get("applied"):
        m = apply_transform(m, -tr_run["angle_deg"], -tr_run["dy"], -tr_run["dx"], cval=0.0)
    out = m[0] > 0.5
    return binary_dilation(out, iterations=1) if out.any() else out


def _movie_traces(movie, px: dict, chunk: int = 32) -> dict:
    """Mean over each footprint per movie frame, in one pass over the movie.
    px: {key: (ys, xs)} in raw-frame pixels."""
    T = movie.shape[0]
    out = {k: np.empty(T, np.float32) for k in px}
    if not px:
        return out
    ys_all = np.concatenate([v[0] for v in px.values()]); xs_all = np.concatenate([v[1] for v in px.values()])
    y0, y1, x0, x1 = ys_all.min(), ys_all.max() + 1, xs_all.min(), xs_all.max() + 1
    idx = {k: (v[0] - y0, v[1] - x0) for k, v in px.items()}
    for s in range(0, T, chunk):
        blk = np.asarray(movie[s:s + chunk, y0:y1, x0:x1], np.float32)
        for k, (yy, xx) in idx.items():
            out[k][s:s + blk.shape[0]] = blk[:, yy, xx].mean(axis=1)
    return out


def _dff(f: np.ndarray, pct: float = 20.0) -> np.ndarray:
    ok = np.isfinite(f)
    if ok.sum() < 5:
        return np.full_like(f, np.nan)
    f0 = float(np.percentile(f[ok], pct))
    if f0 <= 0:
        f0 = float(np.nanmean(f[ok])) or 1.0
    return (f - f0) / f0


def identity_traces(bundles: list[dict], cfu_list: list[dict], ids: list, transforms: list, n_id: int,
                    movies: bool = True, log=print) -> dict:
    """{traces: (n_id, T_total) 0-1 normalised, raw: dF/F, segments: [(start, end, run, condition)],
    detected: (n_id, n_runs) bool, source: 'movie' | 'cfu', onsets: (n_id, T_total)}"""
    from .report import _mean_image, _crop_offset
    from .aqua_io import open_movie, aqua_files
    from .spatial import read_meanimg
    n_runs = len(bundles)
    T_tot = sum(b["T"] for b in bundles)
    segs, start = [], 0
    for b in bundles:
        segs.append((start, start + b["T"], b["run"], b["condition"])); start += b["T"]
    raw = np.full((n_id, T_tot), np.nan)
    onsets = np.zeros((n_id, T_tot))
    detected = np.zeros((n_id, n_runs), bool)
    home = {}
    for i, idv in enumerate(ids):
        for k, cid in enumerate(idv):
            if cid >= 0:
                detected[cid, i] = True
                home.setdefault(int(cid), (i, k))
    # the per-run movie geometry
    use_movie = movies
    geo = []
    if use_movie:
        for b in bundles:
            tif = Path(b["tif"]) if b.get("tif") else None
            g = None
            if tif is not None and (tif.exists() and tif.stat().st_size > 0 or aqua_files(tif)["res"].exists()):
                try:
                    mv, src = open_movie(tif, prefer_registered=True, log=log)
                    if src == "registered":
                        g = {"movie": mv, "offset": (0, 0), "source": src}       # AQuA2's frame = the CFU maps' frame
                    else:
                        raw_mean = _mean_image(tif, mv)
                        exported = read_meanimg(tif)
                        off = _crop_offset(raw_mean, exported) if exported is not None else (0, 0)
                        g = {"movie": mv, "offset": off or (0, 0), "source": src}
                except Exception as e:
                    log(f"    heatmap: movie of {b['run']} not usable ({e})")
            geo.append(g)
        if not all(g is not None for g in geo):
            use_movie = False
            log("    heatmap: not every run has a readable binned movie — rows use AQuA2's CFU dF/F where detected")
    # onsets (and the no-movie fallback) per identity and run
    for cid in range(n_id):
        if cid not in home:
            continue
        for i, (b, (s0, s1, _, _)) in enumerate(zip(bundles, segs)):
            k = int(np.flatnonzero(ids[i] == cid)[0]) if detected[cid, i] else None
            if k is not None and b.get("cfu_onsets") is not None and b["cfu_onsets"].shape[0] > k:
                onsets[cid, s0:s1] = b["cfu_onsets"][k][:b["T"]]
            if not use_movie and k is not None and b.get("cfu_dff") is not None:
                raw[cid, s0:s1] = b["cfu_dff"][k][:b["T"]]
    # movie rows: one pass over each run's movie for every identity's footprint in that run
    if use_movie:
        for i, (b, (s0, s1, _, _)) in enumerate(zip(bundles, segs)):
            g = geo[i]
            mv = g["movie"]
            H, W = mv.shape[1:]
            oy, ox = g["offset"]
            px = {}                                           # identity -> (ys, xs) in the raw frame
            for cid, (hi, hk) in home.items():
                mask_home = cfu_list[hi]["maps"][hk] > 0.1
                m = _footprint_in_run(mask_home, transforms[hi], transforms[i])     # same treatment in every run (dilated 1 px)
                ys, xs = np.nonzero(m)
                ys = ys + oy; xs = xs + ox
                ok = (ys >= 0) & (ys < H) & (xs >= 0) & (xs < W)
                if ok.any():
                    px[cid] = (ys[ok], xs[ok])
            traces = _movie_traces(mv, px)
            for cid, tr in traces.items():
                tr = _dff(tr)
                rec = np.full(b["T"], np.nan)
                dropped = np.asarray(b.get("dropped", np.zeros(b["T"], bool)), bool)
                present = np.flatnonzero(~dropped)
                if present.size == tr.size:
                    rec[present] = tr
                else:
                    n = min(tr.size, b["T"]); rec[:n] = tr[:n]
                raw[cid, s0:s1] = rec
            log(f"    heatmap: {b['run']}: {len(traces)} footprint traces from the {'registered' if g.get('source') == 'registered' else 'binned'} movie")
    for g in geo:                                          # release the HDF5 handles of the registered movies
        if g is not None and hasattr(g.get("movie"), "close"):
            try:
                g["movie"].close()
            except Exception:
                pass
    # excluded frames blanked, 0-1 per CFU over all segments
    for b, (s0, s1, _, _) in zip(bundles, segs):
        ex = np.asarray(b["excluded"], bool)
        raw[:, s0:s1][:, ex[:b["T"]]] = np.nan
    lo = np.nanmin(raw, axis=1, keepdims=True); hi = np.nanmax(raw, axis=1, keepdims=True)
    with np.errstate(all="ignore"):
        norm = (raw - lo) / np.where(hi - lo > 1e-9, hi - lo, np.nan)
    return {"traces": norm, "raw": raw, "segments": segs, "detected": detected, "source": "movie" if use_movie else "cfu",
            "onsets": onsets, "home": home}


# ----------------------------------------------------------------------------- ordering and up/down

def kmeans_order(norm: np.ndarray, segs: list, conditions: list, k: int = 3, seed: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rows ordered by k-means (k groups) on the per-condition means, groups by
    their (last condition - baseline) change, rows within a group by the
    baseline mean. Returns (order, group per row, condition means (n, n_cond))."""
    n = norm.shape[0]
    cm = np.full((n, len(conditions)), np.nan)
    for j, c in enumerate(conditions):
        cols = np.concatenate([np.arange(s0, s1) for s0, s1, _, cc in segs if cc == c] or [np.zeros(0, int)])
        if cols.size:
            with np.errstate(all="ignore"):
                cm[:, j] = np.nanmean(norm[:, cols], axis=1)
    valid = np.isfinite(cm).all(axis=1)
    groups = np.full(n, -1, int)
    if valid.sum() >= k:
        X = cm[valid]
        try:
            from scipy.cluster.vq import kmeans2
            _, lab = kmeans2(X, k, seed=seed, minit="++")
        except Exception:
            lab = np.zeros(X.shape[0], int)
        groups[valid] = lab
    elif valid.any():
        groups[valid] = 0
    change = cm[:, -1] - cm[:, 0] if cm.shape[1] > 1 else cm[:, 0]
    g_change = {g: float(np.nanmean(change[groups == g])) for g in np.unique(groups) if g >= 0}
    g_rank = {g: r for r, g in enumerate(sorted(g_change, key=lambda g: -g_change[g]))}
    key = [((g_rank.get(groups[i], 99)), -np.nan_to_num(cm[i, 0], nan=-1)) for i in range(n)]
    order = np.array(sorted(range(n), key=lambda i: key[i]), int)
    return order, groups, cm


def updown(norm: np.ndarray, segs: list, conditions: list, bin_frames: int = BIN_FRAMES_20S, alpha: float = 0.05) -> list[dict]:
    """Per CFU and condition: 20-s bin means vs the baseline's, Mann-Whitney."""
    from scipy.stats import mannwhitneyu
    base = conditions[0]

    def bins(row, c):
        out = []
        for s0, s1, _, cc in segs:
            if cc != c:
                continue
            seg = row[s0:s1]
            for a in range(0, seg.size, bin_frames):
                v = seg[a:a + bin_frames]
                if np.isfinite(v).sum() >= max(2, bin_frames // 2):
                    out.append(float(np.nanmean(v)))
        return np.array(out)
    rows = []
    for i in range(norm.shape[0]):
        bb = bins(norm[i], base)
        for c in conditions[1:]:
            bc = bins(norm[i], c)
            p = np.nan; d = np.nan
            if bb.size >= 3 and bc.size >= 3:
                d = float(np.mean(bc) - np.mean(bb))
                try:
                    p = float(mannwhitneyu(bc, bb, alternative="two-sided").pvalue)
                except Exception:
                    p = np.nan
            cls = "unchanged"
            if np.isfinite(p) and p < alpha:
                cls = "up" if d > 0 else "down"
            if not (bb.size >= 3 and bc.size >= 3):
                cls = "not in both"
            rows.append({"cfu": i + 1, "condition": c, "baseline": base, "n_bins_baseline": int(bb.size), "n_bins_condition": int(bc.size),
                         "mean_baseline": float(np.mean(bb)) if bb.size else np.nan, "mean_condition": float(np.mean(bc)) if bc.size else np.nan,
                         "difference": d, "p_mannwhitney": p, "class": cls})
    return rows


# ----------------------------------------------------------------------------- figures

def heatmap_figure(exp: str, tr: dict, bundles: list[dict], conditions: list, order: np.ndarray, groups: np.ndarray, out_png: Path,
                   stim_after: str | None = None, title_note: str = ""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    norm = tr["traces"]; segs = tr["segments"]
    n, T = norm.shape
    spf = bundles[0]["spf"]
    fig, axes = plt.subplots(3, 1, figsize=(15, 4.5 + 0.045 * n + 2.6), sharex=True,
                             gridspec_kw={"height_ratios": [max(2.0, 0.045 * n), 1.0, 0.7]})
    ax = axes[0]
    M = norm[order]
    im = ax.imshow(np.where(np.isfinite(M), M, 0), aspect="auto", cmap="viridis", vmin=0, vmax=1, interpolation="nearest",
                   extent=[0, T * spf / 60.0, n, 0])
    nanmask = ~np.isfinite(M)
    if nanmask.any():
        ax.imshow(np.where(nanmask, 1.0, np.nan), aspect="auto", cmap="gray", vmin=0, vmax=1.4, alpha=0.9, interpolation="nearest",
                  extent=[0, T * spf / 60.0, n, 0])
    # group bands on the left
    g_sorted = groups[order]
    for g in np.unique(g_sorted):
        if g < 0:
            continue
        rows = np.flatnonzero(g_sorted == g)
        ax.plot([-0.0, -0.0], [rows.min(), rows.max() + 1], lw=0)
    # separators
    for j, (s0, s1, run, c) in enumerate(segs):
        if j > 0:
            thick = segs[j - 1][3] != c
            for a in axes:
                a.axvline(s0 * spf / 60.0, color="white" if a is axes[0] else "#444", lw=2.2 if thick else 0.9, ls="--")
        ax.text((s0 + s1) / 2 * spf / 60.0, -0.4, f"{c}: {run.split('_')[-3] if run.count('_') >= 3 else run}", ha="center", va="bottom", fontsize=7, color="#333")
    ax.set_ylabel(f"CFU identity ({n}), k-means groups\nsorted by change", fontsize=8)
    ax.set_ylim(n, 0)
    cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.025)
    cb.set_label("dF/F, 0-1 per CFU", fontsize=8); cb.ax.tick_params(labelsize=7)
    # raster
    ax = axes[1]
    on = tr["onsets"][order]
    rr, cc = np.nonzero(on > 0)
    ax.scatter(cc * spf / 60.0, rr, s=4, color="#111", marker="|")
    ax.set_ylim(n, -0.5); ax.set_ylabel("event onsets\n(same rows)", fontsize=8); _style(ax)
    # speed
    ax = axes[2]
    for s0, s1, run, c in segs:
        b = next(x for x in bundles if x["run"] == run)
        t = (np.arange(b["T"]) + s0 + 0.5) * spf / 60.0
        sp = np.asarray(b["speed"], float)
        if np.isfinite(sp).any():
            ax.plot(t, sp, color="#333", lw=0.5)
            ax.fill_between(t, 0, np.where(np.asarray(b["running"], bool), np.nanmax(sp) * 1.05, 0), color="#f2c94c", alpha=0.3, lw=0, step="mid")
    ax.set_ylabel("speed", fontsize=8); ax.set_xlabel("time (min, runs concatenated; white dashed = run boundary, thick = the stimulation / day change)", fontsize=8)
    _style(ax)
    note = ("rows from the movies (AQuA2-registered where available) under each identity's footprint (grey = no frame)" if tr["source"] == "movie"
            else "rows = AQuA2 CFU dF/F where the CFU was detected; grey = not detected in that run (no movie to read)")
    fig.suptitle(f"{exp}: CFU x time, 0-1 per CFU, {note}{title_note}", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)


def updown_figure(exp: str, rows: list[dict], conditions: list, rates: dict, out_png: Path):
    """Pie(s) of up / down / unchanged per condition vs baseline + per-CFU event-rate scatter."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    conds = conditions[1:]
    base = conditions[0]
    fig, axes = plt.subplots(1, 2 * max(len(conds), 1), figsize=(4.4 * 2 * max(len(conds), 1), 4.2), squeeze=False)
    cols = {"up": "#d7263d", "down": "#2a78d6", "unchanged": "#bdbdbd", "not in both": "#f0f0f0"}
    for j, c in enumerate(conds):
        sub = [r for r in rows if r["condition"] == c]
        n = len(sub)
        counts = {k: sum(1 for r in sub if r["class"] == k) for k in ("up", "down", "unchanged", "not in both")}
        ax = axes[0, 2 * j]
        vals = [counts[k] for k in ("up", "down", "unchanged", "not in both") if counts[k]]
        labs = [f"{k} {100 * counts[k] / n:.0f} %" for k in ("up", "down", "unchanged", "not in both") if counts[k]]
        if n:
            ax.pie(vals, labels=labs, colors=[cols[k] for k in ("up", "down", "unchanged", "not in both") if counts[k]],
                   startangle=90, counterclock=False, textprops={"fontsize": 8}, wedgeprops={"lw": 0.6, "edgecolor": "white"})
        ax.set_title(f"{c} vs {base}: Total = {n} CFUs\n(20-s bin means, Mann-Whitney p < 0.05)", fontsize=8)
        ax = axes[0, 2 * j + 1]
        rb = rates.get(base, {}); rc = rates.get(c, {})
        x = np.array([rb.get(r["cfu"], np.nan) for r in sub]); y = np.array([rc.get(r["cfu"], np.nan) for r in sub])
        colr = [cols[r["class"]] for r in sub]
        ax.scatter(x, y, s=16, c=colr, edgecolor="#555", lw=0.3)
        lim = max(np.nanmax(np.r_[x, y, 0.1]) * 1.15, 0.1) if np.isfinite(np.r_[x, y]).any() else 1.0
        ax.plot([0, lim], [0, lim], color="#999", lw=0.8, ls=":")
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel(f"events / min in {base}", fontsize=8); ax.set_ylabel(f"events / min in {c}", fontsize=8)
        ax.set_title("per-CFU event rate (colour = class)", fontsize=8); _style(ax)
    fig.suptitle(f"{exp}: CFUs up / down / unchanged", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)


# ----------------------------------------------------------------------------- entry point

def analyze_experiment(exp: str, bundles: list[dict], conditions: list, cfu_list: list[dict], ids: list, transforms: list,
                       n_id: int, out_dir: Path, movies: bool = True, log=print) -> dict:
    """H1 + H2 of one experiment. bundles and cfu_list in the same run order
    (clock order: the baseline runs first). Writes <exp>_heatmap.png,
    <exp>_updown.png, cfu_updown_<exp>.csv, cfu_traces_<exp>.npz."""
    import pandas as pd
    out_dir = Path(out_dir)
    if n_id == 0:
        return {"lines": [f"  heatmap ({exp}): no CFUs"], "n_id": 0}
    tr = identity_traces(bundles, cfu_list, ids, transforms, n_id, movies=movies, log=log)
    order, groups, cm = kmeans_order(tr["traces"], tr["segments"], conditions)
    rows = updown(tr["traces"], tr["segments"], conditions)
    # per-CFU event rates per condition (identities)
    rates = {}
    for c in conditions:
        minutes = sum(b["minutes_used"] for b in bundles if b["condition"] == c)
        cols = np.concatenate([np.arange(s0, s1) for s0, s1, _, cc in tr["segments"] if cc == c] or [np.zeros(0, int)])
        n_ev = tr["onsets"][:, cols].sum(axis=1) if cols.size else np.zeros(n_id)
        rates[c] = {i + 1: (n_ev[i] / minutes if minutes > 0 else np.nan) for i in range(n_id)}
    for r in rows:
        r["experiment"] = exp
        r["events_per_min_baseline"] = rates[conditions[0]].get(r["cfu"])
        r["events_per_min_condition"] = rates[r["condition"]].get(r["cfu"])
        r["group"] = int(groups[r["cfu"] - 1])
        r["detected_in_runs"] = int(tr["detected"][r["cfu"] - 1].sum())
    pd.DataFrame(rows).to_csv(out_dir / f"cfu_updown_{exp}.csv", index=False)
    np.savez_compressed(out_dir / f"cfu_traces_{exp}.npz", traces=tr["traces"].astype(np.float32), raw=tr["raw"].astype(np.float32),
                        onsets=tr["onsets"].astype(np.uint8), order=order, groups=groups,
                        segments=np.array([(s0, s1) for s0, s1, _, _ in tr["segments"]]),
                        runs=np.array([r for _, _, r, _ in tr["segments"]]), conditions_of_segments=np.array([c for _, _, _, c in tr["segments"]]))
    try:
        heatmap_figure(exp, tr, bundles, conditions, order, groups, out_dir / f"{exp}_heatmap.png")
    except Exception as e:
        log(f"    heatmap figure failed ({e})")
    try:
        updown_figure(exp, rows, conditions, rates, out_dir / f"{exp}_updown.png")
    except Exception as e:
        log(f"    up/down figure failed ({e})")
    lines = [f"  heatmap ({exp}): {n_id} CFU identities, rows from {'the movies (AQuA2-registered where readable)' if tr['source'] == 'movie' else 'AQuA2 CFU dF/F (no movies)'}; "
             f"k-means groups of {len(set(groups[groups >= 0]))}"]
    for c in conditions[1:]:
        sub = [r for r in rows if r["condition"] == c]
        n = len(sub)
        if n:
            up = sum(r["class"] == "up" for r in sub); dn = sum(r["class"] == "down" for r in sub); nb = sum(r["class"] == "not in both" for r in sub)
            lines.append(f"    {c} vs {conditions[0]}: {up} up ({100 * up / n:.0f} %), {dn} down ({100 * dn / n:.0f} %), "
                         f"{n - up - dn - nb} unchanged of {n} CFUs" + (f" ({nb} without data in both)" if nb else "")
                         + " — 20-s bin means of the 0-1 trace, Mann-Whitney p < 0.05 (within one animal: descriptive)")
    summary = {"experiment": exp, "n_id": n_id, "source": tr["source"]}
    for c in conditions[1:]:
        sub = [r for r in rows if r["condition"] == c]
        summary[f"n_up[{c}]"] = sum(r["class"] == "up" for r in sub); summary[f"n_down[{c}]"] = sum(r["class"] == "down" for r in sub)
        summary[f"n_unchanged[{c}]"] = sum(r["class"] == "unchanged" for r in sub)
    return {"lines": lines, "n_id": n_id, "source": tr["source"], "summary": summary}
