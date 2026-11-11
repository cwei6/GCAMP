"""Pupil (arousal) versus afferent activity, per run and pre/post.

Inputs per run (all already produced by the pipeline):
  <run>_..._behavior.csv   pupil_diam_px_mean, pupil_valid_fraction, velocity_mean,
                           still_fraction per RECORDED frame (2-s bins of the eye camera
                           and wheel, same bins as the movie)
  AQuA2 event table        events kept by aqua_analyze (t0 in recorded frames)
  <run>_..._res_cfu.mat    per-CFU dF/F curves and event occurrences (movie frames)
  <run>_..._frames.csv     excluded frames (artefacts + decided state stretches)

Per run it builds three traces on the recorded-frame grid — pupil diameter,
population activity (mean dF/F over the kept CFUs) and event-onset rate
(kept events, smoothed) — masks excluded and invalid frames, and reports:
  * pupil size (median, 10–90 % range, % valid) and its coupling to running
  * pupil ↔ activity: Spearman r at lags −lag_max..+lag_max (positive lag =
    activity follows the pupil), the best lag, and a circular-shift null for
    its p value; the same at lag 0 for the event-onset rate; partial r with
    running speed regressed out; r on still frames only
  * event rate inside low / mid / high pupil terciles (tercile edges pooled
    over the whole experiment, so pre and post are compared at the SAME pupil
    size — an activity change that survives at matched pupil is not just an
    arousal change)
  * per CFU: Spearman r with the pupil and whether it beats its own
    circular-shift null (coupled up / down / not)
Across conditions it summarises pupil, coupling and matched-pupil rates, and
draws <exp>_pupil.png. Statistics are descriptive within a FOV; the unit for a
group test is one value per animal.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np

PUPIL_COL = "pupil_diam_px_mean"
VALID_COL = "pupil_valid_fraction"
SPEED_COL = "velocity_mean"
STILL_COL = "still_fraction"


# ----------------------------------------------------------------------------- helpers

def _rank(x: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    return rankdata(x)


def spearman(x: np.ndarray, y: np.ndarray, min_n: int = 20) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < min_n:
        return float("nan")
    xr, yr = _rank(x[m]), _rank(y[m])
    if xr.std() == 0 or yr.std() == 0:
        return float("nan")
    return float(np.corrcoef(xr, yr)[0, 1])


def lagged_corr(pupil: np.ndarray, act: np.ndarray, lags: np.ndarray, min_n: int = 20) -> np.ndarray:
    """r(lag): act[t] vs pupil[t - lag]; positive lag = activity follows the pupil."""
    out = np.full(lags.size, np.nan)
    T = pupil.size
    for i, lag in enumerate(lags):
        lag = int(lag)
        if lag >= 0:
            p, a = pupil[:T - lag], act[lag:]
        else:
            p, a = pupil[-lag:], act[:T + lag]
        out[i] = spearman(p, a, min_n=min_n)
    return out


def shift_null(pupil: np.ndarray, act: np.ndarray, lags: np.ndarray, n: int = 200, min_shift: int = 15,
               rng=None) -> np.ndarray:
    """max |r| over lags for n circular shifts of the pupil trace (>= min_shift frames)."""
    rng = rng or np.random.default_rng(0)
    T = pupil.size
    if T <= 2 * min_shift + 2:
        return np.array([])
    out = np.empty(n)
    for k in range(n):
        s = int(rng.integers(min_shift, T - min_shift))
        r = lagged_corr(np.roll(pupil, s), act, lags)
        out[k] = np.nanmax(np.abs(r)) if np.isfinite(r).any() else np.nan
    return out


def partial_spearman(x, y, z, min_n: int = 20) -> float:
    """Spearman r between x and y after regressing a linear z out of both."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if m.sum() < min_n or np.nanstd(z[m]) == 0:
        return float("nan")
    xr, yr, zr = _rank(x[m]), _rank(y[m]), _rank(z[m])
    Z = np.column_stack([zr, np.ones(zr.size)])
    rx = xr - Z @ np.linalg.lstsq(Z, xr, rcond=None)[0]
    ry = yr - Z @ np.linalg.lstsq(Z, yr, rcond=None)[0]
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def smooth_box(x: np.ndarray, w: int) -> np.ndarray:
    w = max(1, int(w) | 1)
    k = np.ones(w) / w
    x = np.asarray(x, float)
    m = np.isfinite(x)
    num = np.convolve(np.where(m, x, 0.0), k, mode="same")
    den = np.convolve(m.astype(float), k, mode="same")
    out = num / np.maximum(den, 1e-9)
    out[den < 0.5] = np.nan
    return out


# ----------------------------------------------------------------------------- per run

def run_series(ev, info: dict, beh, frames: dict | None, cfu: dict | None, spf: float,
               smooth_s: float = 10.0, min_valid_fraction: float = 0.5) -> dict | None:
    """Traces on the recorded-frame grid + masks. None without pupil data."""
    if beh is None or PUPIL_COL not in beh.columns:
        return None
    T = int(info.get("frames") or len(beh))
    n = min(T, len(beh))
    pupil = np.full(T, np.nan)
    pupil[:n] = beh[PUPIL_COL].to_numpy(float)[:n]
    if VALID_COL in beh.columns:
        vf = np.zeros(T); vf[:n] = beh[VALID_COL].to_numpy(float)[:n]
        pupil[vf < min_valid_fraction] = np.nan
    speed = np.full(T, np.nan)
    if SPEED_COL in beh.columns:
        speed[:n] = np.abs(beh[SPEED_COL].to_numpy(float)[:n])
    still = np.full(T, np.nan)
    if STILL_COL in beh.columns:
        still[:n] = beh[STILL_COL].to_numpy(float)[:n]
    excluded = np.zeros(T, bool)
    if frames is not None:
        ex = frames["excluded"]
        excluded[:min(T, ex.size)] = ex[:min(T, ex.size)]
    # population dF/F from the CFU curves (movie frames -> recorded frames)
    pop = np.full(T, np.nan)
    cfu_dff_rec = None
    if cfu is not None and cfu.get("dff") is not None and np.asarray(cfu["dff"]).size:
        dff = np.asarray(cfu["dff"], float)                 # nCFU x T_movie
        Tm = dff.shape[1]
        rec_of_movie = None
        if frames is not None and frames["dropped"].any():
            present = np.flatnonzero(~frames["dropped"])
            if present.size == Tm:
                rec_of_movie = frames["frame"][present] - 1
        if rec_of_movie is None:
            rec_of_movie = np.arange(min(Tm, T))
        cfu_dff_rec = np.full((dff.shape[0], T), np.nan)
        cfu_dff_rec[:, rec_of_movie] = dff[:, :rec_of_movie.size]
        pop[:] = np.nanmean(cfu_dff_rec, axis=0) if dff.shape[0] else np.nan
    # event-onset rate (kept events), events per minute, smoothed
    onsets = np.zeros(T)
    if "t0_frame" in ev.columns and len(ev):
        t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, T - 1)
        np.add.at(onsets, t0, 1.0)
    rate = smooth_box(onsets / (spf / 60.0), int(round(smooth_s / spf)))
    good = ~excluded
    for arr in (pupil, speed, still, pop, rate):
        arr[~good] = np.nan
    if cfu_dff_rec is not None:
        cfu_dff_rec[:, ~good] = np.nan
    return {"T": T, "pupil": pupil, "speed": speed, "still": still, "excluded": excluded,
            "pop_dff": pop, "rate": rate, "onsets": onsets, "cfu_dff": cfu_dff_rec, "spf": spf}


def run_metrics(series: dict, lag_max_s: float = 20.0, n_shuffle: int = 200, still_thresh: float = 0.5,
                rng=None) -> dict:
    rng = rng or np.random.default_rng(0)
    spf = series["spf"]
    pupil, pop, rate, speed, still = (series[k] for k in ("pupil", "pop_dff", "rate", "speed", "still"))
    m = {}
    v = pupil[np.isfinite(pupil)]
    m["pupil_frames_valid"] = int(v.size)
    m["pupil_median_px"] = float(np.median(v)) if v.size else np.nan
    m["pupil_p10_px"] = float(np.percentile(v, 10)) if v.size else np.nan
    m["pupil_p90_px"] = float(np.percentile(v, 90)) if v.size else np.nan
    m["pupil_range_pct"] = (100 * (m["pupil_p90_px"] - m["pupil_p10_px"]) / m["pupil_median_px"]
                            if v.size and m["pupil_median_px"] else np.nan)
    m["pupil_speed_r"] = spearman(pupil, speed)
    lag_max = max(1, int(round(lag_max_s / spf)))
    lags = np.arange(-lag_max, lag_max + 1)
    m["lags_s"] = lags * spf
    if np.isfinite(pop).sum() >= 20:
        r = lagged_corr(pupil, pop, lags)
        m["lag_r_dff"] = r
        if np.isfinite(r).any():
            i = int(np.nanargmax(np.abs(r)))
            m["pupil_dff_r_best"] = float(r[i]); m["pupil_dff_lag_s"] = float(lags[i] * spf)
            m["pupil_dff_r0"] = float(r[lag_max])
            null = shift_null(pupil, pop, lags, n=n_shuffle, rng=rng)
            m["pupil_dff_p"] = float((np.sum(null >= abs(r[i])) + 1) / (null.size + 1)) if null.size else np.nan
        m["pupil_dff_partial_r"] = partial_spearman(pupil, pop, speed)
        st = still >= still_thresh
        m["still_frames"] = int(np.nansum(st))
        m["pupil_dff_r_still"] = spearman(np.where(st, pupil, np.nan), pop, min_n=40)
    if np.isfinite(rate).sum() >= 20:
        r = lagged_corr(pupil, rate, lags)
        m["lag_r_rate"] = r
        m["pupil_rate_r0"] = float(r[lag_max])
        if np.isfinite(r).any():
            i = int(np.nanargmax(np.abs(r)))
            m["pupil_rate_r_best"] = float(r[i]); m["pupil_rate_lag_s"] = float(lags[i] * spf)
    # per CFU
    cd = series.get("cfu_dff")
    if cd is not None and cd.shape[0]:
        rs, coupled = [], []
        T = pupil.size
        for j in range(cd.shape[0]):
            rj = spearman(pupil, cd[j])
            rs.append(rj)
            if not np.isfinite(rj):
                coupled.append("n/a"); continue          # too few valid frames or a flat curve
            null = []
            for k in range(100):
                s = int(rng.integers(15, max(16, T - 15)))
                null.append(spearman(np.roll(pupil, s), cd[j]))
            null = np.asarray(null); null = null[np.isfinite(null)]
            if null.size < 50:
                coupled.append("n/a"); continue
            lo, hi = np.percentile(null, [2.5, 97.5])
            coupled.append("up" if rj > hi else "down" if rj < lo else "no")
        m["cfu_r"] = np.asarray(rs); m["cfu_coupled"] = coupled
        m["n_cfu"] = int(cd.shape[0])
        m["n_cfu_pupil_up"] = int(sum(c == "up" for c in coupled))
        m["n_cfu_pupil_down"] = int(sum(c == "down" for c in coupled))
        m["cfu_pupil_r_median"] = float(np.nanmedian(rs)) if np.isfinite(rs).any() else np.nan
    return m


def tercile_rates(series: dict, edges: np.ndarray) -> dict:
    """events/min of the run inside each pupil bin defined by edges (pooled)."""
    pupil, onsets, spf = series["pupil"], series["onsets"], series["spf"]
    out = {}
    names = ["low", "mid", "high"]
    for k in range(3):
        lo = -np.inf if k == 0 else edges[k - 1]
        hi = np.inf if k == 2 else edges[k]
        sel = np.isfinite(pupil) & (pupil >= lo) & (pupil < hi)
        minutes = sel.sum() * spf / 60.0
        n_ev = float(onsets[sel].sum())
        out[f"rate_{names[k]}_pupil"] = n_ev / minutes if minutes > 0.25 else np.nan
        out[f"minutes_{names[k]}_pupil"] = minutes
        out[f"n_events_{names[k]}_pupil"] = n_ev
    return out


# ----------------------------------------------------------------------------- per experiment

def experiment_pupil(exp_name: str, runs_info: list[dict], conditions: list[str], out_dir: Path,
                     lag_max_s: float = 20.0, log=print) -> tuple[list[dict], list[str]]:
    """runs_info: dicts with keys run, condition, series, metrics (filled here).
    Returns (rows for pupil_<exp>.csv, summary lines)."""
    import pandas as pd
    rows, lines = [], []
    have = [r for r in runs_info if r.get("series") is not None]
    if not have:
        return rows, ["  pupil: no behaviour csv with pupil columns (DLC eye file missing at unpack time)"]
    pooled = np.concatenate([r["series"]["pupil"][np.isfinite(r["series"]["pupil"])] for r in have])
    if pooled.size < 30:
        return rows, ["  pupil: too few valid pupil frames"]
    edges = np.percentile(pooled, [100 / 3, 200 / 3])
    for r in have:
        m = run_metrics(r["series"], lag_max_s=lag_max_s)
        m.update(tercile_rates(r["series"], edges))
        r["metrics"] = m
        row = {"experiment": exp_name, "run": r["run"], "condition": r["condition"]}
        for k, v in m.items():
            if isinstance(v, (int, float, np.floating, np.integer)):
                row[k] = round(float(v), 4) if isinstance(v, (float, np.floating)) else int(v)
        row["pupil_tercile_edges_px"] = f"{edges[0]:.1f}/{edges[1]:.1f}"
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / f"pupil_{exp_name}.csv", index=False)
    # per-CFU rows
    cfu_rows = []
    for r in have:
        m = r["metrics"]
        if "cfu_r" in m:
            for j, (rj, c) in enumerate(zip(m["cfu_r"], m["cfu_coupled"])):
                cfu_rows.append({"experiment": exp_name, "run": r["run"], "condition": r["condition"],
                                 "cfu": j + 1, "pupil_r": round(float(rj), 4) if np.isfinite(rj) else np.nan,
                                 "coupled": c})
    if cfu_rows:
        pd.DataFrame(cfu_rows).to_csv(out_dir / f"pupil_cfu_{exp_name}.csv", index=False)
    # summary lines
    lines.append(f"  pupil (terciles pooled over the FOV: < {edges[0]:.1f} / {edges[0]:.1f}-{edges[1]:.1f} / > {edges[1]:.1f} px):")
    for c in conditions:
        sub = df[df["condition"] == c]
        if sub.empty:
            continue
        def med(col):
            return float(np.nanmedian(sub[col])) if col in sub.columns and np.isfinite(sub[col]).any() else np.nan
        pooled_rates = []
        for b in ("low", "mid", "high"):
            n_ev = float(np.nansum(sub[f"n_events_{b}_pupil"])); mins = float(np.nansum(sub[f"minutes_{b}_pupil"]))
            pooled_rates.append(n_ev / mins if mins > 0.25 else np.nan)
        txt = (f"    {c:>10}: pupil median {_f(med('pupil_median_px'), '{:.1f}')} px "
               f"(10-90 % range {_f(med('pupil_range_pct'), '{:.0f}')} %), "
               f"pupil-running r {_f(med('pupil_speed_r'))}; ")
        if "pupil_dff_r_best" in sub.columns:
            txt += (f"pupil-activity r {_f(med('pupil_dff_r_best'), '{:+.2f}')} at lag "
                    f"{_f(med('pupil_dff_lag_s'), '{:+.0f}')} s "
                    f"(partial r vs running {_f(med('pupil_dff_partial_r'), '{:+.2f}')}, "
                    f"still-only {_f(med('pupil_dff_r_still'), '{:+.2f}')}); ")
        if "n_cfu" in sub.columns:
            n_up = int(np.nansum(sub["n_cfu_pupil_up"])); n_dn = int(np.nansum(sub["n_cfu_pupil_down"])); n_all = int(np.nansum(sub["n_cfu"]))
            txt += f"CFUs coupled to the pupil: {n_up} up / {n_dn} down of {n_all}; "
        txt += "events/min at low / mid / high pupil: " + " / ".join(f"{x:.3f}" if np.isfinite(x) else "-" for x in pooled_rates)
        lines.append(txt)
    if len(conditions) > 1:
        base = conditions[0]
        b = df[df["condition"] == base]
        for c in conditions[1:]:
            o = df[df["condition"] == c]
            if b.empty or o.empty:
                continue
            parts = []
            pm_b, pm_o = float(np.nanmedian(b["pupil_median_px"])), float(np.nanmedian(o["pupil_median_px"]))
            if np.isfinite(pm_b) and np.isfinite(pm_o) and pm_b:
                parts.append(f"pupil {pm_b:.1f} -> {pm_o:.1f} px ({100 * (pm_o - pm_b) / pm_b:+.0f} %)")
            if "pupil_dff_r_best" in df.columns and np.isfinite(b["pupil_dff_r_best"]).any() \
                    and np.isfinite(o["pupil_dff_r_best"]).any():
                parts.append(f"pupil-activity r {float(np.nanmedian(b['pupil_dff_r_best'])):+.2f} -> "
                             f"{float(np.nanmedian(o['pupil_dff_r_best'])):+.2f}")
            matched = []
            for bn in ("low", "mid", "high"):
                nb, mb = float(np.nansum(b[f"n_events_{bn}_pupil"])), float(np.nansum(b[f"minutes_{bn}_pupil"]))
                no, mo = float(np.nansum(o[f"n_events_{bn}_pupil"])), float(np.nansum(o[f"minutes_{bn}_pupil"]))
                if mb > 0.5 and mo > 0.5:
                    rb, ro = nb / mb, no / mo
                    matched.append(f"{bn} x{ro / rb:.2f}" if rb > 0 else f"{bn} {ro:.3f}/min vs 0")
            if matched:
                parts.append("rate change at matched pupil: " + ", ".join(matched))
            lines.append(f"    {c} vs {base}: " + "; ".join(parts))
        lines.append("    (a rate change that holds within every pupil tercile is not explained by arousal; "
                     "one that vanishes there is)")
    try:
        pupil_figure(exp_name, have, conditions, edges, out_dir / f"{exp_name}_pupil.png")
    except Exception as e:
        log(f"  pupil figure failed: {e}")
    return rows, lines


def _f(x, fmt: str = "{:.2f}") -> str:
    """Format a number; 'n/a' when it is missing (nan/None)."""
    try:
        return fmt.format(x) if x is not None and np.isfinite(x) else "n/a"
    except (TypeError, ValueError):
        return "n/a"


# ----------------------------------------------------------------------------- figure

COND_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def pupil_figure(exp_name: str, runs: list[dict], conditions: list[str], edges: np.ndarray, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    n = len(runs)
    colour = {c: COND_COLOURS[i % len(COND_COLOURS)] for i, c in enumerate(conditions)}
    row_h = 2.6
    fig_h = row_h * n + 5.4
    fig = plt.figure(figsize=(14, fig_h))
    gs = GridSpec(2 * n + 2, 3, figure=fig, height_ratios=[0.55, 0.55] * n + [1.7, 1.7], hspace=0.5, wspace=0.3,
                  left=0.06, right=0.98, top=1 - 0.75 / fig_h, bottom=0.4 / fig_h)
    # per run: pupil trace (top) and activity trace (bottom), same time axis; lag curve on the right
    for i, r in enumerate(runs):
        s = r["series"]; t = np.arange(s["T"]) * s["spf"] / 60.0
        c = colour.get(r["condition"], "#333333")
        m = r.get("metrics", {})
        ax = fig.add_subplot(gs[2 * i, :2])
        ax.plot(t, s["pupil"], color=c, lw=1.0)
        ax.set_ylabel("pupil px", fontsize=8)
        ax.tick_params(labelsize=7, labelbottom=False)
        for a, b in _runs_of(s["excluded"]):
            ax.axvspan(t[a], t[b] + s["spf"] / 60, color="#888888", alpha=0.25, lw=0)
        ax.set_title(f"{r['run']}  [{r['condition']}]   pupil-activity r {_f(m.get('pupil_dff_r_best'), '{:+.2f}')} "
                     f"at {_f(m.get('pupil_dff_lag_s'), '{:+.0f}')} s (shuffle p {_f(m.get('pupil_dff_p'))}), "
                     f"partial r vs running {_f(m.get('pupil_dff_partial_r'), '{:+.2f}')}; grey = excluded frames",
                     fontsize=8, loc="left")
        ax2 = fig.add_subplot(gs[2 * i + 1, :2], sharex=ax)
        if np.isfinite(s["pop_dff"]).any():
            ax2.plot(t, s["pop_dff"], color="#2a2a2a", lw=0.8)
            ax2.set_ylabel("mean CFU\ndF/F", fontsize=8)
        else:
            ax2.plot(t, s["rate"], color="#2a2a2a", lw=0.8)
            ax2.set_ylabel("events/min", fontsize=8)
        for a, b in _runs_of(s["excluded"]):
            ax2.axvspan(t[a], t[b] + s["spf"] / 60, color="#888888", alpha=0.25, lw=0)
        ax2.tick_params(labelsize=7)
        if i == n - 1:
            ax2.set_xlabel("time (min)", fontsize=8)
        axl = fig.add_subplot(gs[2 * i:2 * i + 2, 2])
        if "lag_r_dff" in m:
            axl.plot(m["lags_s"], m["lag_r_dff"], color=c, lw=1.2)
            axl.axvline(0, color="#999999", lw=0.6); axl.axhline(0, color="#999999", lw=0.6)
            axl.set_ylim(-1, 1)
        axl.set_title("r(pupil, activity) by lag; + = activity follows pupil", fontsize=7)
        axl.tick_params(labelsize=7)
        if i == n - 1:
            axl.set_xlabel("lag (s)", fontsize=8)
    n = 2 * n                                        # grid row where the summary panels start
    # rates by pupil tercile per condition
    ax = fig.add_subplot(gs[n, 0:2])
    width = 0.8 / max(1, len(conditions))
    for ci, cnd in enumerate(conditions):
        sub = [r for r in runs if r["condition"] == cnd and r.get("metrics")]
        for k, bn in enumerate(("low", "mid", "high")):
            n_ev = sum(r["metrics"].get(f"n_events_{bn}_pupil", 0) for r in sub)
            mins = sum(r["metrics"].get(f"minutes_{bn}_pupil", 0) for r in sub)
            pooled = n_ev / mins if mins > 0.25 else np.nan
            x = k + (ci - (len(conditions) - 1) / 2) * width
            ax.bar(x, pooled, width=width * 0.92, color=colour[cnd], label=cnd if k == 0 else None)
            per_run = [r["metrics"].get(f"rate_{bn}_pupil", np.nan) for r in sub]
            ax.plot([x] * len(per_run), per_run, "o", color="#222222", ms=3, alpha=0.7)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([f"low\n< {edges[0]:.0f} px", f"mid", f"high\n> {edges[1]:.0f} px"], fontsize=8)
    ax.set_ylabel("events / min", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title("event rate at matched pupil size (bars = pooled, dots = runs)", fontsize=8, loc="left")
    ax.legend(fontsize=7, frameon=False)
    # per-CFU pupil correlation per condition
    ax = fig.add_subplot(gs[n, 2])
    for ci, cnd in enumerate(conditions):
        vals = np.concatenate([np.asarray(r["metrics"]["cfu_r"], float) for r in runs
                               if r["condition"] == cnd and r.get("metrics") and "cfu_r" in r["metrics"]] or [np.array([])])
        vals = vals[np.isfinite(vals)]
        if vals.size:
            jitter = (np.random.default_rng(ci).random(vals.size) - 0.5) * 0.3
            ax.plot(ci + jitter, vals, "o", color=colour[cnd], ms=3, alpha=0.6)
            ax.plot([ci - 0.25, ci + 0.25], [np.median(vals)] * 2, color="#222222", lw=1.5)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.set_xticks(range(len(conditions))); ax.set_xticklabels(conditions, fontsize=8)
    ax.set_ylabel("per-CFU r (dF/F vs pupil)", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title("CFU coupling to the pupil", fontsize=8, loc="left")
    # pupil vs running scatter per condition
    ax = fig.add_subplot(gs[n + 1, 0])
    for cnd in conditions:
        for r in runs:
            if r["condition"] != cnd:
                continue
            s = r["series"]
            ax.plot(s["speed"], s["pupil"], ".", color=colour[cnd], ms=2.5, alpha=0.35, label=cnd if r is [rr for rr in runs if rr["condition"] == cnd][0] else None)
    ax.set_xlabel("running speed", fontsize=8); ax.set_ylabel("pupil px", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title("pupil vs running (all frames)", fontsize=8, loc="left"); ax.legend(fontsize=7, frameon=False)
    # pupil distribution per condition
    ax = fig.add_subplot(gs[n + 1, 1])
    for ci, cnd in enumerate(conditions):
        vals = np.concatenate([r["series"]["pupil"][np.isfinite(r["series"]["pupil"])] for r in runs if r["condition"] == cnd] or [np.array([])])
        if vals.size:
            ax.hist(vals, bins=30, histtype="step", color=colour[cnd], lw=1.2, label=cnd, density=True)
    for e in edges:
        ax.axvline(e, color="#999999", lw=0.6, ls="--")
    ax.set_xlabel("pupil px", fontsize=8); ax.set_ylabel("density", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title("pupil size per condition (dashed = tercile edges)", fontsize=8, loc="left"); ax.legend(fontsize=7, frameon=False)
    # activity vs pupil scatter (binned) per condition
    ax = fig.add_subplot(gs[n + 1, 2])
    use_dff = any(np.isfinite(r["series"]["pop_dff"]).any() for r in runs)
    for cnd in conditions:
        xs, ys = [], []
        for r in runs:
            if r["condition"] != cnd:
                continue
            s = r["series"]
            y = s["pop_dff"] if np.isfinite(s["pop_dff"]).any() else s["rate"]
            xs.append(s["pupil"]); ys.append(y)
        if not xs:
            continue
        x = np.concatenate(xs); y = np.concatenate(ys)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 10:
            continue
        qs = np.percentile(x[m], np.linspace(0, 100, 9))
        cx, cy, ce = [], [], []
        for a, b in zip(qs[:-1], qs[1:]):
            sel = m & (x >= a) & (x <= b)
            if sel.sum() >= 3:
                cx.append(x[sel].mean()); cy.append(np.nanmean(y[sel])); ce.append(np.nanstd(y[sel]) / np.sqrt(sel.sum()))
        ax.errorbar(cx, cy, yerr=ce, color=colour[cnd], marker="o", ms=3, lw=1.2, capsize=2, label=cnd)
    ax.set_xlabel("pupil px (octiles)", fontsize=8)
    ax.set_ylabel("mean CFU dF/F" if use_dff else "events / min", fontsize=8)
    ax.tick_params(labelsize=7); ax.set_title("activity vs pupil size", fontsize=8, loc="left"); ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{exp_name} — pupil and afferent activity", fontsize=10)
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def _runs_of(flag: np.ndarray):
    idx = np.flatnonzero(flag)
    if idx.size == 0:
        return []
    runs = []
    s = p = int(idx[0])
    for i in idx[1:]:
        i = int(i)
        if i == p + 1:
            p = i
            continue
        runs.append((s, p)); s = p = i
    runs.append((s, p))
    return runs
