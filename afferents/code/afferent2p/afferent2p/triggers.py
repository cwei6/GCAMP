"""triggers — the three triggers of afferent activity, compared pre vs post
(optostim) or across days (Treg) on the same scale-aware statistics as the
rate comparison of the report:

  running onsets   the locomotion trials of locotrials.py (Yu's windows: 10 s
                   still before the onset, 30 s response window);
  pupil dilations  the same windows around dilation onsets of the pupil
                   (DLC diameter per 2-s frame, scaled 0-1 between the run's
                   1st and 99th percentile): a stable baseline window (range
                   <= 0.15, no running) followed by a rise >= 0.2 within 10 s;
                   a dilation with no running in the whole response window is
                   a STILL dilation — arousal without movement, the set the
                   comparison uses, so that a "pupil response" is not a running
                   response in disguise (dilations during running are kept
                   as 'moving' for the figure); the pupil also remains a
                   control elsewhere (matched-pupil terciles, GLM predictor);
  resting activity events in still frames >= 20 s away from any running
                   (rest_only of locotrials.py) — the whole run when the
                   animal never runs.

Per trigger and condition: the evoked event rate (kept events starting in
the response windows per minute of response window; rest: events per minute
of rest) with the rate difference, the exact Poisson interval of the ratio,
the pseudo-count log ratio and the modulation index (afferent2p.stats, as the
main comparison), the responsive-trial fraction (z rule) with a risk ratio
and Fisher's p, and the population response amplitude (peak - baseline
dF/F) with a bootstrap interval and the Mann-Whitney p. A trigger is
TESTABLE for a comparison when both conditions have >= MIN_TRIALS trials
(rest: >= 1 min each); otherwise the row says why (the animal never ran, no
still dilations, ...) and the report reads the other triggers — a session in
which the animal does not run is compared on its pupil dilations and its
resting activity instead of on nothing.

Outputs (per experiment, in aqua_analysis\\): pupil_trials_<exp>.csv,
<exp>_triggers.png; for the mouse: trigger_comparison.csv (one row per
experiment x condition x trigger).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .locotrials import PRE_FRAMES, RUN_FRAMES, POST_FRAMES, PALETTE, trial_responses, _style
from .stats import change_measures, fraction_change, amplitude_change

PUPIL_RISE = 0.20          # rise (fraction of the run's 1-99 % pupil range) within RISE_FRAMES that counts as a dilation onset
PUPIL_STABLE = 0.15        # the baseline window may vary at most this much
RISE_FRAMES = 5            # frames (10 s) in which the rise must happen
MIN_TRIALS = 3             # trials per condition for a comparison to be testable
MIN_REST_MIN = 1.0         # minutes of rest per condition
TRIGGERS = ("running", "pupil", "rest")
TRIGGER_LABEL = {"running": "running onsets", "pupil": "pupil dilations (still)", "rest": "resting activity"}


# ----------------------------------------------------------------------------- pupil-dilation trials

def pupil01(b: dict) -> np.ndarray | None:
    """The run's pupil diameter scaled 0-1 between its 1st and 99th percentile (NaN kept)."""
    p = np.asarray(b.get("pupil", []), float)
    if p.size == 0 or np.isfinite(p).sum() < 20:
        return None
    lo, hi = np.nanpercentile(p, 1), np.nanpercentile(p, 99)
    if not (hi > lo):
        return None
    return np.clip((p - lo) / (hi - lo), 0, 1)


def pupil_trials(b: dict, pre: int = PRE_FRAMES, run: int = RUN_FRAMES, post: int = POST_FRAMES,
                 rise: float = PUPIL_RISE, stable: float = PUPIL_STABLE, rise_frames: int = RISE_FRAMES) -> list[dict]:
    """Dilation-onset trials: onset t = first frame of a rise >= `rise` (of the
    run's range) within rise_frames after a stable (range <= `stable`), still,
    fully valid baseline window of `pre` frames; one onset per response window
    (the next can start `run` frames later); trials touching an excluded frame
    are dropped (counted in b['pupil_trials_dropped']). kind = 'still' (no
    running from the baseline to the end of the response window) or 'moving'."""
    p = pupil01(b)
    T = b["T"]
    b["pupil_trials_dropped"] = 0
    if p is None:
        return []
    running = np.asarray(b.get("running_merged", b["running"]), bool)
    excluded = np.asarray(b["excluded"], bool)
    out, dropped, last = [], 0, -10 ** 9
    for t in range(pre, T - run + 1):
        if t - last < run:
            continue
        base = p[t - pre:t]
        if not np.isfinite(base).all() or not np.isfinite(p[t]):
            continue
        if running[t - pre:t].any():
            continue
        if np.ptp(base) > stable:
            continue
        seg = p[t:min(T, t + rise_frames + 1)]
        if not np.isfinite(seg).any():
            continue
        amp = float(np.nanmax(seg) - base.mean())
        if amp < rise or not (p[t] > base[-1]):
            continue
        if excluded[t - pre:min(T, t + run)].any():
            dropped += 1; last = t
            continue
        still = not running[t:t + run].any()
        sp = np.asarray(b["speed"], float)[t:t + run]
        out.append({"trial": len(out) + 1, "onset": int(t), "onset_s": float(t * b["spf"]), "kind": "still" if still else "moving",
                    "pupil_rise": amp, "rise_frames_to_peak": int(np.nanargmax(seg)),
                    "run_fraction": float(running[t:t + run].mean()),
                    "speed_mean_run": float(np.nanmean(sp)) if np.isfinite(sp).any() else np.nan,
                    "speed_peak_run": float(np.nanmax(sp)) if np.isfinite(sp).any() else np.nan})
        last = t
    b["pupil_trials_dropped"] = dropped
    return out


# ----------------------------------------------------------------------------- the comparison

def _sel(tdf, condition: str, kind: str | None = None):
    if tdf is None or len(tdf) == 0:
        return tdf.iloc[0:0] if tdf is not None else None
    s = tdf[tdf["condition"] == condition]
    if kind is not None and "kind" in s.columns:
        s = s[s["kind"] == kind]
    return s


def _why_no_trials(trigger: str, run_rows: list[dict], condition: str) -> str:
    rr = [r for r in run_rows if r.get("condition") == condition]
    if not rr:
        return "no runs"
    if all(r.get("note") == "no locomotion data" for r in rr):
        return "no behaviour data"
    pct = [r.get("pct_running", np.nan) for r in rr]
    pct = [v for v in pct if np.isfinite(v)]
    if trigger == "running":
        if pct and max(pct) == 0:
            return "the animal never ran"
        if pct and min(pct) > 90:
            return f"the animal ran {np.mean(pct):.0f} % of the time: no still 10 s before a bout"
        return "no running bout with a still 10 s before it"
    if trigger == "pupil":
        if not any(r.get("has_pupil") for r in rr):
            return "no pupil data"
        return "no still pupil dilation (stable baseline, rise >= 20 % of the range)"
    return "less than 1 min of rest"


def compare(exp: str, conditions: list[str], loco_tdf, pupil_tdf, cond_rows: list[dict], run_rows: list[dict], spf: float) -> list[dict]:
    """One row per condition (vs the first = baseline) and trigger."""
    base = conditions[0]
    rows = []
    for c in conditions[1:]:
        for trig in TRIGGERS:
            row = {"experiment": exp, "condition": c, "baseline": base, "trigger": trig, "label": TRIGGER_LABEL[trig]}
            if trig in ("running", "pupil"):
                tdf = loco_tdf if trig == "running" else pupil_tdf
                kind = "still" if trig == "pupil" else None
                tb, tc = _sel(tdf, base, kind), _sel(tdf, c, kind)
                nb, nc = (len(tb) if tb is not None else 0), (len(tc) if tc is not None else 0)
                row.update({"n_base": nb, "n_cond": nc, "window_s": RUN_FRAMES * spf})
                mb, mc = nb * RUN_FRAMES * spf / 60.0, nc * RUN_FRAMES * spf / 60.0
                eb = int(tb["n_events_run"].sum()) if nb else 0; ec = int(tc["n_events_run"].sum()) if nc else 0
                row.update({"n_events_base": eb, "n_events_cond": ec, "minutes_base": mb, "minutes_cond": mc,
                            "rate_base": eb / mb if mb > 0 else np.nan, "rate_cond": ec / mc if mc > 0 else np.nan})
                if nb and nc:
                    cm = change_measures(eb, mb, ec, mc)
                    row.update({"rate_delta": cm.get("delta_rate"), "rate_ratio": (ec / mc) / (eb / mb) if eb > 0 else np.nan,
                                "rate_ratio_ci_lo": cm.get("ratio_ci_lo"), "rate_ratio_ci_hi": cm.get("ratio_ci_hi"),
                                "log2_ratio_pc": cm.get("log2_ratio_pc"), "modulation_index": cm.get("modulation_index"),
                                "low_counts": cm.get("low_counts")})
                    row.update(fraction_change(int(tb["responsive_z"].sum()), nb, int(tc["responsive_z"].sum()), nc))
                    row.update(amplitude_change(tb["pop_response"], tc["pop_response"]))
                testable = nb >= MIN_TRIALS and nc >= MIN_TRIALS
                why = ""
                if not testable:
                    short, reasons = [], []
                    for cond, n in ((base, nb), (c, nc)):
                        if n < MIN_TRIALS:
                            short.append(f"{cond} {n}")
                            if n == 0:
                                reasons.append(_why_no_trials(trig, run_rows, cond))
                    why = f"fewer than {MIN_TRIALS} trials ({', '.join(short)})"
                    if reasons:
                        why += " — " + "; ".join(dict.fromkeys(reasons))
            else:
                rb = next((x for x in cond_rows if x["condition"] == base), {}); rc = next((x for x in cond_rows if x["condition"] == c), {})
                eb, mb = rb.get("rest_n_events", np.nan), rb.get("rest_minutes", np.nan)
                ec, mc = rc.get("rest_n_events", np.nan), rc.get("rest_minutes", np.nan)
                row.update({"n_base": int(rb.get("n_runs", 0)), "n_cond": int(rc.get("n_runs", 0)), "n_events_base": eb, "n_events_cond": ec,
                            "minutes_base": mb, "minutes_cond": mc,
                            "rate_base": eb / mb if np.isfinite(mb) and mb > 0 else np.nan, "rate_cond": ec / mc if np.isfinite(mc) and mc > 0 else np.nan})
                testable = bool(np.isfinite(mb) and np.isfinite(mc) and mb >= MIN_REST_MIN and mc >= MIN_REST_MIN)
                why = "" if testable else f"less than {MIN_REST_MIN:.0f} min of rest ({base}: {mb if np.isfinite(mb) else 0:.1f} min, {c}: {mc if np.isfinite(mc) else 0:.1f} min)"
                if testable:
                    cm = change_measures(eb, mb, ec, mc)
                    row.update({"rate_delta": cm.get("delta_rate"), "rate_ratio": (ec / mc) / (eb / mb) if eb > 0 else np.nan,
                                "rate_ratio_ci_lo": cm.get("ratio_ci_lo"), "rate_ratio_ci_hi": cm.get("ratio_ci_hi"),
                                "log2_ratio_pc": cm.get("log2_ratio_pc"), "modulation_index": cm.get("modulation_index"),
                                "low_counts": cm.get("low_counts")})
            row["testable"] = bool(testable)
            row["why_not"] = why
            row["verdict"] = verdict(row)
            rows.append(row)
    return rows


def verdict(row: dict) -> str:
    """One phrase: what the evoked / resting rate says, scale-free; when the
    rate cannot speak (no or too few events), what the response amplitude
    of the trials says."""
    if not row.get("testable"):
        return "not testable: " + str(row.get("why_not", ""))
    lo, hi = row.get("rate_ratio_ci_lo", np.nan), row.get("rate_ratio_ci_hi", np.nan)
    nb, nc = row.get("n_events_base", 0), row.get("n_events_cond", 0)
    amp = ""
    pa = row.get("amp_p_mwu", np.nan)
    if np.isfinite(row.get("amp_delta", np.nan)):
        amp = (f"response amplitude {row.get('amp_baseline', np.nan):.3f} -> {row.get('amp_condition', np.nan):.3f} dF/F"
               + (f" (MWU p = {pa:.2g}{', differs' if pa < 0.05 else ''})" if np.isfinite(pa) else ""))
    if not np.isfinite(nb) or not np.isfinite(nc) or nb + nc == 0:
        return "no evoked events in either condition" + (f"; {amp}" if amp else "")
    up = np.isfinite(lo) and lo > 1
    down = np.isfinite(hi) and hi < 1
    if row.get("low_counts"):
        if up or down:
            out = f"{'increase' if up else 'decrease'} of the rate (interval excludes 1) but on {int(nb)} vs {int(nc)} events — direction only"
        else:
            out = f"too few events ({int(nb)} vs {int(nc)}) for a conclusion on the rate"
        return out + (f"; {amp}" if amp else "")
    if up:
        out = "an increase of the rate the interval supports"
    elif down:
        out = "a decrease of the rate the interval supports"
    else:
        out = "rate not distinguishable from no change"
    return out + (f"; {amp}" if amp and np.isfinite(pa) and pa < 0.05 else "")


def lines_for(exp: str, rows: list[dict]) -> list[str]:
    out = [f"  triggers compared ({exp}): running onsets / still pupil dilations / rest, same scale-aware statistics as the rate comparison"]
    for r in rows:
        if r["trigger"] == "rest":
            head = (f"    {r['condition']} vs {r['baseline']} — {r['label']}: {r.get('rate_base', np.nan):.3f} -> {r.get('rate_cond', np.nan):.3f} events/min "
                    f"({r.get('n_events_base', 0):.0f} in {r.get('minutes_base', 0):.1f} min -> {r.get('n_events_cond', 0):.0f} in {r.get('minutes_cond', 0):.1f} min)")
        else:
            head = (f"    {r['condition']} vs {r['baseline']} — {r['label']}: {r['n_base']} -> {r['n_cond']} trials; evoked "
                    f"{r.get('rate_base', np.nan):.3f} -> {r.get('rate_cond', np.nan):.3f} events/min in the response window")
            if np.isfinite(r.get("frac_baseline", np.nan)):
                head += (f"; responsive {100 * r['frac_baseline']:.0f} % -> {100 * r['frac_condition']:.0f} % (Fisher p = {r.get('frac_p_fisher', np.nan):.2g})"
                         f"; response {r.get('amp_baseline', np.nan):.3f} -> {r.get('amp_condition', np.nan):.3f} dF/F (MWU p = {r.get('amp_p_mwu', np.nan):.2g})")
        if r.get("testable") and np.isfinite(r.get("rate_ratio_ci_lo", np.nan)):
            hi = r.get("rate_ratio_ci_hi", np.nan)
            head += f"; ratio CI {r['rate_ratio_ci_lo']:.2f}-{'inf' if hi == np.inf else f'{hi:.2f}'}"
        out.append(head + f" -> {r['verdict']}")
    n_nt = sum(1 for r in rows if not r.get("testable"))
    if rows and n_nt == len(rows):
        out.append("    no trigger is testable in this experiment (see why above)")
    return out


# ----------------------------------------------------------------------------- figure

def _trace_panel(ax, traces: dict, conditions: list, cols: dict, kind: str | None, title: str, xlabel: str):
    drawn = False
    for c in conditions:
        tr = traces.get(c)
        if not tr or tr["pop"].shape[0] == 0:
            continue
        sel = np.array([k == kind for k in tr["kinds"]]) if kind else np.ones(len(tr["kinds"]), bool)
        M = tr["pop"][sel]
        if M.shape[0] == 0:
            continue
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                m = np.nanmean(M, axis=0)
                s = np.nanstd(M, axis=0, ddof=1) / np.sqrt(M.shape[0]) if M.shape[0] > 1 else np.zeros_like(m)
        ax.plot(tr["t_s"], m, color=cols[c], label=f"{c} (n={M.shape[0]})")
        ax.fill_between(tr["t_s"], m - s, m + s, color=cols[c], alpha=0.2, lw=0)
        drawn = True
    if drawn:
        t0 = next(iter(traces.values()))["t_s"]
        ax.axvspan(0, RUN_FRAMES * (t0[1] - t0[0]), color="#f2c94c", alpha=0.15, lw=0)
        ax.axvline(0, color="#999", lw=0.8, ls=":")
        ax.legend(fontsize=7, frameon=False)
    else:
        ax.text(0.5, 0.5, "no trials", ha="center", va="center", transform=ax.transAxes, fontsize=9, color="#888")
    ax.set_xlabel(xlabel, fontsize=8); ax.set_ylabel("population dF/F", fontsize=8); ax.set_title(title, fontsize=9); _style(ax)


def _rate_panel(ax, rows: list[dict], conditions: list, cols: dict, trig: str, base: str, ylabel: str, title: str):
    """Evoked (or resting) rate per condition with its exact Poisson interval; the comparison rows give the condition values."""
    from .stats import poisson_rate_ci
    xs, labels, vals = [], [], []
    for i, c in enumerate(conditions):
        r = next((x for x in rows if x["trigger"] == trig and (x["condition"] == c or (c == base and x["baseline"] == base))), None)
        if r is None:
            continue
        if c == base:
            n, m, rate = r.get("n_events_base", np.nan), r.get("minutes_base", np.nan), r.get("rate_base", np.nan)
        else:
            n, m, rate = r.get("n_events_cond", np.nan), r.get("minutes_cond", np.nan), r.get("rate_cond", np.nan)
        lo, hi = poisson_rate_ci(n, m) if np.isfinite(rate) else (np.nan, np.nan)
        vals.append((i, c, n, m, rate, lo, hi))
        xs.append(i); labels.append(c)
    rates = [v[4] for v in vals if np.isfinite(v[4])]
    top = max(max(rates) * 1.6, 0.2) if rates else 1.0
    for i, c, n, m, rate, lo, hi in vals:
        if not np.isfinite(rate):
            ax.text(i, 0, "n/a", ha="center", va="bottom", fontsize=7, color="#888")
            continue
        ax.bar(i, rate, color=cols[c], alpha=0.5, width=0.7)
        if np.isfinite(lo):
            clipped = hi > top
            ax.errorbar(i, rate, yerr=[[rate - lo], [min(hi, top) - rate]], color=cols[c], capsize=0 if clipped else 3, lw=1.2)
            if clipped:
                ax.annotate("", xy=(i, top), xytext=(i, top * 0.93), arrowprops=dict(arrowstyle="->", color=cols[c], lw=1.2))
        ax.text(i, min(rate, top * 0.9), f"{int(n)} ev / {m:.1f} min", ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, top)
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8); ax.set_title(title, fontsize=9); _style(ax)


def figure(exp: str, conditions: list, loco_traces: dict, pupil_traces: dict, rows: list[dict], out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(conditions)}
    base = conditions[0]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.2))
    _trace_panel(axes[0, 0], loco_traces, conditions, cols, None, "running onsets: onset-aligned population dF/F (mean ± sem)", "time from running onset (s)")
    _trace_panel(axes[0, 1], pupil_traces, conditions, cols, "still", "pupil dilations (still): onset-aligned population dF/F", "time from dilation onset (s)")
    _rate_panel(axes[0, 2], rows, conditions, cols, "rest", base, "events / min at rest", "resting activity (>= 20 s from running), exact 95 % CI")
    _rate_panel(axes[1, 0], rows, conditions, cols, "running", base, "kept events / min in the response window", "running onsets: evoked event rate, exact 95 % CI")
    _rate_panel(axes[1, 1], rows, conditions, cols, "pupil", base, "kept events / min in the response window", "pupil dilations (still): evoked event rate, exact 95 % CI")
    ax = axes[1, 2]; ax.axis("off")
    txt = []
    for r in rows:
        head = f"{r['condition']} vs {r['baseline']} — {r['label']}"
        if r["trigger"] != "rest":
            head += f" (n {r['n_base']} / {r['n_cond']})"
        txt.append(head)
        if r.get("testable") and np.isfinite(r.get("rate_ratio_ci_lo", np.nan)):
            hi = r.get("rate_ratio_ci_hi", np.nan)
            txt.append(f"   rate x{r.get('rate_ratio', np.nan):.2f} [{r['rate_ratio_ci_lo']:.2f}-{'inf' if hi == np.inf else f'{hi:.2f}'}], "
                       f"MI {r.get('modulation_index', np.nan):+.2f}" + (f", responsive {100 * r['frac_baseline']:.0f} -> {100 * r['frac_condition']:.0f} %"
                                                                        if np.isfinite(r.get("frac_baseline", np.nan)) else ""))
        txt.append("   -> " + r["verdict"])
    ax.text(0.0, 1.0, "\n".join(txt) if txt else "no comparison (one condition)", va="top", ha="left", fontsize=7.5, family="monospace",
            transform=ax.transAxes, wrap=True)
    ax.set_title("verdicts (exact Poisson interval of the rate ratio)", fontsize=9)
    fig.suptitle(f"{exp}: the three triggers compared — running onsets, still pupil dilations, resting activity "
                 "(a session without running is read on the pupil and the rest)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)


# ----------------------------------------------------------------------------- per experiment

def analyze_experiment(exp: str, bundles: list[dict], conditions: list[str], loco_result: dict | None, out_dir: Path, log=print) -> dict:
    """Pupil-dilation trials of every run (+ responses), the three-trigger
    comparison and its figure. loco_result = locotrials.analyze_experiment's
    return value (trial rows, cond rows, run rows, traces)."""
    import pandas as pd
    out_dir = Path(out_dir)
    pupil_rows, pupil_traces, run_notes = [], {}, []
    for b in bundles:
        if not b.get("has_pupil"):
            run_notes.append({"run": b["run"], "condition": b["condition"], "n_pupil_trials": 0, "note": "no pupil data"})
            continue
        if "running_merged" not in b:
            from .locotrials import bouts
            bouts(b)
        trials = pupil_trials(b)
        if not trials:
            run_notes.append({"run": b["run"], "condition": b["condition"], "n_pupil_trials": 0, "note": "no dilation onset found"})
            continue
        resp = trial_responses(b, trials)
        for r in resp["rows"]:
            r["experiment"] = exp
        pupil_rows += resp["rows"]
        run_notes.append({"run": b["run"], "condition": b["condition"], "n_pupil_trials": len(trials),
                          "n_pupil_trials_still": sum(1 for t in trials if t["kind"] == "still"),
                          "n_pupil_trials_dropped_excluded": int(b.get("pupil_trials_dropped", 0))})
        tr = pupil_traces.setdefault(b["condition"], {"pop": [], "kinds": [], "t_s": resp["t_s"]})
        tr["pop"].append(resp["pop_traces"]); tr["kinds"] += [r["kind"] for r in resp["rows"]]
    for c, tr in pupil_traces.items():
        tr["pop"] = np.vstack([x for x in tr["pop"] if x.size]) if any(x.size for x in tr["pop"]) else np.zeros((0, tr["t_s"].size))
    pdf = pd.DataFrame(pupil_rows) if pupil_rows else pd.DataFrame(columns=["run", "condition", "kind", "n_events_run", "responsive_z", "pop_response"])
    if pupil_rows:
        pdf.to_csv(out_dir / f"pupil_trials_{exp}.csv", index=False)
    lr = loco_result or {}
    ltdf = pd.DataFrame(lr.get("trial_rows") or []) if lr.get("trial_rows") else pd.DataFrame(columns=["run", "condition", "kind", "n_events_run", "responsive_z", "pop_response"])
    run_rows = list(lr.get("run_rows") or [])
    for n in run_notes:                                   # has_pupil for the why-not text
        for r in run_rows:
            if r.get("run") == n["run"]:
                r["has_pupil"] = n.get("note") != "no pupil data"
    spf = float(bundles[0]["spf"]) if bundles else 2.0
    rows = compare(exp, conditions, ltdf, pdf, lr.get("cond_rows") or [], run_rows, spf)
    lines = lines_for(exp, rows)
    n_still = int((pdf["kind"] == "still").sum()) if len(pdf) else 0
    lines.insert(1, f"    pupil dilations: {len(pdf)} onsets in {sum(1 for n in run_notes if n.get('n_pupil_trials'))} run(s), {n_still} still "
                    f"(responsive by z rule: {int(pdf.loc[pdf['kind'] == 'still', 'responsive_z'].sum()) if n_still else 0})")
    try:
        figure(exp, conditions, lr.get("traces") or {}, pupil_traces, rows, out_dir / f"{exp}_triggers.png")
    except Exception as e:
        log(f"    trigger figure failed ({e})")
    return {"rows": rows, "pupil_trial_rows": pupil_rows, "run_notes": run_notes, "lines": lines, "n_pupil_trials": len(pdf), "n_still": n_still}
