"""Time course of the change after the stimulation (optostim paradigm).

In the RearChr protocol the two fields of view alternate (FOV1, FOV2, FOV1,
FOV2 before the stimulation and again after it), so the post runs of a day sit
at four different times after the stimulation (roughly 0-16, 17-33, 34-50 and
51-67 min) and the pre runs at four times before it. Putting every run on the
day's clock answers three things the pooled pre/post ratio cannot:

  * does the change decay, build up or stay flat over the hour after the
    stimulation (slope of the rate over time since the first post run, with
    the behaviour taken out),
  * is the step at the stimulation bigger than the drift between consecutive
    runs of the same field within a condition (otherwise a "change" may be the
    session's trend rather than the stimulation),
  * how long the gap between the last pre run and the first post run was.

Clock: the end of every recording is the modification time of its .sbx (or
_quadrature.mat / .mat) next to the movie; the start is that minus the
recorded duration. When no such file is next to the movie, a run_times.csv
(columns run, end — ISO time or epoch seconds — and optionally start) in the
mouse folder or in aqua_analysis is used; otherwise the runs are placed in
order with assumed gaps and the result is marked "estimated".
"""
from __future__ import annotations

import datetime as _dt
import json
import re
from pathlib import Path

import numpy as np

from .adjust import bin_run, poisson_fit

BIN_TAG = re.compile(r"_bin\d+f_[\d.p]+s$")


def run_base(run_name: str) -> str:
    return BIN_TAG.sub("", run_name)


def day_of(run_name: str) -> str | None:
    m = re.search(r"_(\d{6})_\d{3}", run_name)
    return m.group(1) if m else None


def cond_day(cond: str) -> str | None:
    """The date in a day condition: 'd261010' and 'd261010_post' (v44 run_labels.csv) -> '261010'; None otherwise."""
    m = re.match(r"^d(\d{6})(?:$|_)", str(cond))
    return m.group(1) if m else None


def cond_label(cond: str, day_labels: dict) -> str:
    """The day label of a condition: day_labels['261010_post'] when the file names the exact
    condition, else day_labels['261010'] + the condition's suffix ('D14 post'); '' when none."""
    c = str(cond)
    key = c[1:] if c.startswith("d") else c
    if day_labels.get(key):
        return day_labels[key]
    d = cond_day(c)
    if d and day_labels.get(d):
        suffix = key[len(d):].strip("_").replace("_", " ")
        return day_labels[d] + (f" {suffix}" if suffix else "")
    return ""


def _read_run_times(mouse_dir: Path) -> dict:
    """{run base: (start_s or None, end_s)} from run_times.csv if present."""
    import pandas as pd
    out = {}
    for p in (mouse_dir / "run_times.csv", mouse_dir / "aqua_analysis" / "run_times.csv"):
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p)
        except Exception:
            continue

        def to_s(v):
            if v is None or (isinstance(v, float) and np.isnan(v)):
                return None
            try:
                return float(v)
            except (TypeError, ValueError):
                t = pd.Timestamp(str(v))
                if t.tzinfo is None:
                    t = t.tz_localize("UTC")
                return float(t.timestamp())
        for _, r in df.iterrows():
            end = to_s(r.get("end")); start = to_s(r.get("start")) if "start" in df.columns else None
            if end is not None:
                out[run_base(str(r["run"]))] = (start, end)
        if out:
            break
    return out


def run_clock(tif: Path, info: dict, mouse_dir: Path, cache: dict | None = None) -> dict | None:
    """{'start_s', 'end_s', 'duration_s', 'source'} for one run, or None."""
    tif = Path(tif)
    base = run_base(tif.stem)
    duration = None
    jp = tif.with_suffix(".json")
    if jp.exists():
        try:
            j = json.loads(jp.read_text())
            if j.get("frames_in") and j.get("fs_in_hz"):
                duration = float(j["frames_in"]) / float(j["fs_in_hz"])
        except Exception:
            pass
    if duration is None:
        duration = float(info.get("frames", 0)) * float(info.get("seconds_per_frame", 2.0))
    for suffix, src in ((".sbx", "sbx file time"), ("_quadrature.mat", "quadrature file time"), (".mat", "mat file time")):
        p = tif.parent / (base + suffix)
        if p.exists():
            end = p.stat().st_mtime
            return {"start_s": end - duration, "end_s": end, "duration_s": duration, "source": src}
    if cache is None:
        cache = _read_run_times(mouse_dir)
    if base in cache:
        start, end = cache[base]
        if start is None:
            start = end - duration
        return {"start_s": start, "end_s": end, "duration_s": duration, "source": "run_times.csv"}
    return None


# ----------------------------------------------------------------------------- per day

def read_day_labels(mouse_dir: Path, extra: dict | None = None) -> dict:
    """{day: label} from <mouse>/day_labels.csv (columns day, label) plus any given on the command line."""
    out = {}
    p = Path(mouse_dir) / "day_labels.csv"
    if p.exists():
        try:
            import pandas as pd
            df = pd.read_csv(p, dtype=str)
            for _, r in df.iterrows():
                if str(r.get("day", "")).strip():
                    out[str(r["day"]).strip()] = str(r.get("label", "")).strip()
        except Exception:
            pass
    out.update(extra or {})
    return out


def day_timecourse(day: str, exps: list, out_dir: Path, mouse_dir: Path, bin_min: float = 2.0,
                   assumed_gap_min: float = 0.8, assumed_stim_gap_min: float = 73.0,
                   n_boot: int = 300, label: str | None = None, log=print) -> dict | None:
    """exps: [{experiment, fov, conditions, base, spf, runs: [{run, condition, tif, info,
    onsets, excluded, speed, still, pupil, n_events, minutes_used, props}]}] of one day.
    Writes nothing itself; returns {summary, run_rows, bin_rows, lines, fig}."""
    cond_post = None
    runs = []
    for e in exps:
        for r in e["runs"]:
            rr = dict(r); rr["fov"] = e["fov"]; rr["experiment"] = e["experiment"]; rr["spf"] = e["spf"]
            rr["is_base"] = r["condition"] == e["base"]
            runs.append(rr)
            if not rr["is_base"]:
                cond_post = r["condition"]
    if not runs or cond_post is None:
        return None
    runs.sort(key=lambda r: r["run"])
    # ---- clock
    cache = _read_run_times(mouse_dir)
    clocks = [run_clock(Path(r["tif"]), r["info"], mouse_dir, cache) for r in runs]
    if all(c is not None for c in clocks):
        starts = np.array([c["start_s"] for c in clocks]); ends = np.array([c["end_s"] for c in clocks])
        src = sorted({c["source"] for c in clocks})
        estimated = False
    else:
        estimated = True
        src = ["ESTIMATED (no .sbx / run_times.csv next to the movies)"]
        starts = np.zeros(len(runs)); ends = np.zeros(len(runs))
        t = 0.0
        prev_base = None
        for i, r in enumerate(runs):
            dur = r["info"].get("frames", 0) * r["spf"]
            if prev_base is not None and prev_base and not r["is_base"]:
                t += assumed_stim_gap_min * 60.0
            starts[i] = t; ends[i] = t + dur
            t = ends[i] + assumed_gap_min * 60.0
            prev_base = r["is_base"]
    t_first = starts.min()
    post_idx = [i for i, r in enumerate(runs) if not r["is_base"]]
    pre_idx = [i for i, r in enumerate(runs) if r["is_base"]]
    if not post_idx or not pre_idx:
        return None
    t0_post = starts[post_idx].min()
    gap_min = (t0_post - ends[pre_idx].max()) / 60.0
    # ---- per-FOV baselines
    fovs = sorted({r["fov"] for r in runs}, key=lambda v: (v is None, v))
    pre_rate = {}
    for f in fovs:
        ev = sum(r["n_events"] for r in runs if r["fov"] == f and r["is_base"])
        mn = sum(r["minutes_used"] for r in runs if r["fov"] == f and r["is_base"])
        pre_rate[f] = ev / mn if mn > 0 else np.nan
    # ---- run rows
    run_rows = []
    for i, r in enumerate(runs):
        rate = r["n_events"] / r["minutes_used"] if r["minutes_used"] > 0 else np.nan
        pr = pre_rate.get(r["fov"], np.nan)
        row = {"day": day, "run": r["run"], "fov": r["fov"], "condition": r["condition"],
               "clock_source": "estimated" if estimated else ", ".join(src),
               "start_utc": (None if estimated else _dt.datetime.fromtimestamp(starts[i], _dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")),
               "t_session_min": round((starts[i] - t_first) / 60.0, 1),
               "t_since_first_post_min": round((starts[i] - t0_post) / 60.0, 1),
               "minutes_used": round(r["minutes_used"], 2), "n_events": int(r["n_events"]),
               "events_per_min": round(rate, 4) if np.isfinite(rate) else np.nan,
               "pre_rate_fov": round(pr, 4) if np.isfinite(pr) else np.nan,
               "ratio_to_pre_fov": round(rate / pr, 3) if (np.isfinite(rate) and np.isfinite(pr) and pr > 0) else np.nan}
        for k in ("still_fraction", "velocity_median", "pupil_median", "dff_max_median", "duration_s_median", "area_um2_median"):
            if k in r.get("props", {}):
                row[k] = r["props"][k]
        run_rows.append(row)
    # ---- bins on the day's clock
    bin_rows = []
    for i, r in enumerate(runs):
        t = bin_run(r["onsets"], r["excluded"], r["speed"], r["still"], r["pupil"], r["spf"], bin_s=bin_min * 60.0)
        n = t["events"].size
        for b in range(n):
            if not np.isfinite(t["minutes"][b]) or t["minutes"][b] <= 0.25 * bin_min:
                continue
            tmid = starts[i] + (b + 0.5) * bin_min * 60.0
            pr = pre_rate.get(r["fov"], np.nan)
            rate = t["events"][b] / t["minutes"][b]
            bin_rows.append({"day": day, "run": r["run"], "fov": r["fov"], "condition": r["condition"], "run_index": i,
                             "t_session_min": (tmid - t_first) / 60.0, "t_since_first_post_min": (tmid - t0_post) / 60.0,
                             "events": float(t["events"][b]), "minutes": float(t["minutes"][b]), "events_per_min": rate,
                             "ratio_to_pre_fov": rate / pr if (np.isfinite(pr) and pr > 0) else np.nan,
                             "speed": float(t["speed"][b]) if np.isfinite(t["speed"][b]) else np.nan,
                             "running": float(t["running"][b]) if np.isfinite(t["running"][b]) else np.nan,
                             "pupil": float(t["pupil"][b]) if np.isfinite(t["pupil"][b]) else np.nan})
    # ---- post slope (per-FOV intercepts, time since the first post run, behaviour)
    summary = {"day": day, "clock": "estimated" if estimated else ", ".join(src), "n_runs": len(runs),
               "gap_pre_to_post_min": round(gap_min, 1),
               "post_span_min": round((ends[post_idx].max() - t0_post) / 60.0, 1),
               "n_events_pre": int(sum(runs[i]["n_events"] for i in pre_idx)),
               "n_events_post": int(sum(runs[i]["n_events"] for i in post_idx))}
    post_bins = [b for b in bin_rows if b["condition"] == cond_post]
    slope = _slope(post_bins, fovs, n_boot=n_boot)
    summary.update({k: v for k, v in slope.items() if k != "slope_ci"})
    summary["slope_ci"] = slope.get("slope_ci")
    # ---- early vs late post run per FOV, and the step vs the drift
    steps = []
    for f in fovs:
        rf = [i for i, r in enumerate(runs) if r["fov"] == f]
        for a, b in zip(rf[:-1], rf[1:]):
            ra = run_rows[a]["events_per_min"]; rb = run_rows[b]["events_per_min"]
            kind = ("pre->pre" if runs[a]["is_base"] and runs[b]["is_base"] else
                    "post->post" if not runs[a]["is_base"] and not runs[b]["is_base"] else "pre->post (stimulation)")
            steps.append({"day": day, "fov": f, "from_run": runs[a]["run"], "to_run": runs[b]["run"], "kind": kind,
                          "rate_from": ra, "rate_to": rb,
                          "ratio": (rb / ra) if (np.isfinite(ra) and np.isfinite(rb) and ra > 0) else np.nan,
                          "gap_min": round((starts[b] - ends[a]) / 60.0, 1),
                          "events_from": run_rows[a]["n_events"], "events_to": run_rows[b]["n_events"]})
    early = [i for f in fovs for i in [j for j in post_idx if runs[j]["fov"] == f][:1]]
    late = [i for f in fovs for i in [j for j in post_idx if runs[j]["fov"] == f][1:2]]
    def pooled(ix):
        ev = sum(runs[i]["n_events"] for i in ix); mn = sum(runs[i]["minutes_used"] for i in ix)
        return ev / mn if mn > 0 else np.nan
    re_, rl_ = pooled(early), pooled(late)
    summary["rate_post_early"] = round(re_, 4) if np.isfinite(re_) else np.nan
    summary["rate_post_late"] = round(rl_, 4) if np.isfinite(rl_) else np.nan
    summary["ratio_late_vs_early"] = round(rl_ / re_, 3) if (np.isfinite(re_) and re_ > 0 and np.isfinite(rl_)) else np.nan
    summary["t_early_min"] = round(float(np.mean([(starts[i] + ends[i]) / 2 - t0_post for i in early]) / 60.0), 1) if early else np.nan
    summary["t_late_min"] = round(float(np.mean([(starts[i] + ends[i]) / 2 - t0_post for i in late]) / 60.0), 1) if late else np.nan
    def defined(s):          # a ratio with events on both sides (log-comparable)
        return np.isfinite(s["ratio"]) and s["ratio"] > 0 and s["events_from"] > 0 and s["events_to"] > 0
    stim_steps = [s for s in steps if s["kind"].startswith("pre->post") and defined(s)]
    drift_steps = [s for s in steps if not s["kind"].startswith("pre->post") and defined(s)]
    if stim_steps:
        summary["stim_step_log2_mean"] = round(float(np.mean([np.log2(s["ratio"]) for s in stim_steps])), 3)
    if drift_steps:
        summary["drift_log2_abs_max"] = round(float(np.max([abs(np.log2(s["ratio"])) for s in drift_steps])), 3)
    # ---- text
    lines = []
    summary["label"] = label or ""
    lines.append(f"time course, day {day}{' (' + label + ')' if label else ''} ({'clock ' + ', '.join(src) if not estimated else 'clock ESTIMATED: runs placed in order with assumed gaps'}): "
                 f"first post run started {gap_min:.0f} min after the last pre run ended; the post runs span "
                 f"{summary['post_span_min']:.0f} min")
    if slope.get("slope_per_10min") is not None:
        ci = slope.get("slope_ci")
        adj = slope.get("slope_adj_per_10min")
        trend = ("builds up" if slope["slope_per_10min"] > 1.1 else "decays" if slope["slope_per_10min"] < 0.9 else "is flat")
        lines.append(f"  post period ({summary['n_events_post']} events): the rate changes x{slope['slope_per_10min']:.2f} per 10 min since the first post run"
                     + (f" (95 % CI {ci[0]:.2f}-{ci[1]:.2f})" if ci else "")
                     + (f"; with the behaviour taken out x{adj:.2f}" if adj is not None else "")
                     + f" — over the hour the effect {trend}"
                     + (" (interval includes 1)" if ci and ci[0] <= 1 <= ci[1] else "")
                     + (f"; only {summary['n_events_post']} post events — a weak estimate" if summary["n_events_post"] < 30 else ""))
    else:
        lines.append(f"  post period: slope not fitted ({slope.get('why_not', 'too few events')})")
    if np.isfinite(summary.get("ratio_late_vs_early", np.nan)):
        per = "; ".join(f"FOV{f} x{run_rows[l]['events_per_min'] / run_rows[e]['events_per_min']:.2f}"
                        f" ({run_rows[e]['n_events']}->{run_rows[l]['n_events']} events)"
                        for f, e, l in [(runs[e]["fov"], e, l) for e, l in zip(early, late)]
                        if run_rows[e]["events_per_min"] and run_rows[e]["events_per_min"] > 0)
        lines.append(f"  second post run vs first (about {summary['t_late_min']:.0f} vs {summary['t_early_min']:.0f} min after the first post run): "
                     f"pooled x{summary['ratio_late_vs_early']:.2f}" + (f"; {per}" if per else ""))
    for f in fovs:
        sf = [s for s in steps if s["fov"] == f]
        if not sf:
            continue
        txt = ", ".join(f"{s['kind'].split(' ')[0]} {s['rate_from']:.2f}->{s['rate_to']:.2f}/min"
                        + (f" (x{s['ratio']:.2f})" if np.isfinite(s["ratio"]) else "") for s in sf)
        lines.append(f"  FOV{f} run by run: {txt}")
    if stim_steps and drift_steps:
        st = float(np.mean([np.log2(s["ratio"]) for s in stim_steps]))
        dr = float(np.max([abs(np.log2(s["ratio"])) for s in drift_steps]))
        worst = max(drift_steps, key=lambda s: abs(np.log2(s["ratio"])))
        lines.append(f"  step at the stimulation (same field, last pre run -> first post run; mean over "
                     f"{len(stim_steps)} field(s) with events on both sides) x{2 ** st:.2f} vs the largest change between "
                     f"consecutive runs of one field within a condition x{worst['ratio']:.2f} ({worst['kind']}, FOV{worst['fov']}): "
                     + ("the step is NOT larger than the within-condition drift — a pre/post difference of this size can be the session's trend"
                        if abs(st) <= dr else "the step exceeds the within-condition drift"))
    elif stim_steps:
        lines.append("  step at the stimulation: no within-condition drift to compare with (consecutive runs with events on both sides)")
    else:
        lines.append("  step at the stimulation: undefined (no field has events in both its last pre and its first post run)")
    return {"summary": summary, "run_rows": run_rows, "bin_rows": bin_rows, "steps": steps, "lines": lines,
            "fovs": fovs, "pre_rate": pre_rate, "t0_post_rel_min": (t0_post - t_first) / 60.0,
            "pre_end_rel_min": (ends[pre_idx].max() - t_first) / 60.0, "estimated": estimated, "cond_post": cond_post}


def _slope(post_bins: list, fovs: list, n_boot: int = 300, min_events: int = 10, rng=None) -> dict:
    """Poisson regression of the post bins: per-FOV intercepts + time since the
    first post run (per 10 min) [+ behaviour]. Returns the multiplicative change
    per 10 min, its block-bootstrap CI and the behaviour-adjusted version."""
    rng = rng or np.random.default_rng(0)
    out = {}
    if not post_bins:
        out["why_not"] = "no post bins"; return out
    y = np.array([b["events"] for b in post_bins]); mn = np.array([b["minutes"] for b in post_bins])
    if y.sum() < min_events:
        out["why_not"] = f"only {int(y.sum())} post events"; return out
    t10 = np.array([b["t_since_first_post_min"] for b in post_bins]) / 10.0
    fov = np.array([b["fov"] for b in post_bins], dtype=object)
    run_ix = np.array([b["run_index"] for b in post_bins])
    cols = [(fov == f).astype(float) for f in fovs if (fov == f).any()]
    X0 = np.column_stack(cols + [t10])
    off = np.log(mn)
    beta, _, _ = poisson_fit(X0, y, off, ridge=0.01)
    out["slope_per_10min"] = float(np.exp(beta[-1]))
    # behaviour-adjusted
    sp = np.array([b["speed"] for b in post_bins]); ru = np.array([b["running"] for b in post_bins])
    pu = np.array([b["pupil"] for b in post_bins])
    have_sp = np.isfinite(sp).mean() > 0.8; have_pu = np.isfinite(pu).mean() > 0.8
    bcols = []
    if have_sp:
        lsp = np.log1p(np.maximum(np.nan_to_num(sp, nan=0.0), 0)); bcols.append((lsp - lsp.mean()) / (lsp.std() + 1e-9))
        bcols.append(np.nan_to_num(ru, nan=0.0))
    if have_pu:
        pz = np.nan_to_num(pu, nan=np.nanmean(pu)); bcols.append((pz - pz.mean()) / (pz.std() + 1e-9))
    if bcols:
        X1 = np.column_stack(cols + [t10] + bcols)
        b1, _, _ = poisson_fit(X1, y, off, ridge=0.01)
        out["slope_adj_per_10min"] = float(np.exp(b1[len(cols)]))
    # block bootstrap: blocks of 3 bins within runs
    blk = run_ix * 1000 + np.arange(len(post_bins)) // 3
    ids = np.unique(blk)
    by_run = {}
    for b in ids:
        by_run.setdefault(int(b // 1000), []).append(b)
    boots = []
    for _ in range(n_boot):
        pick = np.concatenate([rng.choice(np.asarray(bl), len(bl), replace=True) for bl in by_run.values()])
        sel = np.concatenate([np.flatnonzero(blk == b) for b in pick])
        Xb = X0[sel]
        keep = [j for j in range(Xb.shape[1] - 1) if Xb[:, j].any()] + [Xb.shape[1] - 1]
        try:
            bb, _, _ = poisson_fit(Xb[:, keep], y[sel], off[sel], max_iter=30, ridge=0.01)
            boots.append(float(np.exp(bb[-1])))
        except Exception:
            continue
    if len(boots) >= 20:
        out["slope_ci"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
        out["slope_ci_lo"], out["slope_ci_hi"] = out["slope_ci"]
    out["n_post_bins"] = int(len(post_bins))
    return out


# ----------------------------------------------------------------------------- figure

def timecourse_figure(days: list, mouse: str, out_png: Path, labels: dict | None = None):
    """days: [day_timecourse results]. One column per day: absolute rates on the
    session clock (top), ratio to the field's pre rate against the time since
    the first post run (bottom, log scale) with the fitted slope."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    days = [d for d in days if d]
    if not days:
        return None
    n = len(days)
    fig, axes = plt.subplots(2, n, figsize=(5.2 * n, 7.2), squeeze=False)
    fcol = {1: "#3b7dd8", 2: "#e8702a", 3: "#2ca02c", 4: "#9467bd"}
    for k, d in enumerate(days):
        ax = axes[0, k]; axb = axes[1, k]
        day = d["summary"]["day"]
        lab = (labels or {}).get(day) or d["summary"].get("label")
        title = f"{day}  —  {lab}" if lab else day
        rb = d["bin_rows"]; rr = d["run_rows"]
        ymax = max([b["events_per_min"] for b in rb] + [0.1])
        ax.axvspan(0, d["pre_end_rel_min"], color="#dfe8f7", alpha=0.6, lw=0)
        ax.axvspan(d["t0_post_rel_min"], max(b["t_session_min"] for b in rb) + 1, color="#fbe3d4", alpha=0.6, lw=0)
        ax.text((d["pre_end_rel_min"] + d["t0_post_rel_min"]) / 2, ymax * 1.08,
                f"stimulation\n(gap {d['summary']['gap_pre_to_post_min']:.0f} min)", ha="center", va="bottom", fontsize=7.5, color="#555")
        for f in d["fovs"]:
            c = fcol.get(f, "k")
            bf = [b for b in rb if b["fov"] == f]
            ax.plot([b["t_session_min"] for b in bf], [b["events_per_min"] for b in bf], ".", color=c, ms=3, alpha=0.5)
            runs_f = [r for r in rr if r["fov"] == f]
            xs = [r["t_session_min"] + r["minutes_used"] / 2 for r in runs_f]; ys = [r["events_per_min"] for r in runs_f]
            ax.plot(xs, ys, "o-", color=c, ms=7, lw=1.4, label=f"FOV{f}", mec="black", mew=0.5)
            for x, yv, r in zip(xs, ys, runs_f):
                ax.text(x, yv, f" {r['n_events']}", fontsize=6.5, color=c, va="bottom")
        # running fraction per bin on a second axis
        ax2 = ax.twinx()
        for ri in sorted({b["run_index"] for b in rb}):
            br_ = [b for b in rb if b["run_index"] == ri]
            ax2.plot([b["t_session_min"] for b in br_], [b["running"] for b in br_], "-", color="0.6", lw=0.8, alpha=0.7)
        ax2.set_ylim(0, 1.05); ax2.set_ylabel("running fraction (grey)", fontsize=7, color="0.5"); ax2.tick_params(labelsize=7, colors="0.5")
        ax.set_ylim(0, ymax * 1.35 + 1e-3)
        ax.set_xlabel("minutes since the first pre recording", fontsize=8)
        ax.set_ylabel("events / min (dots = 2-min bins, circles = runs)", fontsize=8)
        ax.set_title(f"{title}\nclock: {d['summary']['clock']}", fontsize=9, loc="left")
        ax.tick_params(labelsize=8); ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        ax.spines[["top"]].set_visible(False)
        # bottom: ratio vs time since the first post run
        pb = [b for b in rb if b["condition"] == d["cond_post"] and np.isfinite(b["ratio_to_pre_fov"])]
        for f in d["fovs"]:
            c = fcol.get(f, "k")
            bf = [b for b in pb if b["fov"] == f]
            if bf:
                axb.plot([b["t_since_first_post_min"] for b in bf], [max(b["ratio_to_pre_fov"], 0.05) for b in bf],
                         ".", color=c, ms=3, alpha=0.45)
            runs_f = [r for r in rr if r["fov"] == f and r["condition"] == d["cond_post"] and np.isfinite(r["ratio_to_pre_fov"])]
            if runs_f:
                axb.plot([r["t_since_first_post_min"] + r["minutes_used"] / 2 for r in runs_f],
                         [max(r["ratio_to_pre_fov"], 0.05) for r in runs_f], "o", color=c, ms=7, mec="black", mew=0.5)
        s = d["summary"]
        if s.get("slope_per_10min") and pb:
            t = np.linspace(min(b["t_since_first_post_min"] for b in pb), max(b["t_since_first_post_min"] for b in pb), 20)
            # anchor the line at the geometric mean of the run ratios
            rr_post = [r["ratio_to_pre_fov"] for r in rr if r["condition"] == d["cond_post"] and np.isfinite(r["ratio_to_pre_fov"]) and r["ratio_to_pre_fov"] > 0]
            if rr_post:
                g = np.exp(np.mean(np.log(rr_post))); tm = np.mean([r["t_since_first_post_min"] + r["minutes_used"] / 2 for r in rr if r["condition"] == d["cond_post"]])
                axb.plot(t, g * s["slope_per_10min"] ** ((t - tm) / 10.0), "-", color="k", lw=1.2)
            ci = s.get("slope_ci")
            axb.set_title(f"slope x{s['slope_per_10min']:.2f} per 10 min" + (f" (CI {ci[0]:.2f}-{ci[1]:.2f})" if ci else "")
                          + (f", behaviour-adjusted x{s['slope_adj_per_10min']:.2f}" if s.get("slope_adj_per_10min") else "")
                          + (f"; late vs early run x{s['ratio_late_vs_early']:.2f}" if np.isfinite(s.get("ratio_late_vs_early", np.nan)) else "")
                          + (f"\n{s['n_events_post']} post events" + (" — weak" if s["n_events_post"] < 30 else "")),
                          fontsize=8, loc="left")
        else:
            axb.set_title(s.get("why_not", "no slope (too few post events)"), fontsize=8, loc="left")
        axb.axhline(1, color="#444", lw=1)
        axb.set_yscale("log")
        axb.set_xlabel("minutes since the first post recording", fontsize=8)
        axb.set_ylabel("rate / the field's pre rate (log; empty bins at the floor)", fontsize=8)
        axb.tick_params(labelsize=8); axb.grid(axis="y", alpha=0.3)
        axb.spines[["top", "right"]].set_visible(False)
        pre_missing = [f for f in d["fovs"] if not (np.isfinite(d["pre_rate"].get(f, np.nan)) and d["pre_rate"][f] > 0)]
        if pre_missing:
            axb.text(0.98, 0.03, "no pre events in FOV " + ", ".join(str(f) for f in pre_missing) + " — ratio undefined",
                     transform=axb.transAxes, ha="right", fontsize=7, color="#a00")
    fig.suptitle(f"{mouse}: time course of the event rate over the day — does the change after the stimulation decay, build up or continue a trend?",
                 fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png
