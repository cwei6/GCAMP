"""The general representation of the data — the measurements of Yu Yong's
lab-meeting slides (LabMeeting_261012) computed per run / condition so the
report can open with them:

  G1 session overview: event-count histogram on the day's clock with the
     running speed (state shading, bouts, locomotion trials) and the pupil
     under it (slides 19/20/28);
  G2 rates by behavioural state: events/min and event max dF/F in running vs
     resting frames, % running, event overlays running vs resting (s4/s19);
  G3 ECDFs of the normalised event amplitude, duration and area, running vs
     resting x condition (s4);
  G6 pupil measures: pupil 0-1 per FOV and day, variation (SD), resting
     pupil (valid when >= 10 % resting frames), spontaneous dilations at rest
     (s13-s15);
  G8 acute (post / pre, RearChr) and chronic (pre Day# / pre Day1, both
     paradigms) summary of every measure vs day, per FOV.

Everything works on afferent2p.runbundle dicts. The csvs written here are what
the report reads; nothing is recomputed there.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np

from .locotrials import REST_MARGIN_FRAMES, _sem, _style
from .timecourse import cond_label

BASE_COL, COND_COLS = "#8d8d8d", ["#d7263d", "#f46036", "#1b998b", "#7b2cbf", "#2a78d6", "#c08400"]
RUN_COL, REST_COL = "#d7263d", "#2a78d6"


def cond_colour(conditions: list, c: str) -> str:
    i = conditions.index(c) if c in conditions else 0
    return BASE_COL if i == 0 else COND_COLS[(i - 1) % len(COND_COLS)]


# ----------------------------------------------------------------------------- clock

def place_runs(bundles: list[dict], mouse_dir: Path, stim_gap_min: float = 73.0, gap_min: float = 0.8,
               block_gap_min: float = 6.0, blocks: list | None = None) -> dict:
    """Start time (s) of every run on one axis. Real clock from
    timecourse.run_clock when every run has one; else runs in order with
    assumed gaps (optostim: stim_gap_min between the last baseline run and the
    first other run; blocks: block_gap_min between the blocks, used for
    spontaneous experiments where a block = a day)."""
    from .timecourse import run_clock, _read_run_times
    cache = _read_run_times(mouse_dir)
    order = sorted(range(len(bundles)), key=lambda i: bundles[i]["run"])
    clocks = [run_clock(Path(bundles[i]["tif"]), bundles[i]["info"], mouse_dir, cache) if bundles[i].get("tif") else None for i in order]
    starts = np.zeros(len(bundles))
    estimated = not all(c is not None for c in clocks) or blocks is not None
    if not estimated:
        t0 = min(c["start_s"] for c in clocks)
        for i, c in zip(order, clocks):
            starts[i] = c["start_s"] - t0
        return {"starts": starts, "estimated": False, "source": sorted({c["source"] for c in clocks})}
    t = 0.0
    prev = None
    for i in order:
        b = bundles[i]
        key = (b["condition"] if blocks is None else b.get("block", b["condition"]))
        if prev is not None and key != prev:
            t += (block_gap_min if blocks is not None else stim_gap_min) * 60.0
        starts[i] = t
        t += b["T"] * b["spf"] + gap_min * 60.0
        prev = key
    return {"starts": starts, "estimated": True, "source": ["estimated"]}


# ----------------------------------------------------------------------------- G1 session overview

def session_figure(title: str, rows: list[dict], starts_by_run: dict, out_png: Path, conditions: list, base: str,
                   trials_by_run: dict | None = None, stim_span: tuple | None = None, bin_min: float = 1.0,
                   day_labels: dict | None = None, estimated: bool = False, block_labels: dict | None = None):
    """rows: [{label, bundles}] — one histogram row per field of view (optostim: both
    FOVs of the day) or one per experiment; the speed and pupil under them come
    from every run (the same animal)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .locotrials import bouts as _bouts
    all_b = [b for r in rows for b in r["bundles"]]
    if not all_b:
        return
    t_end = max(starts_by_run[b["run"]] + b["T"] * b["spf"] for b in all_b) / 60.0
    has_pupil = any(b["has_pupil"] for b in all_b)
    n_hist = len(rows)
    heights = [1.0] * n_hist + [0.8] + ([0.6] if has_pupil else [])
    fig, axes = plt.subplots(len(heights), 1, figsize=(15, 1.9 * len(heights) + 1.2), sharex=True,
                             gridspec_kw={"height_ratios": heights})
    axes = np.atleast_1d(axes)
    edges = np.arange(0, np.ceil(t_end / bin_min) * bin_min + bin_min, bin_min)
    ymax = 1
    for ax, r in zip(axes[:n_hist], rows):
        for b in r["bundles"]:
            on = np.asarray(b["onsets"], float)
            t = starts_by_run[b["run"]] / 60.0 + (np.flatnonzero(on > 0) + 0.5) * b["spf"] / 60.0
            w = on[on > 0]
            col = cond_colour(conditions, b["condition"])
            if t.size:
                h, _ = np.histogram(t, bins=edges, weights=w)
                ymax = max(ymax, int(h.max()))
                ax.bar(edges[:-1], h, width=bin_min, align="edge", color=col, lw=0)
            s0 = starts_by_run[b["run"]] / 60.0
            ax.axvspan(s0, s0 + b["T"] * b["spf"] / 60.0, color=col, alpha=0.06, lw=0)
            ex = np.asarray(b["excluded"], bool)
            if ex.any():
                for t_ex in np.flatnonzero(ex):
                    ax.axvspan(s0 + t_ex * b["spf"] / 60.0, s0 + (t_ex + 1) * b["spf"] / 60.0, color="#ff8c00", alpha=0.25, lw=0)
            ax.text(s0 + 0.2, 0.97, b["run"].split("_")[-3] if b["run"].count("_") >= 3 else b["run"], transform=ax.get_xaxis_transform(),
                    fontsize=6.5, va="top", color="#555")
        ax.set_ylabel(f"{r['label']}\nevents / {bin_min:g} min", fontsize=8)
        _style(ax)
    for ax in axes[:n_hist]:
        ax.set_ylim(0, ymax * 1.15)
    # speed
    ax = axes[n_hist]
    sp_max = 1.0
    for b in all_b:
        s0 = starts_by_run[b["run"]] / 60.0
        tt = s0 + (np.arange(b["T"]) + 0.5) * b["spf"] / 60.0
        sp = np.asarray(b["speed"], float)
        if np.isfinite(sp).any():
            sp_max = max(sp_max, float(np.nanpercentile(sp, 99.5)))
            ax.plot(tt, sp, color="#333", lw=0.5)
            run = np.asarray(b["running"], bool)
            ax.fill_between(tt, 0, np.where(run, sp_max * 2, 0), color="#f2c94c", alpha=0.25, lw=0, step="mid")
            bts = _bouts(b)
            for bt in bts:
                ax.plot([s0 + bt["start"] * b["spf"] / 60.0, s0 + bt["end"] * b["spf"] / 60.0], [sp_max * 1.05] * 2, color="#b8860b", lw=2, solid_capstyle="butt")
        if trials_by_run and b["run"] in trials_by_run:
            for onset, kind in trials_by_run[b["run"]]:
                ax.plot(s0 + onset * b["spf"] / 60.0, sp_max * 1.15, marker="v", ms=5 if kind == "strict" else 3.5,
                        color="#111" if kind == "strict" else "#888", lw=0)
    ax.set_ylim(0, sp_max * 1.25)
    ax.set_ylabel("speed\n(encoder /s)", fontsize=8)
    ax.set_title("shading = running; bars = bouts; ▼ = locomotion trial (black strict, grey relaxed); orange = excluded frames", fontsize=6.5, loc="left", color="#555", pad=2)
    _style(ax)
    if has_pupil:
        ax = axes[n_hist + 1]
        for b in all_b:
            s0 = starts_by_run[b["run"]] / 60.0
            tt = s0 + (np.arange(b["T"]) + 0.5) * b["spf"] / 60.0
            ax.plot(tt, b["pupil"], color="#5a3e85", lw=0.6)
        ax.set_ylabel("pupil\n(px)", fontsize=8); _style(ax)
    if stim_span is not None:
        for ax in axes:
            ax.axvspan(stim_span[0] / 60.0, stim_span[1] / 60.0, color="#bbb", alpha=0.35, hatch="//", lw=0)
        axes[0].text((stim_span[0] + stim_span[1]) / 120.0, 0.5, "stimulation\n(between the recordings)", transform=axes[0].get_xaxis_transform(),
                     ha="center", va="center", fontsize=7, color="#444")
    if block_labels:
        for x, lab in block_labels.items():
            axes[0].text(x / 60.0, 1.02, lab, transform=axes[0].get_xaxis_transform(), fontsize=8, va="bottom", color="#333")
    axes[-1].set_xlabel("time in the session (min)" + ("  — clock ESTIMATED from the run order" if estimated else ""), fontsize=8)
    axes[-1].set_xlim(0, t_end + 0.5)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)


def session_figures(exp_bundles: list[dict], mouse_dir: Path, out_dir: Path, day_labels: dict, trials_by_run: dict | None,
                    stim_gap_min: float = 73.0, log=print) -> list[dict]:
    """One session figure per optostim day (both FOVs) and per spontaneous
    experiment. exp_bundles: [{experiment, fov, para, conditions, base, bundles}].
    Returns rows describing each figure (for the report)."""
    from .timecourse import day_of
    out = []
    by_day = {}
    for e in exp_bundles:
        if e["para"] == "optostim" and e["bundles"]:
            by_day.setdefault(day_of(e["bundles"][0]["run"]) or e["experiment"], []).append(e)
    mouse = Path(mouse_dir).name
    for day, es in sorted(by_day.items()):
        bl = [b for e in es for b in e["bundles"]]
        base = es[0]["base"]; conditions = es[0]["conditions"]
        pl = place_runs(bl, mouse_dir, stim_gap_min=stim_gap_min)
        starts = {b["run"]: s for b, s in zip(bl, pl["starts"])}
        pre_end = max((starts[b["run"]] + b["T"] * b["spf"] for b in bl if b["condition"] == base), default=None)
        post_start = min((starts[b["run"]] for b in bl if b["condition"] != base), default=None)
        span = (pre_end, post_start) if pre_end is not None and post_start is not None and post_start > pre_end else None
        rows = [{"label": f"FOV{e['fov']}" if e.get("fov") is not None else e["experiment"], "bundles": e["bundles"]} for e in es]
        lab = day_labels.get(day, "")
        png = out_dir / f"{mouse}_session_{day}.png"
        try:
            session_figure(f"{mouse} day {day}{' (' + lab + ')' if lab else ''}: events per minute on the day's clock, "
                           f"{base} grey, {', '.join(c for c in conditions if c != base)} coloured", rows, starts, png, conditions, base,
                           trials_by_run=trials_by_run, stim_span=span, day_labels=day_labels, estimated=pl["estimated"])
            out.append({"key": day, "kind": "day", "png": png.name, "experiments": [e["experiment"] for e in es],
                        "estimated_clock": pl["estimated"], "gap_min": (span[1] - span[0]) / 60.0 if span else np.nan})
        except Exception as ex:
            log(f"  session figure of day {day} failed ({ex})")
    for e in exp_bundles:
        if e["para"] == "optostim" or not e["bundles"]:
            continue
        bl = e["bundles"]
        for b in bl:
            b["block"] = b["condition"]
        pl = place_runs(bl, mouse_dir, blocks=True)
        starts = {b["run"]: s for b, s in zip(bl, pl["starts"])}
        labels = {}
        for c in e["conditions"]:
            first = min((starts[b["run"]] for b in bl if b["condition"] == c), default=None)
            if first is not None:
                lab = cond_label(c, day_labels)
                labels[first] = f"{c}" + (f" = {lab}" if lab else "")
        png = out_dir / f"{e['experiment']}_session.png"
        try:
            session_figure(f"{e['experiment']}: events per minute, one block per day ({e['base']} = baseline day, grey)",
                           [{"label": f"FOV{e['fov']}" if e.get("fov") is not None else e["experiment"], "bundles": bl}], starts, png,
                           e["conditions"], e["base"], trials_by_run=trials_by_run, day_labels=day_labels, estimated=True, block_labels=labels)
            out.append({"key": e["experiment"], "kind": "experiment", "png": png.name, "experiments": [e["experiment"]],
                        "estimated_clock": True, "gap_min": np.nan})
        except Exception as ex:
            log(f"  session figure of {e['experiment']} failed ({ex})")
    return out


# ----------------------------------------------------------------------------- G2 / G3 rates by state, ECDFs

def state_rows(exp: str, bundles: list[dict], still_thresh: float = 0.5) -> list[dict]:
    """Per run: minutes, events, rate and median event max dF/F in running and resting frames, % running."""
    rows = []
    for b in bundles:
        ok = ~np.asarray(b["excluded"], bool) & np.isfinite(b["still"])
        run = np.asarray(b["running"], bool)
        ev = b["events"]
        row = {"experiment": exp, "run": b["run"], "condition": b["condition"], "has_loco": b["has_loco"]}
        for state, m in (("running", run & ok), ("resting", ~run & ok)):
            minutes = float(m.sum()) * b["spf"] / 60.0
            n = int(np.asarray(b["onsets"], float)[m].sum())
            row[f"minutes_{state}"] = minutes; row[f"n_events_{state}"] = n
            row[f"events_per_min_{state}"] = n / minutes if minutes >= 0.25 else np.nan
            if ev is not None and len(ev) and "t0_frame" in ev.columns and "dff_max" in ev.columns:
                t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, b["T"] - 1)
                v = ev["dff_max"].to_numpy(float)[m[t0]]; v = v[np.isfinite(v)]
                row[f"dff_max_median_{state}"] = float(np.median(v)) if v.size else np.nan
                row[f"dff_max_mean_{state}"] = float(np.mean(v)) if v.size else np.nan
            else:
                row[f"dff_max_median_{state}"] = np.nan; row[f"dff_max_mean_{state}"] = np.nan
        row["pct_running"] = 100.0 * float(run[ok].mean()) if ok.any() else np.nan
        rows.append(row)
    return rows


def _overlay(ax, img, boxes, title: str, colours=None):
    from matplotlib.patches import Rectangle
    if img is not None:
        lo, hi = np.percentile(img, [1, 99.7])
        ax.imshow(np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1), cmap="gray", interpolation="nearest")
    rng = np.random.default_rng(3)
    for i, (y0, y1, x0, x1) in enumerate(boxes):
        col = colours[i] if colours is not None else rng.uniform(0.25, 1.0, 3)
        ax.add_patch(Rectangle((x0 - 1, y0 - 1), x1 - x0 + 1, y1 - y0 + 1, fill=True, facecolor=(*col, 0.35), edgecolor=(*col, 0.9), lw=0.8))
    ax.set_title(title, fontsize=8); ax.set_xticks([]); ax.set_yticks([])


def state_figure(exp: str, rows: list[dict], bundles: list[dict], conditions: list, out_png: Path, log=print):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .spatial import read_events, read_meanimg
    import pandas as pd
    df = pd.DataFrame(rows)
    n_c = len(conditions)
    fig = plt.figure(figsize=(15, 4.2 + 3.6 * n_c))
    gs = fig.add_gridspec(1 + n_c, 4, height_ratios=[1.0] + [1.0] * n_c)
    rng = np.random.default_rng(2)
    # bars: rate by state, dF/F by state, % running
    for j, (key, ylab, ttl) in enumerate((("events_per_min", "events / min", "event rate by state (dot = run)"),
                                           ("dff_max_median", "median event max dF/F", "event amplitude by state (dot = run)"))):
        ax = fig.add_subplot(gs[0, j])
        for i, c in enumerate(conditions):
            sub = df[df["condition"] == c]
            for k, (state, col) in enumerate((("running", RUN_COL), ("resting", REST_COL))):
                v = sub[f"{key}_{state}"].to_numpy(float); v = v[np.isfinite(v)]
                x = i * 3 + k
                if v.size:
                    ax.bar(x, v.mean(), color=col, alpha=0.35, width=0.8)
                    if v.size > 1:
                        ax.errorbar(x, v.mean(), yerr=_sem(v), color=col, capsize=3, lw=1)
                    ax.scatter(x + rng.uniform(-0.2, 0.2, v.size), v, s=14, color=col, zorder=3, edgecolor="white", lw=0.4)
                else:
                    ax.text(x, 0, "n/a", ha="center", va="bottom", fontsize=7, color="#888")
        ax.set_xticks([i * 3 + 0.5 for i in range(n_c)]); ax.set_xticklabels(conditions, fontsize=8)
        ax.set_ylabel(ylab, fontsize=8); ax.set_title(ttl, fontsize=9); _style(ax)
        if j == 0:
            ax.scatter([], [], color=RUN_COL, label="running"); ax.scatter([], [], color=REST_COL, label="resting"); ax.legend(fontsize=7, frameon=False)
    ax = fig.add_subplot(gs[0, 2])
    for i, c in enumerate(conditions):
        v = df.loc[df["condition"] == c, "pct_running"].to_numpy(float); v = v[np.isfinite(v)]
        col = cond_colour(conditions, c)
        if v.size:
            ax.bar(i, v.mean(), color=col, alpha=0.4, width=0.7)
            ax.scatter(i + rng.uniform(-0.15, 0.15, v.size), v, s=14, color=col, zorder=3, edgecolor="white", lw=0.4)
    ax.set_xticks(range(n_c)); ax.set_xticklabels(conditions, fontsize=8); ax.set_ylabel("% of frames running", fontsize=8)
    ax.set_title("running (dot = run)", fontsize=9); _style(ax)
    ax = fig.add_subplot(gs[0, 3])
    for i, c in enumerate(conditions):
        sub = df[df["condition"] == c]
        for k, (state, col) in enumerate((("running", RUN_COL), ("resting", REST_COL))):
            v = sub[f"minutes_{state}"].to_numpy(float)
            ax.bar(i * 3 + k, v.sum(), color=col, alpha=0.5, width=0.8)
            ax.text(i * 3 + k, v.sum(), f"{int(sub[f'n_events_{state}'].sum())} ev", ha="center", va="bottom", fontsize=7)
    ax.set_xticks([i * 3 + 0.5 for i in range(n_c)]); ax.set_xticklabels(conditions, fontsize=8)
    ax.set_ylabel("minutes", fontsize=8); ax.set_title("time and events in each state (pooled runs)", fontsize=9); _style(ax)
    # overlays per condition x state, on the mean image of the condition's first run
    for i, c in enumerate(conditions):
        bl = [b for b in bundles if b["condition"] == c]
        img = None
        boxes = {"running": [], "resting": []}
        amps = {"running": [], "resting": []}
        for b in bl:
            if b.get("tif") is None:
                continue
            loc = read_events(Path(b["tif"]))
            if img is None:
                img = read_meanimg(Path(b["tif"]))
            ev = b["events"]
            if loc is None or ev is None or not len(ev) or "event" not in ev.columns:
                continue
            loc = loc[loc["event"].isin(ev["event"])]
            t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, b["T"] - 1)
            state_of = {int(e): ("running" if r else "resting") for e, r in zip(ev["event"], np.asarray(b["running"], bool)[t0])}
            for _, r in loc.iterrows():
                st = state_of.get(int(r["event"]))
                if st and all(np.isfinite([r.get("y0", np.nan), r.get("y1", np.nan), r.get("x0", np.nan), r.get("x1", np.nan)])):
                    boxes[st].append((r["y0"], r["y1"], r["x0"], r["x1"])); amps[st].append(float(r.get("dff_max", np.nan)))
        for k, state in enumerate(("running", "resting")):
            ax = fig.add_subplot(gs[1 + i, 2 * k:2 * k + 2])
            if img is None and not boxes[state]:
                ax.axis("off"); ax.set_title(f"{c} — {state}: no event export", fontsize=8); continue
            _overlay(ax, img, boxes[state], f"{c} — events starting while {state} ({len(boxes[state])}; boxes = AQuA2 event extent, all runs of the condition on the first run's mean image)")
    fig.suptitle(f"{exp}: activity by behavioural state (running = still fraction < 0.5 in the frame)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=100)
    plt.close(fig)


def ecdf_figure(exp: str, bundles: list[dict], conditions: list, out_png: Path, norm_pct: float = 99.0) -> dict:
    """ECDFs of event max dF/F (normalised by the experiment's 99th percentile),
    duration and area, by condition (colour) and state at onset (solid running,
    dashed resting). Returns the normalisation value."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    parts = []
    for b in bundles:
        ev = b["events"]
        if ev is None or not len(ev) or "t0_frame" not in ev.columns:
            continue
        t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, b["T"] - 1)
        d = ev.copy()
        d["state"] = np.where(np.asarray(b["running"], bool)[t0], "running", "resting") if b["has_loco"] else "unknown"
        d["condition"] = b["condition"]
        parts.append(d)
    if not parts:
        return {}
    ev = pd.concat(parts, ignore_index=True)
    p99 = float(np.nanpercentile(ev["dff_max"], norm_pct)) if "dff_max" in ev.columns and ev["dff_max"].notna().any() else np.nan
    ev["dff_norm"] = ev["dff_max"] / p99 if np.isfinite(p99) and p99 > 0 else np.nan
    panels = [("dff_norm", f"event max dF/F / {norm_pct:.0f}th percentile of the FOV", False), ("duration_s", "duration (s)", False), ("area_um2", "area (um2)", True)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    for ax, (col, lab, logx) in zip(axes, panels):
        if col not in ev.columns:
            ax.axis("off"); continue
        for c in conditions:
            for state, ls in (("running", "-"), ("resting", "--")):
                v = ev.loc[(ev["condition"] == c) & (ev["state"] == state), col].to_numpy(float)
                v = np.sort(v[np.isfinite(v)])
                if v.size:
                    ax.plot(v, np.arange(1, v.size + 1) / v.size, ls=ls, color=cond_colour(conditions, c), label=f"{c} {state} (n={v.size})", lw=1.3)
        ax.set_xlabel(lab, fontsize=8); ax.set_ylabel("cumulative fraction", fontsize=8); ax.legend(fontsize=6.5, frameon=False)
        if logx:
            ax.set_xscale("log")
        _style(ax)
    fig.suptitle(f"{exp}: event properties by state at onset and condition (events pooled; the unit for statistics is the run / animal)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)
    return {"dff_norm_p99": p99, "n_events": int(len(ev))}


# ----------------------------------------------------------------------------- G6 pupil measures

def pupil_rows(exp: str, bundles: list[dict], para: str, trial_rows: list[dict] | None = None, min_rest_frac: float = 0.10,
               spont_rise: float = 0.10) -> list[dict]:
    """Per run: pupil normalised 0-1 per FOV and day (common 1st/99th percentile
    over that day's runs of the FOV), pupil variation (SD of the 0-1 trace),
    resting pupil (mean of 0-1 over resting frames; valid when >= 10 % resting),
    locomotion-related dilation (per trial, from locotrials), spontaneous
    dilations at rest (rises >= spont_rise of the range, >= 20 s from running)."""
    from .timecourse import day_of
    rows = []
    groups = {}
    for b in bundles:
        key = (exp, b["condition"]) if para != "optostim" else (exp, day_of(b["run"]) or "")
        groups.setdefault(key, []).append(b)
    lims = {}
    for key, bl in groups.items():
        v = np.concatenate([np.asarray(b["pupil"], float)[np.isfinite(b["pupil"])] for b in bl]) if bl else np.zeros(0)
        lims[key] = (np.nanpercentile(v, 1), np.nanpercentile(v, 99)) if v.size > 20 else (np.nan, np.nan)
    dil_by_run = {}
    for r in trial_rows or []:
        if np.isfinite(r.get("pupil_dilation", np.nan)):
            dil_by_run.setdefault(r["run"], []).append(r["pupil_dilation"])
    for b in bundles:
        key = (exp, b["condition"]) if para != "optostim" else (exp, day_of(b["run"]) or "")
        lo, hi = lims[key]
        row = {"experiment": exp, "run": b["run"], "condition": b["condition"], "has_pupil": b["has_pupil"],
               "norm_group": "+".join(str(k) for k in key), "pupil_norm_lo_px": lo, "pupil_norm_hi_px": hi}
        if not b["has_pupil"] or not np.isfinite(lo) or hi <= lo:
            rows.append(row); continue
        p01 = np.clip((np.asarray(b["pupil"], float) - lo) / (hi - lo), 0, 1)
        ok = ~np.asarray(b["excluded"], bool) & np.isfinite(p01)
        run = np.asarray(b["running"], bool)
        row["pupil_variation"] = float(np.nanstd(p01[ok])) if ok.sum() > 10 else np.nan
        row["pupil_mean_norm"] = float(np.nanmean(p01[ok])) if ok.any() else np.nan
        rest = ~run & ok
        row["resting_fraction"] = float(rest.sum() / max(ok.sum(), 1))
        row["resting_pupil_valid"] = int(row["resting_fraction"] >= min_rest_frac)
        row["resting_pupil"] = float(np.nanmean(p01[rest])) if row["resting_pupil_valid"] else np.nan
        row["resting_pupil_px"] = float(np.nanmean(np.asarray(b["pupil"], float)[rest])) if row["resting_pupil_valid"] else np.nan
        row["running_pupil"] = float(np.nanmean(p01[run & ok])) if (run & ok).sum() > 5 else np.nan
        d = dil_by_run.get(b["run"], [])
        row["loco_dilation_px_mean"] = float(np.mean(d)) if d else np.nan
        row["loco_dilation_norm_mean"] = float(np.mean(d) / (hi - lo)) if d else np.nan
        row["n_loco_trials_pupil"] = len(d)
        # spontaneous dilations at rest
        T = b["T"]
        near = run.copy()
        idx = np.flatnonzero(run)
        for m in range(1, REST_MARGIN_FRAMES + 1):
            near[np.clip(idx + m, 0, T - 1)] = True; near[np.clip(idx - m, 0, T - 1)] = True
        quiet = ~near & ok
        n_rise, amps = 0, []
        try:
            from scipy.signal import find_peaks
            x = np.where(quiet, p01, np.nan)
            seg_starts = np.flatnonzero(np.diff(np.r_[0, quiet.astype(int), 0]) == 1)
            seg_ends = np.flatnonzero(np.diff(np.r_[0, quiet.astype(int), 0]) == -1)
            for s, e in zip(seg_starts, seg_ends):
                if e - s < 6:
                    continue
                seg = x[s:e]
                pk, props = find_peaks(seg, prominence=spont_rise)
                n_rise += pk.size; amps += list(props["prominences"])
        except Exception:
            pass
        quiet_min = float(quiet.sum()) * b["spf"] / 60.0
        row["rest_minutes_pupil"] = quiet_min
        row["spont_dilations_per_min"] = n_rise / quiet_min if quiet_min >= 0.5 else np.nan
        row["spont_dilation_amp_norm"] = float(np.mean(amps)) if amps else np.nan
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------- G8 acute & chronic

MEASURES = [
    # (key, label, source) — source: cond = condition_summary, loco = loco_conditions, pupil = pupil_measures (mean over runs), state = state rows
    ("events_per_min", "event rate (all)", "cond"),
    ("events_per_min_running", "event rate, running", "cond"),
    ("events_per_min_still", "event rate, resting", "cond"),
    ("rest_events_per_min", "event rate, rest-only (≥ 20 s from running)", "loco"),
    ("dff_max_median", "event max dF/F (median)", "cond"),
    ("pct_responsive_z_relaxed", "responsive loco trials (%, z rule)", "loco"),
    ("peak_responsive_mean_relaxed", "responsive-trial peak dF/F", "loco"),
    ("events_per_trial_relaxed", "events per loco trial", "loco"),
    ("pct_cfu_loco_responsive", "loco-responsive CFUs (%)", "loco"),
    ("pct_running_mean", "% running", "loco"),
    ("n_bouts_mean", "running bouts per run", "loco"),
    ("pupil_variation", "pupil variation (SD of 0-1)", "pupil"),
    ("resting_pupil", "resting pupil (0-1)", "pupil"),
    ("loco_dilation_norm_mean", "loco-related pupil dilation (0-1)", "pupil"),
    ("spont_dilations_per_min", "spontaneous pupil dilations / min at rest", "pupil"),
    ("dev_l2f_population_mean", "GLM speed→dF/F deviance explained", "loco"),
    ("dev_lp2f_population_mean", "GLM speed+pupil→dF/F deviance explained", "loco"),
]


def acute_chronic(mouse: str, exp_meta: list[dict], cond_rows: list[dict], loco_cond: list[dict], pupil_rows_: list[dict],
                  day_labels: dict, out_dir: Path, log=print) -> list[dict]:
    """exp_meta: [{experiment, fov, para, conditions, base, day}] . Writes
    acute_chronic.csv (one row per experiment x condition x measure with the
    value, the baseline value, acute = value/baseline (optostim) and chronic =
    baseline(day)/baseline(first day of the FOV)) and <mouse>_acute_chronic.png."""
    import pandas as pd
    cond = pd.DataFrame(cond_rows) if cond_rows else pd.DataFrame()
    loco = pd.DataFrame(loco_cond) if loco_cond else pd.DataFrame()
    pup = pd.DataFrame(pupil_rows_) if pupil_rows_ else pd.DataFrame()

    def value(src, exp, c, key):
        if src == "cond" and len(cond) and key in cond.columns:
            v = cond.loc[(cond["experiment"] == exp) & (cond["condition"] == c), key]
        elif src == "loco" and len(loco) and key in loco.columns:
            v = loco.loc[(loco["experiment"] == exp) & (loco["condition"] == c), key]
        elif src == "pupil" and len(pup) and key in pup.columns:
            v = pup.loc[(pup["experiment"] == exp) & (pup["condition"] == c), key]
            v = v[np.isfinite(v.to_numpy(float))]
            return float(v.mean()) if len(v) else np.nan
        else:
            return np.nan
        v = v.to_numpy(float)
        return float(v[0]) if v.size and np.isfinite(v[0]) else np.nan
    rows = []
    for e in exp_meta:
        for c in e["conditions"]:
            for key, label, src in MEASURES:
                v = value(src, e["experiment"], c, key)
                vb = value(src, e["experiment"], e["base"], key)
                rows.append({"experiment": e["experiment"], "fov": e.get("fov"), "paradigm": e["para"], "condition": c,
                             "day": (e.get("day") if e["para"] == "optostim" else (c[1:] if c.startswith("d") else c)),
                             "measure": key, "label": label, "value": v, "baseline_value": vb,
                             "acute_ratio": (v / vb if (e["para"] == "optostim" and c != e["base"] and np.isfinite(v) and np.isfinite(vb) and vb > 0) else np.nan)})
    df = pd.DataFrame(rows)
    # chronic: the baseline value of each FOV on each day vs its first day
    df["chronic_ratio"] = np.nan
    base_of = {e["experiment"]: e["base"] for e in exp_meta}
    if len(df):
        for (fov, key), sub in df.groupby(["fov", "measure"], dropna=False):
            opt = sub[(sub["paradigm"] == "optostim") & (sub["condition"] == sub["experiment"].map(base_of))].sort_values("day")
            if len(opt):
                first = float(opt.iloc[0]["value"])
                if np.isfinite(first) and first > 0:
                    df.loc[opt.index, "chronic_ratio"] = opt["value"] / first
            for exp, se in sub[sub["paradigm"] != "optostim"].groupby("experiment"):
                bv = se.loc[se["condition"] == base_of.get(exp), "value"]
                first = float(bv.iloc[0]) if len(bv) else np.nan
                if np.isfinite(first) and first > 0:
                    df.loc[se.index, "chronic_ratio"] = se["value"] / first
    df["day_label"] = df["day"].map(lambda d: (day_labels.get(str(d), "") or cond_label("d" + str(d), day_labels)) if d is not None else "")
    df.to_csv(out_dir / "acute_chronic.csv", index=False)
    try:
        acute_chronic_figure(mouse, df, out_dir / f"{mouse}_acute_chronic.png")
    except Exception as e:
        log(f"  acute/chronic figure failed ({e})")
    return rows


def acute_chronic_figure(mouse: str, df, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = [k for k, _, _ in MEASURES if np.isfinite(df.loc[df["measure"] == k, "acute_ratio"]).any() or np.isfinite(df.loc[df["measure"] == k, "chronic_ratio"]).any()]
    if not keys:
        return
    n = len(keys); ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.9 * nrow + 0.8), squeeze=False)
    fovs = sorted(df["fov"].dropna().unique().tolist(), key=lambda v: str(v))
    fcol = {f: ["#2a78d6", "#d7263d", "#1b998b", "#7b2cbf"][i % 4] for i, f in enumerate(fovs)}
    days = sorted(df["day"].dropna().astype(str).unique().tolist())
    xpos = {d: i for i, d in enumerate(days)}
    labels = {d: (d + ("\n" + df.loc[df["day"].astype(str) == d, "day_label"].iloc[0] if df.loc[df["day"].astype(str) == d, "day_label"].iloc[0] else "")) for d in days}
    for ax, key in zip(axes.ravel(), keys):
        sub = df[df["measure"] == key]
        lab = sub["label"].iloc[0]
        for f in fovs:
            s = sub[sub["fov"] == f]
            for ratio, ls, mk, what in (("acute_ratio", "-", "o", "acute post/pre"), ("chronic_ratio", "--", "s", "chronic vs day 1")):
                pts = s[np.isfinite(s[ratio].to_numpy(float))]
                if not len(pts):
                    continue
                pts = pts.sort_values("day")
                x = [xpos[str(d)] for d in pts["day"]]; y = np.clip(pts[ratio].to_numpy(float), 0.05, 20)
                ax.plot(x, y, ls=ls, marker=mk, ms=4, color=fcol[f], lw=1.2, label=f"FOV{f} {what}")
        ax.axhline(1, color="#999", lw=0.8, ls=":")
        ax.set_yscale("log"); ax.set_yticks([0.1, 0.25, 0.5, 1, 2, 4, 10]); ax.set_yticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "10"], fontsize=7)
        ax.set_ylim(0.04, 25)
        ax.set_xticks(range(len(days))); ax.set_xticklabels([labels[d] for d in days], fontsize=7)
        ax.set_title(lab, fontsize=8); _style(ax)
    for ax in axes.ravel()[len(keys):]:
        ax.axis("off")
    h, l = axes.ravel()[0].get_legend_handles_labels()
    if h:
        fig.legend(h, l, loc="lower right", fontsize=7, frameon=False, ncol=2)
    fig.suptitle(f"{mouse}: acute (post / pre, solid) and chronic (pre of day # / pre of day 1, dashed) change of every measure, per FOV "
                 "(1 = no change; values clipped to 0.05-20)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)


def utc(ts: float) -> str:
    return _dt.datetime.fromtimestamp(ts, _dt.timezone.utc).strftime("%Y-%m-%d %H:%M")
