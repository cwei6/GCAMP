#!/usr/bin/env python3
"""aqua_conditions — the stimulation protocols side by side, across animals.

Quick comparison of the pre -> post change of afferent activity by protocol
(2 Hz, 5 Hz, 5 Hz + 40 ms pulses, sham, ...) over every RearChr animal that has
been through tab 7 (aqua_analyze.py). Nothing is recomputed from the movies:
the tool reads each animal's aqua_analysis\\experiments_summary.csv,
condition_summary.csv, heatmap_summary.csv and day_labels.csv, labels every
(day, FOV) experiment with the protocol of its day, and puts the experiments
of one protocol together.

Protocol = read from the day label in <mouse>\\day_labels.csv:
  'sham'                        -> sham
  '... 40 ms ...'               -> <n> Hz + 40 ms   (the pulse-width protocol)
  '<n> Hz ...'                  -> <n> Hz           (power kept separately: '3.3 mW' / 'max power')
  'test run (no stimulation)'   -> none (listed, not compared)
Anything else is its own group under the label as written.

Per experiment (one day x FOV of one animal; pre = the 2 runs before, post =
the 2 runs after the stimulation):
  * event rate pre and post (events/min, all frames) and its change: ratio
    with the exact interval, log2 ratio with pseudo-count, modulation index,
    low-count flag — the numbers of the rate section of the report;
  * the same change on still frames only and on running frames only (the
    locomotion-matched numbers), and the speed-matched ratio when the two
    conditions overlap in behaviour;
  * event amplitude (median peak dF/F) and duration pre and post;
  * CFU fate: how many CFUs went up / down / unchanged (heatmap step);
  * the behavioural state (time still) pre and post — the confound to read
    the raw ratio against.
Per protocol: the experiments' changes averaged (mean +- SEM of the log2
ratios, n experiments, n animals), the pooled rate ratio (all events and
minutes of the group together, exact interval) and, when there are enough
experiments, a Mann-Whitney test of each protocol's log2 ratios against the
sham days — descriptive with these n.

No pupil step is used here (the DLC data may still be in the making); the
locomotion matching comes from the analysis as it was run.

Outputs, in <parent>\\cross_animal\\ (parent = the folder that holds the mouse
folders): conditions_experiments.csv (one row per experiment),
conditions_groups.csv (one row per protocol and measure), conditions.png,
conditions.html (figure + tables + what is missing), conditions.txt.

Usage:
  python aqua_conditions.py "D:\\2photon\\Chao\\GCAMP_affrerents"
  python aqua_conditions.py <parent> --animals NavRearChrTIGRE-2608-1,NavRearChrTIGRE-2608-2
  python aqua_conditions.py <parent> --out <dir>
"""
from __future__ import annotations

import argparse
import html
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from afferent2p import stats as S                          # noqa: E402
from afferent2p.timecourse import read_day_labels          # noqa: E402

GROUP_ORDER = ["sham", "2 Hz", "5 Hz", "5 Hz + 40 ms"]
GROUP_COLOR = {"sham": "#8a8a8a", "2 Hz": "#1f77b4", "5 Hz": "#d62728", "5 Hz + 40 ms": "#9467bd", "none": "#cccccc"}
ANIMAL_MARKERS = ["o", "s", "D", "^", "v", "P", "X", "*"]
MEASURES = [
    ("log2_rate_all", "event rate, all frames (log2 post/pre)"),
    ("log2_rate_still", "event rate, still frames (log2 post/pre)"),
    ("log2_rate_running", "event rate, running frames (log2 post/pre)"),
    ("log2_amp", "event amplitude, median peak dF/F (log2 post/pre)"),
    ("log2_duration", "event duration, median (log2 post/pre)"),
    ("cfu_up_share", "share of CFUs up (of those tested)"),
    ("cfu_down_share", "share of CFUs down (of those tested)"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------------- labels -> protocol

def protocol_of(label: str) -> tuple[str, str, bool]:
    """(group, power, new_fov) from a day label."""
    s = (label or "").strip()
    low = s.lower()
    new_fov = "new fov" in low
    power = "3.3 mW" if "3.3" in low else ("max power" if "max" in low else "")
    if not s:
        return "unlabelled", power, new_fov
    if "sham" in low:
        return "sham", power, new_fov
    if "no stim" in low or "test run" in low:
        return "none", power, new_fov
    m = re.search(r"(\d+(?:\.\d+)?)\s*hz", low)
    if m:
        hz = m.group(1).rstrip("0").rstrip(".") if "." in m.group(1) else m.group(1)
        if "40 ms" in low or "40ms" in low:
            return f"{hz} Hz + 40 ms", power, new_fov
        return f"{hz} Hz", power, new_fov
    return s, power, new_fov


def _f(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def _log2_ratio(nb, tb, nc, tc, pseudo=0.5):
    nb, tb, nc, tc = (_f(v) for v in (nb, tb, nc, tc))
    if not all(np.isfinite(v) for v in (nb, tb, nc, tc)) or tb <= 0 or tc <= 0:
        return np.nan
    return float(np.log2(((nc + pseudo) / tc) / ((nb + pseudo) / tb)))


# ----------------------------------------------------------------------------- gather

def gather(parent: Path, animals: list[str] | None = None) -> tuple[pd.DataFrame, list[str]]:
    """One row per (animal, day, FOV) experiment with the measures; and notes on what is missing."""
    rows, notes = [], []
    dirs = sorted(d for d in parent.iterdir() if d.is_dir() and "rearchr" in d.name.lower())
    if animals:
        want = {a.strip() for a in animals}
        dirs = [d for d in dirs if d.name in want]
        missing = want - {d.name for d in dirs}
        for m in sorted(missing):
            notes.append(f"{m}: folder not found under {parent}")
    if not dirs:
        notes.append(f"no RearChr animal folder under {parent}")
    for d in dirs:
        an = d / "aqua_analysis"
        es = an / "experiments_summary.csv"
        if not es.exists():
            notes.append(f"{d.name}: no aqua_analysis\\experiments_summary.csv — run tab 7 (analyse) for this animal")
            continue
        E = pd.read_csv(es)
        C = pd.read_csv(an / "condition_summary.csv") if (an / "condition_summary.csv").exists() else None
        Hm = pd.read_csv(an / "heatmap_summary.csv") if (an / "heatmap_summary.csv").exists() else None
        labels = read_day_labels(d)
        days_done = set()
        for _, r in E.iterrows():
            exp = str(r["experiment"])
            m = re.search(r"_(\d{6})_FOV(\d+)$", exp)
            if not m:
                continue
            day, fov = m.group(1), int(m.group(2))
            days_done.add(day)
            label = labels.get(day, "")
            group, power, new_fov = protocol_of(label)
            if not label:
                notes.append(f"{d.name} {day}: no label in day_labels.csv — listed as 'unlabelled'")
            row = {"animal": d.name, "day": day, "fov": fov, "experiment": exp, "label": label, "group": group,
                   "power": power, "new_fov": new_fov}
            for k in ("n_runs_baseline", "n_runs_condition", "minutes_baseline", "minutes_condition", "rate_baseline",
                      "rate_condition", "n_events_baseline", "n_events_condition", "ratio_raw", "ratio_ci_lo", "ratio_ci_hi",
                      "log2_ratio_pc", "modulation_index", "low_counts", "still_frac_baseline", "still_frac_condition",
                      "ratio_still", "ratio_running", "ratio_speed_matched", "behaviour_comparable", "ratio_adjusted",
                      "ratio_adjusted_ci_lo", "ratio_adjusted_ci_hi"):
                row[k] = r.get(k, np.nan)
            row["log2_rate_all"] = _f(r.get("log2_ratio_pc"))
            # still / running from the condition table (counts and minutes -> pseudo-count log2 ratios)
            if C is not None:
                cb = C[(C["experiment"] == exp) & (C["condition"] == r.get("baseline", "pre"))]
                cc = C[(C["experiment"] == exp) & (C["condition"] == r.get("condition", "post"))]
                if len(cb) and len(cc):
                    b, c = cb.iloc[0], cc.iloc[0]
                    row["n_events_still_baseline"], row["minutes_still_baseline"] = _f(b.get("n_events_still")), _f(b.get("minutes_still"))
                    row["n_events_still_condition"], row["minutes_still_condition"] = _f(c.get("n_events_still")), _f(c.get("minutes_still"))
                    row["n_events_running_baseline"], row["minutes_running_baseline"] = _f(b.get("n_events_running")), _f(b.get("minutes_running"))
                    row["n_events_running_condition"], row["minutes_running_condition"] = _f(c.get("n_events_running")), _f(c.get("minutes_running"))
                    row["rate_still_baseline"], row["rate_still_condition"] = _f(b.get("events_per_min_still")), _f(c.get("events_per_min_still"))
                    row["rate_running_baseline"], row["rate_running_condition"] = _f(b.get("events_per_min_running")), _f(c.get("events_per_min_running"))
                    row["log2_rate_still"] = _log2_ratio(b.get("n_events_still"), b.get("minutes_still"), c.get("n_events_still"), c.get("minutes_still")) \
                        if min(_f(b.get("minutes_still")), _f(c.get("minutes_still"))) >= 1.0 else np.nan
                    row["log2_rate_running"] = _log2_ratio(b.get("n_events_running"), b.get("minutes_running"), c.get("n_events_running"), c.get("minutes_running")) \
                        if min(_f(b.get("minutes_running")), _f(c.get("minutes_running"))) >= 1.0 else np.nan
                    row["amp_baseline"], row["amp_condition"] = _f(b.get("dff_max_median")), _f(c.get("dff_max_median"))
                    row["duration_baseline"], row["duration_condition"] = _f(b.get("duration_s_median")), _f(c.get("duration_s_median"))
                    row["area_baseline"], row["area_condition"] = _f(b.get("area_um2_median")), _f(c.get("area_um2_median"))
                    ab, ac = row["amp_baseline"], row["amp_condition"]
                    row["log2_amp"] = float(np.log2(ac / ab)) if ab > 0 and ac > 0 else np.nan
                    db, dc = row["duration_baseline"], row["duration_condition"]
                    row["log2_duration"] = float(np.log2(dc / db)) if db > 0 and dc > 0 else np.nan
            if Hm is not None:
                h = Hm[Hm["experiment"] == exp]
                if len(h):
                    h = h.iloc[0]
                    up = next((_f(h[c]) for c in h.index if str(c).startswith("n_up")), np.nan)
                    dn = next((_f(h[c]) for c in h.index if str(c).startswith("n_down")), np.nan)
                    un = next((_f(h[c]) for c in h.index if str(c).startswith("n_unchanged")), np.nan)
                    tested = np.nansum([up, dn, un])
                    row["cfu_up"], row["cfu_down"], row["cfu_unchanged"], row["cfu_tested"] = up, dn, un, tested
                    row["cfu_up_share"] = up / tested if tested > 0 else np.nan
                    row["cfu_down_share"] = dn / tested if tested > 0 else np.nan
            rows.append(row)
        # days present in the folder but not analysed
        for p in sorted(d.glob(f"{d.name}_*_0*")):
            m = re.match(rf"{re.escape(d.name)}_(\d{{6}})_\d{{3}}$", p.name)
            if m and m.group(1) not in days_done:
                tif = list(p.glob("*_bin*s.tif"))
                state = "unpacked, no analysis yet" if tif else "raw only (not unpacked)"
                lab = labels.get(m.group(1), "")
                note = f"{d.name} {m.group(1)} ({lab or 'no label'}): {state}"
                if note not in notes:
                    notes.append(note)
    T = pd.DataFrame(rows)
    return T, notes


# ----------------------------------------------------------------------------- group statistics

def summarise(T: pd.DataFrame) -> pd.DataFrame:
    from scipy import stats as sps
    rows = []
    if T.empty:
        return pd.DataFrame(rows)
    groups = [g for g in GROUP_ORDER if g in set(T["group"])] + sorted(g for g in set(T["group"]) if g not in GROUP_ORDER and g != "none")
    sham = T[T["group"] == "sham"]
    for g in groups:
        G = T[T["group"] == g]
        base = {"group": g, "n_experiments": int(len(G)), "n_animals": int(G["animal"].nunique()),
                "n_days": int(G[["animal", "day"]].drop_duplicates().shape[0]),
                "animals": ", ".join(sorted(G["animal"].unique())), "powers": ", ".join(sorted(p for p in G["power"].unique() if p))}
        for key, label in MEASURES:
            if key not in G.columns:
                continue
            v = G[key].astype(float).values
            v = v[np.isfinite(v)]
            r = dict(base, measure=key, measure_label=label, n=int(len(v)))
            if len(v):
                r["mean"] = float(np.mean(v)); r["sem"] = float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else np.nan
                r["median"] = float(np.median(v)); r["min"] = float(np.min(v)); r["max"] = float(np.max(v))
                if len(v) >= 3 and key.startswith("log2"):
                    try:
                        r["p_wilcoxon_vs_0"] = float(sps.wilcoxon(v).pvalue)
                    except ValueError:
                        r["p_wilcoxon_vs_0"] = np.nan
                if g != "sham" and key in sham.columns:
                    s = sham[key].astype(float).values; s = s[np.isfinite(s)]
                    if len(v) >= 2 and len(s) >= 2:
                        try:
                            r["p_mannwhitney_vs_sham"] = float(sps.mannwhitneyu(v, s, alternative="two-sided").pvalue)
                        except ValueError:
                            r["p_mannwhitney_vs_sham"] = np.nan
            if key == "log2_rate_all":
                nb, tb = G["n_events_baseline"].astype(float).sum(), G["minutes_baseline"].astype(float).sum()
                nc, tc = G["n_events_condition"].astype(float).sum(), G["minutes_condition"].astype(float).sum()
                cm = S.change_measures(nb, tb, nc, tc)
                r.update({"pooled_n_events_pre": nb, "pooled_minutes_pre": tb, "pooled_n_events_post": nc, "pooled_minutes_post": tc,
                          "pooled_rate_pre": nb / tb if tb else np.nan, "pooled_rate_post": nc / tc if tc else np.nan,
                          "pooled_ratio": (nc / tc) / (nb / tb) if tb and tc and nb else np.nan,
                          "pooled_ratio_ci_lo": cm.get("ratio_ci_lo"), "pooled_ratio_ci_hi": cm.get("ratio_ci_hi"),
                          "pooled_log2_ratio_pc": cm.get("log2_ratio_pc")})
            if key == "log2_rate_still" and "n_events_still_baseline" in G.columns:
                nb, tb = G["n_events_still_baseline"].astype(float).sum(), G["minutes_still_baseline"].astype(float).sum()
                nc, tc = G["n_events_still_condition"].astype(float).sum(), G["minutes_still_condition"].astype(float).sum()
                cm = S.change_measures(nb, tb, nc, tc)
                r.update({"pooled_n_events_pre": nb, "pooled_minutes_pre": tb, "pooled_n_events_post": nc, "pooled_minutes_post": tc,
                          "pooled_rate_pre": nb / tb if tb else np.nan, "pooled_rate_post": nc / tc if tc else np.nan,
                          "pooled_ratio": (nc / tc) / (nb / tb) if tb and tc and nb else np.nan,
                          "pooled_ratio_ci_lo": cm.get("ratio_ci_lo"), "pooled_ratio_ci_hi": cm.get("ratio_ci_hi"),
                          "pooled_log2_ratio_pc": cm.get("log2_ratio_pc")})
            rows.append(r)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- figure

def figure(T: pd.DataFrame, G: pd.DataFrame, out_png: Path, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size": 8.5})
    groups = [g for g in GROUP_ORDER if g in set(T["group"])] + sorted(g for g in set(T["group"]) if g not in GROUP_ORDER and g != "none")
    animals = sorted(T["animal"].unique())
    mk = {a: ANIMAL_MARKERS[i % len(ANIMAL_MARKERS)] for i, a in enumerate(animals)}
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle(title, fontsize=12)

    def xpos(gi, n, k):
        return gi + (np.linspace(-0.22, 0.22, n)[k] if n > 1 else 0.0)

    def scatter_measure(ax, key, ylabel, ttl, zero_line=True, ci=None):
        for gi, g in enumerate(groups):
            Gg = T[T["group"] == g].reset_index(drop=True)
            vals = []
            for k, (_, r) in enumerate(Gg.iterrows()):
                v = _f(r.get(key))
                if not np.isfinite(v):
                    continue
                x = xpos(gi, len(Gg), k)
                filled = r.get("power") != "3.3 mW"
                ax.plot(x, v, marker=mk[r["animal"]], ms=7, color=GROUP_COLOR.get(g, "#333"), mfc=GROUP_COLOR.get(g, "#333") if filled else "white",
                        mew=1.3, ls="none", zorder=3)
                if ci is not None:
                    lo, hi = ci(r)
                    if np.isfinite(lo) and np.isfinite(hi):
                        ax.plot([x, x], [max(lo, -6), min(hi, 6)], color=GROUP_COLOR.get(g, "#333"), lw=0.8, alpha=0.6, zorder=2)
                if r.get("low_counts") in (True, "True") and key == "log2_rate_all":
                    ax.plot(x, v, marker="x", ms=4, color="k", ls="none", zorder=4)
                vals.append(v)
            if vals:
                m = float(np.mean(vals)); sem = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
                ax.plot([gi - 0.3, gi + 0.3], [m, m], color=GROUP_COLOR.get(g, "#333"), lw=2.2, zorder=5)
                ax.fill_between([gi - 0.3, gi + 0.3], m - sem, m + sem, color=GROUP_COLOR.get(g, "#333"), alpha=0.15, lw=0, zorder=1)
                ax.text(gi, ax.get_ylim()[1] if False else np.nan, "", fontsize=7)
        if zero_line:
            ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(len(groups))); ax.set_xticklabels([f"{g}\n(n = {int((T['group'] == g).sum())})" for g in groups])
        ax.set_ylabel(ylabel); ax.set_title(ttl, fontsize=9.5)
        ax.grid(axis="y", alpha=0.25)

    # A: pre -> post rates (paired)
    ax = axes[0, 0]
    for gi, g in enumerate(groups):
        Gg = T[T["group"] == g].reset_index(drop=True)
        for k, (_, r) in enumerate(Gg.iterrows()):
            x = xpos(gi, len(Gg), k)
            rb, rc = _f(r.get("rate_baseline")), _f(r.get("rate_condition"))
            if not (np.isfinite(rb) and np.isfinite(rc)):
                continue
            ax.plot([x - 0.08, x + 0.08], [rb, rc], color=GROUP_COLOR.get(g, "#333"), lw=1.2, alpha=0.85)
            ax.plot(x - 0.08, rb, marker=mk[r["animal"]], ms=5, color=GROUP_COLOR.get(g, "#333"), mfc="white", ls="none")
            ax.plot(x + 0.08, rc, marker=mk[r["animal"]], ms=5, color=GROUP_COLOR.get(g, "#333"), ls="none")
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups)
    ax.set_ylabel("events / min (all frames)"); ax.set_title("A  event rate pre (open) -> post (filled), one line per day x FOV", fontsize=9.5)
    ax.grid(axis="y", alpha=0.25)
    # B: log2 ratio all frames with CI
    def ci_all(r):
        lo, hi = _f(r.get("ratio_ci_lo")), _f(r.get("ratio_ci_hi"))
        return (np.log2(lo) if lo > 0 else -np.inf), (np.log2(hi) if np.isfinite(hi) and hi > 0 else np.inf)
    scatter_measure(axes[0, 1], "log2_rate_all", "log2 (post / pre), pseudo-count", "B  rate change, all frames (bars: exact 95 % interval; x = few events)", ci=ci_all)
    scatter_measure(axes[0, 2], "log2_rate_still", "log2 (post / pre), still frames", "C  rate change on still frames only (locomotion-matched)")
    scatter_measure(axes[1, 0], "log2_rate_running", "log2 (post / pre), running frames", "D  rate change on running frames only")
    scatter_measure(axes[1, 1], "log2_amp", "log2 (post / pre) of the median peak dF/F", "E  event amplitude change")
    # F: CFU fate + behaviour
    ax = axes[1, 2]
    for gi, g in enumerate(groups):
        Gg = T[T["group"] == g].reset_index(drop=True)
        for k, (_, r) in enumerate(Gg.iterrows()):
            x = xpos(gi, len(Gg), k)
            up, dn, un = _f(r.get("cfu_up")), _f(r.get("cfu_down")), _f(r.get("cfu_unchanged"))
            tot = np.nansum([up, dn, un])
            if tot > 0:
                ax.bar(x, up / tot, 0.12, color="#d62728", alpha=0.85)
                ax.bar(x, -dn / tot, 0.12, color="#1f77b4", alpha=0.85)
                ax.text(x, (up / tot) + 0.03, f"{int(tot)}", ha="center", fontsize=6.5, color="#333")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylim(-1.05, 1.15)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups)
    ax.set_ylabel("share of CFUs: up (red, +) / down (blue, -)"); ax.set_title("F  CFU fate post vs pre (number = CFUs tested)", fontsize=9.5)
    ax.grid(axis="y", alpha=0.25)
    # legend: animals + power
    handles = [Line2D([], [], marker=mk[a], color="#444", ls="none", ms=7, label=a) for a in animals]
    handles += [Line2D([], [], marker="o", color="#444", mfc="#444", ls="none", ms=7, label="max power (filled)"),
                Line2D([], [], marker="o", color="#444", mfc="white", ls="none", ms=7, label="3.3 mW (open)")]
    fig.legend(handles=handles, loc="lower center", ncol=min(6, len(handles)), fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- text + html

def _fmt(v, nd=2):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return ""
    return "" if not np.isfinite(v) else f"{v:.{nd}f}"


def write_text(T: pd.DataFrame, G: pd.DataFrame, notes: list[str], out_txt: Path, parent: Path) -> str:
    lines = [f"Stimulation protocols compared across animals — {parent}", f"made {time.strftime('%Y-%m-%d %H:%M')} by aqua_conditions.py", ""]
    if T.empty:
        lines.append("no analysed RearChr experiment found")
    else:
        groups = [g for g in GROUP_ORDER if g in set(T["group"])] + sorted(g for g in set(T["group"]) if g not in GROUP_ORDER)
        for g in groups:
            Gg = T[T["group"] == g]
            lines.append(f"{g}: {len(Gg)} experiment(s) (day x FOV) from {Gg['animal'].nunique()} animal(s): "
                         + "; ".join(f"{r['animal'].replace('NavRearChrTIGRE-', '')} {r['day']} FOV{r['fov']} [{r['label']}]" for _, r in Gg.iterrows()))
            for _, r in Gg.iterrows():
                lines.append(f"    {r['animal'].replace('NavRearChrTIGRE-', '')} {r['day']} FOV{r['fov']}: rate {_fmt(r.get('rate_baseline'), 3)} -> {_fmt(r.get('rate_condition'), 3)} /min "
                             f"({int(_f(r.get('n_events_baseline')) or 0)} -> {int(_f(r.get('n_events_condition')) or 0)} events), log2 ratio {_fmt(r.get('log2_rate_all'))}"
                             f"{' (few events)' if r.get('low_counts') in (True, 'True') else ''}; still-only log2 {_fmt(r.get('log2_rate_still'))}, "
                             f"running-only {_fmt(r.get('log2_rate_running'))}; peak dF/F {_fmt(r.get('amp_baseline'))} -> {_fmt(r.get('amp_condition'))}; "
                             f"CFUs up/down/same {_fmt(r.get('cfu_up'), 0)}/{_fmt(r.get('cfu_down'), 0)}/{_fmt(r.get('cfu_unchanged'), 0)}; "
                             f"time still {_fmt(100 * _f(r.get('still_frac_baseline')), 0)} % -> {_fmt(100 * _f(r.get('still_frac_condition')), 0)} %")
            if g == "none":
                continue
            gg = G[(G["group"] == g)]
            for _, s in gg.iterrows():
                if s["measure"] in ("log2_rate_all", "log2_rate_still", "log2_amp"):
                    extra = ""
                    if np.isfinite(_f(s.get("pooled_ratio"))):
                        extra = (f"; pooled {_fmt(s.get('pooled_rate_pre'), 3)} -> {_fmt(s.get('pooled_rate_post'), 3)} /min = x{_fmt(s.get('pooled_ratio'))} "
                                 f"(95 % CI {_fmt(s.get('pooled_ratio_ci_lo'))}-{_fmt(s.get('pooled_ratio_ci_hi'))})")
                    p1 = f", Wilcoxon vs 0 p = {_fmt(s.get('p_wilcoxon_vs_0'), 3)}" if np.isfinite(_f(s.get("p_wilcoxon_vs_0"))) else ""
                    p2 = f", vs sham (Mann-Whitney) p = {_fmt(s.get('p_mannwhitney_vs_sham'), 3)}" if np.isfinite(_f(s.get("p_mannwhitney_vs_sham"))) else ""
                    lines.append(f"  {s['measure_label']}: mean {_fmt(s.get('mean'))} +- {_fmt(s.get('sem'))} SEM over {int(s['n'])} experiment(s)"
                                 f"{p1}{p2}{extra}")
            lines.append("")
        lines.append("Reading: a log2 ratio of +1 = the rate doubled, -1 = halved; 'few events' = fewer than 10 events in both conditions together "
                      "(the direction may be clear, the size is not). The still-only change is the one to trust when the animal ran a different "
                      "share of the time pre and post. With 2-6 experiments per protocol the tests are descriptive.")
    if notes:
        lines += ["", "Not included / to do:"] + [f"  - {n}" for n in notes]
    txt = "\n".join(lines)
    out_txt.write_text(txt, encoding="utf-8")
    return txt


def write_html(T: pd.DataFrame, G: pd.DataFrame, notes: list[str], txt: str, png: Path, out_html: Path, parent: Path) -> None:
    cols = ["animal", "day", "fov", "label", "group", "power", "n_events_baseline", "n_events_condition", "minutes_baseline", "minutes_condition",
            "rate_baseline", "rate_condition", "ratio_raw", "ratio_ci_lo", "ratio_ci_hi", "log2_rate_all", "low_counts", "still_frac_baseline",
            "still_frac_condition", "log2_rate_still", "log2_rate_running", "ratio_speed_matched", "behaviour_comparable", "amp_baseline",
            "amp_condition", "log2_amp", "duration_baseline", "duration_condition", "cfu_up", "cfu_down", "cfu_unchanged"]
    cols = [c for c in cols if c in T.columns]

    def table(df, cols):
        h = ["<table><thead><tr>" + "".join(f"<th>{html.escape(str(c))}</th>" for c in cols) + "</tr></thead><tbody>"]
        for _, r in df.iterrows():
            cells = []
            for c in cols:
                v = r.get(c)
                if isinstance(v, (float, np.floating)):
                    cells.append(_fmt(v, 3))
                else:
                    cells.append(html.escape(str(v)))
            h.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
        h.append("</tbody></table>")
        return "\n".join(h)

    gcols = [c for c in ["group", "measure_label", "n", "n_animals", "mean", "sem", "median", "p_wilcoxon_vs_0", "p_mannwhitney_vs_sham",
                         "pooled_rate_pre", "pooled_rate_post", "pooled_ratio", "pooled_ratio_ci_lo", "pooled_ratio_ci_hi"] if c in G.columns]
    doc = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Stimulation protocols across animals</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;margin:24px;color:#222;max-width:1500px}} table{{border-collapse:collapse;font-size:12px;margin:8px 0}}
th,td{{border:1px solid #ccc;padding:2px 6px;text-align:right}} th{{background:#f0f0f0}} td:first-child,th:first-child{{text-align:left}}
pre{{background:#f7f7f7;padding:10px;font-size:12px;white-space:pre-wrap}} img{{max-width:100%;border:1px solid #ddd}}</style></head><body>
<h1>Stimulation protocols compared across animals</h1>
<p>{html.escape(str(parent))} — made {time.strftime('%Y-%m-%d %H:%M')} by aqua_conditions.py from each animal's aqua_analysis tables (tab 7). No pupil step; locomotion matching as in the analysis.</p>
<h2>Figure</h2><p><a href="{png.name}" target="_blank"><img src="{png.name}" alt="conditions figure"></a> (click for full size)</p>
<h2>Summary</h2><pre>{html.escape(txt)}</pre>
<h2>Per protocol</h2>{table(G, gcols) if not G.empty else '<p>nothing</p>'}
<h2>Per experiment (day x FOV)</h2>{table(T, cols) if not T.empty else '<p>nothing</p>'}
</body></html>"""
    out_html.write_text(doc, encoding="utf-8")


# ----------------------------------------------------------------------------- main

def run(parent: Path, animals: list[str] | None = None, out_dir: Path | None = None, log=log) -> dict:
    parent = Path(parent)
    out_dir = Path(out_dir) if out_dir else parent / "cross_animal"
    out_dir.mkdir(parents=True, exist_ok=True)
    T, notes = gather(parent, animals)
    log(f"{len(T)} analysed experiment(s) from {T['animal'].nunique() if len(T) else 0} animal(s)")
    for n in notes:
        log(f"  note: {n}")
    G = summarise(T)
    T.to_csv(out_dir / "conditions_experiments.csv", index=False)
    G.to_csv(out_dir / "conditions_groups.csv", index=False)
    png = out_dir / "conditions.png"
    if len(T):
        figure(T[T["group"] != "none"] if (T["group"] != "none").any() else T, G, png, f"Stimulation protocols across animals — pre -> post change per day x FOV ({parent.name})")
    txt = write_text(T, G, notes, out_dir / "conditions.txt", parent)
    write_html(T, G, notes, txt, png, out_dir / "conditions.html", parent)
    log(txt)
    log(f"wrote {out_dir / 'conditions.html'}, conditions.png, conditions_experiments.csv, conditions_groups.csv, conditions.txt")
    return {"table": T, "groups": G, "notes": notes, "html": out_dir / "conditions.html", "png": png}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parent", help="folder that holds the mouse folders (e.g. D:\\2photon\\Chao\\GCAMP_affrerents)")
    ap.add_argument("--animals", default=None, help="comma list of mouse folder names (default: every *RearChr* folder)")
    ap.add_argument("--out", default=None, help="output folder (default <parent>\\cross_animal)")
    a = ap.parse_args(argv)
    animals = [s for s in (a.animals or "").split(",") if s.strip()] or None
    run(Path(a.parent), animals, Path(a.out) if a.out else None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
