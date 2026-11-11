"""The pre/post change at matched running speed, and how much of the recording
can be matched at all.

The afferents follow locomotion, and a session can have the mouse running the
whole time or sitting the whole time. Two conditions can then only be compared
where their behaviour overlaps, and the honest number is the ratio computed
*within* behaviour bands plus a statement of how much of each condition falls
in bands the other one also has.

Bands (set once per experiment on the pooled 20-s bins of all its runs):
    still              running fraction of the bin < 0.5
    slow / medium / fast   terciles of the mean speed of the running bins

Per condition and band: events, minutes, rate. The STANDARDISED rate of a
condition is the band rates averaged with one fixed set of weights — the
pooled minutes per band over all conditions — restricted to the bands in which
every condition has at least `min_minutes`; the standardised ratio is the
ratio of these (direct standardisation, as in epidemiology). The OVERLAP of a
condition is the share of its minutes that falls in such shared bands: at
100 % the whole condition is comparable, at 10 % the ratio rests on a sliver
and is reported as not comparable (below `min_overlap`).

The regression in adjust.py uses the same information with a model; this is
the model-free version that also says when no model can help, because the
conditions never visited the same state.
"""
from __future__ import annotations

import numpy as np

BANDS = ("still", "slow", "medium", "fast")


def _bins(glm_runs: list, conditions: list, spf: float, bin_s: float):
    from .adjust import bin_run
    tabs = []
    for ri, r in enumerate(glm_runs):
        t = bin_run(r["onsets"], r["excluded"], r["speed"], r["still"], r["pupil"], spf, bin_s)
        n = t["events"].size
        t["cond"] = np.full(n, conditions.index(r["condition"]) if r["condition"] in conditions else -1)
        t["run"] = np.full(n, ri)
        tabs.append(t)
    return {k: np.concatenate([t[k] for t in tabs]) for k in tabs[0]}


def speed_standardised(glm_runs: list, conditions: list, spf: float, bin_s: float = 20.0,
                       min_minutes: float = 1.0, min_overlap: float = 0.5) -> dict:
    """Band table + standardised ratios of every condition vs the first."""
    d = _bins(glm_runs, conditions, spf, bin_s)
    use = (d["minutes"] > 0.25 * bin_s / 60.0) & (d["cond"] >= 0) & np.isfinite(d["running"])
    out = {"bands": {}, "edges": None, "n_bins": int(use.sum())}
    if use.sum() < 6:
        out["why_not"] = "no behaviour data"
        return out
    running = d["running"] >= 0.5
    sp = np.where(np.isfinite(d["speed"]), d["speed"], 0.0)
    run_sp = sp[use & running]
    if run_sp.size >= 6:
        edges = np.percentile(run_sp, [100 / 3, 200 / 3])
    else:
        edges = np.array([np.inf, np.inf])
    out["edges"] = [float(e) for e in edges]
    band = np.full(d["events"].size, -1)
    band[use & ~running] = 0
    band[use & running & (sp < edges[0])] = 1
    band[use & running & (sp >= edges[0]) & (sp < edges[1])] = 2
    band[use & running & (sp >= edges[1])] = 3
    table = {}
    for ci, c in enumerate(conditions):
        for bi, bname in enumerate(BANDS):
            m = (d["cond"] == ci) & (band == bi)
            ev = float(d["events"][m].sum()); mins = float(d["minutes"][m].sum())
            table[(c, bname)] = {"events": ev, "minutes": mins, "rate": ev / mins if mins > 0 else np.nan}
    out["bands"] = table
    base = conditions[0]
    tot = {c: sum(table[(c, b)]["minutes"] for b in BANDS) for c in conditions}
    out["minutes_total"] = tot
    for c in conditions[1:]:
        shared = [b for b in BANDS if table[(base, b)]["minutes"] >= min_minutes and table[(c, b)]["minutes"] >= min_minutes]
        ov_b = sum(table[(base, b)]["minutes"] for b in shared) / tot[base] if tot[base] > 0 else 0.0
        ov_c = sum(table[(c, b)]["minutes"] for b in shared) / tot[c] if tot[c] > 0 else 0.0
        res = {"shared_bands": shared, "overlap_baseline": float(ov_b), "overlap_condition": float(ov_c),
               "comparable": bool(shared) and min(ov_b, ov_c) >= min_overlap}
        if shared:
            w = np.array([table[(base, b)]["minutes"] + table[(c, b)]["minutes"] for b in shared])
            rb = np.array([table[(base, b)]["rate"] for b in shared]); rc = np.array([table[(c, b)]["rate"] for b in shared])
            sb = float((w * rb).sum() / w.sum()); sc = float((w * rc).sum() / w.sum())
            res["rate_std_baseline"] = sb; res["rate_std_condition"] = sc
            res["ratio_std"] = sc / sb if sb > 0 else np.nan
            res["band_ratios"] = {b: (table[(c, b)]["rate"] / table[(base, b)]["rate"]
                                      if table[(base, b)]["rate"] and table[(base, b)]["rate"] > 0 else np.nan) for b in shared}
        else:
            res["ratio_std"] = np.nan
        out[c] = res
    return out


def describe(res: dict, base: str, c: str) -> str:
    """One summary line for the text report."""
    if "why_not" in res or c not in res:
        return f"  at matched speed: not computed ({res.get('why_not', 'no bins')})"
    r = res[c]
    tab = res["bands"]
    occ = ", ".join(f"{b} {tab[(base, b)]['minutes']:.0f}/{tab[(c, b)]['minutes']:.0f} min"
                    for b in BANDS)
    if not r["comparable"]:
        return (f"  at matched speed: NOT COMPARABLE — only {100 * r['overlap_baseline']:.0f} % of the {base} minutes and "
                f"{100 * r['overlap_condition']:.0f} % of the {c} minutes fall in speed bands both conditions visited "
                f"({base}/{c} minutes per band: {occ}); no statistic can take the behaviour out of this comparison "
                f"— the two conditions have to be recorded in overlapping states")
    br = ", ".join(f"{b} x{v:.2f}" for b, v in r["band_ratios"].items() if np.isfinite(v))
    if np.isfinite(r["ratio_std"]):
        head = f"x{r['ratio_std']:.2f} ({r['rate_std_condition']:.3f} vs {r['rate_std_baseline']:.3f}/min); per band {br or 'n/a'}"
    elif r["rate_std_condition"] > 0:
        head = f"{r['rate_std_condition']:.3f}/min vs 0 in {base} (no {base} events in the shared bands — no ratio)"
    else:
        head = "no events in the shared bands of either condition"
    return (f"  at matched speed (standardised over {', '.join(r['shared_bands'])}; weights = pooled minutes): {head}; "
            f"behaviour overlap {100 * r['overlap_baseline']:.0f} % ({base}) / {100 * r['overlap_condition']:.0f} % ({c}) "
            f"of the minutes; {base}/{c} minutes per band: {occ}")


def band_rows(exp_name: str, res: dict, conditions: list) -> list[dict]:
    rows = []
    for c in conditions:
        for b in BANDS:
            t = res["bands"].get((c, b))
            if t is None:
                continue
            rows.append({"experiment": exp_name, "condition": c, "band": b,
                         "speed_lo": 0 if b == "still" else (0 if b == "slow" else res["edges"][0] if b == "medium" else res["edges"][1]),
                         "events": t["events"], "minutes": round(t["minutes"], 2),
                         "events_per_min": round(t["rate"], 4) if np.isfinite(t["rate"]) else np.nan})
    return rows


def tuning_figure(results: list, out_png, log=print):
    """results: [{experiment, conditions, res, label}] -> one panel per experiment:
    rate per speed band, one bar per condition, minutes annotated; title with
    the standardised ratio and the overlap."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    results = [r for r in results if r["res"].get("bands")]
    if not results:
        return None
    n = len(results)
    ncol = min(3, n); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow), squeeze=False)
    cols = ["#3b7dd8", "#e8702a", "#2ca02c", "#9467bd"]
    for ax, r in zip(axes.ravel(), results):
        res, conds = r["res"], r["conditions"]
        tab = res["bands"]
        x = np.arange(len(BANDS)); w = 0.8 / len(conds)
        ymax = 0.0
        for ci, c in enumerate(conds):
            rates = [tab[(c, b)]["rate"] if np.isfinite(tab[(c, b)]["rate"]) else 0.0 for b in BANDS]
            mins = [tab[(c, b)]["minutes"] for b in BANDS]
            ax.bar(x + (ci - (len(conds) - 1) / 2) * w, rates, w, color=cols[ci % len(cols)], alpha=0.85, label=c)
            for xi, (rt, mn) in enumerate(zip(rates, mins)):
                if mn > 0:
                    ax.text(xi + (ci - (len(conds) - 1) / 2) * w, rt, f"{mn:.0f}'", ha="center", va="bottom", fontsize=6.5, color="#333")
                    ymax = max(ymax, rt)
        e = res.get("edges") or [np.nan, np.nan]
        ax.set_xticks(x)
        ax.set_xticklabels(["still", f"slow\n<{e[0]:.0f}", f"medium\n{e[0]:.0f}-{e[1]:.0f}", f"fast\n>{e[1]:.0f}"], fontsize=8)
        ax.set_ylabel("events / min", fontsize=8)
        ax.set_ylim(0, ymax * 1.3 + 1e-3)
        ax.tick_params(labelsize=8)
        c1 = conds[1] if len(conds) > 1 else None
        rr = res.get(c1, {}) if c1 else {}
        if rr.get("comparable") and np.isfinite(rr.get("ratio_std", np.nan)):
            t = f"matched-speed ratio x{rr['ratio_std']:.2f}  (overlap {100 * rr['overlap_baseline']:.0f} % / {100 * rr['overlap_condition']:.0f} %)"
        elif rr.get("comparable"):
            t = f"no ratio: no {conds[0]} events in the shared bands  (overlap {100 * rr['overlap_baseline']:.0f} % / {100 * rr['overlap_condition']:.0f} %)"
        elif c1:
            t = f"NOT comparable: overlap {100 * rr.get('overlap_baseline', 0):.0f} % / {100 * rr.get('overlap_condition', 0):.0f} %"
        else:
            t = ""
        ax.set_title(f"{r.get('label', r['experiment'])}\n{t}", fontsize=8.5, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("Event rate per running-speed band (numbers = minutes in the band); a change that holds in the shared bands "
                 "is not a change of state", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png
