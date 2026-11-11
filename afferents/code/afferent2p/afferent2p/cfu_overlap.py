"""cfu_overlap — the overlap-threshold check for AQuA2's CFU rule (v42).

The pipeline defines a CFU with AQuA2's published rule: events whose duration-
weighted footprints overlap (weighted intersection over union) by at least the
overlap threshold join one unit, single linkage. The threshold is one number
for the dataset (0.3, chosen on all recordings of the three RearChr animals on
2 Nov 2026). This module re-asks the question on ANY animal from its own
results and says whether that number still fits:

  * the CFU batch (v42 scripts) writes next to each result
    <run>_AQuA2_cfu_pre.mat = AQuA2's own pairwise overlap table
    (cfu.CFU_tmp_function's s_t0) plus, per event, the duration-weighted
    centre, the footprint size and the first/last frame;
  * every overlapping pair of events is labelled WITHOUT using the overlap:
    same place  — the two centres lie within POS_RADIUS x the smaller event's
                  equivalent radius of each other (AQuA2's "consistent spatial
                  foundation": one place, events of different extent);
    two places  — each centre lies outside the other event's extent
                  (distance > the larger radius); the rest is not scored;
  * for every threshold the events of each recording are clustered as AQuA2
    does it (single linkage, so chains count) and the agreement of the units
    with the labels is measured (F1 of "same unit"); per experiment
    (day x FOV, pre and post together), per animal, and pooled.

Verdict per experiment: OK when the chosen threshold reaches at least
OK_RATIO of that experiment's best F1 (or the experiment has too few labelled
pairs to say anything: < MIN_POS same-place pairs), otherwise CHECK with the
experiment's own optimum. The outputs land in <mouse>\\aqua_batch\\:
cfu_overlap_check.csv (one row per experiment + pooled), cfu_overlap_check.txt
(plain words, read by the report) and cfu_overlap_check.png.

It changes nothing by itself: the threshold stays the one in tab 4 unless you
change it there (which redoes every CFU file and the curation).
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np

from .aqua_io import aqua_files, read_mat_fields, relocate

THRESHOLDS = tuple(round(x, 3) for x in np.arange(0.10, 0.5001, 0.025))
POS_RADIUS = 0.5          # same place: centre distance <= POS_RADIUS x the smaller equivalent radius
MIN_POS = 30              # fewer same-place pairs than this: the experiment cannot judge the threshold
OK_RATIO = 0.90           # chosen F1 / best F1 at or above this: the chosen threshold fits the experiment


# ----------------------------------------------------------------------------- reading

def read_cfu_pre(path: Path) -> dict | None:
    """<run>_AQuA2_cfu_pre.mat -> {'pairs': (n,3) [i, j, overlap] 1-based, 'centre': (nE,2) [y, x],
    'n_pixels': (nE,), 't0_t1': (nE,2), 'overlap_thr': float, 'min_num_evt': int}; None if unreadable."""
    try:
        d = read_mat_fields(Path(path), ["pairs", "centre_yx", "n_pixels", "t0_t1", "overlap_thr", "min_num_evt"])
    except Exception:
        return None
    if d.get("pairs") is None or d.get("centre_yx") is None:
        return None
    pairs = np.asarray(d["pairs"], float).reshape(-1, 3) if np.size(d["pairs"]) else np.zeros((0, 3))
    centre = np.asarray(d["centre_yx"], float).reshape(-1, 2) if np.size(d["centre_yx"]) else np.zeros((0, 2))
    npx = np.asarray(d["n_pixels"], float).ravel() if d.get("n_pixels") is not None else np.zeros(len(centre))
    tt = np.asarray(d["t0_t1"], float).reshape(-1, 2) if d.get("t0_t1") is not None and np.size(d["t0_t1"]) else np.zeros((len(centre), 2))
    thr = float(np.asarray(d["overlap_thr"]).ravel()[0]) if d.get("overlap_thr") is not None else float("nan")
    mne = int(np.asarray(d["min_num_evt"]).ravel()[0]) if d.get("min_num_evt") is not None else -1
    return {"pairs": pairs, "centre": centre, "n_pixels": npx, "t0_t1": tt, "overlap_thr": thr, "min_num_evt": mne}


# ----------------------------------------------------------------------------- the measurement

def label_pairs(pre: dict, pos_radius: float = POS_RADIUS):
    """For every overlapping pair: (i, j, overlap, is_same_place, is_two_places, centre distance px)."""
    P = pre["pairs"]
    if len(P) == 0:
        return np.zeros((0, 3)), np.zeros(0, bool), np.zeros(0, bool), np.zeros(0)
    i = P[:, 0].astype(int) - 1
    j = P[:, 1].astype(int) - 1
    v = P[:, 2]
    c = pre["centre"]
    r = np.sqrt(np.maximum(pre["n_pixels"], 1) / np.pi)
    d = np.hypot(c[i, 0] - c[j, 0], c[i, 1] - c[j, 1])
    rmin = np.minimum(r[i], r[j])
    rmax = np.maximum(r[i], r[j])
    pos = d <= pos_radius * rmin
    neg = d > rmax
    return np.column_stack([i, j, v]), pos, neg, d


def cluster_labels(n: int, pairs: np.ndarray, thr: float) -> np.ndarray:
    """Single linkage exactly as cfu.hierarchicalClusteringMinMeasure: pairs with overlap > thr join."""
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j, v in pairs:
        if v > thr:
            a, b = find(int(i)), find(int(j))
            if a != b:
                parent[b] = a
    return np.array([find(k) for k in range(n)])


def same_unit_matrix(n_events: int, pairs: np.ndarray, thresholds=THRESHOLDS) -> np.ndarray:
    """(n_thresholds, n_pairs) bool: are the two events of each pair in one unit at that threshold."""
    out = np.zeros((len(thresholds), len(pairs)), bool)
    if len(pairs) == 0:
        return out
    I = pairs[:, 0].astype(int)
    J = pairs[:, 1].astype(int)
    for t_idx, t in enumerate(thresholds):
        lab = cluster_labels(n_events, pairs, t)
        out[t_idx] = lab[I] == lab[J]
    return out


def f1_curve(same: np.ndarray, pos: np.ndarray, neg: np.ndarray) -> np.ndarray:
    """F1 of 'same unit' against the labels, per threshold (NaN when there is no same-place pair)."""
    if same.shape[1] == 0 or pos.sum() == 0:
        return np.full(same.shape[0], np.nan)
    tp = same[:, pos].sum(1)
    fp = same[:, neg].sum(1)
    fn = (~same[:, pos]).sum(1)
    return 2 * tp / np.maximum(2 * tp + fp + fn, 1)


def tpr_fpr(same: np.ndarray, pos: np.ndarray, neg: np.ndarray):
    tpr = same[:, pos].mean(1) if pos.sum() else np.full(same.shape[0], np.nan)
    fpr = same[:, neg].mean(1) if neg.sum() else np.full(same.shape[0], np.nan)
    return tpr, fpr


def nearest_threshold(t: float, thresholds=THRESHOLDS) -> int:
    return int(np.argmin(np.abs(np.asarray(thresholds) - t)))


# ----------------------------------------------------------------------------- the check

def check_mouse(mouse_dir: Path, chosen: float | None = None, pos_radius: float = POS_RADIUS,
                min_pos: int = MIN_POS, ok_ratio: float = OK_RATIO, log=print) -> dict:
    """Run the check on one animal (reads aqua_batch\\manifest.json for the experiments and runs).
    Returns the result dict (also written to aqua_batch\\cfu_overlap_check.csv/.txt/.png)."""
    mouse_dir = Path(mouse_dir)
    bdir = mouse_dir / "aqua_batch"
    mp = bdir / "manifest.json"
    if not mp.exists():
        raise FileNotFoundError(f"{mp} not found: build the batch first (tab 4)")
    man = json.loads(mp.read_text())
    if chosen is None:
        chosen = float((man.get("cfu") or {}).get("overlap", 0.3))
    thresholds = THRESHOLDS
    t_chosen = nearest_threshold(chosen, thresholds)
    runs_out, missing, stale = [], [], []
    per_exp: dict[str, dict] = {}
    for e in man.get("experiments", []):
        for r in e.get("runs", []):
            tif = Path(relocate(r["tif"], mouse_dir))     # the manifest's path, re-rooted when the folder moved
            f = aqua_files(tif)
            if not f["cfu"].exists():
                continue                      # no CFU result yet: nothing to check
            if not f["cfu_pre"].exists():
                missing.append(r["run"]); continue
            if f["cfu_pre"].stat().st_mtime + 2 < f["cfu"].stat().st_mtime - 3600:
                stale.append(r["run"])
            pre = read_cfu_pre(f["cfu_pre"])
            if pre is None:
                missing.append(r["run"]); continue
            n_ev = len(pre["centre"])
            pairs, pos, neg, dist = label_pairs(pre, pos_radius)
            same = same_unit_matrix(n_ev, pairs, thresholds)
            rec = {"run": r["run"], "experiment": e["name"], "condition": r.get("condition", ""), "n_events": n_ev,
                   "n_pairs": int(len(pairs)), "n_pos": int(pos.sum()), "n_neg": int(neg.sum()),
                   "same": same, "pos": pos, "neg": neg, "overlap_thr_in_file": pre["overlap_thr"],
                   "n_units_chosen": int(len(np.unique(cluster_labels(n_ev, pairs, thresholds[t_chosen])))) if n_ev else 0}
            runs_out.append(rec)
            x = per_exp.setdefault(e["name"], {"experiment": e["name"], "runs": [], "same": [], "pos": [], "neg": []})
            x["runs"].append(r["run"]); x["same"].append(same); x["pos"].append(pos); x["neg"].append(neg)

    def summarise(name, sames, poss, negs):
        same = np.concatenate(sames, axis=1) if sames else np.zeros((len(thresholds), 0), bool)
        pos = np.concatenate(poss) if poss else np.zeros(0, bool)
        neg = np.concatenate(negs) if negs else np.zeros(0, bool)
        f1 = f1_curve(same, pos, neg)
        tpr, fpr = tpr_fpr(same, pos, neg)
        n_pos, n_neg = int(pos.sum()), int(neg.sum())
        if n_pos < min_pos or np.all(np.isnan(f1)):
            best_t, best_f1, verdict = float("nan"), float("nan"), "too few pairs"
            ratio = float("nan")
        else:
            b = int(np.nanargmax(f1)); best_t, best_f1 = float(thresholds[b]), float(f1[b])
            ratio = float(f1[t_chosen] / best_f1) if best_f1 > 0 else float("nan")
            verdict = "OK" if ratio >= ok_ratio else "CHECK"
        return {"name": name, "n_pos": n_pos, "n_neg": n_neg, "f1": f1, "tpr": tpr, "fpr": fpr,
                "best_t": best_t, "best_f1": best_f1, "f1_chosen": float(f1[t_chosen]) if not np.all(np.isnan(f1)) else float("nan"),
                "tpr_chosen": float(tpr[t_chosen]) if not np.all(np.isnan(tpr)) else float("nan"),
                "fpr_chosen": float(fpr[t_chosen]) if not np.all(np.isnan(fpr)) else float("nan"),
                "ratio": ratio, "verdict": verdict}

    exps = [summarise(k, x["same"], x["pos"], x["neg"]) | {"n_runs": len(x["runs"])} for k, x in per_exp.items()]
    pooled = summarise("all experiments", [r["same"] for r in runs_out], [r["pos"] for r in runs_out], [r["neg"] for r in runs_out])
    pooled["n_runs"] = len(runs_out)
    res = {"mouse": mouse_dir.name, "chosen": float(thresholds[t_chosen]), "chosen_requested": float(chosen),
           "thresholds": list(thresholds), "pos_radius": pos_radius, "min_pos": min_pos, "ok_ratio": ok_ratio,
           "experiments": exps, "pooled": pooled, "runs": runs_out, "missing": missing, "stale": stale,
           "n_units_chosen": int(sum(r["n_units_chosen"] for r in runs_out)),
           "n_events": int(sum(r["n_events"] for r in runs_out)), "made": time.strftime("%Y-%m-%d %H:%M")}
    bdir.mkdir(exist_ok=True)
    write_outputs(res, bdir, log=log)
    return res


# ----------------------------------------------------------------------------- outputs

def verdict_words(res: dict) -> str:
    """The plain-words summary (also the .txt)."""
    t = res["chosen"]
    p = res["pooled"]
    lines = [f"CFU overlap-threshold check — {res['mouse']} — {res['made']}",
             f"chosen threshold {t:g} (tab 4 / the batch); {p['n_runs']} recording(s) with results, {res['n_events']} events, "
             f"{res['n_units_chosen']} units at {t:g}; labels: same place = centres within {res['pos_radius']:g} x the smaller radius, "
             f"two places = each centre outside the other's extent"]
    if np.isnan(p["best_t"]):
        lines.append(f"POOLED: too few labelled pairs ({p['n_pos']} same-place, {p['n_neg']} two-place) to judge the threshold on this animal.")
    else:
        lines.append(f"POOLED: best agreement at {p['best_t']:g} (F1 {p['best_f1']:.2f}); at {t:g}: F1 {p['f1_chosen']:.2f} "
                     f"= {100 * p['ratio']:.0f} % of the best, {100 * p['tpr_chosen']:.0f} % of same-place pairs in one unit, "
                     f"{100 * p['fpr_chosen']:.0f} % of two-place pairs in one unit "
                     f"({p['n_pos']} same-place / {p['n_neg']} two-place pairs) -> {p['verdict']}")
    n_check = sum(1 for e in res["experiments"] if e["verdict"] == "CHECK")
    n_ok = sum(1 for e in res["experiments"] if e["verdict"] == "OK")
    n_few = sum(1 for e in res["experiments"] if e["verdict"] == "too few pairs")
    lines.append(f"per experiment (day x FOV, pre and post together): {n_ok} OK, {n_check} CHECK, {n_few} with too few pairs to judge")
    for e in sorted(res["experiments"], key=lambda x: x["name"]):
        if e["verdict"] == "too few pairs":
            lines.append(f"  {e['name']}: too few pairs ({e['n_pos']} same-place / {e['n_neg']} two-place, {e['n_runs']} runs)")
        else:
            lines.append(f"  {e['name']}: {e['verdict']} — best {e['best_t']:g} (F1 {e['best_f1']:.2f}), at {t:g} F1 {e['f1_chosen']:.2f} "
                         f"= {100 * e['ratio']:.0f} % ({e['n_pos']} / {e['n_neg']} pairs, {e['n_runs']} runs)")
    if n_check:
        lines.append(f"CHECK means: that experiment's own best threshold gives a clearly better unit/label agreement than {t:g} "
                     f"(below {100 * res['ok_ratio']:.0f} % of its best). One threshold per dataset is still the rule — change it in tab 4 "
                     "only if several experiments say so, and know that every CFU file and the curation are then redone.")
    else:
        lines.append(f"Verdict: {t:g} fits every experiment that can judge it; keep it.")
    if res["missing"]:
        lines.append(f"{len(res['missing'])} result(s) without an overlap table (_cfu_pre.mat, written by the v42 CFU batch) — not included: "
                     + ", ".join(res["missing"][:6]) + (" …" if len(res["missing"]) > 6 else "")
                     + ". Re-run the CFU batch (tab 5) to get them.")
    if res["stale"]:
        lines.append(f"{len(res['stale'])} overlap table(s) are much older than their CFU result (made by different runs?): " + ", ".join(res["stale"][:6]))
    return "\n".join(lines)


def write_outputs(res: dict, bdir: Path, log=print):
    import csv
    thresholds = res["thresholds"]
    rows = res["experiments"] + [res["pooled"]]
    with open(bdir / "cfu_overlap_check.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["experiment", "n_runs", "n_same_place_pairs", "n_two_place_pairs", "best_threshold", "best_f1",
                    "chosen_threshold", "f1_at_chosen", "ratio", "same_place_merged_at_chosen", "two_place_merged_at_chosen", "verdict"]
                   + [f"f1_at_{t:g}" for t in thresholds])
        for e in rows:
            w.writerow([e["name"], e["n_runs"], e["n_pos"], e["n_neg"], e["best_t"] if not np.isnan(e["best_t"]) else "",
                        round(e["best_f1"], 4) if not np.isnan(e["best_f1"]) else "",
                        res["chosen"], round(e["f1_chosen"], 4) if not np.isnan(e["f1_chosen"]) else "",
                        round(e["ratio"], 3) if not np.isnan(e["ratio"]) else "",
                        round(e["tpr_chosen"], 3) if not np.isnan(e["tpr_chosen"]) else "",
                        round(e["fpr_chosen"], 3) if not np.isnan(e["fpr_chosen"]) else "", e["verdict"]]
                       + [round(float(x), 4) if not np.isnan(x) else "" for x in e["f1"]])
    txt = verdict_words(res)
    (bdir / "cfu_overlap_check.txt").write_text(txt, encoding="utf-8")
    try:
        figure(res, bdir / "cfu_overlap_check.png")
    except Exception as e:                    # the picture is a convenience; the csv/txt are the result
        log(f"(picture not drawn: {e})")
    log(txt)


def figure(res: dict, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    TXT, TXT2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
    thr = np.asarray(res["thresholds"])
    t = res["chosen"]
    p = res["pooled"]
    fig = plt.figure(figsize=(13, 4.8), facecolor=SURF)
    gs = fig.add_gridspec(1, 3, wspace=0.3, left=0.05, right=0.99, top=0.76, bottom=0.16)
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(thr, p["tpr"], color=PAL[0], lw=2, marker="o", ms=4, label="same-place pairs in one unit (want high)")
    ax.plot(thr, p["fpr"], color=PAL[1], lw=2, marker="s", ms=4, label="two-place pairs in one unit (want low)")
    ax.axvline(t, color=TXT2, lw=1, ls=(0, (3, 3)))
    ax.set_xlim(0.52, 0.08); ax.set_ylim(0, 1.02)
    ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2); ax.set_ylabel("fraction of labelled pairs", fontsize=8.5, color=TXT2)
    ax.legend(fontsize=7.3, frameon=False, loc="upper left")
    ax.set_title(f"A  Units after single linkage, all recordings\n{p['n_pos']:,} same-place pairs, {p['n_neg']:,} two-place pairs",
                 fontsize=9.5, color=TXT, loc="left")
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(thr, p["f1"], color=PAL[0], lw=2.4, marker="o", ms=5)
    ax.axvline(t, color=TXT2, lw=1, ls=(0, (3, 3)))
    ax.set_xlim(0.52, 0.08); ax.set_ylim(0, 1.0)
    ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2); ax.set_ylabel("F1 of 'same unit' vs the labels", fontsize=8.5, color=TXT2)
    sub = (f"best at {p['best_t']:g}, at {t:g} {100 * p['ratio']:.0f} % of it -> {p['verdict']}" if not np.isnan(p["best_t"])
           else "too few labelled pairs to judge")
    ax.set_title(f"B  Agreement with the labels, pooled\n{sub}", fontsize=9.5, color=TXT, loc="left")
    ax = fig.add_subplot(gs[0, 2])
    exps = [e for e in sorted(res["experiments"], key=lambda x: x["name"]) if e["verdict"] != "too few pairs"]
    for k, e in enumerate(exps):
        short = e["name"].split("_", 1)[1] if "_" in e["name"] else e["name"]
        ax.plot(thr, e["f1"], color=PAL[k % 8], lw=1.3 if e["verdict"] == "OK" else 2.2,
                ls="-" if e["verdict"] == "OK" else (0, (4, 2)), marker="o", ms=3,
                label=f"{short}: best {e['best_t']:g} ({e['verdict']})")
    ax.axvline(t, color=TXT2, lw=1, ls=(0, (3, 3)))
    ax.set_xlim(0.52, 0.08); ax.set_ylim(0, 1.0)
    ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2); ax.set_ylabel("F1 per experiment", fontsize=8.5, color=TXT2)
    if exps:
        ax.legend(fontsize=6.3, frameon=False, loc="lower center", ncol=1 if len(exps) <= 8 else 2)
    n_few = sum(1 for e in res["experiments"] if e["verdict"] == "too few pairs")
    ax.set_title(f"C  Per day x FOV (pre + post): {len(exps)} judged, {n_few} too few pairs\n"
                 f"dashed = CHECK (chosen < {100 * res['ok_ratio']:.0f} % of its own best)", fontsize=9.5, color=TXT, loc="left")
    for a in fig.axes:
        a.set_facecolor(SURF)
        for s in a.spines.values():
            s.set_edgecolor(GRID)
        a.tick_params(colors=TXT2, labelsize=8); a.grid(axis="y", color=GRID, lw=0.6)
    fig.suptitle(f"CFU overlap-threshold check — {res['mouse']} — chosen threshold {t:g} (dashed line)\n"
                 f"labels without the overlap: same place = centres within {res['pos_radius']:g} x the smaller radius; "
                 "two places = each centre outside the other's extent", fontsize=10, color=TXT, x=0.05, ha="left", y=0.985, va="top")
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)
