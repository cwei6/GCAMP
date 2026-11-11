"""Optimise AQuA2's overlap threshold across ALL runs of the three RearChr animals with labels that do not use the overlap:
  NEG (two places for sure): the two units were active at the same time (one place cannot fire twice at once)
  POS (one place): not co-active, and their duration-weighted centres lie within R px of each other
Everything else is ambiguous and not scored. For each threshold t: TPR = POS pairs with overlap >= t, FPR = NEG pairs with
overlap >= t; the optimum maximises Youden's J = TPR - FPR (also reported: F1 on the labelled pairs)."""
import glob, json, sys
import numpy as np
sys.path.insert(0, "/home/claude/afferent2p")
from afferent2p.aqua_io import read_mat_fields, _cfu_records

SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
OUT = SP + "/opt_threshold_pairs.npz"

rows = []        # run_idx, i, j, overlap, coactive, centre distance (px), n_ev_i, n_ev_j, area_i, area_j, same-day gap (frames)
runs = []
for k, f in enumerate(sorted(glob.glob(SP + "/overlap_all/*.json"))):
    r = json.load(open(f))
    if r["n_cfu"] < 2 or not r["pairs"]:
        continue
    recs = _cfu_records(read_mat_fields(r["path"], ["cfuInfo1"])["cfuInfo1"])
    assert len(recs) == r["n_cfu"], (r["run"], len(recs), r["n_cfu"])
    cents, tws, areas = [], [], []
    for rec in recs:
        m = np.asarray(rec["weightMap"], np.float32)
        if m.ndim == 3:
            m = m[:, :, 0] if m.shape[2] == 1 else m.max(axis=2)
        ys, xs = np.nonzero(m > 0)
        w = m[ys, xs]
        cents.append((float((ys * w).sum() / w.sum()), float((xs * w).sum() / w.sum())) if w.sum() > 0 else (np.nan, np.nan))
        tw = np.asarray(rec["timeWindow"]).ravel().astype(bool)
        tws.append(tw)
        areas.append(int(len(ys)))
    runs.append(r["run"])
    for (i, j, v) in r["pairs"]:
        co = bool(np.any(tws[i] & tws[j]))
        # gap between the two activity windows (frames; 0 when they touch or overlap)
        ti = np.flatnonzero(tws[i]); tj = np.flatnonzero(tws[j])
        gap = 0 if co else int(max(tj.min() - ti.max(), ti.min() - tj.max()))
        d = float(np.hypot(cents[i][0] - cents[j][0], cents[i][1] - cents[j][1]))
        rows.append((len(runs) - 1, i, j, v, co, d, r["n_events_per_cfu"][i], r["n_events_per_cfu"][j], areas[i], areas[j], gap))
A = np.array(rows, dtype=float)
np.savez(OUT, pairs=A, runs=np.array(runs))
print("runs with pairs:", len(runs), "pairs:", len(A))
v, co, d = A[:, 3], A[:, 4].astype(bool), A[:, 5]
print(f"co-active pairs: {co.sum()} ({100 * co.mean():.1f} %)")
thrs = np.round(np.arange(0.10, 0.51, 0.025), 3)
for R in (3, 4, 6):
    pos = (~co) & (d <= R)
    neg = co
    print(f"\n--- POS = not co-active & centres <= {R} px ({R * 0.792:.1f} um): {pos.sum()} pairs; NEG = co-active: {neg.sum()} pairs")
    best = None
    for t in thrs:
        tpr = (v[pos] >= t).mean(); fpr = (v[neg] >= t).mean()
        tp = (v[pos] >= t).sum(); fp = (v[neg] >= t).sum(); fn = (v[pos] < t).sum()
        f1 = 2 * tp / max(2 * tp + fp + fn, 1)
        J = tpr - fpr
        if best is None or J > best[1]:
            best = (t, J)
        print(f"  t={t:.3f}  TPR={tpr:.3f}  FPR={fpr:.3f}  J={J:.3f}  F1={f1:.3f}  merged pairs={int((v >= t).sum())}")
    print("  best J at t =", best)
# per animal (R = 4)
pos = (~co) & (d <= 4); neg = co
for a in ("2608-1", "2608-2", "2609-1"):
    sel = np.array([a in runs[int(k)] for k in A[:, 0]])
    Js = [((v[pos & sel] >= t).mean() - (v[neg & sel] >= t).mean(), t) for t in thrs]
    print(a, "pairs", sel.sum(), "POS", (pos & sel).sum(), "NEG", (neg & sel).sum(), "best t:", max(Js)[1], "J=%.3f" % max(Js)[0])
# what do co-active pairs look like: overlap distribution
print("co-active pairs overlap quantiles:", np.round(np.quantile(v[co], [.5, .75, .9, .95, .99]), 2))
print("POS(4px) pairs overlap quantiles:", np.round(np.quantile(v[pos], [.05, .1, .25, .5, .75]), 2))
