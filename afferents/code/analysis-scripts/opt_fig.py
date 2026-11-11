"""Unit-level optimisation of the overlap threshold across all runs (chains included): figure + numbers."""
import glob, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
OUT = "/mnt/user-data/outputs/cfu_rules/overlap_threshold_optimisation.png"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TXT, TXT2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"

z = np.load(SP + "/opt_threshold_pairs.npz", allow_pickle=True); A = z["pairs"]; runs = list(z["runs"])
v, d, ai, aj = A[:, 3], A[:, 5], A[:, 8], A[:, 9]
ri, rj = np.sqrt(ai / np.pi), np.sqrt(aj / np.pi); rmin, rmax = np.minimum(ri, rj), np.maximum(ri, rj)
pos = d <= 0.5 * rmin           # the two cores coincide: same place by AQuA2's definition (consistent spatial foundation)
neg = d > rmax                  # each centre lies outside the other's extent: two places
R = {json.load(open(f))["run"]: json.load(open(f)) for f in glob.glob(SP + "/overlap_all/*.json")}
thrs = np.round(np.arange(0.10, 0.51, 0.025), 3)
run_idx = A[:, 0].astype(int); I = A[:, 1].astype(int); J = A[:, 2].astype(int)


def uf_labels(n, pairs, t):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j, vv in pairs:
        if vv > t:
            a, b = find(i), find(j)
            if a != b:
                parent[b] = a
    return [find(i) for i in range(n)]


animals = ("2608-1", "2608-2", "2609-1")
res = {"all": [], **{a: [] for a in animals}}
for t in thrs:
    same = np.zeros(len(A), bool)
    for k, run in enumerate(runs):
        r = R[run]; lab = uf_labels(r["n_cfu"], r["pairs"], t)
        sel = run_idx == k
        same[sel] = [lab[i] == lab[j] for i, j in zip(I[sel], J[sel])]
    for key in res:
        sel = np.ones(len(A), bool) if key == "all" else np.array([key in runs[k] for k in run_idx])
        p, n = pos & sel, neg & sel
        tp = same[p].sum(); fp = same[n].sum(); fn = (~same[p]).sum()
        res[key].append((same[p].mean(), same[n].mean(), 2 * tp / max(2 * tp + fp + fn, 1)))
for key in res:
    res[key] = np.array(res[key])
    f1 = res[key][:, 2]; jj = res[key][:, 0] - res[key][:, 1]
    print(f"{key}: POS {int((pos & (np.ones(len(A), bool) if key == 'all' else np.array([key in runs[k] for k in run_idx]))).sum())} "
          f"NEG {int((neg & (np.ones(len(A), bool) if key == 'all' else np.array([key in runs[k] for k in run_idx]))).sum())} "
          f"| best F1 {f1.max():.3f} at t={thrs[f1.argmax()]} | best J {jj.max():.3f} at t={thrs[jj.argmax()]}")

fig = plt.figure(figsize=(13, 4.8), facecolor=SURF)
gs = fig.add_gridspec(1, 3, wspace=0.3, left=0.05, right=0.99, top=0.78, bottom=0.16)
ax = fig.add_subplot(gs[0, 0])
ax.plot(thrs, res["all"][:, 0], color=PAL[0], lw=2, marker="o", ms=4, label="same-place pairs that end up in ONE unit (want high)")
ax.plot(thrs, res["all"][:, 1], color=PAL[1], lw=2, marker="s", ms=4, label="two-place pairs that end up in one unit (want low)")
ax.axvline(0.3, color=TXT2, lw=1, ls=(0, (3, 3)))
ax.set_xlim(0.52, 0.08); ax.set_ylim(0, 1.02)
ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2); ax.set_ylabel("fraction of labelled pairs", fontsize=8.5, color=TXT2)
ax.legend(fontsize=7.3, frameon=False, loc="upper left")
ax.set_title(f"A  Units after single linkage (chains included), {len(runs)} runs\n"
             f"{int(pos.sum()):,} same-place pairs, {int(neg.sum()):,} two-place pairs", fontsize=9.5, color=TXT, loc="left")
ax = fig.add_subplot(gs[0, 1])
ax.plot(thrs, res["all"][:, 2], color=PAL[0], lw=2.4, marker="o", ms=5, label="all three animals")
for a, c in zip(animals, (PAL[2], PAL[3], PAL[4])):
    ax.plot(thrs, res[a][:, 2], color=c, lw=1.2, marker="o", ms=3, label=a)
ax.axvline(0.3, color=TXT2, lw=1, ls=(0, (3, 3)))
ax.set_xlim(0.52, 0.08); ax.set_ylim(0, 0.85)
ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2); ax.set_ylabel("F1 of 'same unit' vs the spatial labels", fontsize=8.5, color=TXT2)
ax.legend(fontsize=7.5, frameon=False, loc="lower center")
bt = thrs[res["all"][:, 2].argmax()]
ax.set_title(f"B  Agreement with the labels (F1): best at {bt} pooled\nper animal: " +
             ", ".join(f"{a} {thrs[res[a][:, 2].argmax()]}" for a in animals), fontsize=9.5, color=TXT, loc="left")
ax = fig.add_subplot(gs[0, 2])
jj = res["all"][:, 0] - res["all"][:, 1]
ax.plot(thrs, jj, color=PAL[0], lw=2.4, marker="o", ms=5, label="all three animals")
for a, c in zip(animals, (PAL[2], PAL[3], PAL[4])):
    ax.plot(thrs, res[a][:, 0] - res[a][:, 1], color=c, lw=1.2, marker="o", ms=3, label=a)
ax.axvline(0.3, color=TXT2, lw=1, ls=(0, (3, 3)))
ax.set_xlim(0.52, 0.08); ax.set_ylim(0, 0.8)
ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2); ax.set_ylabel("merged same-place − merged two-place (Youden J)", fontsize=8.5, color=TXT2)
ax.legend(fontsize=7.5, frameon=False, loc="lower center")
ax.set_title(f"C  Balance (Youden J): best at {thrs[jj.argmax()]} pooled\nper animal: " +
             ", ".join(f"{a} {thrs[(res[a][:, 0] - res[a][:, 1]).argmax()]}" for a in animals), fontsize=9.5, color=TXT, loc="left")
for a in fig.axes:
    a.set_facecolor(SURF)
    for s in a.spines.values():
        s.set_edgecolor(GRID)
    a.tick_params(colors=TXT2, labelsize=8); a.grid(axis="y", color=GRID, lw=0.6)
fig.suptitle("Overlap threshold optimised on all runs of the three RearChr animals — labels that do not use the overlap:\n"
             "same place = the two events' weighted centres lie within half the smaller event's radius of each other; two places = each centre lies outside the other event's extent",
             fontsize=10.5, color=TXT, x=0.05, ha="left", y=0.985, va="top")
fig.savefig(OUT, dpi=150, facecolor=SURF)
print("saved", OUT)
