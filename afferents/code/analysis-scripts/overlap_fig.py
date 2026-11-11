"""Pooled overlap distribution and threshold effects across the 88 runs of the three RearChr animals."""
import glob, json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/claude/afferent2p")
from afferent2p.aqua_io import read_cfu

SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
OUT = "/mnt/user-data/outputs/cfu_rules/overlap_threshold_all_runs.png"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TXT, TXT2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
UM = 0.792
THRS = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]

R = {}
for f in sorted(glob.glob(SP + "/overlap_all/*.json")):
    r = json.load(open(f)); R[r["run"]] = r
vals = np.array([v for r in R.values() for (_, _, v) in r["pairs"]])
tot = sum(r["n_cfu"] for r in R.values())
n_events = sum(sum(r["n_events_per_cfu"]) for r in R.values())

# per-threshold: count, and spread of the merged clusters (needs the maps once per run)
spread = {t: [] for t in THRS}
counts = {t: 0 for t in THRS}
per_animal = {a: {t: 0 for t in THRS} for a in ("2608-1", "2608-2", "2609-1")}
for run, r in R.items():
    animal = [a for a in per_animal if a in run][0]
    for t in THRS:
        cl = r["clusters"][str(t)]
        counts[t] += len(cl)
        per_animal[animal][t] += len(cl)
    if r["n_cfu"] < 2:
        continue
    maps = np.asarray(read_cfu(r["path"])["maps"]) > 0
    cents = [np.argwhere(m).mean(0) if m.any() else np.array([np.nan, np.nan]) for m in maps]
    for t in THRS:
        for g in r["clusters"][str(t)]:
            if len(g) < 2:
                continue
            c = np.array([cents[i] for i in g])
            dmax = max(np.linalg.norm(c[i] - c[j]) for i in range(len(g)) for j in range(i + 1, len(g)))
            spread[t].append(dmax * UM)

fig = plt.figure(figsize=(16, 4.6), facecolor=SURF)
gs = fig.add_gridspec(1, 4, wspace=0.36, left=0.045, right=0.99, top=0.80, bottom=0.17)

ax = fig.add_subplot(gs[0, 0])
edges = np.arange(0, 0.6001, 0.025)
ax.hist(vals, bins=edges, color=PAL[0], edgecolor=SURF, linewidth=0.6)
for t, c in ((0.5, TXT2), (0.3, PAL[1])):
    ax.axvline(t, color=c, lw=1.2, ls=(0, (3, 3)))
    ax.text(t + 0.006, ax.get_ylim()[1] * 0.95, f"{t}", color=c, fontsize=8.5, va="top")
ax.set_xlabel("AQuA2 overlap between two units (weighted Jaccard)", fontsize=8.5, color=TXT2)
ax.set_ylabel("pairs of units", fontsize=8.5, color=TXT2)
ax.set_title(f"A  All overlapping unit pairs, {len(R)} runs ({len(vals):,} pairs)\nno gap: the threshold is a choice, not a break in the data",
             fontsize=9.5, color=TXT, loc="left")

ax = fig.add_subplot(gs[0, 1])
xs = THRS
ax.plot(xs, [counts[t] / tot * 100 for t in xs], color=PAL[0], lw=2, marker="o", ms=5, label="all three animals")
for k, (a, col) in enumerate(zip(("2608-1", "2608-2", "2609-1"), (PAL[2], PAL[3], PAL[4]))):
    n0 = per_animal[a][0.5]
    ax.plot(xs, [per_animal[a][t] / n0 * 100 for t in xs], color=col, lw=1.2, marker="o", ms=3, label=a)
ax.axvline(0.3, color=PAL[1], lw=1.2, ls=(0, (3, 3)))
ax.set_xlim(0.52, 0.08)
ax.set_ylim(30, 105)
ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2)
ax.set_ylabel("units kept (% of the count at 0.5)", fontsize=8.5, color=TXT2)
ax.legend(fontsize=7.5, frameon=False, loc="lower left")
ax.set_title(f"B  Unit count vs threshold ({tot:,} units / {n_events:,} events at 0.5)\n0.3 → {counts[0.3]:,} units ({counts[0.3]/tot*100:.0f} %), 0.2 → {counts[0.2]:,} ({counts[0.2]/tot*100:.0f} %)",
             fontsize=9.5, color=TXT, loc="left")

ax = fig.add_subplot(gs[0, 2])
med = [np.median(spread[t]) if spread[t] else np.nan for t in xs]
p90 = [np.percentile(spread[t], 90) if spread[t] else np.nan for t in xs]
frac = [100 * np.mean(np.array(spread[t]) > 12) if spread[t] else np.nan for t in xs]
ax.plot(xs, med, color=PAL[0], lw=2, marker="o", ms=5, label="median")
ax.plot(xs, p90, color=PAL[0], lw=1.2, ls=(0, (3, 3)), marker="o", ms=3, label="90th percentile")
ax.set_xlim(0.52, 0.08)
ax.set_ylim(0, max(np.nanmax(p90) * 1.15, 10))
ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2)
ax.set_ylabel("largest centre-to-centre distance\nwithin a merged group (µm)", fontsize=8.5, color=TXT2)
ax.axvline(0.3, color=PAL[1], lw=1.2, ls=(0, (3, 3)))
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
ax.set_title(f"C  How far apart are the units that merge?\nat 0.3 median {med[xs.index(0.3)]:.1f} µm, at 0.2 {med[xs.index(0.2)]:.1f} µm", fontsize=9.5, color=TXT, loc="left")

ax = fig.add_subplot(gs[0, 3])
ax.plot(xs, frac, color=PAL[1], lw=2, marker="s", ms=4)
ax.axvline(0.3, color=PAL[1], lw=1.2, ls=(0, (3, 3)))
ax.set_xlim(0.52, 0.08)
ax.set_ylim(0, max(np.nanmax(frac) * 1.2, 10))
ax.set_xlabel("overlap threshold", fontsize=8.5, color=TXT2)
ax.set_ylabel("merged groups spanning > 12 µm (%)", fontsize=8.5, color=TXT2)
ax.set_title(f"D  Chains along fibres (> 12 µm ≈ 15 px)\nat 0.3: {frac[xs.index(0.3)]:.0f} % of groups, at 0.2: {frac[xs.index(0.2)]:.0f} %", fontsize=9.5, color=TXT, loc="left")
for a in fig.axes:
    a.set_facecolor(SURF)
    for s in a.spines.values():
        s.set_edgecolor(GRID)
    a.tick_params(colors=TXT2, labelsize=8)
    a.grid(axis="y", color=GRID, lw=0.6)
fig.suptitle(f"AQuA2 published CFU rule on the three RearChr animals ({len(R)} runs, 2608-1 / 2608-2 / 2609-1): what the overlap threshold does",
             fontsize=11.5, color=TXT, x=0.05, ha="left", y=0.97)
fig.savefig(OUT, dpi=150, facecolor=SURF)
print("saved", OUT)
print({t: (counts[t], round(med[i], 1), round(p90[i], 1), round(frac[i], 1)) for i, t in enumerate(xs)})
