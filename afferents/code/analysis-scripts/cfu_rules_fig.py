"""Compare the three AQuA2 CFU rule sets on 260911_001 (figure for the user)."""
import sys
import numpy as np
import tifffile
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from skimage import measure

sys.path.insert(0, "/home/claude/afferent2p")
from afferent2p.aqua_io import read_cfu

SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
UP = "/mnt/user-data/uploads/GCAMP_affrerents/NavRearChrTIGRE-2608-1/NavRearChrTIGRE-2608-1_260911_001/"
STEM = "NavRearChrTIGRE-2608-1_260911_001_bin31f_1p98s_AQuA2"
OUT = "/mnt/user-data/outputs/cfu_rules/260911_001_cfu_rules_compared.png"

PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TXT, TXT2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"

old = read_cfu(SP + "/cfu_old/260911_001/" + STEM + "_res_cfu.mat")
new = read_cfu(UP + STEM + "_res_cfu.mat")
old_maps = np.asarray(old["maps"], dtype=np.float32)          # one event per CFU, index = event id - 1
new_maps = np.asarray(new["maps"], dtype=np.float32)
new_events = [list(np.atleast_1d(e).astype(int)) for e in new["events"]]
ov = np.load(SP + "/cfu_old/260911_001/overlap_matrix.npy")
n_ev = old_maps.shape[0]

mean_img = tifffile.imread(UP + STEM + "_meanimg.tif").astype(float)
lo, hi = np.percentile(mean_img, (1, 99.7))
mean_img = np.clip((mean_img - lo) / (hi - lo), 0, 1)


def single_linkage(thr):
    parent = list(range(n_ev))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n_ev):
        for j in range(i + 1, n_ev):
            if ov[i, j] > thr:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    groups = {}
    for i in range(n_ev):
        groups.setdefault(find(i), []).append(i + 1)
    return sorted(groups.values(), key=len, reverse=True)


groups03 = single_linkage(0.3)
n03 = len(groups03)
n_shared = sum(1 for e in range(1, n_ev + 1) if sum(e in ev for ev in new_events) > 1)
max_split = max(sum(e in ev for ev in new_events) for e in range(1, n_ev + 1))


def show_bg(ax, y0, y1, x0, x1):
    ax.imshow(mean_img[y0:y1, x0:x1], cmap="gray", vmin=0, vmax=1, extent=(x0, x1, y1, y0), interpolation="nearest")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor("#c9c8c3")


def outline(ax, mask, color, lw=1.4, label=None, fill=False, alpha=0.35):
    mask = mask > 0
    if fill:
        rgba = np.zeros(mask.shape + (4,))
        rgb = matplotlib.colors.to_rgb(color)
        rgba[mask] = (*rgb, alpha)
        ax.imshow(rgba, extent=(0, mask.shape[1], mask.shape[0], 0), interpolation="nearest")
    for c in measure.find_contours(mask.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color=color, lw=lw)
    if label is not None:
        ys, xs = np.nonzero(mask)
        ax.text(xs.mean(), ys.mean(), label, color=TXT, fontsize=7.5, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.85))


def point(ax, mask, text, xy_text):
    ys, xs = np.nonzero(mask > 0)
    ax.annotate(text, xy=(xs.mean(), ys.mean()), xytext=xy_text, fontsize=7.5, color=TXT, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none", alpha=0.9),
                arrowprops=dict(arrowstyle="-", color=TXT2, lw=0.8))


circ = [1, 2, 5, 22, 50, 51, 52]
Y0, Y1, X0, X1 = 222, 266, 436, 476          # the circled spot
fig = plt.figure(figsize=(13, 9.2), facecolor=SURF)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.32, wspace=0.12, left=0.03, right=0.99, top=0.86, bottom=0.05)

# --- row 1: the circled spot under the three rule sets -------------------------------------
ax = fig.add_subplot(gs[0, 0])
show_bg(ax, Y0, Y1, X0, X1)
for k, e in enumerate(circ):
    outline(ax, old_maps[e - 1], PAL[k % 8])
ax.set_title("A  Published rule, overlap 0.5 (what you curated)\n7 one-event CFUs: events 1, 2, 5, 22, 50, 51, 52", fontsize=9.5, color=TXT, loc="left")

ax = fig.add_subplot(gs[0, 1])
show_bg(ax, Y0, Y1, X0, X1)
g6 = [g for g in groups03 if 1 in g][0]
union6 = np.max(old_maps[[e - 1 for e in g6]], axis=0)
outline(ax, union6, PAL[0], fill=True)
outline(ax, old_maps[21], PAL[1], fill=True)
point(ax, union6, "1 CFU: events 1, 2, 5, 50, 51, 52", (440, 258))
point(ax, old_maps[21], "1 CFU: event 22", (462, 230))
ax.set_title("B  Published rule, overlap setting 0.3\n2 CFUs: {1,2,5,50,51,52} + {22}", fontsize=9.5, color=TXT, loc="left")

ax = fig.add_subplot(gs[0, 2])
show_bg(ax, Y0, Y1, X0, X1)
outline(ax, new_maps[0], PAL[0], fill=True)
outline(ax, new_maps[75], PAL[1], fill=True)
outline(ax, new_maps[76], PAL[3], fill=True)
point(ax, new_maps[0], "CFU 1: events 1, 2, 5, 50, 51, 52", (440, 258))
point(ax, new_maps[75], "CFU 76: event 22", (440, 230))
point(ax, new_maps[76], "CFU 77: event 22", (462, 230))
ax.set_title("C  Current release (Oct 2026), defaults\n3 CFUs: {1,2,5,50,51,52} + event 22 split in 2 basins", fontsize=9.5, color=TXT, loc="left")

# --- row 2: the side effect of the current release: single events split into basins ----------
Y0b, Y1b, X0b, X1b = 24, 116, 218, 312        # events 7 and 11 (two neighbouring single events)
ax = fig.add_subplot(gs[1, 0])
show_bg(ax, Y0b, Y1b, X0b, X1b)
outline(ax, old_maps[6], PAL[0], fill=True, label="event 7")
outline(ax, old_maps[10], PAL[1], fill=True, label="event 11")
ax.set_title("D  Published rule (0.5 or 0.3): events 7 and 11\n2 CFUs (one event each)", fontsize=9.5, color=TXT, loc="left")

ax = fig.add_subplot(gs[1, 1])
show_bg(ax, Y0b, Y1b, X0b, X1b)
k = 0
for i, ev in enumerate(new_events):
    if ev in ([7], [11]):
        outline(ax, new_maps[i], PAL[k % 8], fill=True, label=str(i + 1), lw=1.0)
        k += 1
ax.set_title(f"E  Current release: the same two events\n{k} CFUs (spatial basins, each holding that one event)", fontsize=9.5, color=TXT, loc="left")

# --- counts panel ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 2])
ax.axis("off")
counts = [("Published rule\noverlap 0.5", 53, "every CFU = one event"),
          ("Published rule\noverlap 0.3", n03, "2 multi-event CFUs\n(6 + 2 events); no event split"),
          ("Current release\ndefaults", len(new_events), f"{n_shared} of 53 events sit in\n2–{max_split} CFUs each")]
ax.set_title("F  Whole recording: 53 events →", fontsize=9.5, color=TXT, loc="left")
ax.set_xlim(-0.6, 2.6)
ax.set_ylim(-75, 175)
for i, (name, n, note) in enumerate(counts):
    ax.bar(i, n, width=0.62, color=PAL[i], edgecolor=SURF, linewidth=1)
    ax.text(i, n + 4, f"{n} CFUs", ha="center", va="bottom", fontsize=9.5, color=TXT)
    ax.text(i, -6, name, ha="center", va="top", fontsize=8.2, color=TXT)
    ax.text(i, -34, note, ha="center", va="top", fontsize=7.2, color=TXT2)
ax.axhline(53, color=TXT2, lw=0.8, ls=(0, (3, 3)))
ax.text(2.55, 55, "53 events", ha="right", va="bottom", fontsize=7.5, color=TXT2)
fig.suptitle("NavRearChrTIGRE-2608-1  260911_001 — the same 53 AQuA2 events, three AQuA2 CFU rule sets",
             fontsize=12, color=TXT, x=0.03, ha="left", y=0.985, va="top")
fig.text(0.03, 0.955, "Top row: the blue-circled spot.  Bottom row: what the current release's basin splitting does to a single event.  "
         "Published rule = single-linkage clustering of events by weighted spatial overlap (the AQuA2 paper); "
         "current release adds watershed splitting of each cluster into basins and a post-merge of correlated neighbours.",
         fontsize=8.2, color=TXT2, wrap=True, va="top")
import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150, facecolor=SURF)
print("saved", OUT, "| groups at 0.3:", [g for g in groups03 if len(g) > 1], "| n03", n03, "| shared", n_shared, "max_split", max_split)
