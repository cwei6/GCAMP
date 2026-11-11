"""Curation preview of NavRearChrTIGRE-2608-1 at the new CFU rule (overlap >= 0.3, min events 1):
one panel per run, drawn as the curation window draws it (mean image + one outline per unit, overlapping
units never share a colour; dotted = the window's 'looks non-fibre' suggestion), grouped by day x FOV.
Second figure: 260911_001 at 0.5 (before) vs 0.3 (now) with the formerly circled group as an inset."""
import sys, glob, json, re
sys.path.insert(0, "/home/claude/afferent2p")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import tifffile as tf
from afferent2p.aqua_io import read_cfu
from afferent2p.cfu_shape import cfu_features, fibre_mask_for
from cfu_curate import cfu_masks, mask_outline, assign_colours, PALETTE

SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
UP = "/mnt/user-data/uploads/GCAMP_affrerents/NavRearChrTIGRE-2608-1"
OUT = "/mnt/user-data/outputs/curation_preview"
import os; os.makedirs(OUT, exist_ok=True)
TXT, TXT2, SURF, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e6e5e1"
UM = 0.792

man = json.load(open(UP + "/aqua_batch/manifest.json"))
old_n = {json.load(open(f))["run"]: json.load(open(f))["n_cfu"] for f in glob.glob(SP + "/overlap_all/*.json")}


def load(run_dir, run):
    f = f"{UP}/{run_dir}/{run}_AQuA2_res_cfu.mat"
    d = read_cfu(f)
    mean = np.asarray(tf.imread(f"{UP}/{run_dir}/{run}_AQuA2_meanimg.tif"), np.float32)
    masks = cfu_masks(d["maps"])
    cols = [PALETTE[c] for c in assign_colours(masks)] if masks else []
    feats = []
    if masks:
        try:
            fmask, _ = fibre_mask_for(mean, UM)
            feats = cfu_features(masks, mean, UM, fmask=fmask)
        except Exception as e:
            print("  no suggestions for", run, e)
    n_ev = [int(np.size(e)) for e in d["events"]]
    return {"cfu": d, "mean": mean, "masks": masks, "cols": cols, "feats": feats, "n_events": n_ev, "rules": d.get("rules", "")}


def draw(ax, L, lw=0.8):
    img = L["mean"]
    lo, hi = np.percentile(img, [1, 99.7])
    ax.imshow(np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1), cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    for k, (m, c) in enumerate(zip(L["masks"], L["cols"])):
        sugg = k < len(L["feats"]) and L["feats"][k]["suggest"] != "fibre"
        for o in mask_outline(m):
            ax.plot(o[:, 0], o[:, 1], color=c, lw=lw, ls=(0, (1.2, 1.2)) if sugg else "-")
    ax.set_xlim(-0.5, img.shape[1] - 0.5); ax.set_ylim(img.shape[0] - 0.5, -0.5)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor(GRID)


# ------------------------------------------------------------------ figure 1: every run, grouped by day x FOV
rows = []
for e in man["experiments"]:
    day = re.search(r"_(\d{6})_FOV(\d)", e["name"])
    rows.append((f"{day.group(1)}  FOV{day.group(2)}", e))
H, W = 499, 782
pw, ph = 4.1, 4.1 * H / W + 0.55
fig = plt.figure(figsize=(0.55 + 4 * pw, 1.35 + len(rows) * ph), facecolor=SURF)
gs = fig.add_gridspec(len(rows), 4, left=0.03, right=0.995, top=1 - 0.95 / fig.get_figheight(), bottom=0.45 / fig.get_figheight(),
                      wspace=0.04, hspace=0.34)
tot_new = tot_old = tot_sugg = tot_ev = 0
per_exp = []
for i, (label, e) in enumerate(rows):
    nu = no = ns = ne = 0
    for j, r in enumerate(e["runs"]):
        run = r["run"]; run_dir = run.replace("_bin31f_1p98s", "")
        L = load(run_dir, run)
        ax = fig.add_subplot(gs[i, j]); draw(ax, L)
        n = L["cfu"]["n"]; o = old_n.get(run_dir, None); sugg = sum(1 for f_ in L["feats"] if f_["suggest"] != "fibre")
        nu += n; no += (o or 0); ns += sugg; ne += sum(L["n_events"])
        short = run_dir.split("_")[-2] + "_" + run_dir.split("_")[-1]
        ax.set_title(f"[{r['condition']}] {short}\n{n} units" + (f" (was {o} at 0.5)" if o is not None and o != n else (" (same at 0.5)" if o is not None else ""))
                     + f" · {sum(L['n_events'])} events" + (f" · {sugg} look non-fibre" if sugg else ""),
                     fontsize=8.6, color=TXT, loc="left", pad=3)
        if j == 0:
            ax.text(-0.012, 0.5, label, transform=ax.transAxes, rotation=90, ha="right", va="center", fontsize=10, color=TXT, fontweight="bold")
    per_exp.append((label, nu, no, ne, ns))
    tot_new += nu; tot_old += no; tot_sugg += ns; tot_ev += ne
fig.suptitle(f"NavRearChrTIGRE-2608-1 — what the curation window will show (v43): AQuA2 published CFU rule, overlap ≥ 0.3, min events 1\n"
             f"40 runs, 0 curated under these rules · {tot_new} units in all (was {tot_old} at overlap 0.5) · {tot_ev} events · "
             f"{tot_sugg} units look non-fibre by the window's shape rule",
             fontsize=13, color=TXT, x=0.03, ha="left", y=1 - 0.12 / fig.get_figheight(), va="top")
fig.text(0.03, 0.1 / fig.get_figheight(),
         "Each panel = one run's AQuA2 mean image with one outline per unit, coloured as in the curation window (units that overlap or touch never share a colour; "
         "the colour itself means nothing). Dotted outline = the window's 'looks' column calls the unit off-fibre or a blob — what 'Select suggested' (u) picks; "
         "a suggestion to look at, not a decision. 'was N at 0.5' = the same run's unit count under the previous threshold.",
         fontsize=9, color=TXT2, ha="left", va="bottom", wrap=True)
fig.savefig(OUT + "/NavRearChrTIGRE-2608-1_curation_preview_0p3.png", dpi=150, facecolor=SURF)
print("saved figure 1")
for row in per_exp:
    print("%-16s units %4d (was %4d)  events %4d  non-fibre %3d" % row)
print("total", tot_new, tot_old, tot_ev, tot_sugg)

# ------------------------------------------------------------------ figure 2: 260911_001 before / after, with the circled group
run_dir = "NavRearChrTIGRE-2608-1_260911_001"; run = run_dir + "_bin31f_1p98s"
new = load(run_dir, run)
oldd = read_cfu(glob.glob(SP + "/cfu_old/260911_001/*_res_cfu.mat")[0])
old = {"cfu": oldd, "mean": new["mean"], "masks": cfu_masks(oldd["maps"]), "feats": [], "n_events": [int(np.size(e)) for e in oldd["events"]]}
old["cols"] = [PALETTE[c] for c in assign_colours(old["masks"])]
# the group: 0.5-units 1,2,5,50,51,52 (1-based) became one 0.3-unit; its bounding box with a margin
grp = [0, 1, 4, 49, 50, 51]
ys, xs = np.nonzero(np.any([old["masks"][k] for k in grp], axis=0))
y0, y1, x0, x1 = ys.min() - 30, ys.max() + 30, xs.min() - 45, xs.max() + 45
new_unit = [k for k, m in enumerate(new["masks"]) if m[ys, xs].any()]
old_in_box = [k for k, m in enumerate(old["masks"]) if m[y0:y1, x0:x1].any()]
new_in_box = [k for k, m in enumerate(new["masks"]) if m[y0:y1, x0:x1].any()]
fig = plt.figure(figsize=(17, 13.2), facecolor=SURF)
gs = fig.add_gridspec(2, 2, left=0.03, right=0.99, top=0.90, bottom=0.05, wspace=0.05, hspace=0.22, height_ratios=[0.64, 0.80])
for j, (L, ttl) in enumerate(((old, f"BEFORE — overlap 0.5: {old['cfu']['n']} units, {sum(old['n_events'])} events"),
                              (new, f"NOW — overlap 0.3: {new['cfu']['n']} units, {sum(new['n_events'])} events"))):
    ax = fig.add_subplot(gs[0, j]); draw(ax, L, lw=0.9)
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="#2a78d6", lw=1.6))
    ax.set_title(ttl, fontsize=11.5, color=TXT, loc="left")
    ax2 = fig.add_subplot(gs[1, j]); draw(ax2, L, lw=1.8)
    ax2.set_xlim(x0, x1); ax2.set_ylim(y1, y0)
    inbox = old_in_box if j == 0 else new_in_box
    ev = [L["n_events"][k] for k in inbox]
    ax2.set_title(f"inside the blue box: {len(inbox)} units ({', '.join(str(v) for v in ev)} events each)\n"
                  + ("the same spot firing repeatedly came out as separate one-event units" if j == 0
                     else "the repeated events of that place are now one unit (6 events); the rest unchanged"),
                  fontsize=10.5, color=TXT, loc="left")
    for k in inbox:
        m = L["masks"][k]; yy, xx = np.nonzero(m)
        if y0 < yy.mean() < y1 and x0 < xx.mean() < x1:
            ax2.text(xx.mean(), yy.min() - 2, f"{k + 1}", color=L["cols"][k], fontsize=9, ha="center", va="bottom", clip_on=True)
fig.suptitle("NavRearChrTIGRE-2608-1  260911_001 (pre) — the run of the overlap question, before and after the threshold change\n"
             "same AQuA2 result and events; only the published CFU rule's overlap threshold changed (0.5 → 0.3, min events 1 in both)",
             fontsize=13, color=TXT, x=0.03, ha="left", y=0.975, va="top")
fig.text(0.03, 0.02, "Outlines as in the curation window (overlapping units never share a colour); numbers in the insets are the unit numbers of each result.",
         fontsize=9, color=TXT2)
fig.savefig(OUT + "/NavRearChrTIGRE-2608-1_260911_001_before_after.png", dpi=150, facecolor=SURF)
print("saved figure 2; box", (y0, y1, x0, x1), "old in box", [k + 1 for k in old_in_box], "new in box", [k + 1 for k in new_in_box],
      "new unit(s) at the group", [k + 1 for k in new_unit])
