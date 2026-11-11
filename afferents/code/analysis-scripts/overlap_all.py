"""Pairwise AQuA2 CFU overlap (weighted Jaccard, as cfu.CFU_tmp_function) for every run of the three
RearChr animals, from the CFU weight maps of the current (overlap 0.5, min events 1) results.
For a one-event CFU the map equals the event's weighted footprint exactly."""
import glob, json, os, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, "/home/claude/afferent2p")
UP = "/mnt/user-data/uploads/GCAMP_affrerents"
SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
OUT = SP + "/overlap_all"
os.makedirs(OUT, exist_ok=True)
THRS = [0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]


def single_linkage(n, pairs, thr):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j, v in pairs:
        if v > thr:
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[rj] = ri
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def one(path):
    from afferent2p.aqua_io import read_cfu
    d = read_cfu(path)
    maps = np.asarray(d["maps"], dtype=np.float32)
    n = int(maps.shape[0]) if maps.ndim == 3 else 0
    if n == 0:
        return {"path": path, "run": os.path.basename(path).split("_bin31f")[0], "n_cfu": 0,
                "n_events_per_cfu": [], "n_px": [], "pairs": [], "clusters": {str(t): [] for t in THRS}}
    flat = maps.reshape(n, -1)
    nz = [np.flatnonzero(f) for f in flat]
    npx = [int(len(z)) for z in nz]
    nev = [int(np.size(e)) for e in d["events"]]
    # bounding boxes for a cheap pre-check
    H, W = maps.shape[1:]
    boxes = []
    for z in nz:
        if len(z) == 0:
            boxes.append((0, -1, 0, -1)); continue
        ys, xs = np.divmod(z, W)
        boxes.append((ys.min(), ys.max(), xs.min(), xs.max()))
    pairs = []
    for i in range(n):
        y0, y1, x0, x1 = boxes[i]
        for j in range(i + 1, n):
            b = boxes[j]
            if b[0] > y1 or b[1] < y0 or b[2] > x1 or b[3] < x0:
                continue
            inter = np.intersect1d(nz[i], nz[j], assume_unique=True)
            if inter.size == 0:
                continue
            idx = np.union1d(nz[i], nz[j])
            a = flat[i, idx]; c = flat[j, idx]
            v = float(np.minimum(a, c).sum() / np.maximum(a, c).sum())
            pairs.append((i, j, v))
    res = {"path": path, "run": os.path.basename(path).split("_bin31f")[0], "n_cfu": n,
           "n_events_per_cfu": nev, "n_px": npx, "pairs": pairs,
           "clusters": {str(t): [[int(i) for i in g] for g in single_linkage(n, pairs, t)] for t in THRS}}
    return res


if __name__ == "__main__":
    files = sorted(glob.glob(UP + "/NavRearChrTIGRE-26*/NavRearChrTIGRE-26*_*/*_AQuA2_res_cfu.mat"))
    files = [f for f in files if "260911_001" not in f]
    files.append(SP + "/cfu_old/260911_001/NavRearChrTIGRE-2608-1_260911_001_bin31f_1p98s_AQuA2_res_cfu.mat")
    print(len(files), "files")
    with ProcessPoolExecutor(max_workers=2) as ex:
        for res in ex.map(one, files):
            json.dump(res, open(OUT + "/" + res["run"] + ".json", "w"))
            nmulti = sum(1 for k in res["n_events_per_cfu"] if k > 1)
            print(res["run"], "CFUs", res["n_cfu"], "multi-event", nmulti, "pairs>0", len(res["pairs"]),
                  "n@0.3", len(res["clusters"]["0.3"]), "n@0.2", len(res["clusters"]["0.2"]), flush=True)
