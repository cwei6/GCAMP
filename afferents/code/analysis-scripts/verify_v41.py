"""After the v41 CFU re-run: compare the new (overlap 0.3) results with the prediction from the 0.5 maps.
usage: python3 verify_v41.py <staged root of re-staged _res_cfu.mat files>"""
import glob, json, os, sys
import numpy as np
sys.path.insert(0, "/home/claude/afferent2p")
from afferent2p.aqua_io import read_cfu, cfu_rules_of

SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
root = sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/uploads/GCAMP_affrerents"
pred = {}
for f in glob.glob(SP + "/overlap_all/*.json"):
    r = json.load(open(f)); pred[r["run"]] = (r["n_cfu"], len(r["clusters"]["0.3"]), sum(r["n_events_per_cfu"]))
rows = []
for f in sorted(glob.glob(root + "/NavRearChrTIGRE-26*/NavRearChrTIGRE-26*_*/*_AQuA2_res_cfu.mat")):
    run = os.path.basename(f).split("_bin31f")[0]
    rules = cfu_rules_of(f)
    d = read_cfu(f)
    n_ev = sum(int(np.size(e)) for e in d["events"])
    p = pred.get(run)
    rows.append((run, rules, d["n"], n_ev, p))
bad = 0
tot_new = tot_old = tot_pred = 0
for run, rules, n, n_ev, p in rows:
    flag = ""
    if "0.3" not in rules:
        flag = "  <-- NOT 0.3"; bad += 1
    elif p and (n != p[1] or n_ev != p[2]):
        flag = f"  (predicted {p[1]} units / {p[2]} events)"
    print(f"{run:36s} {rules:38s} units {n:4d} (was {p[0] if p else '?':>4}) events {n_ev:4d}{flag}")
    if p:
        tot_new += n; tot_old += p[0]; tot_pred += p[1]
print(f"\n{len(rows)} runs; units {tot_old} at 0.5 -> {tot_new} now (predicted {tot_pred}); {bad} run(s) not at 0.3")
