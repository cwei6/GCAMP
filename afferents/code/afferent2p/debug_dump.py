#!/usr/bin/env python3
"""Write a small, shareable dump of one experiment's ROI traces for tuning the
afferent filter offline.

    python debug_dump.py <experiment output folder> [--n 150] [--out file.npz]

Picks the N ROIs with the highest suite2p 'snr' (lag-1 autocorrelation), all
currently kept ROIs, and 50 random others; saves their neuropil-corrected
dF/F quantized to int8 (clipped to [-1, 5] dF/F, 0.025 steps; ~8 MB for 200
ROIs x 60k frames), the curation features, ROI outlines, the green and red
mean images and run boundaries. Nothing else from the recording.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.analyze import load_rois, dff_and_z          # noqa: E402
from afferent2p.curate import roi_features, curate_rois      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", type=Path, help="experiment output folder (contains suite2p/, manifest.json)")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    plane0 = a.out_dir / "suite2p" / "plane0"
    man = json.loads((a.out_dir / "manifest.json").read_text())
    fs = float(man["fs"])
    cfg = yaml.safe_load((a.out_dir / "config_used.yaml").read_text()) if (a.out_dir / "config_used.yaml").exists() else {}
    an = cfg.get("analysis", {})
    F, Fneu, spks, stat, iscell, names = load_rois(plane0)
    n = F.shape[0]
    _, _, dff, _ = dff_and_z(F - float(an.get("neuropil_coeff", 0.7)) * Fneu, fs,
                             float(an.get("baseline_win_s", 32.0)), float(an.get("baseline_percentile", 10)))
    ops = np.load(plane0 / "ops.npy", allow_pickle=True).item()
    red = ops.get("meanImg_chan2")
    feats, _ = roi_features(dff, fs, stat=stat, F=F, Fneu=Fneu,
                            red_img=None if red is None else np.asarray(red, np.float32))
    keep, score, crit, p = curate_rois(feats, an)
    rng = np.random.default_rng(0)
    s2p_snr = np.array([float(s.get("snr", 0)) for s in stat])
    sel = set(np.argsort(-s2p_snr)[: a.n].tolist()) | set(np.flatnonzero(keep).tolist()) \
        | set(rng.choice(n, min(50, n), replace=False).tolist())
    sel = np.array(sorted(sel))
    out = a.out or (a.out_dir / "curation_debug.npz")
    q = np.clip(np.round((dff[sel] + 1.0) / 0.025), 0, 240).astype(np.uint8)   # dff = q*0.025 - 1
    np.savez_compressed(
        out, roi_index=sel, dff_q=q, dff_scale=0.025, dff_offset=-1.0, fs=fs,
        run_nframes=np.array([r["nframes"] for r in man["runs"]]),
        run_conditions=np.array([r["condition"] for r in man["runs"]]),
        keep=keep[sel], score=score[sel],
        features=json.dumps({k: np.asarray(v)[sel].tolist() for k, v in feats.items()}),
        ypix=np.array([stat[i]["ypix"] for i in sel], dtype=object),
        xpix=np.array([stat[i]["xpix"] for i in sel], dtype=object),
        meanImg=ops.get("meanImg"), max_proj=ops.get("max_proj"),
        meanImg_chan2=ops.get("meanImg_chan2"),
        Ly=int(ops["Ly"]), Lx=int(ops["Lx"]), n_total=n,
        thresholds=json.dumps(p),
    )
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB): {len(sel)} of {n} ROIs, "
          f"{int(keep.sum())} currently kept")


if __name__ == "__main__":
    main()
