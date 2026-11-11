"""Frame -> run/condition labels for the concatenated suite2p output.

suite2p registers and detects ROIs on ALL runs of an experiment concatenated
(pre + post, or the same FOV across days), so one ROI has one identity across
conditions.  Everything downstream needs to know which frame belongs to which
run.  This writes that mapping in three redundant, plain formats next to the
suite2p output (suite2p/plane0/):

  runs.json            list of runs with start/stop frame (0-based, stop
                       exclusive), condition, sbx path, nframes, fs, prealign
                       transform
  frame_labels.csv     one row per frame: frame, run_index, run, condition,
                       time_in_run_s
  afferent2p_runs.mat  the same for MATLAB (run names, conditions, start/stop
                       1-based inclusive, condition index per frame)

suite2p's own ops['frames_per_folder'] / ops['nframes'] carry the counts too,
but not the condition names.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

import numpy as np


def write_frame_labels(plane0: Path, run_names, run_conditions, run_nframes, fs,
                       sbx_paths=None, transforms=None) -> None:
    plane0 = Path(plane0)
    plane0.mkdir(parents=True, exist_ok=True)
    bounds = np.r_[0, np.cumsum(run_nframes)]
    runs = []
    for i, (n, c, f) in enumerate(zip(run_names, run_conditions, run_nframes)):
        r = {"index": i, "run": n, "condition": c, "start_frame": int(bounds[i]),
             "stop_frame": int(bounds[i + 1]), "nframes": int(f)}
        if sbx_paths is not None:
            r["sbx"] = str(sbx_paths[i])
        if transforms is not None and transforms[i]:
            t = transforms[i]
            r["prealign"] = {k: t.get(k) for k in ("angle_deg", "dy", "dx", "apply")}
        runs.append(r)
    conditions = list(dict.fromkeys(run_conditions))
    meta = {"fs": float(fs), "total_frames": int(bounds[-1]), "conditions": conditions,
            "note": "suite2p ran on all runs concatenated in this order; frames are "
                    "0-based, stop exclusive", "runs": runs}
    (plane0 / "runs.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    with open(plane0 / "frame_labels.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "run_index", "run", "condition", "time_in_run_s"])
        for i, (n, c, f) in enumerate(zip(run_names, run_conditions, run_nframes)):
            t = np.arange(f) / float(fs)
            w.writerows(zip(range(int(bounds[i]), int(bounds[i + 1])), [i] * f, [n] * f, [c] * f,
                            np.round(t, 4)))
    try:
        from scipy.io import savemat
        cond_idx = np.concatenate([np.full(f, conditions.index(c) + 1, np.int16)
                                   for c, f in zip(run_conditions, run_nframes)])
        run_idx = np.concatenate([np.full(f, i + 1, np.int16) for i, f in enumerate(run_nframes)])
        savemat(str(plane0 / "afferent2p_runs.mat"), {
            "fs": float(fs), "run_names": np.array(list(run_names), dtype=object),
            "run_conditions": np.array(list(run_conditions), dtype=object),
            "conditions": np.array(conditions, dtype=object),
            "run_start": (bounds[:-1] + 1).astype(np.int64), "run_stop": bounds[1:].astype(np.int64),
            "run_nframes": np.asarray(run_nframes, np.int64),
            "frame_run": run_idx, "frame_condition": cond_idx,
            "note": "1-based inclusive frame ranges; frame_condition indexes 'conditions'",
        }, do_compression=True)
    except Exception as e:  # pragma: no cover
        print(f"  (could not write afferent2p_runs.mat: {e})")
