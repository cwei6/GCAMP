#!/usr/bin/env python3
"""plan — group a mouse folder's recordings into FOV-level experiments and run them.

    python plan.py <mouse folder> [<mouse folder> ...]            # show the plan only
    python plan.py <mouse folder> --run                            # write configs + run
    python plan.py <mouse folder> --paradigm optostim --n-fov 2 --run
    python plan.py <mouse folder> --paradigm spontaneous --n-fov 3 --run
    python plan.py <mouse folder> --stim-after-run 4 --run         # optostim: runs 1-4 pre
    python plan.py <mouse folder> --run-labels                     # v44: create / complete <mouse>\run_labels.csv

Paradigm is guessed from the animal name ('RearChr' -> optostim, otherwise
spontaneous); see afferent2p/plan.py for the grouping rules. Settings other
than the experiment block come from --template (default: config_example.yaml
next to this script) or an existing config yaml.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from afferent2p.plan import plan_mouse, describe, experiment_config, write_run_labels  # noqa: E402


def load_template(path: Path | None) -> dict:
    p = path or (HERE / "config_example.yaml")
    cfg = yaml.safe_load(Path(p).read_text())
    cfg.pop("experiment", None)
    return cfg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+", help="mouse folder(s)")
    ap.add_argument("--paradigm", default="auto", choices=["auto", "optostim", "spontaneous"])
    ap.add_argument("--n-fov", type=int, default=None, help="FOVs per day (optostim default 2; spontaneous = max runs/day)")
    ap.add_argument("--stim-after-run", type=int, default=None, help="optostim: last PRE run number (default: half the runs)")
    ap.add_argument("--run-labels", action="store_true",
                    help="create <mouse>\\run_labels.csv with every run and the rule's FOV / condition, or add the runs an "
                         "existing file does not list (edit fov / condition / use / note, then look at the plan again)")
    ap.add_argument("--template", type=Path, default=None, help="yaml whose non-experiment settings to use")
    ap.add_argument("--steps", default=None, help="pipeline steps (default: all)")
    ap.add_argument("--run", action="store_true", help="write configs and run the pipeline for each experiment")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    template = load_template(a.template)
    all_exps = []
    for f in a.folders:
        if a.run_labels:
            write_run_labels(Path(f), paradigm=a.paradigm, n_fov=a.n_fov, stim_after_run=a.stim_after_run)
        exps = plan_mouse(Path(f), paradigm=a.paradigm, n_fov=a.n_fov, stim_after_run=a.stim_after_run)
        if not exps:
            print(f"{f}: no recordings recognised"); continue
        print(f"\n{f}  ->  {len(exps)} experiment(s)\n" + describe(exps))
        all_exps += exps
    if not a.run:
        print("\n(dry run — add --run to write configs and process)")
        return
    for e in all_exps:
        cfg = experiment_config(e, template)
        e.output_dir.mkdir(parents=True, exist_ok=True)
        cfg_path = e.output_dir / "afferent2p_config.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        cmd = [sys.executable, str(HERE / "run_pipeline.py"), str(cfg_path)]
        if a.steps:
            cmd += ["--steps", a.steps]
        if a.overwrite:
            cmd.append("--overwrite")
        print(f"\n===== {e.name} =====\nconfig: {cfg_path}")
        rc = subprocess.call(cmd)
        print(f"===== {e.name}: {'done' if rc == 0 else 'FAILED'} =====")


if __name__ == "__main__":
    main()
