#!/usr/bin/env python3
r"""aqua_reset — a fresh start for one animal (v35): park or delete the
pipeline's outputs so that every step runs again from the raw recordings with
the current version (GUI: the "Fresh start…" button next to "Run all").

    python aqua_reset.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --dry-run      # list what would go
    python aqua_reset.py "<mouse>"                                     # park everything in <mouse>\_previous_<date>\
    python aqua_reset.py "<mouse>" --stages aqua,analysis              # only the AQuA2 + CFU results and the analysis
    python aqua_reset.py "<mouse>" --discard-decisions                 # also the frame decisions and the CFU curation
    python aqua_reset.py "<mouse>" --delete                            # delete instead of parking (no way back)

Then: tab 1 … tab 8 (or "Run all (guided)"), or the commands of GUIDE.md.
Stages: unpack (binned movies + sidecars), tune (aqua_tune\ + aqua_opts\),
aqua (AQuA2 + CFU results + aqua_batch\), analysis (aqua_analysis\). Never
touched: the recordings, day_labels.csv, hand-saved *_Opt.csv, hand-drawn
*_fibre_mask.tif, reg_compare\, afferent2p_out\.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from afferent2p import reset as R   # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mouse_dir")
    ap.add_argument("--stages", default=",".join(R.STAGES), help="comma list of " + ", ".join(R.STAGES))
    ap.add_argument("--discard-decisions", action="store_true", help="also move the frame decisions and the CFU curation")
    ap.add_argument("--delete", action="store_true", help="delete instead of parking in _previous_<date>")
    ap.add_argument("--dry-run", action="store_true", help="only list what would go")
    a = ap.parse_args(argv)
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    bad = [s for s in stages if s not in R.STAGES]
    if bad:
        raise SystemExit(f"unknown stage(s) {bad}; choose from {R.STAGES}")
    m = Path(a.mouse_dir)
    if not m.is_dir():
        raise SystemExit(f"not a folder: {m}")
    R.reset_mouse(m, stages=stages, discard_decisions=a.discard_decisions, delete=a.delete, dry_run=a.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
