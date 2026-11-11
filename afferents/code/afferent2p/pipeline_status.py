#!/usr/bin/env python3
r"""pipeline_status — where an animal is in the pipeline, from the files on disk.

    python pipeline_status.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" [more mouse folders]

One line per step (1 Unpack … 8 Report: done / partial / stale / pending / missing)
and the next thing to do — the same strip the GUI shows above its tabs.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.pipeline_state import mouse_state, format_state     # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+", help="mouse folders")
    ap.add_argument("--sec", type=float, default=2.0, help="which binned tiff (s per frame)")
    a = ap.parse_args()
    for f in a.folders:
        print(format_state(mouse_state(Path(f), sec=a.sec)))
        print()


if __name__ == "__main__":
    main()
