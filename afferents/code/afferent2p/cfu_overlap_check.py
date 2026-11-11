#!/usr/bin/env python3
"""cfu_overlap_check — does the CFU overlap threshold fit this animal? (v42)

    python cfu_overlap_check.py <mouse folder> [--thr 0.3] [--pos-radius 0.5] [--min-pairs 30] [--ok-ratio 0.9]

Reads the overlap tables the CFU batch writes next to every result
(<run>_AQuA2_cfu_pre.mat: AQuA2's own pairwise event overlaps + per-event
centre / size / frames), labels every overlapping pair without the overlap
(same place / two places, see afferent2p/cfu_overlap.py), clusters the events
of each recording as AQuA2 does at every threshold 0.10 … 0.50, and reports per
experiment (day x FOV) and pooled how well the units agree with the labels at
the chosen threshold compared with that experiment's own best. Writes
<mouse>\\aqua_batch\\cfu_overlap_check.csv / .txt / .png. Changes nothing.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p import cfu_overlap as C                                  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mouse", type=Path, help="mouse folder (with aqua_batch\\manifest.json)")
    ap.add_argument("--thr", type=float, default=None,
                    help="the threshold to judge (default: the one recorded in aqua_batch\\manifest.json, else 0.3)")
    ap.add_argument("--pos-radius", type=float, default=C.POS_RADIUS,
                    help=f"same place = centres within this x the smaller event's radius (default {C.POS_RADIUS})")
    ap.add_argument("--min-pairs", type=int, default=C.MIN_POS,
                    help=f"same-place pairs an experiment needs before it can judge (default {C.MIN_POS})")
    ap.add_argument("--ok-ratio", type=float, default=C.OK_RATIO,
                    help=f"OK when the chosen threshold reaches this share of the experiment's best F1 (default {C.OK_RATIO})")
    a = ap.parse_args()
    C.check_mouse(a.mouse, chosen=a.thr, pos_radius=a.pos_radius, min_pos=a.min_pairs, ok_ratio=a.ok_ratio)


if __name__ == "__main__":
    main()
