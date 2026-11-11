#!/usr/bin/env python3
r"""aqua_report — the HTML reports of one animal (GUI tab 8).

    python aqua_report.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1"
        -> <mouse>\aqua_analysis\<mouse>_report.html   the MAIN report (fast: figures and tables of the analysis)
    python aqua_report.py <mouse> --qc
        -> <mouse>\aqua_analysis\<mouse>_qc.html       the QC report (reads the movies: detection maps, rasters,
                                                       clips, a montage of EVERY event, CFU table; minutes)
    python aqua_report.py <mouse> --all               both

Main report: QC paragraph, 1 general representation (sessions on the day's
clock, activity by running / resting state, ECDFs, acute / chronic summary),
2 heatmaps, 3 locomotion-evoked activity, 4 pupil, 5 coupling GLMs, 6 the
behaviour-controlled change with a verdict per experiment, 7 time course and
matched speed, 8 per experiment, 9 settings and tuning, 10 methods.

Needs the batch manifest (tab 4), the MATLAB results + per-event exports and
the analysis (tab 7); the movies for the QC montages.
Options: --max-events N (QC montages for the first N events per recording;
0 = all), --clips N (clips per recording, default 6; --no-clips), --out <file>.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.report import build_report          # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mouse", type=Path, help="mouse folder")
    ap.add_argument("--out", type=Path, default=None, help="output html (default aqua_analysis\\<mouse>_report.html / _qc.html)")
    ap.add_argument("--qc", action="store_true", help="build the QC report (per-recording evidence) instead of the main report")
    ap.add_argument("--all", action="store_true", help="build both reports")
    ap.add_argument("--max-events", type=int, default=0, help="QC: montages for at most N events per recording (0 = every event)")
    ap.add_argument("--clips", type=int, default=6, help="QC: animated clips per recording (default 6)")
    ap.add_argument("--no-clips", action="store_true", help="QC: no animated clips")
    ap.add_argument("--um-per-px", type=float, default=0.792, help="pixel size for the CFU 'looks' rule in the QC report")
    ap.add_argument("--no-qc", action="store_true", help="(kept for old scripts: the main report never reads the movies now)")
    a = ap.parse_args()
    if a.qc or a.all:
        from afferent2p.qc_report import build_qc_report
        out = build_qc_report(a.mouse, a.out if a.qc else None, max_events=a.max_events, clips_per_run=a.clips, no_clips=a.no_clips,
                              um_per_px=a.um_per_px)
        print(f"open {out}")
    if not a.qc:
        out = build_report(a.mouse, a.out, max_events=a.max_events, clips_per_run=a.clips, no_clips=a.no_clips, no_qc=a.no_qc)
        print(f"open {out}")


if __name__ == "__main__":
    main()
