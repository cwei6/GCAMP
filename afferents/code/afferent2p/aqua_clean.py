#!/usr/bin/env python3
"""aqua_clean — the brightness step on binned tiffs: flatten whole-field
brightness, flag bad frames, check them frame by frame, decide what the
statistics do with them. New unpacks get it automatically from unpack_bin.py;
run this for tiffs unpacked earlier, to re-flag with the current rules
(--force), to preview (--report-only), to change decisions, or to rebuild the
ImageJ check stacks (--stacks-only).

    python aqua_clean.py <mouse folder> [<mouse folder> ...]   # every run's channel-1 tiff
    python aqua_clean.py <one binned tiff>
    python aqua_clean.py <mouse folder> --report-only          # flags + QC + check files, movie untouched
    python aqua_clean.py <mouse folder> --force                # re-flag corrected tiffs with the current rules
    python aqua_clean.py <one binned tiff> --exclude 63-87 --keep 13-19          # decisions per stretch
    python aqua_clean.py <one binned tiff> --exclude-all 63-87                   # the whole stretch
    python aqua_clean.py <one binned tiff> --exclude-frames "63-87:68,76-77"     # exactly these frames
    python aqua_clean.py <mouse folder> --stacks-only          # (re)write the two check stacks only
    python aqua_clean.py <mouse folder> --drop --force         # leave excluded frames out of the movies
    python aqua_clean.py <mouse folder> --restore-raw          # put <run>_raw.tif back, then re-run

Per tiff (afferent2p/globalfix.py):
  * the level (spatial median per frame) is flattened by a per-frame gain;
    dF/F is scale-free, so events are unchanged, bleaching and steps vanish
  * two classes of flags:
      artefact — fast changes: a jump between consecutive frames beyond
                 max(4 robust SD, 3 %), a spike (<= 10 s, returning), or the
                 content correlation collapsing (flash, shutter, brief shift);
                 tails while the level is still moving, +1 frame each side
      state    — slow changes: a level plateau longer than a spike or not
                 returning, or the p10/median ratio (which a pure gain change
                 cannot move) more than 5 % away from its clean-frame
                 baseline; a laser/focus change or widespread activity — the
                 trace cannot tell which
  * the frame-by-frame check: every frame of a flagged stretch is compared
    with the kept frames around it. The FLAGGED FRAMES of an artefact stretch
    are the frames that differ from those neighbours by >= 10 % of the frame
    and >= 3 x the clean floor (--core-min-pct / --core-ratio) — the flash
    itself, not the tails and margins around it.
  * decisions per stretch: keep | exclude (= its flagged frames only) |
    exclude-all (the whole stretch) | frames:68,76-77 (exactly these).
    Excluded frames: events starting in them (+/- a margin) leave the
    statistics, their time leaves the denominator. Defaults: artefact with
    flagged frames -> exclude; artefact without -> keep; state -> keep. Your
    own decisions go in <base>_frames_decisions.json (the flags above, or the
    GUI) and survive re-runs.
  * the frames STAY IN THE MOVIE (default). --drop leaves the excluded frames
    out (the csv maps movie frames to recorded frames); --fill replaces them
    by a frozen filler. Both are explicit choices, never the default.
  * files next to the tiff: <base>_frames.csv (one row per recorded frame:
    level, gain, class, stretch, decision, flagged_frame, excluded,
    tiff_frame, reason, the per-frame check), <base>_globalfix.png (level and
    shape traces; red = artefact, amber = state, black bar = excluded), the
    check files — <base>_check_frames.tif (16-bit ImageJ stack: every flagged
    stretch with 6 frames of context on each side, each slice labelled with
    its recorded frame number, stretch, numbers, flagged / EXCLUDED — scroll
    through it in ImageJ), <base>_check_movie.tif (the same slices as a
    composite: frame, % difference from the neighbours, traces with the
    cursor), <base>_flagged_frames.png/.csv (overview sheet + table) — and a
    "global_correction" block in the .json. A corrected tiff is skipped unless
    --force.

If AQuA2 has already been run on a tiff, decisions are applied by
aqua_analyze.py without re-running AQuA2; only --drop/--fill (or --force after
--drop) change the movie and make the AQuA2 result stale.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

try:
    import tifffile  # noqa: F401
    import scipy  # noqa: F401
    import matplotlib  # noqa: F401
except ImportError as e:
    raise SystemExit(f"\nMissing Python package: {e.name}\n"
                     f'Install it: "{sys.executable}" -m pip install tifffile scipy matplotlib\n')

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.aqua_io import find_binned_tifs, aqua_files      # noqa: E402
from afferent2p.globalfix import apply, apply_decisions, restore_raw, make_check_stacks  # noqa: E402
from afferent2p.plan import skip_parts                           # noqa: E402


def collect(targets) -> list[Path]:
    out = []
    for t in targets:
        t = Path(t)
        if t.is_file():
            out.append(t)
            continue
        for d in sorted({p.parent for p in t.rglob("*.sbx")} | {p.parent for p in t.rglob("*_bin*s.tif")} | {t}):
            if skip_parts(d.relative_to(t).parts):
                continue
            out += find_binned_tifs(d)
    seen, uniq = set(), []
    for p in out:
        if p not in seen:
            seen.add(p); uniq.append(p)
    return uniq


def parse_ids(txt: str | None) -> list[str]:
    if not txt:
        return []
    return [x.strip() for x in txt.replace(";", ",").split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets", nargs="+", help="mouse folders, run folders or binned tiffs")
    ap.add_argument("--report-only", action="store_true",
                    help="flags, frames.csv, QC plot and the check sheet only; the movie is not touched")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--drop", action="store_true", help="leave the excluded frames out of the movie (explicit choice)")
    g.add_argument("--fill", action="store_true", help="replace the excluded frames by a frozen filler")
    ap.add_argument("--exclude", default=None,
                    help="stretch ids whose FLAGGED frames are excluded, e.g. 63-87,94-102 (decisions, one tiff)")
    ap.add_argument("--exclude-all", default=None, help="stretch ids excluded whole, e.g. 63-87 (one tiff)")
    ap.add_argument("--exclude-frames", default=None,
                    help='exactly these frames of a stretch: "63-87:68,76-77;94-102:99" (one tiff)')
    ap.add_argument("--keep", default=None, help="stretch ids to keep entirely, e.g. 13-19 (decisions, one tiff)")
    ap.add_argument("--decisions-only", action="store_true",
                    help="only apply the decisions to the existing flags (no re-flagging, movie untouched)")
    ap.add_argument("--core-ratio", type=float, default=3.0,
                    help="a frame is FLAGGED when its difference from the clean neighbours is at least this many "
                         "times the clean floor (default 3)")
    ap.add_argument("--core-min-pct", type=float, default=10.0,
                    help="... and at least this %% of the frame (smoothed RMS; default 10)")
    ap.add_argument("--no-stacks", action="store_true", help="do not write the two ImageJ check stacks")
    ap.add_argument("--stacks-only", action="store_true",
                    help="(re)write the check stacks from the existing flags and decisions; nothing else")
    ap.add_argument("--stack-context", type=int, default=6, help="frames of context on each side of a stretch (default 6)")
    ap.add_argument("--stack-ds", type=int, default=1,
                    help="block-average the 16-bit stack by this factor (default 1 = full resolution)")
    ap.add_argument("--restore-raw", action="store_true",
                    help="put <base>_raw.tif back as the movie and forget the correction, then run the step again")
    ap.add_argument("--no-normalise", action="store_true", help="do not flatten the global level")
    ap.add_argument("--k", type=float, default=4.0, help="flag threshold in robust SDs (default 4)")
    ap.add_argument("--min-rel", type=float, default=0.03, help="level floor for jumps/spikes (default 0.03 = 3%%)")
    ap.add_argument("--shape-rel", type=float, default=0.05,
                    help="floor for the p10/median ratio (state class; default 0.05 = 5%%)")
    ap.add_argument("--shape-win", type=float, default=300.0,
                    help="window (s) of the clean-frame moving median behind the shape rule (default 300)")
    ap.add_argument("--max-spike", type=float, default=10.0,
                    help="an excursion up to this long (s) that returns counts as a spike = artefact (default 10)")
    ap.add_argument("--min-clean", type=float, default=10.0,
                    help="clean stretches shorter than this (s) between two artefact stretches are flagged too")
    ap.add_argument("--win", type=float, default=60.0, help="moving-median window in seconds (default 60)")
    ap.add_argument("--force", action="store_true",
                    help="re-flag an already corrected tiff with the current rules (from its recorded raw levels)")
    a = ap.parse_args()
    mode = "drop" if a.drop else "fill" if a.fill else "keep"
    decisions = {**{i: "exclude" for i in parse_ids(a.exclude)}, **{i: "exclude-all" for i in parse_ids(a.exclude_all)},
                 **{i: "keep" for i in parse_ids(a.keep)}}
    if a.exclude_frames:
        for part in a.exclude_frames.split(";"):
            if not part.strip():
                continue
            if ":" not in part:
                raise SystemExit(f'--exclude-frames wants "<stretch id>:<frames>" pairs, e.g. "63-87:68,76-77" (got {part!r})')
            sid, frames = part.split(":", 1)
            decisions[sid.strip()] = "frames:" + frames.strip()
    tifs = collect(a.targets)
    if not tifs:
        raise SystemExit("no binned tiffs found (run unpack_bin.py first)")
    if decisions and len(tifs) > 1:
        raise SystemExit("decisions name stretches of ONE tiff: give that tiff (or its run folder) as the target")
    print(f"{len(tifs)} tiff(s)")
    if a.decisions_only:
        for t in tifs:
            apply_decisions(t, decisions)
        return
    if a.stacks_only:
        for t in tifs:
            if decisions:
                apply_decisions(t, decisions)
            make_check_stacks(t, context=a.stack_context, ds=a.stack_ds)
        return
    n_done = 0
    for t in tifs:
        if a.restore_raw:
            if not restore_raw(t):
                continue
        rec = apply(t, win_s=a.win, k=a.k, min_rel=a.min_rel, shape_rel=a.shape_rel, shape_win_s=a.shape_win,
                    min_clean_s=a.min_clean, max_spike_s=a.max_spike, mode=mode, normalise=not a.no_normalise,
                    force=a.force or a.restore_raw, decisions=decisions or None, report_only=a.report_only,
                    core_ratio=a.core_ratio, core_min_pct=a.core_min_pct, stacks=not a.no_stacks,
                    stack_context=a.stack_context, stack_ds=a.stack_ds)
        if rec and not a.report_only:
            n_done += 1
            if mode != "keep":
                stale = [p for k, p in aqua_files(t).items() if k in ("res", "features", "cfu") and p.exists()]
                if stale:
                    print(f"    note: AQuA2 results exist for this tiff and were computed on the PREVIOUS movie "
                          f"— rebuild the batch with --overwrite and re-run MATLAB for it")
    if not a.report_only:
        print(f"done: {n_done} tiff(s)")


if __name__ == "__main__":
    main()
