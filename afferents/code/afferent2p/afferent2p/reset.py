"""reset — a fresh start for one animal: park (or delete) the pipeline's
outputs so that every step runs again from the raw recordings with the
current version, instead of a mixture of results made by different versions.

What is an OUTPUT (reset by stage) and what is never touched:

  stage 'unpack'    in every run folder: the binned tiffs and their sidecars
                    (<run>_bin*s.tif / _ch2.tif, .json, _meanproj.tif,
                    _behavior.csv, _check_frames.tif, _check_movie.tif,
                    _flagged_frames.csv/.png, _frames.csv, _globalfix.png,
                    _raw.tif, the older dotted-name twins, _fullreg* / _legacy*
                    sidecars of the registration experiments, _motion.csv/.json)
  stage 'tune'      <mouse>\\aqua_tune\\ (the sweep: movies, jobs, reports) and
                    <mouse>\\aqua_opts\\ (the settings the sweep wrote)
  stage 'aqua'      in every run folder: AQuA2's results (<run>_..._AQuA2.mat,
                    _AQuA2_Ch1.csv, _AQuA2_Ch1_curves.xlsx, _AQuA2_events.csv,
                    _AQuA2_eventmap.tif, _AQuA2_meanimg.tif, _AQuA2_res_cfu.mat,
                    the _risingMaps_* folders) and <mouse>\\aqua_batch\\
  stage 'analysis'  <mouse>\\aqua_analysis\\ (tables, figures, both reports)

  never touched     the recordings (.sbx, the info .mat, _quadrature.mat, the
                    eye videos and DLC files), day_labels.csv, settings files
                    saved by hand (*_Opt.csv in the mouse folder or a run
                    folder), hand-drawn fibre masks (*_fibre_mask.tif),
                    reg_compare\\, afferent2p_out\\ (the old suite2p route),
                    anything else that is not listed above
  decisions         kept unless discard_decisions: the frame decisions of tab 2
                    (<run>_..._frames_decisions.json) and the CFU curation of
                    tab 6 (<run>_..._AQuA2_res_cfu_curated.json) — they are
                    yours, and they still apply to the re-made results

The files are MOVED into <mouse>\\_previous_<yymmdd_HHMM>\\ with their folder
structure (a rename on the same drive: instant, and reversible by moving them
back) — delete that folder when the new run is through and you are happy
with it; with delete=True they are deleted at once (no way back). A dry run
only lists what would go.
"""
from __future__ import annotations

import re
import shutil
import time
from pathlib import Path

STAGES = ("unpack", "tune", "aqua", "analysis")
STAGE_LABEL = {"unpack": "1 unpacked movies + sidecars", "tune": "3 tuning sweep + its settings",
               "aqua": "4-5 AQuA2 + CFU results + batch", "analysis": "7-8 analysis + reports"}
RUN_DIR_RE = re.compile(r"^(?P<mouse>.+?)_(?P<date>\d{6})_(?P<run>\d{3})$")
SKIP_DIRS = {"afferent2p_out", "aqua_batch", "aqua_analysis", "aqua_opts", "aqua_tune", "reg_compare",
             "qc", "qc_raw", "analysis", "analysis_manual", "behavior", "aqua_input", "aqua_analysis_manual", "suite2p"}

# per run folder, by stage: glob patterns relative to the run folder (a trailing / marks a folder pattern)
RUN_PATTERNS = {
    "unpack": ["*_bin*s.tif", "*_bin*s_ch2.tif", "*_bin*s.json", "*_bin*s_ch2.json", "*_bin*s_meanproj.tif", "*_bin*s_ch2_meanproj.tif",
               "*_bin*s_behavior.csv", "*_bin*s_check_frames.tif", "*_bin*s_check_movie.tif", "*_bin*s_flagged_frames.csv",
               "*_bin*s_flagged_frames.png", "*_bin*s_frames.csv", "*_bin*s_globalfix.png", "*_raw.tif", "*_bin*s_raw.tif",
               "*_bin*s_fullreg*", "*_bin*s_legacy*", "*_bin*s_motion.csv", "*_bin*s_motion.json", "*_bin*s_prereg*"],
    "aqua": ["*_AQuA2.mat", "*_AQuA2_Ch1.csv", "*_AQuA2_Ch1_curves.xlsx", "*_AQuA2_events.csv", "*_AQuA2_eventmap.tif",
             "*_AQuA2_meanimg.tif", "*_AQuA2_res_cfu.mat", "*_risingMaps_*/", "*_AQuA2_curves*", "*_AQuA2_overlay*",
             "*_bin*s_motion.csv", "*_bin*s_motion.json"],
}
DECISION_PATTERNS = {"unpack": ["*_frames_decisions.json"], "aqua": ["*_AQuA2_res_cfu_curated.json"]}
MOUSE_DIRS = {"tune": ["aqua_tune", "aqua_opts"], "aqua": ["aqua_batch"], "analysis": ["aqua_analysis"]}
KEEP_NAMES = {"day_labels.csv", "run_labels.csv"}


def run_folders(mouse_dir: Path) -> list[Path]:
    """The run folders (<mouse>_<yymmdd>_<nnn>) of the animal, pipeline folders excluded."""
    out = []
    for p in sorted(Path(mouse_dir).iterdir()):
        if p.is_dir() and RUN_DIR_RE.match(p.name) and p.name.lower() not in SKIP_DIRS:
            out.append(p)
    return out


def _is_kept(p: Path) -> bool:
    n = p.name
    return n in KEEP_NAMES or n.endswith("_Opt.csv") or n.endswith("_fibre_mask.tif")


def collect(mouse_dir: Path, stages=STAGES, discard_decisions: bool = False) -> dict:
    """{stage: [paths]} of what a reset of these stages would move, plus 'decisions'."""
    mouse_dir = Path(mouse_dir)
    found: dict = {s: [] for s in STAGES}
    found["decisions"] = []
    runs = run_folders(mouse_dir)
    for stage in stages:
        for rf in runs:
            for pat in RUN_PATTERNS.get(stage, []):
                folder = pat.endswith("/")
                for p in sorted(rf.glob(pat.rstrip("/"))):
                    if folder and not p.is_dir():
                        continue
                    if not folder and not p.is_file():
                        continue
                    if _is_kept(p) or p.name.endswith("_frames_decisions.json") or p.name.endswith("_res_cfu_curated.json"):
                        continue
                    if p not in found[stage]:
                        found[stage].append(p)
            if discard_decisions:
                for pat in DECISION_PATTERNS.get(stage, []):
                    for p in sorted(rf.glob(pat)):
                        if p.is_file() and p not in found["decisions"]:
                            found["decisions"].append(p)
        for name in MOUSE_DIRS.get(stage, []):
            d = mouse_dir / name
            if d.is_dir():
                found[stage].append(d)
    return found


def _size(p: Path) -> int:
    try:
        if p.is_file():
            return p.stat().st_size
        return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
    except OSError:
        return 0


def describe(found: dict, mouse_dir: Path) -> list[str]:
    lines = []
    total = 0
    for stage in list(STAGES) + ["decisions"]:
        items = found.get(stage, [])
        if not items:
            continue
        gb = sum(_size(p) for p in items) / 1e9
        total += gb
        label = STAGE_LABEL.get(stage, "your decisions (frame decisions, CFU curation)")
        n_files = sum(1 for p in items if p.is_file())
        n_dirs = sum(1 for p in items if p.is_dir())
        lines.append(f"  {label}: {n_files} file(s)" + (f" + {n_dirs} folder(s)" if n_dirs else "") + f", {gb:.1f} GB")
    lines.append(f"  total {total:.1f} GB")
    return lines


def reset_mouse(mouse_dir: Path, stages=STAGES, discard_decisions: bool = False, delete: bool = False,
                dry_run: bool = False, log=print) -> dict:
    """Park (or delete) the outputs of the given stages. Returns {park_dir, moved, deleted, bytes, found}."""
    mouse_dir = Path(mouse_dir)
    stages = [s for s in STAGES if s in set(stages)]
    found = collect(mouse_dir, stages, discard_decisions)
    items = []
    for p in [q for s in stages for q in found[s]] + list(found["decisions"]):
        if p not in items:                                   # (the motion files belong to two stages)
            items.append(p)
    log(f"fresh start for {mouse_dir.name}: stages {', '.join(stages)}" + (", decisions discarded too" if discard_decisions else ", decisions kept")
        + (" — DRY RUN, nothing is moved" if dry_run else (" — files are DELETED" if delete else " — files are parked")))
    for ln in describe(found, mouse_dir):
        log(ln)
    if not items:
        log("  nothing to reset")
        return {"park_dir": None, "moved": 0, "deleted": 0, "bytes": 0, "found": found}
    nbytes = sum(_size(p) for p in items)
    if dry_run:
        for p in items[:400]:
            log(f"    would {'delete' if delete else 'park'} {p.relative_to(mouse_dir)}")
        if len(items) > 400:
            log(f"    ... and {len(items) - 400} more")
        return {"park_dir": None, "moved": 0, "deleted": 0, "bytes": nbytes, "found": found}
    park = None
    moved = deleted = 0
    if not delete:
        park = mouse_dir / f"_previous_{time.strftime('%y%m%d_%H%M')}"
        park.mkdir(exist_ok=True)
    for p in items:
        rel = p.relative_to(mouse_dir)
        try:
            if delete:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
                deleted += 1
            else:
                dst = park / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.exists():
                    dst = dst.with_name(dst.name + f".{int(time.time())}")
                try:
                    p.rename(dst)
                except OSError:
                    shutil.move(str(p), str(dst))
                moved += 1
        except Exception as e:                                           # noqa: BLE001
            log(f"  could not {'delete' if delete else 'move'} {rel}: {e}")
    if delete:
        log(f"  deleted {deleted} item(s), {nbytes / 1e9:.1f} GB")
    else:
        log(f"  parked {moved} item(s), {nbytes / 1e9:.1f} GB, in {park}\\ — delete that folder when the new run is through "
            "(or move things back from it if you change your mind)")
    return {"park_dir": park, "moved": moved, "deleted": deleted, "bytes": nbytes, "found": found}


def estimate_hours(mouse_dir: Path, stages=STAGES) -> dict:
    """A rough idea of the time the re-run costs, from the run count."""
    n_runs = len(run_folders(Path(mouse_dir)))
    days = {RUN_DIR_RE.match(p.name)["date"] for p in run_folders(Path(mouse_dir))}
    est = {"n_runs": n_runs, "n_days": len(days)}
    if "unpack" in stages:
        est["unpack_h"] = n_runs * 3 / 60
    if "tune" in stages:
        est["tune_h"] = len(days) * 4 * 48 * 1.4 / 60          # 4 movies per day (pre + post x 2 FOVs), 16 settings x 3 variants
    if "aqua" in stages:
        est["aqua_h"] = n_runs * 6 / 60
    if "analysis" in stages:
        est["analysis_h"] = 0.5
    return est
