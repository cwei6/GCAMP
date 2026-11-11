#!/usr/bin/env python3
"""aqua_analyze — read AQuA2 batch results and compare event properties across
conditions (days for spontaneous mice, pre vs post for RearChr mice).

    python aqua_analyze.py <mouse folder> [<mouse folder> ...] [options]

Reads, for every run of every FOV experiment planned by aqua_batch.py:
    <run>_bin31f_1p98s_AQuA2_Ch1.csv     event feature table (AQuA2)
    <run>_bin31f_1p98s.json              s per frame, frames, crop (unpack_bin)
    <run>_bin31f_1p98s_behavior.csv      velocity / still fraction / pupil per frame
    <run>_bin31f_1p98s_AQuA2_res_cfu.mat CFUs (optional; from run_cfu_batch.m)

Writes into <mouse>\\aqua_analysis\\:
    events_<experiment>.csv       every event, with run / condition / time / state
    run_summary.csv               one row per run   (the unit you can average)
    condition_summary.csv         one row per experiment x condition
    cfu_matched_<experiment>.csv  CFUs matched across runs by map overlap (if CFU files exist)
    <experiment>_overview.png     rate + property distributions per condition
    <experiment>_cfu_map.png      matched CFU maps (if CFU files exist)
    summary.txt                   the numbers in words, with the statistics caveats
  v32 — the general representation (afferent2p.general), the locomotion-evoked
  activity (afferent2p.locotrials) and the heatmaps (afferent2p.heatmaps):
    <mouse>_session_<day>.png / <experiment>_session.png   the sessions on the day's clock
    <experiment>_state.png, state_rates.csv                 activity by running / resting state, event overlays
    <experiment>_ecdf.png                                   ECDFs by condition x state
    <experiment>_heatmap.png, cfu_traces_<experiment>.npz   CFU identities x time (rows from the movies)
    <experiment>_updown.png, cfu_updown_<experiment>.csv    CFUs up / down / unchanged
    <experiment>_loco.png, _loco_cfu.png, loco_trials_*.csv, loco_cfu_*.csv, loco_glm_*.csv,
    loco_summary.csv, loco_conditions.csv                   locomotion trials, responses, GLMs
    motion_summary.csv                                      per-run motion / registration check (v33)
    pupil_measures.csv                                      pupil 0-1, variation, resting pupil, dilations
    acute_chronic.csv, <mouse>_acute_chronic.png            post/pre and day#/day 1 of every measure
  (--no-loco / --no-glm / --no-heatmaps / --no-movies / --no-motion switch parts off)

Event rate is events per minute (and per 1000 µm² of imaged field, since a
larger field sees more events). Amplitude/size/duration are summarised as
medians over events, and also as the median over runs of each run's median, so
the comparison can be made on the run level — which is what you can carry to a
group statistic.
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

try:
    import numpy  # noqa: F401
    import pandas  # noqa: F401
    import tifffile  # noqa: F401
    import scipy  # noqa: F401
except ImportError as e:
    raise SystemExit(
        f"\nMissing Python package: {e.name}\nInstall it into THIS Python "
        f"({sys.executable}):\n\n  \"{sys.executable}\" -m pip install -r requirements.txt\n"
        "(pandas tifffile scipy matplotlib h5py scikit-image; h5py reads AQuA2's v7.3 .mat files, "
        "scikit-image segments the fibres for the maps)\n")

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.aqua_io import (aqua_files, find_binned_tifs, read_feature_table, relocate, AQUA_SUFFIX,   # noqa: E402
                                read_cfu, read_res_summary, read_curation,
                                pick_binned_tif, curation_path)
from afferent2p.plan import plan_mouse                                              # noqa: E402
from afferent2p.globalfix import read_frames_csv, frame_map                         # noqa: E402

import matplotlib                                                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                     # noqa: E402

# event-property columns summarised per run / condition (short name -> label)
PROPS = {
    "dff_max": "peak dF/F",
    "df_max": "peak dF",
    "area_um2": "area (um2)",
    "duration_s": "duration (s)",
    "width55_s": "width 50-50 (s)",
    "rise19_s": "rise 10-90 (s)",
    "fall91_s": "decay 90-10 (s)",
    "dff_auc": "dF/F AUC",
    "circularity": "circularity",
}


# ----------------------------------------------------------------------------- per-run reading

def read_behavior(path: Path, log=print):
    """<movie>_behavior.csv as a DataFrame (None when there is none). When the
    csv lacks the pupil columns but a DLC eye csv sits next to the .sbx — the
    DLC file arrived after unpacking, or the movie predates the pupil columns —
    the csv is written again from it first (the movie is not touched)."""
    path = Path(path)
    try:
        from afferent2p.behavior import behavior_files, has_pupil, rebuild
        tif = Path(str(path)[:-len("_behavior.csv")] + ".tif") if path.name.endswith("_behavior.csv") else None
        if tif is not None and tif.exists() and not has_pupil(path) and behavior_files(tif)["dlc"] is not None:
            log(f"  {tif.name}: {'no behaviour csv' if not path.exists() else 'behaviour csv without pupil columns'} "
                f"but a DLC eye csv exists — writing {path.name} from it")
            rebuild(tif, log=log)
    except Exception as e:                     # never let this stop the analysis
        log(f"  behaviour csv rebuild skipped ({e})")
    if not path.exists():
        return None
    import pandas as pd
    return pd.read_csv(path)


def rejected_event_ids(tif: Path, cur: dict | None):
    """AQuA2 event ids that belong to CFUs marked as noise in cfu_curate.py."""
    if not cur or not cur.get("rejected"):
        return set(), set()
    f = aqua_files(tif)
    if not f["cfu"].exists():
        return set(), set()
    try:
        c = read_cfu(f["cfu"])
    except Exception:
        return set(), set()
    bad, good = set(), set()
    for k in range(c["n"]):
        ids = set(np.asarray(c["events"][k]).ravel().astype(int).tolist())
        (bad if (k + 1) in set(cur["rejected"]) else good).update(ids)
    return bad - good, good            # an event shared with a kept CFU is kept


def run_events(tif: Path, min_p: float = 0.0, gate: str = "none",
               still_thresh: float = 0.5, event_source: str = "drop-rejected-cfu",
               flag_margin: int = 1, keep_flagged: bool = False, log=print):
    """Event table for one run + run-level info. Returns (DataFrame, info) or None."""
    import pandas as pd
    f = aqua_files(tif)
    if not f["features"].exists():
        return None
    ev = read_feature_table(f["features"])
    info = {"run": tif.stem, "tif": str(tif)}
    spf = frames = None
    crop = None
    if f["json"].exists():
        try:
            j = json.loads(f["json"].read_text())
            spf = float(j.get("seconds_per_frame"))
            # frames_recorded is set when bad frames were dropped from the movie
            frames = int(j.get("frames_recorded") or j.get("frames_out"))
            info["um_per_px"] = None
            info["height"], info["width"] = int(j.get("height", 0)), int(j.get("width", 0))
            crop = j.get("crop_LRTB_px")
        except Exception:
            pass
    if spf is None or frames is None:                 # fall back to the AQuA2 result
        rs = read_res_summary(f["res"]) if f["res"].exists() else {}
        spf = spf or rs.get("frameRate") or 2.0
        sz = rs.get("sz")
        if frames is None and sz is not None and np.size(sz) >= 4:
            frames = int(np.asarray(sz).ravel()[3])
        info["height"] = info.get("height") or (int(np.asarray(sz).ravel()[0]) if sz is not None else 0)
        info["width"] = info.get("width") or (int(np.asarray(sz).ravel()[1]) if sz is not None else 0)
    # AQuA2 analyses the frame cropped by regMaskGap (burst.prep1), so the field
    # that produced these events is smaller than the tiff — use the analysed size
    # for the per-area rate when we can see it
    analysed = None
    n_result_frames = None                     # movie length the AQuA2 result was computed on
    if f["cfu"].exists():
        try:
            c = read_cfu(f["cfu"])
            analysed = c["maps"].shape[1:] if c["maps"].ndim == 3 and min(c["maps"].shape[1:]) > 0 else None
            if c["occurrence"].ndim == 2 and c["occurrence"].shape[1] > 1:
                n_result_frames = int(c["occurrence"].shape[1])
        except Exception:
            analysed = None
    if (analysed is None or n_result_frames is None) and f["res"].exists():
        try:
            sz = read_res_summary(f["res"]).get("sz")
        except Exception as e:                 # a truncated / half-written result (MATLAB stopped) must not end the analysis
            log(f"  {f['res'].name}: unreadable ({type(e).__name__}) — run AQuA2 on this tiff again (tab 4 with overwrite, then tab 5)")
            sz = None
        if sz is not None and np.size(sz) >= 2 and analysed is None:
            analysed = tuple(int(v) for v in np.asarray(sz).ravel()[:2])
        if sz is not None and np.size(sz) >= 4 and n_result_frames is None:
            n_result_frames = int(np.asarray(sz).ravel()[3])
    if analysed and all(analysed):
        info["crop_gap_px"] = int(round((info.get("height", analysed[0]) - analysed[0]) / 2))
        info["height"], info["width"] = int(analysed[0]), int(analysed[1])
    info["seconds_per_frame"] = float(spf)
    # frames flagged by the global-brightness step (unpack_bin / aqua_clean).
    # Since v23 they are dropped from the movie, so AQuA2's frame numbers are
    # movie frames: map them back to recorded frames (behaviour csv, time axis)
    fr = read_frames_csv(Path(str(tif.with_suffix("")) + "_frames.csv"))
    if fr is not None and fr["frame"].size > (frames or 0):
        frames = int(fr["frame"].size)
    info["frames"] = int(frames or 0)
    info["minutes"] = info["frames"] * info["seconds_per_frame"] / 60.0
    info["crop_LRTB_px"] = crop
    fmap = frame_map(fr)
    info["frames_dropped"] = int(fr["dropped"].sum()) if fr is not None else 0
    if fmap is not None and n_result_frames is not None and n_result_frames != fmap.size:
        if n_result_frames == fr["frame"].size:
            log(f"      WARNING: the AQuA2 result of {tif.name} has {n_result_frames} frames but the "
                f"movie now has {fmap.size}: it was computed before bad frames were dropped — "
                f"re-run the batch with --overwrite; treating its frame numbers as recorded frames")
            fmap = None
        else:
            log(f"      WARNING: the AQuA2 result of {tif.name} has {n_result_frames} frames, the movie "
                f"{fmap.size} — frame mapping may be wrong; re-run the batch with --overwrite")
        info["stale_result"] = True
    if "t0_frame" in ev.columns:
        t0_movie = ev["t0_frame"].to_numpy(int)
        ev["t0_frame_movie"] = t0_movie
        ev["t0_frame"] = fmap[np.clip(t0_movie - 1, 0, fmap.size - 1)] if fmap is not None else t0_movie
        ev["t0_s"] = (ev["t0_frame"] - 1) * info["seconds_per_frame"]
    n_all = len(ev)
    if min_p > 0 and "neg_log10_p" in ev.columns:
        ev = ev[ev["neg_log10_p"] >= min_p]
    # frames excluded by the brightness step (artefacts, plus state stretches
    # decided 'exclude'): events starting in them or next to them (a cut /
    # filler seam) leave the table, and their time leaves the denominator
    info["flagged_frames"] = int(fr["flagged"].sum()) if fr is not None else 0
    info["state_frames"] = int(fr["state"].sum()) if fr is not None else 0
    info["excluded_frames"] = int(fr["excluded"].sum()) if fr is not None else 0
    info["global_qc"] = fr is not None
    bad_frames = np.zeros(max(info["frames"], 1), bool)
    absent = np.zeros(max(info["frames"], 1), bool)
    if fr is not None and fr["excluded"].any():
        idx = fr["frame"][fr["excluded"]] - 1
        for m in range(-flag_margin, flag_margin + 1):
            j = np.clip(idx + m, 0, bad_frames.size - 1)
            bad_frames[j] = True
    if fr is not None and fr["dropped"].any():
        d = fr["frame"][fr["dropped"]] - 1
        absent[np.clip(d, 0, absent.size - 1)] = True
    if fr is not None and bad_frames.any() and not keep_flagged and "t0_frame" in ev.columns:
        before = len(ev)
        t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, bad_frames.size - 1)
        ev = ev[~bad_frames[t0]]
        info["n_events_dropped_flagged"] = before - len(ev)
    # state stretches the user kept: their events count — say how many, so the
    # decision can be revisited with numbers
    if fr is not None and fr["state"].any() and "t0_frame" in ev.columns:
        kept_state = fr["state"] & ~fr["excluded"]
        if kept_state.any():
            ks = np.zeros(max(info["frames"], 1), bool)
            ks[np.clip(fr["frame"][kept_state] - 1, 0, ks.size - 1)] = True
            t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, ks.size - 1)
            info["n_events_in_kept_state"] = int(ks[t0].sum())
            info["kept_state_frames"] = int(kept_state.sum())
    excluded = absent | (bad_frames if not keep_flagged else np.zeros_like(bad_frames))
    info["flagged_minutes"] = float(bad_frames.sum()) * info["seconds_per_frame"] / 60.0
    info["excluded_minutes"] = float(excluded.sum()) * info["seconds_per_frame"] / 60.0
    # manual CFU curation (cfu_curate.py)
    cur = read_curation(tif)
    info["curated"] = cur is not None
    info["n_cfu_rejected"] = len(cur["rejected"]) if cur else 0
    info["n_cfu_added"] = len(cur.get("added", [])) if cur else 0
    if cur and event_source != "all" and "event" in ev.columns:
        bad, good = rejected_event_ids(tif, cur)
        before = len(ev)
        if event_source == "kept-cfu-only":
            ev = ev[ev["event"].isin(good)]
        else:                                       # drop-rejected-cfu (default)
            ev = ev[~ev["event"].isin(bad)]
        info["n_events_dropped_by_curation"] = before - len(ev)
    beh = read_behavior(f["behavior"], log=log)
    info["has_behavior"] = beh is not None
    if beh is not None and "t0_frame" in ev.columns:
        idx = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, len(beh) - 1)
        for col, name in (("still_fraction", "still_fraction"), ("velocity_mean", "velocity"),
                          ("pupil_diam_px_mean", "pupil")):
            if col in beh.columns:
                ev[name] = beh[col].to_numpy()[idx]
        if "still_fraction" in beh.columns:
            info["still_frames"] = int((beh["still_fraction"].to_numpy() >= still_thresh).sum())
            info["still_minutes"] = info["still_frames"] * info["seconds_per_frame"] / 60.0
            if "velocity_mean" in beh.columns:
                info["velocity_median"] = float(np.nanmedian(beh["velocity_mean"]))
        if "pupil_diam_px_mean" in beh.columns:
            info["pupil_median"] = float(np.nanmedian(beh["pupil_diam_px_mean"]))
    # events and time split by locomotion state (still = still_fraction >= thresh at the
    # onset frame / in the frame), excluded frames left out — the locomotion-matched
    # comparison of the conditions; independent of --gate
    if beh is not None and "still_fraction" in beh.columns and "still_fraction" in ev.columns:
        sf = beh["still_fraction"].to_numpy(float)
        n = min(sf.size, excluded.size)
        ok = ~excluded[:n]
        st_frame = sf[:n] >= still_thresh
        info["still_minutes_used"] = float((st_frame & ok).sum()) * info["seconds_per_frame"] / 60.0
        info["running_minutes_used"] = float((~st_frame & np.isfinite(sf[:n]) & ok).sum()) * info["seconds_per_frame"] / 60.0
        evs = ev["still_fraction"].to_numpy(float)
        info["n_events_still"] = int(np.sum(evs >= still_thresh))
        info["n_events_running"] = int(np.sum(evs < still_thresh))
    if gate == "still":
        if "still_fraction" in ev.columns:
            ev = ev[ev["still_fraction"] >= still_thresh]
            info["gated"] = "still"
            still_ok = (beh["still_fraction"].to_numpy() >= still_thresh) if beh is not None else None
            if still_ok is not None:
                n = min(still_ok.size, excluded.size)
                info["minutes_used"] = float((still_ok[:n] & ~excluded[:n]).sum()) * info["seconds_per_frame"] / 60.0
            else:
                info["minutes_used"] = info.get("still_minutes", info["minutes"])
        else:
            info["gated"] = "still requested, no behaviour data"
            info["minutes_used"] = info["minutes"] - info["excluded_minutes"]
    else:
        info["gated"] = "none"
        info["minutes_used"] = info["minutes"] - info["excluded_minutes"]
    info["n_events_detected"] = n_all
    info["n_events_used"] = len(ev)
    ev = ev.reset_index(drop=True)
    ev.insert(0, "run", tif.stem)
    return ev, info


def run_metrics(ev, info: dict, um_per_px: float) -> dict:
    """One row of run_summary.csv."""
    px = info.get("height", 0) * info.get("width", 0)
    area_mm = px * (um_per_px ** 2) / 1e6                   # field area in mm2
    minutes = max(info.get("minutes_used", 0.0), 1e-9)
    m = {"run": info["run"], "minutes": round(info["minutes"], 2),
         "minutes_used": round(info.get("minutes_used", 0.0), 2),
         "seconds_per_frame": info["seconds_per_frame"],
         "field_area_mm2": round(area_mm, 4),
         "n_events": int(info["n_events_used"]),
         "n_events_detected": int(info["n_events_detected"]),
         "events_per_min": round(info["n_events_used"] / minutes, 4),
         "events_per_min_per_mm2": round(info["n_events_used"] / minutes / area_mm, 3) if area_mm else np.nan,
         "gated": info.get("gated", "none")}
    for k in ("flagged_frames", "state_frames", "excluded_frames", "frames_dropped", "n_events_dropped_flagged",
              "n_events_in_kept_state", "kept_state_frames", "n_cfu_rejected", "n_cfu_added",
              "n_events_dropped_by_curation", "crop_gap_px"):
        if k in info:
            m[k] = int(info[k])
    m["curated"] = int(bool(info.get("curated")))
    m["global_qc"] = int(bool(info.get("global_qc")))
    m["flagged_minutes"] = round(float(info.get("flagged_minutes", 0.0)), 3)
    for k in ("still_minutes", "velocity_median", "pupil_median"):
        if k in info:
            m[k] = round(float(info[k]), 3)
    if "still_minutes_used" in info:
        m["minutes_still"] = round(float(info["still_minutes_used"]), 3)
        m["minutes_running"] = round(float(info["running_minutes_used"]), 3)
        m["n_events_still"] = int(info["n_events_still"])
        m["n_events_running"] = int(info["n_events_running"])
        m["events_per_min_still"] = (round(info["n_events_still"] / info["still_minutes_used"], 4)
                                     if info["still_minutes_used"] >= 0.25 else np.nan)
        m["events_per_min_running"] = (round(info["n_events_running"] / info["running_minutes_used"], 4)
                                       if info["running_minutes_used"] >= 0.25 else np.nan)
    for c in PROPS:
        if c in ev.columns and len(ev):
            v = ev[c].to_numpy(float)
            v = v[np.isfinite(v)]
            m[f"{c}_median"] = round(float(np.median(v)), 4) if v.size else np.nan
        else:
            m[f"{c}_median"] = np.nan
    if "active_area_um2" not in m and "area_um2" in ev.columns and len(ev):
        m["area_um2_total"] = round(float(ev["area_um2"].sum()), 1)
    return m


# ----------------------------------------------------------------------------- CFU matching

def match_cfus(cfu_list, region_thr: float = 0.1, iou_thr: float = 0.3,
               align: bool = True, log=print):
    """Match CFUs across runs of one FOV.

    Runs are separate AQuA2 results (each registered internally), so the field
    can be shifted/rotated between them — the mean images (datPro) are aligned
    to the first run first (rotation + translation, afferent2p.prealign), then
    CFUs are matched greedily by intersection-over-union of their masks, like
    AQuA2's tools/multi_session_CFU_comparison.m does for two sessions.
    Returns (labels per run, n_matched, transforms).
    """
    from afferent2p.prealign import estimate_transform, apply_transform
    # AQuA2's registration crops every movie by its motion, so maps and mean
    # images differ by a few pixels between runs: pad to one canvas (top-left)
    shapes = [tuple(c["maps"].shape[1:]) for c in cfu_list if c["maps"].ndim == 3 and min(c["maps"].shape[1:]) > 0]
    shapes += [tuple(c["datPro"].shape) for c in cfu_list if c.get("datPro") is not None and np.ndim(c["datPro"]) == 2]
    if not shapes:
        return [np.full(c["n"], -1, int) for c in cfu_list], 0, [{"applied": False}] * len(cfu_list), 0
    H = max(sh[0] for sh in shapes); W = max(sh[1] for sh in shapes)
    if len(set(shapes)) > 1:
        log(f"    CFU maps / mean images differ in size between runs — padded to {H}x{W}")
    for c in cfu_list:
        if c["maps"].ndim == 3 and c["maps"].shape[1:] != (H, W):
            m = np.zeros((c["maps"].shape[0], H, W), np.float32)
            if c["maps"].size:
                m[:, :c["maps"].shape[1], :c["maps"].shape[2]] = c["maps"]
            c["maps"] = m
        if c.get("datPro") is not None and np.ndim(c["datPro"]) == 2 and c["datPro"].shape != (H, W):
            d = np.full((H, W), float(np.median(c["datPro"])), np.float32)
            d[:c["datPro"].shape[0], :c["datPro"].shape[1]] = c["datPro"]
            c["datPro"] = d
    ref = cfu_list[0]
    masks, transforms = [], []
    for i, c in enumerate(cfu_list):
        maps = c["maps"]
        tr = {"angle_deg": 0.0, "dy": 0, "dx": 0, "ncc": np.nan, "applied": False}
        if i > 0 and align and c["datPro"] is not None and ref["datPro"] is not None \
                and c["datPro"].shape == ref["datPro"].shape:
            t = estimate_transform(ref["datPro"], c["datPro"], max_angle_deg=4.0)
            before, after = float(t.get("ncc_before", 0.0)), float(t["ncc"])
            moved = abs(t["angle_deg"]) >= 0.5 or abs(t["dy"]) >= 2 or abs(t["dx"]) >= 2
            # only trust an alignment that clearly improves the match of the two
            # mean images; otherwise the fields are too dissimilar (different
            # plane, or no common structure) and warping would invent overlap
            if moved and after > max(0.3, 1.1 * before):
                tr = {**t, "applied": True}
                if maps.size:
                    maps = apply_transform(maps.astype(np.float32), t["angle_deg"],
                                           t["dy"], t["dx"], cval=0.0)
                log(f"    run {i + 1}: aligned to run 1 by {t['angle_deg']:+.2f} deg, "
                    f"(dy={t['dy']}, dx={t['dx']}); image match {before:.2f} -> {after:.2f}")
            else:
                tr = {**t, "applied": False}
                log(f"    run {i + 1}: NOT aligned (best transform {t['angle_deg']:+.2f} deg, "
                    f"dy={t['dy']}, dx={t['dx']} would give match {after:.2f} vs {before:.2f}) — "
                    "matching CFUs on the unshifted maps")
        transforms.append(tr)
        if maps.size == 0:
            masks.append(np.zeros((0, H, W), bool))
        else:
            masks.append(maps > region_thr)
    # greedy matching onto run 1's CFUs, then unmatched CFUs start new identities
    ids = [np.full(len(m), -1, int) for m in masks]
    ids[0] = np.arange(len(masks[0]))
    n_id = len(masks[0])
    for i in range(1, len(masks)):
        taken = set()
        for j in range(len(masks[i])):
            best, best_iou = -1, 0.0
            a = masks[i][j]
            if not a.any():
                continue
            for k in range(len(masks[0])):
                if k in taken:
                    continue
                b = masks[0][k]
                inter = np.logical_and(a, b).sum()
                if not inter:
                    continue
                iou = inter / np.logical_or(a, b).sum()
                if iou > best_iou:
                    best, best_iou = k, iou
            if best >= 0 and best_iou >= iou_thr:
                ids[i][j] = best
                taken.add(best)
            else:
                ids[i][j] = n_id
                n_id += 1
    n_matched = int(sum(1 for k in range(len(masks[0]))
                        if all((ids[i] == k).any() for i in range(len(masks)))))
    return ids, n_matched, transforms, n_id


def cfu_table(exp_name, runs, cfu_list, ids, n_id, um_per_px, conditions,
              kept_event_ids=None, log=print):
    """Per-CFU table: one row per matched CFU identity, its event rate in every
    run and pooled per condition, plus a paired comparison over the CFUs that
    are present in all runs (this is the readout CFU matching buys you: the same
    spatial unit before and after, so the comparison is paired)."""
    keys, seen = [], {}
    for r in runs:
        seen[r["condition"]] = seen.get(r["condition"], 0) + 1
        keys.append(f"{r['condition']}_r{seen[r['condition']]}")
    rows = []
    for cid in range(n_id):
        row = {"experiment": exp_name, "cfu": cid + 1}
        pool = {c: [0, 0.0] for c in conditions}          # events, minutes
        present = 0
        for i_run, (key, rinfo, c, idv) in enumerate(zip(keys, runs, cfu_list, ids)):
            cond = rinfo["condition"]
            minutes = max(rinfo.get("minutes_used", rinfo["minutes"]), 1e-9)
            pool[cond][1] += minutes
            hit = np.flatnonzero(idv == cid)
            if hit.size == 0:
                row[f"n_events[{key}]"] = 0
                row[f"events_per_min[{key}]"] = 0.0
                continue
            present += 1
            k = int(hit[0])
            ev_ids = np.asarray(c["events"][k]).ravel().astype(int)
            if kept_event_ids is not None:              # same gating as the event table
                keep_set = kept_event_ids[i_run]
                ev_ids = np.array([e for e in ev_ids if e in keep_set], int)
            n_ev = int(ev_ids.size)
            pool[cond][0] += n_ev
            row[f"n_events[{key}]"] = n_ev
            row[f"events_per_min[{key}]"] = round(n_ev / minutes, 4)
            row[f"area_um2[{key}]"] = round(float((c["maps"][k] > 0.1).sum()) * um_per_px ** 2, 1)
            dff = c["dff"][k]
            row[f"dff_p99[{key}]"] = round(float(np.nanpercentile(dff, 99)), 4) if dff.size else np.nan
        for cond in conditions:
            n_ev, minutes = pool[cond]
            row[f"events_per_min[{cond}]"] = round(n_ev / max(minutes, 1e-9), 4)
            row[f"n_events[{cond}]"] = n_ev
        base = conditions[0]
        for cond in conditions[1:]:
            b, o = row[f"events_per_min[{base}]"], row[f"events_per_min[{cond}]"]
            row[f"rate_change_per_min[{cond}_vs_{base}]"] = round(o - b, 4)
            row[f"rate_ratio[{cond}_vs_{base}]"] = round(o / b, 3) if b > 0 else np.nan
        row["runs_present"] = present
        row["in_all_runs"] = int(present == len(runs))
        rows.append(row)
    # paired comparison over CFUs seen in every run
    notes = []
    keep = [r for r in rows if r["in_all_runs"]]
    base = conditions[0]
    for cond in conditions[1:]:
        b = np.array([r[f"events_per_min[{base}]"] for r in keep], float)
        o = np.array([r[f"events_per_min[{cond}]"] for r in keep], float)
        if b.size < 3:
            notes.append(f"  CFU-level {cond} vs {base}: only {b.size} CFU(s) in every run — no test")
            continue
        up, down = int((o > b).sum()), int((o < b).sum())
        txt = (f"  CFU-level {cond} vs {base} ({b.size} CFUs in every run): "
               f"{up} up / {down} down, median {np.median(b):.3f} -> {np.median(o):.3f} events/min")
        try:
            from scipy.stats import wilcoxon
            if np.any(o != b):
                _, pv = wilcoxon(o, b)
                txt += f", Wilcoxon p = {pv:.3g} (within-FOV, exploratory)"
        except Exception:
            pass
        notes.append(txt)
    return rows, notes


def manual_cfu_rows(exp_name, infos, runs, conditions, use_curation=True):
    """Rows for regions added by hand in cfu_curate.py. They are matched across
    runs by their order in each run's sidecar (a region drawn on run 1 and then
    on run 2 is region 1 in both), and flagged source = manual."""
    if not use_curation:
        return []
    per_run = []
    for r in runs:
        cur = read_curation(Path(r["tif"]))
        per_run.append((cur or {}).get("added", []))
    n = max((len(a) for a in per_run), default=0)
    if n == 0:
        return []
    keys, seen = [], {}
    for info in infos:
        seen[info["condition"]] = seen.get(info["condition"], 0) + 1
        keys.append(f"{info['condition']}_r{seen[info['condition']]}")
    rows = []
    for j in range(n):
        row = {"experiment": exp_name, "cfu": f"manual{j + 1}", "source": "manual"}
        pool = {c: [0, 0.0] for c in conditions}
        present = 0
        for key, info, added in zip(keys, infos, per_run):
            minutes = max(info.get("minutes_used", info["minutes"]), 1e-9)
            pool[info["condition"]][1] += minutes
            if j >= len(added):
                row[f"n_events[{key}]"] = 0
                row[f"events_per_min[{key}]"] = 0.0
                continue
            present += 1
            rec = added[j]
            pool[info["condition"]][0] += int(rec["n_events"])
            row[f"n_events[{key}]"] = int(rec["n_events"])
            row[f"events_per_min[{key}]"] = round(rec["n_events"] / minutes, 4)
            row[f"area_um2[{key}]"] = None
            row[f"dff_p99[{key}]"] = rec.get("peak_dff")
        for cond in conditions:
            n_ev, minutes = pool[cond]
            row[f"events_per_min[{cond}]"] = round(n_ev / max(minutes, 1e-9), 4)
            row[f"n_events[{cond}]"] = n_ev
        base = conditions[0]
        for cond in conditions[1:]:
            b, o = row[f"events_per_min[{base}]"], row[f"events_per_min[{cond}]"]
            row[f"rate_change_per_min[{cond}_vs_{base}]"] = round(o - b, 4)
            row[f"rate_ratio[{cond}_vs_{base}]"] = round(o / b, 3) if b > 0 else None
        row["runs_present"] = present
        row["in_all_runs"] = int(present == len(infos))
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------- figures

def overview_figure(exp_name, events, run_rows, conditions, out_png, log=print, rate_ymax=None):
    import pandas as pd
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    rs = pd.DataFrame(run_rows)
    # 1. event rate per run, grouped by condition, with the exact Poisson interval of each run;
    #    the y axis is shared by every experiment of the animal (rate_ymax) so days are comparable at a glance
    ax = axes[0, 0]
    x = 0
    ticks, labels = [], []
    for c in conditions:
        sub = rs[rs["condition"] == c]
        xs = np.arange(len(sub)) + x
        ax.bar(xs, sub["events_per_min"], color="#2a78d6", width=0.8)
        for xi, (n, mins) in zip(xs, zip(sub["n_events"], sub["minutes_used"])):
            lo, hi = poisson_rate_ci(float(n), float(mins))
            if np.isfinite(lo):
                ax.plot([xi, xi], [lo, hi], color="#0b0b0b", lw=1.2)
            ax.text(xi, (hi if np.isfinite(hi) else 0) + 0.01 * (rate_ymax or 1), str(int(n)), fontsize=7, ha="center", va="bottom", color="#52514e")
        ticks.append(x + (len(sub) - 1) / 2); labels.append(c)
        x += len(sub) + 1
    ax.set_xticks(ticks); ax.set_xticklabels(labels)
    if rate_ymax:
        ax.set_ylim(0, rate_ymax)
    ax.set_ylabel("events / min"); ax.set_title("events/min per run (95 % Poisson bars" + ("; y axis shared by all experiments)" if rate_ymax else ")"), fontsize=9)
    # 2-6. property distributions (ECDF) per condition
    panels = [("dff_max", axes[0, 1]), ("area_um2", axes[0, 2]), ("duration_s", axes[1, 0]),
              ("rise19_s", axes[1, 1]), ("fall91_s", axes[1, 2])]
    for col, ax in panels:
        if col not in events.columns:
            ax.axis("off"); continue
        for c in conditions:
            v = events.loc[events["condition"] == c, col].to_numpy(float)
            v = np.sort(v[np.isfinite(v)])
            if v.size:
                ax.plot(v, np.arange(1, v.size + 1) / v.size, label=f"{c} (n={v.size})")
        ax.set_xlabel(PROPS.get(col, col)); ax.set_ylabel("cumulative fraction")
        ax.legend(fontsize=8)
        if col in ("area_um2", "dff_auc"):
            ax.set_xscale("log")
    fig.suptitle(f"{exp_name}: AQuA2 events by condition "
                 "(ECDFs pool events; statistics belong on run/animal level)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def cfu_figure(exp_name, runs, cfu_list, ids, out_png):
    n = len(cfu_list)
    fig, axes = plt.subplots(1, n + 1, figsize=(4.2 * (n + 1), 4.4))
    axes = np.atleast_1d(axes)
    for i, (rinfo, c) in enumerate(zip(runs, cfu_list)):
        ax = axes[i]
        if c["datPro"] is not None:
            ax.imshow(c["datPro"], cmap="gray")
        lab = np.zeros(c["maps"].shape[1:], np.float32) if c["maps"].size else None
        if lab is not None:
            for k in range(c["maps"].shape[0]):
                lab[c["maps"][k] > 0.1] = ids[i][k] + 1
            ax.imshow(np.ma.masked_where(lab == 0, lab), cmap="tab20", alpha=0.6, interpolation="nearest")
        ax.set_title(f"{rinfo['condition']}  {rinfo['run'][-18:]}\n{c['n']} CFUs", fontsize=9)
        ax.axis("off")
    ax = axes[-1]
    shared = None
    for i, c in enumerate(cfu_list):
        if c["maps"].size == 0:
            continue
        m = np.zeros(c["maps"].shape[1:], bool)
        for k in range(c["maps"].shape[0]):
            m |= c["maps"][k] > 0.1
        shared = m if shared is None else (shared & m)
    if shared is not None:
        ax.imshow(shared, cmap="gray")
    ax.set_title("area covered by CFUs in ALL runs", fontsize=9); ax.axis("off")
    fig.suptitle(f"{exp_name}: CFU maps, colour = matched identity")
    fig.tight_layout(); fig.savefig(out_png, dpi=110); plt.close(fig)


# ----------------------------------------------------------------------------- figures: behaviour adjustment, effects

COND_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"


def _log_ticks(ax, axis: str):
    """Plain-number ticks on a log ratio axis, no minor labels."""
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    lo, hi = (ax.get_xlim() if axis == "x" else ax.get_ylim())
    cands = [0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 6, 8, 10, 15, 20]
    ticks = [t for t in cands if lo * 0.999 <= t <= hi * 1.001] or [1]
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks)); a.set_minor_locator(NullLocator())
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def _style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)


def adjusted_figure(exp_name: str, adj: dict, conditions: list, glm_runs: list, spf: float, out_png: Path):
    """Bins: events/min against running speed per condition with the model's
    prediction; the locomotion-onset-triggered rate; the rate ratios."""
    from afferent2p.adjust import bin_run
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw=dict(wspace=0.35, left=0.06, right=0.98, top=0.82, bottom=0.16))
    colour = {c: COND_COLOURS[i % len(COND_COLOURS)] for i, c in enumerate(conditions)}
    ax = axes[0]; _style(ax)
    coef = adj.get("coef", {})
    for c in conditions:
        xs, ys = [], []
        for r in glm_runs:
            if r["condition"] != c:
                continue
            t = bin_run(r["onsets"], r["excluded"], r["speed"], r["still"], r["pupil"], spf, adj.get("bin_s", 20.0))
            ok = t["minutes"] > 0
            xs.append(t["speed"][ok]); ys.append(t["events"][ok] / t["minutes"][ok])
        if xs:
            x = np.concatenate(xs); y = np.concatenate(ys)
            ax.plot(x, y, "o", ms=4, color=colour[c], alpha=0.45, mec="white", mew=0.4, label=c)
    ax.set_xlabel("running speed in the bin (encoder units)", fontsize=8)
    ax.set_ylabel(f"events / min per {adj.get('bin_s', 20):.0f}-s bin", fontsize=8)
    ax.set_title("every bin: rate vs running, per condition", fontsize=8, loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]; _style(ax)
    trig = adj.get("onset_triggered", {})
    for c, t in trig.items():
        ax.plot(t["t_s"], t["rate"], color=colour.get(c, "#333"), lw=2, label=f"{c} (n = {t['n_onsets']} onsets)")
    ax.axvline(0, color="#999999", lw=0.8)
    ax.set_xlabel("time from running onset (s)", fontsize=8); ax.set_ylabel("events / min", fontsize=8)
    ax.set_title("event rate around running onsets (after >= 20 s still)", fontsize=8, loc="left")
    if trig:
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.text(0.5, 0.5, "fewer than 3 running onsets", ha="center", va="center", transform=ax.transAxes, color=INK2, fontsize=8)
    ax = axes[2]; _style(ax)
    base = conditions[0]
    labels, vals, los, his = [], [], [], []
    for c in conditions[1:]:
        for name, key in (("raw", f"raw_ratio_{c}"), ("adjusted", f"adjusted_ratio_{c}")):
            v = adj.get(key, np.nan)
            ci = adj.get(f"adjusted_ratio_ci_{c}") if name == "adjusted" else None
            labels.append(f"{c} vs {base}\n{name}"); vals.append(v)
            los.append(max(v - ci[0], 0) if ci else 0); his.append(max(ci[1] - v, 0) if ci else 0)
    xpos = np.arange(len(labels))
    ax.axhline(1, color="#999999", lw=0.8)
    ax.errorbar(xpos, vals, yerr=[los, his], fmt="o", color=INK, ms=7, capsize=4, lw=1.5)
    ax.set_xticks(xpos); ax.set_xticklabels(labels, fontsize=8); ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylabel("rate ratio", fontsize=8); ax.set_yscale("log")
    _log_ticks(ax, "y")
    ax.set_title(f"rate ratio: raw vs at the same behaviour (95 % CI, block bootstrap)\n"
                 f"behaviour explains {100 * adj.get('deviance_explained_behaviour', np.nan):.0f} % of the rate variation, "
                 f"the condition {100 * adj.get('deviance_explained_condition', np.nan):.0f} %", fontsize=8, loc="left")
    fig.suptitle(f"{exp_name} — activity and running", fontsize=10, x=0.06, ha="left")
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


from afferent2p.stats import poisson_rate_ci, rate_ratio_ci, change_measures   # noqa: E402,F401  (moved to the package in v34)


def read_tuning_json(mouse_dir, exp_name: str) -> dict | None:
    """aqua_opts\\<experiment>_tuning.json written by aqua_tune --apply (v36): which rule chose the day's settings
    (matched sensitivity / one common setting) and, per condition, the measured recall and false-alarm rate."""
    p = Path(mouse_dir) / "aqua_opts" / f"{exp_name}_tuning.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def detection_fields(tj: dict | None, base: str, cond: str, rate_b, rate_c) -> dict:
    """The detection's own numbers for one comparison: settings per condition, sensitivity (recall of the injected
    transients), false-alarm rate (null movie), and the rates with the expected false alarms subtracted."""
    out: dict = {}
    if not tj:
        return out
    cb, cc = (tj.get("conditions") or {}).get(base), (tj.get("conditions") or {}).get(cond)
    out["detection_rule"] = tj.get("rule")
    out["detection_comparable"] = bool(tj.get("comparable", True))
    if not (cb and cc):
        return out

    def lab(d):
        st = d.get("setting") or {}
        return " ".join(f"{k}={v}" for k, v in st.items() if k in ("thrARScl", "minSize", "maxSize", "smoXY", "minDur"))
    out["setting_baseline"], out["setting_condition"] = lab(cb), lab(cc)
    out["recall_baseline"], out["recall_condition"] = cb.get("recall"), cc.get("recall")
    out["fp_rate_baseline"], out["fp_rate_condition"] = cb.get("fp_rate_all_per_min", cb.get("fp_rate_per_min")), cc.get("fp_rate_all_per_min", cc.get("fp_rate_per_min"))
    rb_, rc_ = out["recall_baseline"], out["recall_condition"]
    out["sensitivity_gap"] = abs(float(rb_) - float(rc_)) if (rb_ is not None and rc_ is not None) else np.nan
    out["sensitivity_matched"] = bool(tj.get("matched")) or (np.isfinite(out["sensitivity_gap"]) and out["sensitivity_gap"] <= 0.05)
    fb, fc = out["fp_rate_baseline"], out["fp_rate_condition"]
    try:
        if rate_b is not None and rate_c is not None and fb is not None and fc is not None:
            rbc, rcc = max(float(rate_b) - float(fb), 0.0), max(float(rate_c) - float(fc), 0.0)
            out["rate_baseline_fpcorr"], out["rate_condition_fpcorr"] = rbc, rcc
            out["ratio_fpcorr"] = (rcc / rbc) if rbc > 0 else np.nan
    except (TypeError, ValueError):
        pass
    return out


def describe_detection(d: dict, base: str, cond: str) -> str:
    """One line for summary.txt: how pre and post were detected and what the false-alarm subtraction gives."""
    if not d or not d.get("detection_rule"):
        return ""
    rule = d["detection_rule"]
    if "setting_baseline" not in d:
        return f"  detection: settings by the '{rule}' rule" + ("" if d.get("detection_comparable", True) else " — NOT COMPARABLE")
    head = ("matched sensitivity" if rule == "matched" else "one common setting")
    txt = (f"  detection ({head}): {base} {d['setting_baseline']} (recall {d.get('recall_baseline', float('nan')):.2f}, "
           f"false alarms {d.get('fp_rate_baseline', float('nan')):.2f}/min) / {cond} {d['setting_condition']} "
           f"(recall {d.get('recall_condition', float('nan')):.2f}, false alarms {d.get('fp_rate_condition', float('nan')):.2f}/min)")
    if "ratio_fpcorr" in d:
        txt += (f"; with the expected false alarms subtracted: {d['rate_baseline_fpcorr']:.3f} -> {d['rate_condition_fpcorr']:.3f}/min"
                + (f" = x{d['ratio_fpcorr']:.2f}" if np.isfinite(d.get("ratio_fpcorr", np.nan)) else ""))
    if not d.get("detection_comparable", True):
        txt += ("  -> NOT COMPARABLE: the two conditions could not be detected with the same sensitivity "
                f"(recall gap {d.get('sensitivity_gap', float('nan')):.2f}); read the change with that in mind")
    elif np.isfinite(d.get("sensitivity_gap", np.nan)) and d["sensitivity_gap"] > 0.05:
        txt += f"  (sensitivity gap {d['sensitivity_gap']:.2f}: the raw ratio carries part of it)"
    return txt


def describe_change(e: dict) -> str:
    """One line for summary.txt: the change with its uncertainty, in words that do not depend on the scale."""
    nb, nc = e.get("n_events_baseline"), e.get("n_events_condition")
    rb, rc = e.get("rate_baseline", np.nan), e.get("rate_condition", np.nan)
    lo, hi = e.get("ratio_ci_lo", np.nan), e.get("ratio_ci_hi", np.nan)
    d = e.get("delta_rate", np.nan)
    if nb is None or nc is None:
        return ""
    ci = ""
    if np.isfinite(lo) and (np.isfinite(hi) or hi == np.inf):
        ci = f" (95 % CI {lo:.2f}-{'inf' if hi == np.inf else f'{hi:.2f}'})"
    if nb + nc == 0:
        return "    scale: no events in either condition"
    up = np.isfinite(lo) and lo > 1
    down = np.isfinite(hi) and hi < 1
    if e.get("low_counts"):
        if up or down:
            verdict = (f"the direction ({'increase' if up else 'decrease'}) is unlikely to be chance (interval excludes 1), but with "
                       f"{nb} vs {nc} events the SIZE cannot be put on the same scale as a well-sampled experiment — compare the "
                       f"rate difference ({d:+.3f}/min), not the ratio")
        else:
            verdict = f"too few events ({nb} vs {nc}) for a conclusion"
    elif up:
        verdict = "an increase the interval supports"
    elif down:
        verdict = "a decrease the interval supports"
    else:
        verdict = "not distinguishable from no change (the interval includes 1)"
    ratio = f"ratio x{rc / rb:.2f}" if rb > 0 else "ratio undefined (no baseline events)"
    return (f"    scale: {nb} -> {nc} events, {rb:.3f} -> {rc:.3f} events/min, difference {d:+.3f}/min, "
            f"{ratio}{ci}, log2 ratio with pseudo-count {e.get('log2_ratio_pc', np.nan):+.2f}, "
            f"modulation index {e.get('modulation_index', np.nan):+.2f} -> {verdict}")


def effects_figure(effects: list, out_png: Path, day_labels: dict | None = None):
    """One row per experiment: the rate ratio vs baseline by every control. The raw
    ratio carries its exact Poisson interval; an experiment without baseline
    events is drawn at its pseudo-count ratio (hollow) so it does not vanish;
    rows with fewer than 10 events in total are labelled as such."""
    import pandas as pd
    df = pd.DataFrame(effects)
    day_labels = day_labels or {}

    def cond_name(c):
        c = str(c)
        lab = day_labels.get(c[1:]) if c.startswith("d") and c[1:].isdigit() else None
        return f"{c} ({lab})" if lab else c
    measures = [("ratio_raw", "raw (95 % interval)", "#8a8a8a"), ("ratio_still", "still frames only", "#1baf7a"),
                ("ratio_running", "running frames only", "#eda100"), ("ratio_pupil_matched", "matched pupil size", "#e87ba4"),
                ("ratio_adjusted", "behaviour-adjusted (regression, 95 % CI)", "#2a78d6")]
    measures = [m for m in measures if m[0] in df.columns and np.isfinite(df[m[0]].astype(float)).any()]
    n = len(df)
    fig, ax = plt.subplots(figsize=(10, 1.0 + 0.55 * n), gridspec_kw=dict(left=0.3, right=0.97, top=0.82, bottom=0.12))
    _style(ax); ax.grid(axis="x", color=GRID, lw=0.8); ax.grid(axis="y", visible=False)
    ax.axvline(1, color="#999999", lw=1.0)
    ys = np.arange(n)[::-1]
    k = len(measures)
    xmin, xmax = 0.08, 12.0
    for j, (key, label, col) in enumerate(measures):
        off = (j - (k - 1) / 2) * (0.7 / max(k, 1))
        v = df[key].astype(float).to_numpy()
        if key == "ratio_adjusted" and "ratio_adjusted_ci_lo" in df.columns:
            lo = df["ratio_adjusted_ci_lo"].astype(float).to_numpy(); hi = df["ratio_adjusted_ci_hi"].astype(float).to_numpy()
            ok = np.isfinite(v) & np.isfinite(lo) & np.isfinite(hi)
            ax.hlines(ys[ok] + off, lo[ok], hi[ok], color=col, lw=2, alpha=0.6)
        if key == "ratio_raw" and "ratio_ci_lo" in df.columns:
            lo = df["ratio_ci_lo"].astype(float).to_numpy(); hi = df["ratio_ci_hi"].astype(float).to_numpy()
            for i in range(n):
                if not np.isfinite(lo[i]):
                    continue
                l_ = max(lo[i], xmin); h_ = min(hi[i], xmax) if np.isfinite(hi[i]) else xmax
                ax.hlines(ys[i] + off, l_, h_, color=col, lw=2, alpha=0.55)
                if not np.isfinite(hi[i]) or hi[i] > xmax:
                    ax.annotate("", xy=(xmax, ys[i] + off), xytext=(xmax * 0.8, ys[i] + off),
                                arrowprops=dict(arrowstyle="->", color=col, lw=1.2))
            # no baseline events: show the pseudo-count ratio hollow
            if "log2_ratio_pc" in df.columns:
                pc = 2.0 ** df["log2_ratio_pc"].astype(float).to_numpy()
                miss = ~np.isfinite(v) & np.isfinite(pc)
                if miss.any():
                    ax.plot(np.clip(pc[miss], xmin, xmax), ys[miss] + off, "o", mfc="white", mec=col, ms=7, mew=1.4,
                            label="raw, pseudo-count (no baseline events)")
        ax.plot(np.clip(v, xmin, xmax), ys + off, "o", color=col, ms=7, mec="white", mew=0.6, label=label)
    ax.set_yticks(ys)
    labs = []
    for r in df.itertuples():
        nb = getattr(r, "n_events_baseline", None); nc = getattr(r, "n_events_condition", None)
        cnt = f"  [{int(nb)} -> {int(nc)} events" + (", too few" if bool(getattr(r, "low_counts", False)) else "") + "]" \
            if nb is not None and nc is not None and np.isfinite(float(nb)) and np.isfinite(float(nc)) else ""
        labs.append(f"{r.experiment}\n{cond_name(r.condition)} vs {cond_name(r.baseline)}  (still {100 * (r.still_frac_baseline or 0):.0f} % -> "
                    f"{100 * (r.still_frac_condition or 0):.0f} %){cnt}")
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_xscale("log"); ax.set_xlim(xmin, xmax)
    _log_ticks(ax, "x")
    ax.set_xlabel("event-rate ratio vs baseline (log scale; 1 = no change; bars = 95 % intervals, arrows = beyond the axis)", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, 1.2), ncol=3)
    ax.set_title("change in afferent event rate — raw and with each behaviour control",
                 fontsize=9, loc="left", pad=34)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def rates_figure(effects: list, all_cond: list, all_runs: list, out_png: Path, day_labels: dict | None = None):
    """Every condition of every experiment on ONE rate axis (so a change of
    0 -> 0.3 events/min and one of 1.6 -> 4.6 are seen at their true size),
    with exact Poisson intervals and the single runs; below it the change as
    the log2 ratio with a pseudo-count and its exact interval, on one axis."""
    import pandas as pd
    day_labels = day_labels or {}
    if not effects or not all_cond:
        return None
    ef = pd.DataFrame(effects); cd = pd.DataFrame(all_cond); rn = pd.DataFrame(all_runs) if all_runs else None

    def cond_name(c):
        c = str(c)
        lab = day_labels.get(c[1:]) if c.startswith("d") and c[1:].isdigit() else None
        return f"{c}\n{lab}" if lab else c
    exps = list(dict.fromkeys(ef["experiment"]))
    # x layout: one group per experiment, its conditions in order (baseline first)
    groups, x = [], 0
    for e in exps:
        base = ef[ef["experiment"] == e]["baseline"].iloc[0]
        conds = [base] + list(ef[ef["experiment"] == e]["condition"])
        groups.append((e, conds, x))
        x += len(conds) + 1
    width = max(8.0, 0.75 * x)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, 7.6), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 1.0], hspace=0.22))
    _style(ax1); _style(ax2)
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    ticks, tlabels, glabels = [], [], []
    ymax = 0.0
    for e, conds, x0 in groups:
        xs, ys = [], []
        for i, c in enumerate(conds):
            row = cd[(cd["experiment"] == e) & (cd["condition"] == c)]
            if not len(row):
                continue
            row = row.iloc[0]
            n, mins = float(row["n_events"]), float(row["minutes"])
            r = n / mins if mins > 0 else np.nan
            lo, hi = poisson_rate_ci(n, mins)
            xi = x0 + i
            col = cols[min(i, len(cols) - 1)]
            ax1.vlines(xi, lo, hi, color=col, lw=2, alpha=0.6)
            ax1.plot([xi], [r], "o", color=col, ms=7, mec="white", mew=0.6, zorder=3)
            ax1.text(xi, hi, f"{int(n)}", fontsize=7, color=INK2, ha="center", va="bottom")
            if rn is not None and "events_per_min" in rn.columns:
                sub = rn[(rn["experiment"] == e) & (rn["condition"] == c)]
                ax1.plot(np.full(len(sub), xi) + np.linspace(-0.18, 0.18, max(len(sub), 1)), sub["events_per_min"], ".", color=col, ms=5, alpha=0.7, zorder=2)
            xs.append(xi); ys.append(r); ticks.append(xi); tlabels.append(cond_name(c))
            ymax = max(ymax, hi if np.isfinite(hi) else r)
        if len(xs) > 1:
            ax1.plot(xs, ys, "-", color="#9a9aa2", lw=1, zorder=1)
        # the change rows of this experiment
        for i, c in enumerate(conds[1:], start=1):
            er = ef[(ef["experiment"] == e) & (ef["condition"] == c)]
            if not len(er):
                continue
            er = er.iloc[0]
            l2 = float(er.get("log2_ratio_pc", np.nan))
            lo, hi = float(er.get("ratio_ci_lo", np.nan)), float(er.get("ratio_ci_hi", np.nan))
            xi = x0 + i
            col = cols[min(i, len(cols) - 1)]
            if np.isfinite(lo) and lo > 0 and np.isfinite(hi):
                ax2.vlines(xi, np.log2(lo), np.log2(hi), color=col, lw=2, alpha=0.6)
            elif np.isfinite(lo) and lo > 0:                       # open upper end (no baseline events)
                ax2.vlines(xi, np.log2(lo), 6, color=col, lw=2, alpha=0.35, ls="--")
            elif np.isfinite(hi) and hi > 0 and hi < np.inf:       # no condition events: open lower end
                ax2.vlines(xi, -6, np.log2(hi), color=col, lw=2, alpha=0.35, ls="--")
            if np.isfinite(l2):
                low = bool(er.get("low_counts", False))
                nb_, nc_ = er.get("n_events_baseline", np.nan), er.get("n_events_condition", np.nan)
                none = np.isfinite(float(nb_)) and np.isfinite(float(nc_)) and float(nb_) + float(nc_) == 0
                ax2.plot([xi], [l2], "o", color=col, mfc="white" if low else col, mec=col, ms=7, mew=1.4, zorder=3)
                if none:
                    ax2.text(xi, l2, "   no events", fontsize=7, color=INK2, va="center")
                elif low:
                    ax2.text(xi, l2, f"   {int(nb_)} -> {int(nc_)} events", fontsize=7, color=INK2, va="center")
        glabels.append((x0 + (len(conds) - 1) / 2, _exp_short(e, day_labels)))
    ax1.set_ylim(0, ymax * 1.15 if ymax > 0 else 1)
    ax1.set_ylabel("events / min (one axis for every experiment)", fontsize=9)
    ax1.set_title("afferent event rate per condition on ONE axis\n(dots = conditions with exact 95 % Poisson intervals, small dots = single runs, numbers = events)",
                  fontsize=9, loc="left")
    ax2.axhline(0, color="#999999", lw=1)
    ax2.set_ylim(-6, 6)
    ax2.set_yticks([-4, -2, -1, 0, 1, 2, 4]); ax2.set_yticklabels(["÷16", "÷4", "÷2", "same", "×2", "×4", "×16"])
    ax2.set_ylabel("change vs baseline\n(log2 ratio, pseudo-count 0.5)", fontsize=9)
    ax2.set_title("the change vs baseline on ONE axis\n(bars = exact 95 % interval of the rate ratio; dashed = open-ended, a condition without "
                  "events; hollow = too few events to size the change: < 10 in total or < 3 in a condition)", fontsize=9, loc="left")
    ax2.set_xticks(ticks); ax2.set_xticklabels(tlabels, fontsize=8)
    for xg, lab in glabels:
        ax1.text(xg, ax1.get_ylim()[1] * 0.98, lab, fontsize=8, ha="center", va="top", color=INK2,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _exp_short(exp: str, day_labels: dict) -> str:
    """'260917 sham FOV1' for <mouse>_260917_FOV1, 'FOV1 (days)' for <mouse>_FOV1_days."""
    m = re.search(r"_(\d{6})_FOV(\d+)$", exp)
    if m:
        lab = day_labels.get(m.group(1))
        return f"{m.group(1)}{' ' + lab if lab else ''}  FOV{m.group(2)}"
    m = re.search(r"_FOV(\d+)_days$", exp)
    if m:
        return f"FOV{m.group(1)} across days"
    return exp


# ----------------------------------------------------------------------------- main per mouse

def discover(mouse_dir: Path, sec: float | None, paradigm="auto", n_fov=None,
             stim_after_run=None, log=print):
    """[(experiment name, fov, paradigm, [ {run, condition, tif} ])] from the
    aqua_batch manifest if present, else by re-planning the mouse folder."""
    man = mouse_dir / "aqua_batch" / "manifest.json"
    if man.exists():
        m = json.loads(man.read_text())
        out = []
        for e in m["experiments"]:
            runs = [{"run": r["run"], "condition": r["condition"], "tif": relocate(r["tif"], mouse_dir)}
                    for r in e.get("runs", [])]
            if runs:
                out.append((e["name"], e.get("fov"), e.get("paradigm"), runs))
        if out:
            return out
        log(f"  {man} has no runs — re-planning")
    exps = plan_mouse(mouse_dir, paradigm=paradigm, n_fov=n_fov, stim_after_run=stim_after_run)
    out = []
    for e in exps:
        runs = []
        for sbx_name, cond in e.runs:
            hits = list(mouse_dir.rglob(sbx_name))
            if not hits:
                continue
            tif = pick_binned_tif(hits[0].parent, sec, log=log)
            if tif:
                runs.append({"run": tif.stem, "condition": cond, "tif": tif})
        if runs:
            out.append((e.name, e.fov, e.paradigm, runs))
    return out


def _cfu_for(tif: Path, use_curation: bool = True):
    """CFU result of one run with the CFUs marked noise removed."""
    c = read_cfu(aqua_files(tif)["cfu"])
    cu = read_curation(tif) if use_curation else None
    c["orig_index"] = np.arange(c["n"]) + 1                 # AQuA2's CFU numbers (1-based) of the rows kept
    if cu and cu.get("rejected"):
        keep = np.array([i + 1 not in set(cu["rejected"]) for i in range(c["n"])])
        c["maps"] = c["maps"][keep] if c["maps"].size else c["maps"]
        c["events"] = [e for e, k in zip(c["events"], keep) if k]
        c["occurrence"] = c["occurrence"][keep] if c["occurrence"].size else c["occurrence"]
        c["dff"] = c["dff"][keep] if c["dff"].size else c["dff"]
        c["freq"] = [f for f, k in zip(c["freq"], keep) if k]
        c["orig_index"] = c["orig_index"][keep]
        c["n"] = int(keep.sum())
    return c


def analyze_mouse(mouse_dir: Path, um_per_px: float = 0.792, sec: float | None = 2.0,
                  gate: str = "none", still_thresh: float = 0.5, min_p: float = 0.0,
                  baseline: str | None = None, cfu: bool = True, cfu_iou: float = 0.3,
                  paradigm="auto", n_fov=None, stim_after_run=None,
                  event_source: str = "drop-rejected-cfu", use_curation: bool = True,
                  flag_margin: int = 1, keep_flagged: bool = False, pupil_analysis: bool = True,
                  maps: bool = True, tile_um: float = 25.0, timecourse: bool = True, bin_min: float = 2.0,
                  assumed_stim_gap_min: float = 73.0, day_labels: dict | None = None,
                  loco: bool = True, loco_glm: bool = True, heatmaps: bool = True, movies: bool = True,
                  motion: bool = True, log=print):
    import pandas as pd
    mouse_dir = Path(mouse_dir)
    out_dir = mouse_dir / "aqua_analysis"
    out_dir.mkdir(exist_ok=True)
    motion_rows = []                              # v33: per-run motion / registration check (afferent2p.motion)
    effects = []                                  # one row per experiment x condition (vs its baseline)
    speed_rows, tuning_results, tc_store = [], [], []   # matched-speed bands; data for the day time course
    # the general representation / locomotion trials / heatmaps (v32) collect per-run bundles
    exp_bundles, exp_meta = [], []
    loco_run_rows, loco_cond_rows, loco_trial_rows, state_rows_all, pupil_rows_all, heat_rows = [], [], [], [], [], []
    trigger_rows = []
    trials_by_run = {}
    cfu_cache = {}

    def cfu_of(tif: Path):
        key = str(tif)
        if key not in cfu_cache:
            try:
                cfu_cache[key] = _cfu_for(tif, use_curation) if aqua_files(tif)["cfu"].exists() else None
            except Exception as e:
                log(f"  CFU file of {tif.name} unreadable ({e})")
                cfu_cache[key] = None
        return cfu_cache[key]
    from afferent2p.timecourse import read_day_labels, day_of
    day_labels = read_day_labels(mouse_dir, day_labels)   # {yymmdd: 'sham' ...} from day_labels.csv / --day-label

    def exp_label(exp_name: str, runs_of_exp: list) -> str:
        """Experiment name with the day's label (optostim: one day per experiment) or,
        for an experiment that spans days (spontaneous / Treg: conditions = the
        dates), every day's label: d260911 = day -1, d260918 = day 7."""
        days = list(dict.fromkeys(day_of(r["run"]) for r in runs_of_exp if day_of(r["run"])))
        if len(days) <= 1:
            lab = day_labels.get(days[0] if days else "", "")
            return f"{exp_name}  ({lab})" if lab else exp_name
        parts = [f"d{d} = {day_labels[d]}" for d in days if day_labels.get(d)]
        return f"{exp_name}  ({'; '.join(parts)})" if parts else exp_name
    exps = discover(mouse_dir, sec, paradigm, n_fov, stim_after_run, log=log)
    if not exps:
        log(f"{mouse_dir}: nothing to analyse (run aqua_batch.py and the MATLAB scripts first)")
        return None
    all_runs, all_cond, lines = [], [], []
    pending_overviews = []                      # drawn after the loop, on one rate axis for every experiment
    lines.append(f"AQuA2 event comparison — {mouse_dir.name}")
    lines.append("=" * 60)
    lines.append(f"pixel size {um_per_px} um, curation: "
                 + ("honoured, events " + event_source if use_curation else "IGNORED (--no-curation)"))
    lines.append("bad frames: frames dropped from the movie are never counted; "
                 + ("events next to excluded frames kept (--keep-flagged-events)" if keep_flagged else
                    f"events starting in an excluded frame (the flagged frames of artefact stretches, or what was "
                    f"decided in Clean & check) or next to one (+/- {flag_margin}) dropped, that time removed from the "
                    f"denominator")
                 + "; event frames (t0_frame) are RECORDED frames, t0_frame_movie the AQuA2 movie frame")
    lines.append(f"gating: {gate}"
                 + (f" (still_fraction >= {still_thresh})" if gate == "still" else "")
                 + (f", events kept with -log10(p) >= {min_p}" if min_p else ""))
    for k_exp, (exp_name, fov, para, runs) in enumerate(exps, 1):
        log(f"\n=== {exp_name} ({para}, FOV {fov}) ===  [{k_exp}/{len(exps)}]")
        ev_list, infos = [], []
        for r in runs:
            got = run_events(Path(r["tif"]), min_p=min_p, gate=gate,
                             still_thresh=still_thresh,
                             event_source=event_source if use_curation else "all",
                             flag_margin=flag_margin, keep_flagged=keep_flagged, log=log)
            if got is None:
                log(f"  {Path(r['tif']).name}: no AQuA2 feature table yet — skipped")
                continue
            ev, info = got
            info["condition"] = r["condition"]
            ev["condition"] = r["condition"]
            ev["experiment"] = exp_name
            ev_list.append(ev); infos.append(info)
            if info.get("excluded_frames") or info.get("flagged_frames") or info.get("state_frames"):
                log(f"      brightness QC: {info.get('flagged_frames', 0)} artefact frame(s), "
                    f"{info.get('state_frames', 0)} state-change frame(s); {info.get('excluded_frames', 0)} excluded"
                    + (f", {info['frames_dropped']} not in the movie" if info.get("frames_dropped") else "")
                    + (f" -> {info.get('n_events_dropped_flagged', 0)} event(s) starting in/next to them dropped"
                       if not keep_flagged else " (events kept: --keep-flagged-events)")
                    + (f"; {info['n_events_in_kept_state']} event(s) start in KEPT state stretches "
                       f"({info['kept_state_frames']} frames) — review in the GUI if in doubt"
                       if info.get("n_events_in_kept_state") else ""))
            elif not info.get("global_qc"):
                log("      (no *_frames.csv: run aqua_clean.py for global-brightness QC)")
            if info.get("curated"):
                log(f"      curation: {info['n_cfu_rejected']} CFU(s) marked noise, "
                    f"{info['n_cfu_added']} added"
                    + (f", {info.get('n_events_dropped_by_curation', 0)} events dropped"
                       if event_source != "all" else " (events not filtered: --event-source all)"))
            log(f"  [{r['condition']:>7}] {info['run']}: {info['n_events_used']} events"
                + (f" of {info['n_events_detected']}" if info['n_events_used'] != info['n_events_detected'] else "")
                + f" in {info['minutes_used']:.1f} min"
                + ("" if info["has_behavior"] else "  (no behaviour csv)"))
        if not ev_list:
            continue
        events = pd.concat(ev_list, ignore_index=True)
        events.to_csv(out_dir / f"events_{exp_name}.csv", index=False)
        rows = []
        for ev, info in zip(ev_list, infos):
            m = run_metrics(ev, info, um_per_px)
            m.update({"experiment": exp_name, "fov": fov, "paradigm": para,
                      "condition": info["condition"]})
            rows.append(m)
        all_runs.extend(rows)
        conditions = list(dict.fromkeys([r["condition"] for r in rows]))
        base = baseline if baseline in conditions else ("pre" if "pre" in conditions else conditions[0])
        conditions = [base] + [c for c in conditions if c != base]
        # ---- per condition
        lines.append("")
        lines.append(f"{exp_label(exp_name, runs)}  [{para}, FOV {fov}]   baseline = {base}")
        rs = pd.DataFrame(rows)
        for c in conditions:
            sub = rs[rs["condition"] == c]
            evc = events[events["condition"] == c]
            row = {"experiment": exp_name, "fov": fov, "paradigm": para, "condition": c,
                   "n_runs": len(sub), "minutes": round(float(sub["minutes_used"].sum()), 1),
                   "n_events": int(sub["n_events"].sum()),
                   "events_per_min": round(float(sub["n_events"].sum() / max(sub["minutes_used"].sum(), 1e-9)), 4),
                   "events_per_min_median_of_runs": round(float(sub["events_per_min"].median()), 4)}
            for p in PROPS:
                if p in evc.columns and len(evc):
                    v = evc[p].to_numpy(float); v = v[np.isfinite(v)]
                    row[f"{p}_median"] = round(float(np.median(v)), 4) if v.size else np.nan
                    row[f"{p}_median_of_runs"] = round(float(sub[f"{p}_median"].median()), 4)
            if "pupil_median" in sub.columns:
                row["pupil_median"] = round(float(sub["pupil_median"].median()), 2)
            if "velocity_median" in sub.columns:
                row["velocity_median"] = round(float(sub["velocity_median"].median()), 2)
            if "minutes_still" in sub.columns:
                for st in ("still", "running"):
                    mins = float(sub[f"minutes_{st}"].sum()); n_ev = int(sub[f"n_events_{st}"].sum())
                    row[f"minutes_{st}"] = round(mins, 2); row[f"n_events_{st}"] = n_ev
                    row[f"events_per_min_{st}"] = round(n_ev / mins, 4) if mins >= 0.5 else np.nan
                row["still_fraction"] = round(row["minutes_still"] / max(row["minutes_still"] + row["minutes_running"], 1e-9), 3)
            all_cond.append(row)
            lines.append(f"  {c:>10}: {row['n_events']:5d} events in {row['minutes']:6.1f} min"
                         f"  = {row['events_per_min']:7.3f}/min   peak dF/F {row.get('dff_max_median', float('nan')):.3f}"
                         f"   area {row.get('area_um2_median', float('nan')):.1f} um2"
                         f"   duration {row.get('duration_s_median', float('nan')):.1f} s")
        # ---- change vs baseline + event-level test
        try:
            from scipy.stats import mannwhitneyu
        except Exception:
            mannwhitneyu = None
        b = [r for r in all_cond if r["experiment"] == exp_name and r["condition"] == base][0]
        for c in conditions[1:]:
            o = [r for r in all_cond if r["experiment"] == exp_name and r["condition"] == c][0]
            ratio = o["events_per_min"] / b["events_per_min"] if b["events_per_min"] else np.nan
            txt = (f"  {c} vs {base}: rate x{ratio:.2f} "
                   f"({o['events_per_min']:.3f} vs {b['events_per_min']:.3f}/min)")
            if mannwhitneyu is not None and "dff_max" in events.columns:
                for p in ("dff_max", "area_um2", "duration_s"):
                    if p not in events.columns:
                        continue
                    x = events.loc[events["condition"] == c, p].to_numpy(float)
                    y = events.loc[events["condition"] == base, p].to_numpy(float)
                    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
                    if x.size > 5 and y.size > 5:
                        u, pv = mannwhitneyu(x, y, alternative="two-sided")
                        txt += f"; {p} {np.median(x):.3g} vs {np.median(y):.3g} (event-level p={pv:.1e})"
            lines.append(txt)
            sc_line = describe_change({"rate_baseline": b["events_per_min"], "rate_condition": o["events_per_min"],
                                       **change_measures(b.get("n_events"), b.get("minutes"), o.get("n_events"), o.get("minutes"))})
            if sc_line:
                lines.append(sc_line)
            # the same change with the locomotion state held fixed
            if "events_per_min_still" in o and "events_per_min_still" in b:
                parts = []
                for st in ("still", "running"):
                    rb, ro = b.get(f"events_per_min_{st}", np.nan), o.get(f"events_per_min_{st}", np.nan)
                    if np.isfinite(rb) and np.isfinite(ro):
                        parts.append(f"{st} x{ro / rb:.2f} ({ro:.3f} vs {rb:.3f}/min, "
                                     f"{o[f'minutes_{st}']:.1f} vs {b[f'minutes_{st}']:.1f} min)" if rb > 0
                                     else f"{st} {ro:.3f}/min vs 0")
                    else:
                        parts.append(f"{st}: not enough {st} time in both conditions "
                                     f"({o.get(f'minutes_{st}', 0):.1f} vs {b.get(f'minutes_{st}', 0):.1f} min)")
                lines.append(f"    at matched locomotion: " + "; ".join(parts))
                lines.append(f"    time spent still: {base} {100 * b.get('still_fraction', np.nan):.0f} %, "
                             f"{c} {100 * o.get('still_fraction', np.nan):.0f} %"
                             + (" — the conditions were recorded in different behavioural states; trust the matched "
                                "numbers (still / running, and the pupil terciles below) over the raw ratio"
                                if abs(b.get("still_fraction", 0) - o.get("still_fraction", 0)) > 0.25 else ""))
        lines.append("  (event-level p values treat events as independent — descriptive only; "
                     "the statistic that counts is one value per animal/FOV)")
        pending_overviews.append((exp_name, events, rows, conditions))
        hm_data = None                               # (cfu_list, ids, transforms, n_id) for the CFU heatmap
        # ---- CFUs
        if cfu:
            cfu_paths = [aqua_files(Path(r["tif"]))["cfu"] for r in runs]
            have = [p.exists() for p in cfu_paths]
            if all(have) and len(cfu_paths) > 1:
                log("  CFUs: reading and matching across runs")
                try:
                    cfu_list = [read_cfu(p) for p in cfu_paths]
                    curs = [read_curation(Path(r["tif"])) if use_curation else None for r in runs]
                    n_drop = 0
                    for c, cu in zip(cfu_list, curs):
                        if not cu or not cu.get("rejected"):
                            continue
                        keep = np.array([i + 1 not in set(cu["rejected"]) for i in range(c["n"])])
                        n_drop += int((~keep).sum())
                        c["maps"] = c["maps"][keep] if c["maps"].size else c["maps"]
                        c["events"] = [e for e, k in zip(c["events"], keep) if k]
                        c["occurrence"] = c["occurrence"][keep] if c["occurrence"].size else c["occurrence"]
                        c["dff"] = c["dff"][keep] if c["dff"].size else c["dff"]
                        c["freq"] = [f for f, k in zip(c["freq"], keep) if k]
                        c["n"] = int(keep.sum())
                    if n_drop:
                        log(f"  curation: {n_drop} CFU(s) marked noise are excluded from matching")
                    ids, n_matched, transforms, n_id = match_cfus(cfu_list, iou_thr=cfu_iou, log=log)
                    hm_data = (cfu_list, ids, transforms, n_id)
                    kept_ids = [set(np.asarray(e["event"]).astype(int).tolist())
                                if "event" in e.columns else None for e in ev_list]
                    rows_cfu, cfu_notes = cfu_table(exp_name, infos, cfu_list, ids, n_id,
                                                    um_per_px, conditions,
                                                    kept_event_ids=kept_ids if all(
                                                        k is not None for k in kept_ids) else None,
                                                    log=log)
                    manual_rows = manual_cfu_rows(exp_name, infos, runs, conditions,
                                                  use_curation=use_curation)
                    if manual_rows:
                        rows_cfu = rows_cfu + manual_rows
                        lines.append(f"  {len(manual_rows)} hand-added region(s) included "
                                     "(source = manual; events counted by threshold, not by AQuA2)")
                    pd.DataFrame(rows_cfu).to_csv(out_dir / f"cfu_matched_{exp_name}.csv", index=False)
                    cfu_figure(exp_name, infos, cfu_list, ids, out_dir / f"{exp_name}_cfu_map.png")
                    counts = [c["n"] for c in cfu_list]
                    moved = [f"run{i + 1} {t['angle_deg']:+.1f}deg dy{t['dy']:+d} dx{t['dx']:+d}"
                             for i, t in enumerate(transforms) if t.get("applied")]
                    lines.append(f"  CFUs per run: {counts}; matched in every run: {n_matched} "
                                 f"(IoU >= {cfu_iou}); total identities {n_id} "
                                 f"-> cfu_matched_{exp_name}.csv")
                    lines.append("  CFU map alignment: " + (", ".join(moved) if moved else
                                                            "no run needed shifting"))
                    lines += cfu_notes
                except Exception as e:
                    log(f"  CFUs: FAILED ({e})")
                    lines.append(f"  CFUs: could not be read ({e})")
            elif any(have):
                lines.append("  CFUs: only some runs have _res_cfu.mat — run run_cfu_batch.m again")
            else:
                lines.append("  CFUs: no _res_cfu.mat (run run_cfu_batch.m to get per-CFU comparisons)")
        # ---- pupil (arousal) vs activity, and pre/post at matched pupil
        if pupil_analysis:
            try:
                from afferent2p.pupil_coupling import run_series, experiment_pupil
                runs_info = []
                for ev, info in zip(ev_list, infos):
                    tif = Path(info["tif"])
                    fl = aqua_files(tif)
                    beh = read_behavior(fl["behavior"], log=log)
                    fr = read_frames_csv(Path(str(tif.with_suffix("")) + "_frames.csv"))
                    cfu_one = cfu_of(tif)
                    series = run_series(ev, info, beh, fr, cfu_one, info["seconds_per_frame"])
                    runs_info.append({"run": info["run"], "condition": info["condition"], "series": series})
                if any(r["series"] is not None for r in runs_info):
                    log("  pupil: coupling with activity and pre/post at matched pupil size")
                    prow, plines = experiment_pupil(exp_name, runs_info, conditions, out_dir, log=log)
                    by_run = {r["run"]: r for r in prow}
                    for row in rows:
                        pr = by_run.get(row["run"])
                        if pr:
                            for k in ("pupil_median_px", "pupil_range_pct", "pupil_speed_r", "pupil_dff_r_best",
                                      "pupil_dff_lag_s", "pupil_dff_p", "pupil_dff_partial_r", "pupil_dff_r_still",
                                      "pupil_rate_r0", "n_cfu_pupil_up", "n_cfu_pupil_down", "cfu_pupil_r_median",
                                      "rate_low_pupil", "rate_mid_pupil", "rate_high_pupil",
                                      "n_events_low_pupil", "n_events_mid_pupil", "n_events_high_pupil",
                                      "minutes_low_pupil", "minutes_mid_pupil", "minutes_high_pupil"):
                                if k in pr:
                                    row[k] = pr[k]
                    lines += plines
                    lines.append(f"  -> pupil_{exp_name}.csv, pupil_cfu_{exp_name}.csv, {exp_name}_pupil.png")
                else:
                    lines.append("  pupil: no pupil data (no DLC eye csv *eye*DLC*.csv next to the .sbx of any run; "
                                 "put it there and run the analysis again)")
            except Exception as e:
                log(f"  pupil analysis FAILED ({e})")
                lines.append(f"  pupil analysis failed ({e})")
        # ---- behaviour-adjusted rate ratio (all bins, running + pupil as covariates)
        cond_rows = {r["condition"]: r for r in all_cond if r["experiment"] == exp_name}
        for c in conditions[1:]:
            b_, o_ = cond_rows.get(base, {}), cond_rows.get(c, {})
            effects.append({"experiment": exp_name, "fov": fov, "paradigm": para, "baseline": base, "condition": c,
                            "n_runs_baseline": b_.get("n_runs"), "n_runs_condition": o_.get("n_runs"),
                            "minutes_baseline": b_.get("minutes"), "minutes_condition": o_.get("minutes"),
                            "rate_baseline": b_.get("events_per_min"), "rate_condition": o_.get("events_per_min"),
                            "ratio_raw": (o_.get("events_per_min", np.nan) / b_["events_per_min"]
                                          if b_.get("events_per_min") else np.nan),
                            "still_frac_baseline": b_.get("still_fraction"), "still_frac_condition": o_.get("still_fraction"),
                            "ratio_still": (o_.get("events_per_min_still", np.nan) / b_["events_per_min_still"]
                                            if b_.get("events_per_min_still") else np.nan),
                            "ratio_running": (o_.get("events_per_min_running", np.nan) / b_["events_per_min_running"]
                                              if b_.get("events_per_min_running") else np.nan),
                            "pupil_baseline": b_.get("pupil_median"), "pupil_condition": o_.get("pupil_median"),
                            **change_measures(b_.get("n_events"), b_.get("minutes"), o_.get("n_events"), o_.get("minutes"))})
        eff = {e["condition"]: e for e in effects if e["experiment"] == exp_name}
        # how the conditions were detected (aqua_tune's matched-sensitivity / common rule, v36): settings, recall and
        # false-alarm rate per condition, and the rates with the expected false alarms subtracted
        tj = read_tuning_json(mouse_dir, exp_name)
        if tj:
            for c in conditions[1:]:
                if c in eff:
                    d = detection_fields(tj, base, c, eff[c].get("rate_baseline"), eff[c].get("rate_condition"))
                    eff[c].update(d)
                    dl_ = describe_detection(d, base, c)
                    if dl_:
                        lines.append(dl_)
        # pupil-matched ratio: minutes-weighted mean of the log ratios over the terciles with data in both conditions
        rs_now = pd.DataFrame(rows)
        if "n_events_low_pupil" in rs_now.columns:
            for c in conditions[1:]:
                if c not in eff:
                    continue
                num = den = 0.0
                for b in ("low", "mid", "high"):
                    sb = rs_now[rs_now["condition"] == base]; sc = rs_now[rs_now["condition"] == c]
                    nb_, mb_ = float(sb[f"n_events_{b}_pupil"].sum()), float(sb[f"minutes_{b}_pupil"].sum())
                    nc_, mc_ = float(sc[f"n_events_{b}_pupil"].sum()), float(sc[f"minutes_{b}_pupil"].sum())
                    if mb_ > 0.5 and mc_ > 0.5 and nb_ > 0 and nc_ > 0:
                        w = min(mb_, mc_)
                        num += w * np.log((nc_ / mc_) / (nb_ / mb_)); den += w
                    eff[c][f"ratio_pupil_{b}"] = (nc_ / mc_) / (nb_ / mb_) if (mb_ > 0.5 and mc_ > 0.5 and nb_ > 0) else np.nan
                eff[c]["ratio_pupil_matched"] = float(np.exp(num / den)) if den > 0 else np.nan
        try:
            from afferent2p.adjust import adjusted_ratio, onset_triggered
            glm_runs, spf = [], infos[0]["seconds_per_frame"]
            for ev, info in zip(ev_list, infos):
                tif = Path(info["tif"])
                beh = read_behavior(aqua_files(tif)["behavior"])
                fr = read_frames_csv(Path(str(tif.with_suffix("")) + "_frames.csv"))
                T = int(info["frames"])
                onsets = np.zeros(T)
                if "t0_frame" in ev.columns and len(ev):
                    np.add.at(onsets, np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, T - 1), 1.0)
                excluded = np.zeros(T, bool)
                if fr is not None:
                    ex = fr["excluded"].astype(bool); n = min(T, ex.size); excluded[:n] = ex[:n]
                def col(name):
                    v = np.full(T, np.nan)
                    if beh is not None and name in beh.columns:
                        x = beh[name].to_numpy(float); n = min(T, x.size); v[:n] = x[:n]
                    return v
                glm_runs.append({"condition": info["condition"], "onsets": onsets, "excluded": excluded,
                                 "speed": np.abs(col("velocity_mean")), "still": col("still_fraction"),
                                 "pupil": col("pupil_diam_px_mean")})
            adj = adjusted_ratio(glm_runs, conditions, spf)
            if adj.get("covariates"):
                cov = [c for c in adj["covariates"] if not c.startswith("cond") and c != "intercept"]
                lines.append(f"  behaviour-adjusted (Poisson regression on {adj['n_bins']} bins of {adj['bin_s']:.0f} s, "
                             f"covariates {', '.join(cov) or 'none'}; behaviour explains "
                             f"{100 * adj.get('deviance_explained_behaviour', np.nan):.0f} % of the rate variation, "
                             f"the condition {100 * adj.get('deviance_explained_condition', np.nan):.0f} %):")
                for c in conditions[1:]:
                    ci = adj.get(f"adjusted_ratio_ci_{c}")
                    lines.append(f"    {c} vs {base}: adjusted rate ratio x{adj.get(f'adjusted_ratio_{c}', np.nan):.2f}"
                                 + (f" (95 % CI {ci[0]:.2f}-{ci[1]:.2f})" if ci else "")
                                 + f"; raw on the same bins x{adj.get(f'raw_ratio_{c}', np.nan):.2f}"
                                 + (f"; model rate still {adj.get(f'model_rate_still_{base}', np.nan):.2f} -> "
                                    f"{adj.get(f'model_rate_still_{c}', np.nan):.2f}/min, running "
                                    f"{adj.get(f'model_rate_running_{base}', np.nan):.2f} -> "
                                    f"{adj.get(f'model_rate_running_{c}', np.nan):.2f}/min"
                                    if f"model_rate_still_{c}" in adj else ""))
                    if c in eff:
                        eff[c].update({"ratio_adjusted": adj.get(f"adjusted_ratio_{c}"),
                                       "ratio_adjusted_ci_lo": ci[0] if ci else np.nan,
                                       "ratio_adjusted_ci_hi": ci[1] if ci else np.nan,
                                       "deviance_explained_behaviour": adj.get("deviance_explained_behaviour"),
                                       "deviance_explained_condition": adj.get("deviance_explained_condition"),
                                       "glm_bins": adj.get("n_bins")})
                # locomotion-onset-triggered rate per condition (pooled over its runs)
                trig = {}
                for c in conditions:
                    on = np.concatenate([r["onsets"] for r in glm_runs if r["condition"] == c] or [np.zeros(0)])
                    st = np.concatenate([r["still"] for r in glm_runs if r["condition"] == c] or [np.zeros(0)])
                    ex = np.concatenate([r["excluded"] for r in glm_runs if r["condition"] == c] or [np.zeros(0, bool)])
                    if on.size and np.isfinite(st).any():
                        t = onset_triggered(on, np.nan_to_num(st, nan=1.0), ex, spf)
                        if t:
                            trig[c] = t
                if trig:
                    for c, t in trig.items():
                        pre = np.nanmean(t["rate"][t["t_s"] < -5]); post = np.nanmean(t["rate"][(t["t_s"] > 0) & (t["t_s"] < 20)])
                        lines.append(f"    running onsets ({c}, n={t['n_onsets']}): {pre:.2f} events/min before -> "
                                     f"{post:.2f} in the first 20 s of running")
                    adj["onset_triggered"] = trig
                try:
                    adjusted_figure(exp_name, adj, conditions, glm_runs, spf, out_dir / f"{exp_name}_adjusted.png")
                except Exception as e:
                    log(f"  adjusted figure failed ({e})")
            else:
                lines.append("  behaviour-adjusted ratio: not computed — " + adj.get("why_not", "not enough bins with behaviour data"))
            # the same change at matched running speed, model-free, and how much of it can be matched at all
            try:
                from afferent2p.statematch import speed_standardised, describe as sm_describe, band_rows
                sm_res = speed_standardised(glm_runs, conditions, spf)
                for c in conditions[1:]:
                    lines.append(sm_describe(sm_res, base, c))
                    if c in eff and c in sm_res:
                        eff[c].update({"ratio_speed_matched": sm_res[c].get("ratio_std"),
                                       "behaviour_overlap_baseline": sm_res[c].get("overlap_baseline"),
                                       "behaviour_overlap_condition": sm_res[c].get("overlap_condition"),
                                       "behaviour_comparable": sm_res[c].get("comparable")})
                speed_rows += band_rows(exp_name, sm_res, conditions)
                tuning_results.append({"experiment": exp_name, "conditions": conditions, "res": sm_res,
                                       "label": exp_label(exp_name, runs)})
            except Exception as e:
                log(f"  matched-speed comparison FAILED ({e})")
            # keep what the day time course needs (optostim: both FOVs of a day on one clock)
            tc_runs = []
            for g, info, m in zip(glm_runs, infos, rows):
                tc_runs.append({"run": info["run"], "condition": info["condition"], "tif": info["tif"], "info": info,
                                "onsets": g["onsets"], "excluded": g["excluded"], "speed": g["speed"], "still": g["still"],
                                "pupil": g["pupil"], "n_events": int(info.get("n_events_used", 0)),
                                "minutes_used": float(info.get("minutes_used", 0.0)),
                                "props": dict({k: m.get(k) for k in ("velocity_median", "pupil_median", "dff_max_median",
                                                                     "duration_s_median", "area_um2_median") if k in m},
                                              still_fraction=(round(m["minutes_still"] / (m["minutes_still"] + m["minutes_running"]), 3)
                                                              if m.get("minutes_still") is not None and
                                                              (m.get("minutes_still", 0) + m.get("minutes_running", 0)) > 0 else np.nan))})
            tc_store.append({"experiment": exp_name, "fov": fov, "para": para, "conditions": conditions, "base": base,
                             "spf": spf, "runs": tc_runs})
        except Exception as e:
            log(f"  behaviour adjustment FAILED ({e})")
            lines.append(f"  behaviour adjustment failed ({e})")
        # ---- where the activity is: common region, fibre atlas, tile maps, fibre segments
        if maps:
            try:
                from afferent2p.spatial import read_events, read_meanimg, experiment_maps, empty_events
                sp_runs = []
                for ev, info in zip(ev_list, infos):
                    tif = Path(info["tif"])
                    loc = read_events(tif)
                    if loc is None and len(ev) == 0:
                        loc = empty_events()                                   # nothing to locate
                    if loc is not None and "event" in ev.columns and "event" in loc.columns:
                        loc = loc[loc["event"].isin(ev["event"])]           # the kept events only
                    fallback = None
                    if not (Path(str(tif.with_suffix("")) + AQUA_SUFFIX + "_meanimg.tif")).exists() and aqua_files(tif)["cfu"].exists():
                        try:
                            fallback = read_cfu(aqua_files(tif)["cfu"])["datPro"]
                        except Exception:
                            fallback = None
                    sp_runs.append({"run": info["run"], "condition": info["condition"], "tif": tif, "events": loc,
                                    "minutes_used": info.get("minutes_used", 0.0),
                                    "meanimg": read_meanimg(tif, fallback=fallback)})
                if all(r["events"] is not None for r in sp_runs):
                    log("  maps: aligning runs, fibre atlas, tile maps, fibre segments")
                    bfac = {c: (eff[c]["ratio_raw"] / eff[c]["ratio_adjusted"])
                            for c in conditions[1:] if c in eff and eff[c].get("ratio_adjusted") and eff[c].get("ratio_raw")}
                    sm = experiment_maps(exp_name, sp_runs, conditions, um_per_px, out_dir, tile_um=tile_um,
                                         behaviour_factor=bfac, log=log)
                    if sm:
                        lines.append(f"  where: {sm['runs_aligned']} of {sm['n_runs'] - 1} runs match the first "
                                     f"({sm['runs_moved']} needed a shift/rotation"
                                     + (f"; {sm['runs_not_aligned']} do NOT match — check {exp_name}_maps.png" if sm["runs_not_aligned"] else "")
                                     + ")"
                                     + f"; common region {100 * sm['common_region_frac']:.0f} % of the field, fibres "
                                     f"{100 * sm['fibre_frac_of_region']:.0f} % of it ({sm['fibre_area_um2']:.0f} um2, "
                                     f"{sm['n_segments']} segments); {sm['n_valid_tiles']} fibre tiles of {sm['tile_um']:.0f} um")
                        for c in conditions[1:]:
                            cm = sm["comparisons"].get(c, {})
                            sg = sm["segments"]
                            B = cm.get("behaviour_factor", 1.0)
                            beh_txt = (f" (behaviour accounts for x{B:.2f}, so x{cm.get('region_ratio_beyond_behaviour', np.nan):.2f} beyond it)"
                                       if abs(B - 1) > 1e-6 else "")
                            lines.append(f"    {c} vs {base} in the common region: rate x{cm.get('region_ratio', np.nan):.2f} "
                                         f"({sm[f'rate_region_{c}']:.3f} vs {sm[f'rate_region_{base}']:.3f}/min; per mm2 of fibre "
                                         f"{sm[f'rate_per_fibre_mm2_{c}']:.1f} vs {sm[f'rate_per_fibre_mm2_{base}']:.1f}){beh_txt}")
                            lines.append(f"    tile maps ({base} vs {c} with the behaviour factor taken out): r {cm.get('map_r', np.nan):+.2f}; "
                                         f"active tiles {cm.get('active_tiles_a', 0)} -> {cm.get('active_tiles_b', 0)}, "
                                         f"{100 * cm.get('frac_switched', np.nan):.0f} % of the active tiles are active in one condition only "
                                         f"({100 * cm.get('frac_new', np.nan):.0f} % new, {100 * cm.get('frac_lost', np.nan):.0f} % lost); "
                                         f"tile-wise log2 ratio {cm.get('tile_log2_ratio_mean', np.nan):+.2f}"
                                         + (f" (spatial bootstrap CI {cm['tile_log2_ratio_ci'][0]:+.2f}..{cm['tile_log2_ratio_ci'][1]:+.2f})"
                                            if cm.get("tile_log2_ratio_ci") else ""))
                            if sg.get("segments_used") is not None:
                                med = sg.get("segments_log2_ratio_median", np.nan)
                                lines.append(f"    fibre segments: {sg['segments_active_base']} active in {base}, "
                                             f"{sg['segments_active_cond']} in {c}; of {sg['segments_used']} with >= 3 events "
                                             f"{sg['segments_up']} up, {sg['segments_down']} down beyond the behaviour factor "
                                             f"(binomial p < 0.05)"
                                             + (f"; median log2 ratio {med:+.2f}" if np.isfinite(med) else "")
                                             + (f", Wilcoxon p = {sg['segments_wilcoxon_p']:.2g}" if np.isfinite(sg.get("segments_wilcoxon_p", np.nan)) else ""))
                            if c in eff:
                                eff[c].update({"ratio_common_region": cm.get("region_ratio"),
                                               "ratio_common_region_beyond_behaviour": cm.get("region_ratio_beyond_behaviour"),
                                               "map_r": cm.get("map_r"),
                                               "frac_tiles_switched": cm.get("frac_switched"),
                                               "tile_log2_ratio": cm.get("tile_log2_ratio_mean"),
                                               "segments_up": sg.get("segments_up"), "segments_down": sg.get("segments_down"),
                                               "segments_used": sg.get("segments_used"),
                                               "segments_log2_ratio_median": sg.get("segments_log2_ratio_median"),
                                               "runs_not_aligned": sm["runs_not_aligned"],
                                               "common_region_frac": sm["common_region_frac"]})
                        lines.append(f"  -> tiles_{exp_name}.csv, segments_{exp_name}.csv, {exp_name}_maps.png")
                else:
                    lines.append("  where: no per-event locations (run run_cfu_batch.m from v27 — it writes "
                                 "<run>_AQuA2_events.csv / _meanimg.tif)")
            except Exception as e:
                log(f"  maps FAILED ({e})")
                lines.append(f"  maps failed ({e})")
        # ---- v33: how much each run moved, what AQuA2's registration left, do the events follow the motion
        if motion:
            try:
                from afferent2p import motion as MO
                for info in infos:
                    tif = Path(info["tif"])
                    summ = MO.run_motion(tif, log=log)
                    if summ:
                        e_ = summ.get("events", {})
                        motion_rows.append({"experiment": exp_name, "run": summ["run"], "condition": info.get("condition", ""),
                                            "verdict": summ["verdict"][0], "why": summ["verdict"][1],
                                            "shift_median_px": summ["shift_median"], "shift_p95_px": summ["shift_p95"], "shift_max_px": summ["shift_max"],
                                            "jumps_ge2px_per_min": summ["jumps_ge2_per_min"], "jumps_per_min_still": summ.get("jumps_ge2_per_min_still"),
                                            "jumps_per_min_running": summ.get("jumps_ge2_per_min_run"), "running_pct": summ["running_pct"],
                                            "corr_shift_speed": summ.get("corr_shift_speed"), "sharp_drop_running_pct": summ.get("sharp_drop_run_pct"),
                                            "registered_movie": summ.get("registered", {}).get("available", False),
                                            "residual_rigid_p95_px": summ.get("residual_rigid_p95"), "deformation_median_px": summ.get("nonrigid_median"),
                                            "deformation_jumps_per_min": summ.get("nonrigid_jumps_per_min"),
                                            "n_event_onsets": e_.get("n_onsets"), "onsets_near_motion": e_.get("near_motion"),
                                            "onsets_near_motion_expected": e_.get("expected_near_motion"), "onsets_near_motion_ratio": e_.get("ratio_motion")})
                n_meas = sum(1 for r in motion_rows if r["experiment"] == exp_name)
                n_bad = sum(1 for r in motion_rows if r["experiment"] == exp_name and r["verdict"] != "ok")
                lines.append(f"  motion: {n_meas} run(s) measured" + (f", {n_bad} to check (motion_summary.csv, QC report)" if n_bad
                             else (", nothing to flag" if n_meas else " (no readable movie)")))
            except Exception as e:
                log(f"  motion check FAILED ({e})")
                lines.append(f"  motion check failed ({e})")
        # ---- v32: per-run bundles -> locomotion trials, rates by state, ECDFs, pupil measures, CFU heatmap
        try:
            from afferent2p import runbundle
            bundles = []
            for ev, info in zip(ev_list, infos):
                tif = Path(info["tif"])
                beh = read_behavior(aqua_files(tif)["behavior"], log=lambda *a: None)
                fr = read_frames_csv(Path(str(tif.with_suffix("")) + "_frames.csv"))
                b = runbundle.build(ev, info, beh, fr, cfu_of(tif), info["seconds_per_frame"], still_thresh=still_thresh)
                bundles.append(b)
            rank = {c: i for i, c in enumerate(conditions)}
            perm = sorted(range(len(bundles)), key=lambda i: (rank.get(bundles[i]["condition"], 99), bundles[i]["run"]))
            ordered = [bundles[i] for i in perm]
            exp_bundles.append({"experiment": exp_name, "fov": fov, "para": para, "conditions": conditions, "base": base,
                                "bundles": ordered, "day": day_of(runs[0]["run"])})
            exp_meta.append({"experiment": exp_name, "fov": fov, "para": para, "conditions": conditions, "base": base,
                             "day": day_of(runs[0]["run"])})
            lres = None
            if loco:
                try:
                    from afferent2p import locotrials
                    log("  locomotion trials: bouts, trial responses, per-CFU responsiveness, rest-only activity"
                        + (", GLMs" if loco_glm else ""))
                    lres = locotrials.analyze_experiment(exp_name, ordered, conditions, out_dir, glm=loco_glm, log=log)
                    loco_run_rows += lres["run_rows"]; loco_cond_rows += lres["cond_rows"]; loco_trial_rows += lres["trial_rows"]
                    for r in lres["trial_rows"]:
                        trials_by_run.setdefault(r["run"], []).append((int(r["onset"]), r["kind"]))
                    lines += lres["lines"]
                    lines.append(f"  -> loco_trials_{exp_name}.csv, loco_cfu_{exp_name}.csv, loco_glm_{exp_name}.csv, {exp_name}_loco.png, {exp_name}_loco_cfu.png")
                except Exception as e:
                    log(f"  locomotion trials FAILED ({e})")
                    lines.append(f"  locomotion trials failed ({e})")
                # v34: the three triggers (running onsets, still pupil dilations, rest) compared on the same statistics
                try:
                    from afferent2p import triggers
                    log("  triggers: pupil-dilation trials, three-trigger comparison (running / pupil / rest)")
                    tres = triggers.analyze_experiment(exp_name, ordered, conditions, lres, out_dir, log=log)
                    trigger_rows += tres["rows"]
                    lines += tres["lines"]
                    lines.append(f"  -> pupil_trials_{exp_name}.csv, {exp_name}_triggers.png")
                except Exception as e:
                    log(f"  trigger comparison FAILED ({e})")
                    lines.append(f"  trigger comparison failed ({e})")
            try:
                from afferent2p import general
                srows = general.state_rows(exp_name, ordered, still_thresh=still_thresh)
                state_rows_all += srows
                general.state_figure(exp_name, srows, ordered, conditions, out_dir / f"{exp_name}_state.png", log=log)
                general.ecdf_figure(exp_name, ordered, conditions, out_dir / f"{exp_name}_ecdf.png")
                pupil_rows_all += general.pupil_rows(exp_name, ordered, para, trial_rows=(lres or {}).get("trial_rows"))
                lines.append(f"  -> {exp_name}_state.png (rates by running / resting, event overlays), {exp_name}_ecdf.png")
            except Exception as e:
                log(f"  state / ECDF / pupil measures FAILED ({e})")
                lines.append(f"  general representation failed ({e})")
            if heatmaps:
                try:
                    from afferent2p import heatmaps as HM
                    if hm_data is None:
                        one = [(i, cfu_of(Path(info["tif"]))) for i, info in enumerate(infos)]
                        one = [(i, c) for i, c in one if c is not None and c["n"] > 0]
                        if len(one) == 1:
                            i, c = one[0]
                            hm_ids = [np.full(cfu_of(Path(inf["tif"]))["n"] if cfu_of(Path(inf["tif"])) is not None else 0, -1, int) for inf in infos]
                            hm_ids[i] = np.arange(c["n"])
                            hm_data = ([cfu_of(Path(inf["tif"])) or {"maps": np.zeros((0, 0, 0)), "n": 0} for inf in infos], hm_ids,
                                       [{"applied": False}] * len(infos), int(c["n"]))
                    if hm_data is not None:
                        cfu_list_, ids_, transforms_, n_id_ = hm_data
                        log("  heatmap: CFU identities x time" + (" from the movies (AQuA2-registered where readable)" if movies else ""))
                        hres = HM.analyze_experiment(exp_name, ordered, conditions, [cfu_list_[i] for i in perm], [ids_[i] for i in perm],
                                                     [transforms_[i] for i in perm], n_id_, out_dir, movies=movies, log=log)
                        lines += hres["lines"]
                        if hres.get("summary"):
                            heat_rows.append(hres["summary"])
                    else:
                        lines.append("  heatmap: no CFUs (run run_cfu_batch.m)")
                except Exception as e:
                    log(f"  heatmap FAILED ({e})")
                    lines.append(f"  heatmap failed ({e})")
        except Exception as e:
            log(f"  v32 readouts FAILED ({e})")
            lines.append(f"  general representation / locomotion / heatmap readouts failed ({e})")
    # ---- v32: session overviews on the day's clock, acute / chronic summary
    if exp_bundles:
        try:
            from afferent2p import general
            srows = general.session_figures(exp_bundles, mouse_dir, out_dir, day_labels, trials_by_run,
                                            stim_gap_min=assumed_stim_gap_min, log=log)
            if srows:
                pd.DataFrame(srows).to_csv(out_dir / "session_figures.csv", index=False)
        except Exception as e:
            log(f"  session figures FAILED ({e})")
        if loco_run_rows:
            pd.DataFrame(loco_run_rows).to_csv(out_dir / "loco_summary.csv", index=False)
        if loco_cond_rows:
            pd.DataFrame(loco_cond_rows).to_csv(out_dir / "loco_conditions.csv", index=False)
        if trigger_rows:
            pd.DataFrame(trigger_rows).to_csv(out_dir / "trigger_comparison.csv", index=False)
            n_t = sum(1 for r in trigger_rows if r.get("testable")); n_nt = len(trigger_rows) - n_t
            lines += ["", f"Three-trigger comparison (v34): {len(trigger_rows)} comparison(s) in trigger_comparison.csv — {n_t} testable, "
                      f"{n_nt} not (why in the file / report); per experiment <exp>_triggers.png, pupil_trials_<exp>.csv"]
        if state_rows_all:
            pd.DataFrame(state_rows_all).to_csv(out_dir / "state_rates.csv", index=False)
        if pupil_rows_all:
            pd.DataFrame(pupil_rows_all).to_csv(out_dir / "pupil_measures.csv", index=False)
        if heat_rows:
            pd.DataFrame(heat_rows).to_csv(out_dir / "heatmap_summary.csv", index=False)
        try:
            from afferent2p import general
            general.acute_chronic(mouse_dir.name, exp_meta, all_cond, loco_cond_rows, pupil_rows_all, day_labels, out_dir, log=log)
            lines += ["", f"General representation (v32): per experiment <exp>_state.png / _ecdf.png / _loco.png / _heatmap.png / _updown.png, "
                      f"per day {mouse_dir.name}_session_<day>.png; acute (post/pre) and chronic (pre day#/pre day 1) change of every "
                      f"measure in acute_chronic.csv and {mouse_dir.name}_acute_chronic.png"]
        except Exception as e:
            log(f"  acute/chronic summary FAILED ({e})")
    if motion_rows:
        pd.DataFrame(motion_rows).to_csv(out_dir / "motion_summary.csv", index=False)
        n_warn = sum(1 for r in motion_rows if r["verdict"] == "warn"); n_check = sum(1 for r in motion_rows if r["verdict"] == "check")
        lines += ["", f"Motion / registration check (v33): {len(motion_rows)} run(s) in motion_summary.csv — {n_warn} warn, {n_check} check, "
                  f"{len(motion_rows) - n_warn - n_check} ok (per run: <run>_motion.csv / .json next to the tiff; figures in the QC report)"]
    if speed_rows:
        pd.DataFrame(speed_rows).to_csv(out_dir / "speed_bands.csv", index=False)
        try:
            from afferent2p.statematch import tuning_figure
            tuning_figure(tuning_results, out_dir / f"{mouse_dir.name}_speed_tuning.png", log=log)
        except Exception as e:
            log(f"  speed-tuning figure failed ({e})")
    if timecourse and any(e["para"] == "optostim" for e in tc_store):
        try:
            from afferent2p.timecourse import day_timecourse, timecourse_figure
            by_day = {}
            for e in tc_store:
                if e["para"] != "optostim" or not e["runs"]:
                    continue
                by_day.setdefault(day_of(e["runs"][0]["run"]) or e["experiment"], []).append(e)
            results, trows, brows, srows, sums, clock_rows = [], [], [], [], [], []
            lines += ["", "Time since the stimulation (optostim days; both fields of a day on one clock):"]
            for day in sorted(by_day):
                res = day_timecourse(day, by_day[day], out_dir, mouse_dir, bin_min=bin_min,
                                     assumed_stim_gap_min=assumed_stim_gap_min, label=day_labels.get(day), log=log)
                if res is None:
                    lines.append(f"  day {day}: no pre/post pair (nothing to place on a clock)")
                    continue
                results.append(res); trows += res["run_rows"]; brows += res["bin_rows"]; srows += res["steps"]
                sums.append(res["summary"]); lines += res["lines"]
            if results:
                pd.DataFrame(trows).to_csv(out_dir / "timecourse_runs.csv", index=False)
                pd.DataFrame(brows).to_csv(out_dir / "timecourse_bins.csv", index=False)
                pd.DataFrame(srows).to_csv(out_dir / "timecourse_steps.csv", index=False)
                pd.DataFrame([{k: v for k, v in r.items() if k != "slope_ci"} for r in sums]).to_csv(
                    out_dir / "timecourse_summary.csv", index=False)
                if not any(r["estimated"] for r in results):
                    pd.DataFrame([{"run": r["run"], "start": r["start_utc"], "end": None} for r in trows]).to_csv(
                        out_dir / "run_times_found.csv", index=False)
                try:
                    timecourse_figure(results, mouse_dir.name, out_dir / f"{mouse_dir.name}_timecourse.png", labels=day_labels)
                    lines.append(f"  -> timecourse_runs.csv, timecourse_bins.csv, timecourse_steps.csv, timecourse_summary.csv, "
                                 f"{mouse_dir.name}_timecourse.png")
                except Exception as e:
                    log(f"  time-course figure failed ({e})")
                lines.append("  (the clock comes from the recording files; the stimulation itself is not timed — "
                             "'minutes since the first post recording' is the axis, add the stimulation-to-recording delay yourself)")
        except Exception as e:
            log(f"  time course FAILED ({e})")
            lines.append(f"  time course failed ({e})")
    # per-experiment overview figures, their rate panel on one axis for the whole animal
    rate_ymax = None
    if all_runs:
        v = np.array([r.get("events_per_min", np.nan) for r in all_runs], float)
        rate_ymax = float(np.nanmax(v)) * 1.1 if np.isfinite(v).any() and np.nanmax(v) > 0 else None
    for exp_name, events, rows, conditions in pending_overviews:
        try:
            overview_figure(exp_name, events, rows, conditions, out_dir / f"{exp_name}_overview.png", log=log, rate_ymax=rate_ymax)
        except Exception as e:
            log(f"  overview figure of {exp_name} failed ({e})")
    if effects:
        pd.DataFrame(effects).to_csv(out_dir / "experiments_summary.csv", index=False)
        try:
            effects_figure(effects, out_dir / f"{mouse_dir.name}_effects.png", day_labels=day_labels)
            try:
                rates_figure(effects, all_cond, all_runs, out_dir / f"{mouse_dir.name}_rates.png", day_labels=day_labels)
            except Exception as e:
                log(f"  rates figure failed ({e})")
        except Exception as e:
            log(f"  effects figure failed ({e})")
        lines += ["", "Per experiment, the change vs its baseline with the behaviour controls (experiments_summary.csv, "
                  f"{mouse_dir.name}_effects.png): raw, at matched locomotion (still / running frames), "
                  "behaviour-adjusted (regression), in the common region, per fibre segment."]
    if all_runs:
        pd.DataFrame(all_runs).to_csv(out_dir / "run_summary.csv", index=False)
    if all_cond:
        pd.DataFrame(all_cond).to_csv(out_dir / "condition_summary.csv", index=False)
    lines += ["", "How to take this to a group statistic:",
              "  1. run_summary.csv is one row per recording; average the runs of a condition",
              "     within a FOV, then take ONE value per animal (or per FOV if you treat FOVs",
              "     as independent, and say so).",
              "  2. compare the baseline->condition change between treated and control mice",
              "     (unpaired) or pre->post within mice (paired) on those per-animal values.",
              "  3. keep the AQuA2 settings identical within a FOV (aqua_batch.py enforces this)",
              "     and report them; event rates are threshold-dependent.",
              "  4. mCitrine bleaching lowers F0 over a session: prefer event rate and peak dF",
              "     over peak dF/F when comparing pre with post."]
    (out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log("\n".join(lines[-14:]))
    log(f"\nWrote {out_dir}")
    return out_dir


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+", help="mouse folders")
    ap.add_argument("--um-per-px", type=float, default=0.792)
    ap.add_argument("--sec", type=float, default=2.0, help="which binned tiff (s per frame)")
    ap.add_argument("--gate", default="none", choices=["none", "still"],
                    help="still: keep only events starting in still frames (needs the behaviour csv)")
    ap.add_argument("--still-thresh", type=float, default=0.5,
                    help="still_fraction of a binned frame to call it still (default 0.5)")
    ap.add_argument("--min-p", type=float, default=0.0,
                    help="keep events with -log10(p on max dF/F) >= this (0 = keep all)")
    ap.add_argument("--baseline", default=None, help="condition to compare against (default: pre / earliest day)")
    ap.add_argument("--no-cfu", action="store_true")
    ap.add_argument("--event-source", default="drop-rejected-cfu",
                    choices=["all", "drop-rejected-cfu", "kept-cfu-only"],
                    help="how manual CFU curation affects the EVENT statistics: "
                         "all = ignore it, drop-rejected-cfu = drop events of CFUs marked noise "
                         "(default), kept-cfu-only = keep only events inside kept CFUs")
    ap.add_argument("--no-curation", action="store_true",
                    help="ignore the *_res_cfu_curated.json sidecars entirely")
    ap.add_argument("--keep-flagged-events", action="store_true",
                    help="keep events that start in (or right next to) frames the brightness step excluded; "
                         "frames dropped from the movie are never counted")
    ap.add_argument("--flag-margin", type=int, default=1,
                    help="frames on each side of a flagged frame that also count as bad (default 1)")
    ap.add_argument("--no-pupil", action="store_true", help="skip the pupil-versus-activity analysis")
    ap.add_argument("--no-maps", action="store_true", help="skip the common-region / fibre-atlas / tile maps")
    ap.add_argument("--tile-um", type=float, default=25.0, help="tile size of the activity maps (default 25 um)")
    ap.add_argument("--no-timecourse", action="store_true",
                    help="skip the day time course (rate vs time since the first post run, optostim days)")
    ap.add_argument("--bin-min", type=float, default=2.0, help="bin length in minutes for the time course (default 2)")
    ap.add_argument("--day-label", action="append", default=[], metavar="YYMMDD=LABEL",
                    help="name a day's condition, e.g. 260917=sham (repeatable; or put day,label rows in <mouse>\\day_labels.csv)")
    ap.add_argument("--assumed-stim-gap", type=float, default=73.0,
                    help="minutes between the last pre run and the first post run when no recording clock is found "
                         "(.sbx / _quadrature.mat / .mat next to the movies, or run_times.csv); default 73")
    ap.add_argument("--no-loco", action="store_true", help="skip the locomotion-trial readouts (Yu's pipeline)")
    ap.add_argument("--no-glm", action="store_true", help="skip the Loco2Fluor / LocoPupil2Fluor GLMs (the slowest part)")
    ap.add_argument("--no-heatmaps", action="store_true", help="skip the CFU x time heatmaps")
    ap.add_argument("--no-motion", action="store_true", help="skip the per-run motion / registration check (afferent2p.motion)")
    ap.add_argument("--no-movies", action="store_true",
                    help="heatmap rows from AQuA2's CFU dF/F only (do not read the movies under the footprints)")
    ap.add_argument("--cfu-iou", type=float, default=0.3)
    ap.add_argument("--paradigm", default="auto", choices=["auto", "spontaneous", "optostim"])
    ap.add_argument("--n-fov", type=int, default=None)
    ap.add_argument("--stim-after-run", type=int, default=None)
    a = ap.parse_args()
    for f in a.folders:
        analyze_mouse(Path(f), um_per_px=a.um_per_px, sec=a.sec, gate=a.gate,
                      still_thresh=a.still_thresh, min_p=a.min_p, baseline=a.baseline,
                      cfu=not a.no_cfu, cfu_iou=a.cfu_iou, paradigm=a.paradigm,
                      n_fov=a.n_fov, stim_after_run=a.stim_after_run,
                      event_source=a.event_source, use_curation=not a.no_curation,
                      flag_margin=a.flag_margin, keep_flagged=a.keep_flagged_events,
                      pupil_analysis=not a.no_pupil, maps=not a.no_maps, tile_um=a.tile_um,
                      timecourse=not a.no_timecourse, bin_min=a.bin_min, assumed_stim_gap_min=a.assumed_stim_gap,
                      day_labels={k.strip(): v.strip() for k, v in (x.split("=", 1) for x in a.day_label if "=" in x)},
                      loco=not a.no_loco, loco_glm=not a.no_glm, heatmaps=not a.no_heatmaps, movies=not a.no_movies,
                      motion=not a.no_motion)


if __name__ == "__main__":
    main()
