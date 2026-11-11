"""One run's traces on the RECORDED-frame grid, in one dict, for the readouts
that need the behaviour, the events and the CFU curves together (locomotion
trials, the general representation, the heatmaps).

    b = build(ev, info, beh, frames, cfu, spf)
    b["speed"], b["still"], b["running"], b["pupil"], b["excluded"]   per recorded frame
    b["onsets"]        kept-event onsets per frame (count)
    b["events"]        the kept-event table (t0_frame = recorded frame, 1-based)
    b["cfu_dff"]       (nCFU, T) CFU dF/F on the recorded grid (NaN where the movie has no frame)
    b["cfu_onsets"]    (nCFU, T) kept-event onsets per CFU

Frames dropped from the movie (unpack / aqua_clean) are NaN in the curves;
excluded frames (Clean & check, tab 2) are kept in the arrays and flagged in b["excluded"]
so every readout can drop what overlaps them.
"""
from __future__ import annotations

import numpy as np

SPEED_COL, STILL_COL, PUPIL_COL, VALID_COL = "velocity_mean", "still_fraction", "pupil_diam_px_mean", "pupil_valid_fraction"


def build(ev, info: dict, beh, frames: dict | None, cfu: dict | None, spf: float,
          still_thresh: float = 0.5, min_valid_fraction: float = 0.5) -> dict:
    T = int(info.get("frames") or (len(beh) if beh is not None else 0) or 0)
    if T <= 0 and cfu is not None and np.asarray(cfu.get("dff", [])).size:
        T = int(np.asarray(cfu["dff"]).shape[1])
    T = max(T, 1)
    n = min(T, len(beh)) if beh is not None else 0

    def col(name, default=np.nan):
        v = np.full(T, default, float)
        if beh is not None and name in beh.columns:
            v[:n] = beh[name].to_numpy(float)[:n]
        return v
    speed = np.abs(col(SPEED_COL))
    still = col(STILL_COL)
    pupil = col(PUPIL_COL)
    if beh is not None and VALID_COL in beh.columns:
        vf = col(VALID_COL, 0.0)
        pupil[vf < min_valid_fraction] = np.nan
    excluded = np.zeros(T, bool)
    dropped = np.zeros(T, bool)
    if frames is not None:
        ex = np.asarray(frames["excluded"], bool); m = min(T, ex.size); excluded[:m] = ex[:m]
        dr = np.asarray(frames["dropped"], bool); m = min(T, dr.size); dropped[:m] = dr[:m]
    has_loco = bool(np.isfinite(still).any())
    running = np.where(np.isfinite(still), still < still_thresh, False)
    # kept-event onsets
    onsets = np.zeros(T)
    t0 = None
    if ev is not None and "t0_frame" in ev.columns and len(ev):
        t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, T - 1)
        np.add.at(onsets, t0, 1.0)
    # CFU curves: movie frames -> recorded frames
    cfu_dff = cfu_occ = cfu_onsets = None
    n_cfu = 0
    cfu_ids = np.zeros(0, int)
    if cfu is not None and np.asarray(cfu.get("dff", [])).size:
        dff = np.asarray(cfu["dff"], float)
        n_cfu, Tm = dff.shape
        cfu_ids = np.asarray(cfu.get("orig_index", np.arange(n_cfu) + 1), int)
        rec_of_movie = None
        if frames is not None and np.asarray(frames["dropped"]).any():
            present = np.flatnonzero(~np.asarray(frames["dropped"], bool))
            if present.size == Tm:
                rec_of_movie = np.asarray(frames["frame"])[present] - 1
        if rec_of_movie is None:
            rec_of_movie = np.arange(min(Tm, T))
        rec_of_movie = rec_of_movie[rec_of_movie < T]
        cfu_dff = np.full((n_cfu, T), np.nan)
        cfu_dff[:, rec_of_movie] = dff[:, :rec_of_movie.size]
        occ = np.asarray(cfu.get("occurrence", np.zeros((n_cfu, 0))), bool)
        cfu_occ = np.zeros((n_cfu, T), bool)
        if occ.ndim == 2 and occ.shape[0] == n_cfu and occ.shape[1] >= rec_of_movie.size:
            cfu_occ[:, rec_of_movie] = occ[:, :rec_of_movie.size]
        cfu_onsets = np.zeros((n_cfu, T))
        if t0 is not None and "event" in ev.columns:
            ev_ids = ev["event"].to_numpy(int)
            where = {int(e): int(t) for e, t in zip(ev_ids, t0)}
            for k in range(n_cfu):
                for e in np.asarray(cfu["events"][k]).ravel().astype(int):
                    if int(e) in where:
                        cfu_onsets[k, where[int(e)]] += 1.0
    return {"run": info.get("run"), "condition": info.get("condition"), "tif": info.get("tif"), "T": T, "spf": float(spf),
            "speed": speed, "still": still, "running": running, "pupil": pupil, "excluded": excluded, "dropped": dropped,
            "onsets": onsets, "events": ev, "cfu_dff": cfu_dff, "cfu_occ": cfu_occ, "cfu_onsets": cfu_onsets,
            "n_cfu": n_cfu, "cfu_ids": cfu_ids, "has_loco": has_loco, "has_pupil": bool(np.isfinite(pupil).any()),
            "minutes_used": float(info.get("minutes_used", 0.0)), "info": info}


def population_dff(b: dict) -> np.ndarray:
    """Mean dF/F over the kept CFUs (NaN where no CFU / no frame)."""
    if b.get("cfu_dff") is None or b["n_cfu"] == 0:
        return np.full(b["T"], np.nan)
    with np.errstate(all="ignore"):
        return np.nanmean(b["cfu_dff"], axis=0)


def robust_sd(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return np.nan
    return float(1.4826 * np.median(np.abs(x - np.median(x))))
