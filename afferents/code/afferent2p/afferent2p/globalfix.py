"""Global-brightness correction, bad-frame flagging and per-stretch decisions
for the binned AQuA2 input.

Why: a whole-field brightness change (laser/PMT step, a hemodynamic swing, the
brain shifting when the animal brakes, stray light) is brighter than any real
afferent event and, being everywhere at once, is carved by AQuA2 into hundreds
of "events" in a few frames. In one real recording 227 of 291 events started
within 9 frames of such a step. No detection threshold can separate it from
real events, so it is handled before AQuA2 sees the movie, identically for
every recording:

  1. level trace   g(t) = spatial median of frame t. Sparse afferent events do
                   not move the median of ~400 000 pixels, global changes do.
  2. normalise     F'(t) = F(t) * G / g(t)   with G = median over t of g(t).
                   Removes every global multiplicative change, fast and slow;
                   dF/F (and so detection) is scale-free, so events are
                   unchanged and bleaching simply disappears.
  3. flag          two classes (flag_frames):
                   ARTEFACT — fast changes: a jump between consecutive frames
                     beyond max(k robust SD, 3 %), a spike (<= 10 s, returning)
                     or the content correlation collapsing; +tails while the
                     level is still moving, +1 frame each side, islands < 10 s.
                   STATE — slow changes: a level plateau longer than a spike or
                     not returning, or the p10/median ratio (which no gain
                     change can move) > 5 % from its clean-frame baseline. A
                     laser/focus change or widespread activity — the trace
                     cannot tell which, so it is only annotated.
  4. check         every frame of a flagged stretch is compared with the kept
                   frames around it (check_frames): whole-field change and the
                   smoothed RMS % difference, relative to how much two clean
                   frames differ. The FLAGGED FRAMES of an artefact stretch are
                   the frames that differ by >= 10 % of the frame and >= 3 x
                   the clean floor — the flash itself, not the tails and
                   margins around it (frames no longer on disk count too).
  5. decide        per stretch: keep | exclude (= its flagged frames only) |
                   exclude-all (the whole stretch) | frames:68,76-77 (exactly
                   these). Excluded frames: aqua_analyze drops the events
                   starting in them +/- a margin and their time. Defaults:
                   artefact with flagged frames -> exclude; artefact without
                   -> keep; state -> keep (annotated only). Decisions live in
                   <base>_frames_decisions.json (GUI / aqua_clean.py) and
                   survive re-runs.
  6. movie         frames STAY in the movie (mode "keep", default). "drop"
                   leaves the excluded frames out (column tiff_frame maps movie
                   frames to recorded frames), "fill" replaces them by a frozen
                   filler; both explicit.

Everything is written next to the tiff: <base>_frames.csv (per recorded
frame: level, gain, class, stretch, decision, flagged_frame, excluded,
tiff_frame, reason, the per-frame check), <base>_globalfix.png (QC traces;
red = artefact, amber = state, black bar = excluded), the check files —
<base>_check_frames.tif (16-bit ImageJ stack of every flagged stretch with a
few frames of context, each slice labelled with its recorded frame number,
stretch, numbers, flagged / excluded), <base>_check_movie.tif (the same
slices as a composite: frame, % difference from the neighbours, traces with
the cursor), <base>_flagged_frames.png/.csv (overview sheet + table) — and a
"global_correction" block in the .json so the step is never applied twice.
`apply(..., force=True)` re-flags from the raw statistics recorded in
<base>_frames.csv (frames dropped or replaced by an earlier run stay
excluded: their content is gone; restore_raw() puts a <base>_raw.tif back).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np


def global_trace(tif_path: Path, corr_step: int = 7, ref_every: int = 5, log=print) -> dict:
    """Per-frame level statistics, read page by page (memory-flat)."""
    import tifffile
    with tifffile.TiffFile(str(tif_path)) as tf:
        pages = tf.pages
        T = len(pages)
        H, W = pages[0].shape[:2]
        # reference: median of a subsample of frames
        idx = list(range(0, T, max(1, ref_every)))[:200]
        ref = np.median(np.stack([pages[i].asarray().astype(np.float32) for i in idx]), axis=0)
        rflat = ref.ravel()[::corr_step]
        rflat = rflat - rflat.mean()
        rnorm = np.sqrt((rflat ** 2).sum()) + 1e-9
        med = np.zeros(T); mean = np.zeros(T); p10 = np.zeros(T); p90 = np.zeros(T); corr = np.zeros(T)
        for t in range(T):
            a = pages[t].asarray()
            med[t] = np.median(a)
            mean[t] = a.mean()
            p10[t], p90[t] = np.percentile(a, [10, 90])
            f = a.ravel()[::corr_step].astype(np.float32)
            f = f - f.mean()
            corr[t] = float((f * rflat).sum() / (np.sqrt((f ** 2).sum()) * rnorm + 1e-9))
    return {"T": T, "H": H, "W": W, "median": med, "mean": mean, "p10": p10, "p90": p90,
            "corr_ref": corr, "ref": ref}


def _moving_median(x: np.ndarray, win: int) -> np.ndarray:
    from scipy.ndimage import median_filter
    win = max(3, int(win) | 1)
    return median_filter(x, size=win, mode="nearest")


def _masked_moving_median(x: np.ndarray, win: int, valid: np.ndarray,
                          extrapolate_edges: bool = False) -> np.ndarray:
    """Moving median of x using only frames marked valid. A window with fewer
    than half of its frames valid gives no estimate; those positions are
    bridged linearly between the nearest estimated ones (constant at the ends).
    With every frame valid this is the plain moving median.
    extrapolate_edges: where the window is cut by the start/end of the
    recording, continue the trend of the baseline linearly instead of using
    the one-sided median (which lags a drift by half a window)."""
    T = x.size
    win = max(3, int(win) | 1)
    h = win // 2
    out = np.full(T, np.nan)
    for t in range(T):
        lo, hi = max(0, t - h), min(T, t + h + 1)
        v = valid[lo:hi]
        if v.sum() >= max(3, (hi - lo) // 2):
            out[t] = np.median(x[lo:hi][v])
    ok = ~np.isnan(out)
    if not ok.any():
        return np.full(T, float(np.median(x[valid])) if valid.any() else float(np.median(x)))
    if not ok.all():
        out[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), out[ok])
    if extrapolate_edges and T > 3 * h + 2:
        inner = np.arange(h, T - h)              # positions with a full window
        for edge in (slice(0, h), slice(T - h, T)):
            ref = inner[:h] if edge.start == 0 else inner[-h:]
            a, b = np.polyfit(ref, out[ref], 1)
            out[edge] = a * np.arange(edge.start, edge.stop) + b
    return out


def _fill_islands(flag: np.ndarray, min_len: int) -> np.ndarray:
    """Flag clean runs shorter than min_len that sit between flagged runs."""
    if min_len <= 1 or not flag.any():
        return flag
    out = flag.copy()
    T = flag.size
    t = 0
    while t < T:
        if out[t]:
            t += 1
            continue
        s = t
        while t < T and not out[t]:
            t += 1
        # clean run s..t-1 ; bounded by flagged frames on both sides?
        if s > 0 and t < T and (t - s) < min_len:
            out[s:t] = True
    return out


def _robust_sd(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return 1e-9
    return float(1.4826 * np.median(np.abs(x - np.median(x))) + 1e-9)


def _extend(flag: np.ndarray, weak: np.ndarray, max_steps: int):
    """Grow each True run of `flag` over adjacent `weak` frames, at most
    max_steps frames on each side. Returns (extended, added)."""
    T = flag.size
    cur = flag.copy()
    added = np.zeros(T, bool)
    for _ in range(int(max_steps)):
        nb = np.zeros(T, bool)
        nb[1:] |= cur[:-1]
        nb[:-1] |= cur[1:]
        ext = nb & weak & ~cur
        if not ext.any():
            break
        cur |= ext
        added |= ext
    return cur, added


def flag_frames(med: np.ndarray, corr: np.ndarray, spf: float, win_s: float = 60.0,
                k: float = 4.0, min_rel: float = 0.03, corr_k: float = 5.0,
                min_corr_drop: float = 0.05, grow: int = 1, p10: np.ndarray | None = None,
                shape_rel: float = 0.05, shape_win_s: float = 300.0, min_clean_s: float = 10.0,
                max_spike_s: float = 10.0, max_iter: int = 8):
    """Two classes of flags, returned as (artefact, reasons, info); info["state"]
    holds the second class.

    ARTEFACT (excluded from the event statistics by default):
      jump   : the level changes by more than max(k robust SD, min_rel) from one
               frame to the next — both frames of the jump; a flash, a shutter
               frame, the edge of a step
      spike  : the level deviates from a 1-min moving median by more than that
               floor for at most max_spike_s and comes back to where it was
      content: the correlation with the reference image collapses by more
               than max(corr_k robust SD, min_corr_drop)
      A stretch then extends over neighbours still beyond half the floor (the
      decaying tail), grows by `grow` frames and swallows clean islands shorter
      than min_clean_s between two stretches.
    STATE (annotated; excluded only when its normalised background would be
    BRIGHTER than normal, which AQuA2 would call events — see apply()):
      plateau: the level deviates from the moving median for longer than a
               spike, or does not come back: a slower change of state whose
               cause (laser, focus, widespread activity) the trace cannot tell
      shape  : the p10/median ratio, which a pure gain change cannot move,
               deviates by more than max(k robust SD, shape_rel) from a 5-min
               moving median over clean frames (iterated to a fixed point):
               the fluorescence changed but the background did not, or the
               other way round
    All floors are absolute on purpose: a sparse afferent event never moves the
    frame median by 3 % or the p10/median ratio by 5 %.
    """
    T = med.size
    win = max(5, int(round(win_s / max(spf, 1e-6))) | 1)
    swin = max(win, int(round(shape_win_s / max(spf, 1e-6))) | 1)
    min_island = max(1, int(round(min_clean_s / max(spf, 1e-6))))
    max_spike = max(1, int(round(max_spike_s / max(spf, 1e-6))))
    base = _moving_median(med, win)
    rel = med / np.maximum(base, 1e-9) - 1.0
    sd = _robust_sd(rel)
    thr = max(k * sd, min_rel)
    level_flag = np.abs(rel) > thr
    # frame-to-frame jumps
    jump_rel = np.zeros(T)
    jump_rel[1:] = (med[1:] - med[:-1]) / np.maximum(base[1:], 1e-9)
    jump_flag = np.zeros(T, bool)
    big = np.abs(jump_rel) > thr
    jump_flag |= big
    jump_flag[:-1] |= big[1:]                      # the frame before the jump too
    # spikes vs plateaus among the level excursions
    spike_flag = np.zeros(T, bool)
    plateau_flag = np.zeros(T, bool)
    for (a, b) in _runs(level_flag):
        before = med[max(0, a - 5):a]
        after = med[b + 1:b + 6]
        returns = (before.size and after.size
                   and abs(float(np.median(after)) - float(np.median(before))) / max(float(base[a]), 1e-9) < thr)
        if (b - a + 1) <= max_spike and returns:
            spike_flag[a:b + 1] = True
        else:
            plateau_flag[a:b + 1] = True
    csd = _robust_sd(corr)
    corr_thr = float(np.median(corr) - max(corr_k * csd, min_corr_drop))
    corr_flag = corr < corr_thr
    # artefact class: fast changes, with short tails (frames in which the level
    # is still moving fast — a decaying flash; a plateau after a step is not),
    # margins and islands
    art = jump_flag | spike_flag | corr_flag
    moving = np.abs(jump_rel) > thr / 2
    moving[:-1] |= np.abs(jump_rel[1:]) > thr / 2
    cur, hyst = _extend(art, moving, max_spike)
    if grow > 0 and cur.any():
        from scipy.ndimage import binary_dilation
        cur = binary_dilation(cur, iterations=int(grow))
    artefact = _fill_islands(cur, min_island)
    # state class: plateaus + the shape rule (iterated; artefact and state
    # frames are excluded from its baseline)
    shape = (p10 / np.maximum(med, 1e-9)) if p10 is not None else None
    shape_flag = np.zeros(T, bool)
    srel = np.zeros(T); sbase = np.full(T, np.nan); sthr = float("nan")
    state = plateau_flag & ~artefact
    n_iter = 0
    if shape is not None:
        for n_iter in range(1, max_iter + 1):
            prev = state
            valid = ~(artefact | state)
            if valid.sum() < max(3, T // 10):
                valid = np.ones(T, bool)
            sbase = _masked_moving_median(shape, swin, valid, extrapolate_edges=True)
            srel = shape / np.maximum(sbase, 1e-9) - 1.0
            ssd = _robust_sd(srel[valid])
            sthr = max(k * ssd, shape_rel)
            shape_flag = np.abs(srel) > sthr
            core = (plateau_flag | shape_flag) & ~artefact
            # a state stretch is one stretch: extend over neighbours still
            # half-way out (no limit: a plateau is long by definition) and
            # over short clean islands, but never into artefact frames
            weak_state = ((np.abs(srel) > sthr / 2) | (np.abs(rel) > thr / 2)) & ~artefact
            ext, _ = _extend(core, weak_state, T)
            state = _fill_islands(ext | artefact, min_island) & ~artefact
            if np.array_equal(state, prev):
                break
    reasons = np.array(["" for _ in range(T)], dtype=object)
    any_art_rule = jump_flag | spike_flag | corr_flag
    for t in range(T):
        r = []
        if jump_flag[t]:
            j = jump_rel[t] if big[t] else (jump_rel[t + 1] if t + 1 < T else 0.0)
            r.append(f"jump {j * 100:+.1f}% between frames")
        if spike_flag[t]:
            r.append(f"spike {rel[t] * 100:+.1f}%")
        if corr_flag[t]:
            r.append(f"content corr {corr[t]:.2f}")
        if artefact[t] and not r:
            if hyst[t]:
                r.append(f"tail of an artefact stretch (level {rel[t] * 100:+.1f}%)")
            elif any_art_rule[max(0, t - grow):t + grow + 1].any():
                r.append("next to an artefact frame")
            else:
                r.append("short clean stretch between artefact frames")
        if state[t]:
            if plateau_flag[t]:
                r.append(f"state: level {rel[t] * 100:+.1f}% (plateau)")
            if shape_flag[t]:
                r.append(f"state: shape p10/median {srel[t] * 100:+.1f}%")
            if not (plateau_flag[t] or shape_flag[t]):
                r.append(f"state stretch (level {rel[t] * 100:+.1f}%, shape {srel[t] * 100:+.1f}%)")
        reasons[t] = "; ".join(r)
    return artefact, reasons, {"rel": rel, "base": base, "robust_sd": sd, "threshold": thr,
                               "jump_rel": jump_rel, "corr_threshold": corr_thr, "shape": shape,
                               "shape_base": sbase, "shape_rel": srel, "shape_threshold": sthr,
                               "state": state, "n_iter": n_iter,
                               "n_jump": int(jump_flag.sum()), "n_spike": int(spike_flag.sum()),
                               "n_corr": int(corr_flag.sum()), "n_plateau": int((plateau_flag & ~artefact).sum()),
                               "n_shape": int((shape_flag & ~artefact).sum()),
                               "n_level": int(level_flag.sum())}


def correct_tiff(tif_path: Path, gain: np.ndarray, flag: np.ndarray, mode: str = "drop",
                 present: np.ndarray | None = None, out_path: Path | None = None,
                 n_fill: int = 10, log=print) -> dict:
    """Write the corrected movie (uint16, ImageJ layout). In place when
    out_path is None: written to a temp file and swapped in.

    gain[t] multiplies original frame t. `present[t]` says which original
    frames the source tiff still holds (all, unless frames were dropped by an
    earlier run); page i of the source is the i-th present frame.
    mode "drop" : flagged frames are left out — the movie gets shorter and
                  the returned tiff_frame map (1-based, 0 = not in the movie)
                  says where each original frame went.
    mode "fill" : flagged frames are replaced by a linear blend between the
                  median of up to n_fill clean frames before the stretch and
                  after it (frozen, information-free filler, movie length kept).
    mode "keep" : every frame kept (gain only); frames missing from the source
                  are filled as in "fill"."""
    import tifffile
    tif_path = Path(tif_path)
    gain = np.asarray(gain, dtype=np.float64)
    T = gain.size
    flag = np.asarray(flag, bool)
    present = np.ones(T, bool) if present is None else np.asarray(present, bool)
    page_of = np.cumsum(present) - 1                    # source page of original frame t
    if mode == "keep":
        replace = ~present
    else:
        replace = flag | ~present
    good = np.flatnonzero(~replace)                     # usable as neighbours / to keep
    if mode == "drop":
        write_idx = good
    else:
        write_idx = np.arange(T)
    tiff_frame = np.zeros(T, int)
    tiff_frame[write_idx] = np.arange(1, write_idx.size + 1)
    tmp = (out_path or tif_path).with_name((out_path or tif_path).name + ".partial")
    with tifffile.TiffFile(str(tif_path)) as tf:
        pages = tf.pages
        if len(pages) != int(present.sum()):
            raise ValueError(f"{tif_path.name}: {len(pages)} pages but {int(present.sum())} frames expected")
        H, W = pages[0].shape[:2]
        n_out = int(write_idx.size)
        gb = n_out * H * W * 2 / 1e9
        kw = dict(imagej=True, metadata={"axes": "TYX"}) if gb < 3.9 else dict(bigtiff=True)
        dst = tifffile.memmap(str(tmp), shape=(n_out, H, W), dtype=np.uint16, **kw)

        def corrected(t):
            return pages[int(page_of[t])].asarray().astype(np.float32) * np.float32(gain[t])

        def clean_median(idx):
            """Per-pixel median of the clean frames idx (robust to an event in
            one of them), rescaled so its spatial median equals the level of
            those frames: with skewed shot noise a temporal median sits a few %
            below a single frame's spatial median, which would read as a dip."""
            frames = [corrected(int(i)) for i in idx]
            level = float(np.median([np.median(f) for f in frames]))
            fill = frames[0] if len(frames) == 1 else np.median(np.stack(frames), axis=0)
            return fill * np.float32(level / max(float(np.median(fill)), 1e-6))

        fill_cache: dict[tuple, tuple] = {}
        for o, t in enumerate(write_idx):
            t = int(t)
            if replace[t] and good.size:
                lo = good[good < t]
                hi = good[good > t]
                a = int(lo[-1]) if lo.size else None
                b = int(hi[0]) if hi.size else None
                key = (a, b)
                if key not in fill_cache:
                    A = clean_median(lo[-n_fill:]) if lo.size else None
                    B = clean_median(hi[:n_fill]) if hi.size else None
                    fill_cache = {key: (A, B)}          # one stretch at a time
                A, B = fill_cache[key]
                if A is not None and B is not None:
                    w = (t - a) / float(b - a)
                    fr = (1 - w) * A + w * B
                elif A is not None:
                    fr = A
                else:
                    fr = B
            elif replace[t]:
                fr = np.zeros((H, W), np.float32)
            else:
                fr = corrected(t)
            dst[o] = np.clip(np.round(fr), 0, 65535).astype(np.uint16)
        dst.flush(); del dst
    target = out_path or tif_path
    tmp.replace(target)
    return {"gain_min": float(gain.min()), "gain_max": float(gain.max()),
            "n_written": n_out, "n_dropped": int(T - n_out),
            "n_filled": int(replace[write_idx].sum()), "tiff_frame": tiff_frame}


def drop_frames(tif_path: Path, keep: np.ndarray, log=print) -> int:
    """Rewrite a companion tiff (e.g. the red channel) with only the frames
    marked keep (bool over its pages). Returns the number of pages written."""
    import tifffile
    tif_path = Path(tif_path)
    keep = np.asarray(keep, bool)
    idx = np.flatnonzero(keep)
    tmp = tif_path.with_name(tif_path.name + ".partial")
    with tifffile.TiffFile(str(tif_path)) as tf:
        pages = tf.pages
        if len(pages) != keep.size:
            raise ValueError(f"{tif_path.name}: {len(pages)} pages, expected {keep.size}")
        H, W = pages[0].shape[:2]
        gb = idx.size * H * W * 2 / 1e9
        kw = dict(imagej=True, metadata={"axes": "TYX"}) if gb < 3.9 else dict(bigtiff=True)
        dst = tifffile.memmap(str(tmp), shape=(idx.size, H, W), dtype=np.uint16, **kw)
        for o, t in enumerate(idx):
            dst[o] = pages[int(t)].asarray()
        dst.flush(); del dst
    tmp.replace(tif_path)
    return int(idx.size)


def _runs(flag: np.ndarray):
    """[(start, end)] inclusive, 0-based, of True runs."""
    idx = np.flatnonzero(flag)
    if idx.size == 0:
        return []
    runs = []
    s = p = int(idx[0])
    for i in idx[1:]:
        i = int(i)
        if i == p + 1:
            p = i
            continue
        runs.append((s, p)); s = p = i
    runs.append((s, p))
    return runs


def _runs_str(flag: np.ndarray, max_runs: int = 12) -> str:
    """'14-17, 39-41, 52-100' (1-based)."""
    runs = _runs(np.asarray(flag, bool))
    if not runs:
        return ""
    txt = ", ".join(f"{a + 1}" if a == b else f"{a + 1}-{b + 1}" for a, b in runs[:max_runs])
    return txt + (" ..." if len(runs) > max_runs else "")


# ----------------------------------------------------------------------------- stretches + decisions

def build_stretches(artefact: np.ndarray, state: np.ndarray, info: dict, spf: float) -> list[dict]:
    """One dict per contiguous flagged stretch of either class, in frame order.
    id = '<first>-<last>' (1-based) — what the decisions sidecar refers to.
    The default decision is provisional here (artefact -> exclude, state ->
    keep); check_frames() finalises it once it knows which frames of an
    artefact stretch really differ from their neighbours (its FLAGGED frames):
    a stretch without any keeps everything."""
    rel = info.get("rel"); srel = info.get("shape_rel"); jump = info.get("jump_rel")
    sthr = info.get("shape_threshold", float("nan"))
    out = []
    for cls, flag in (("artefact", artefact), ("state", state)):
        for (a, b) in _runs(np.asarray(flag, bool)):
            seg = slice(a, b + 1)
            lvl = rel[seg] * 100 if rel is not None else np.zeros(b - a + 1)
            shp = srel[seg] * 100 if srel is not None else np.zeros(b - a + 1)
            jmp = np.abs(jump[seg]) * 100 if jump is not None else np.zeros(b - a + 1)
            shape_med = float(np.median(shp))
            if cls == "artefact":
                default = "exclude"
                why = "fast change (flash / jump / content): its flagged frames are excluded, the rest kept"
            else:
                bright = bool(np.isfinite(sthr) and shape_med > sthr * 100)
                default = "keep"
                why = ("slow change; background brighter than normal after the correction (p10/median "
                       f"{shape_med:+.0f}%) — kept, check that AQuA2 does not call it events"
                       if bright else "slow change, background not brighter after the correction — kept")
            out.append({"id": f"{a + 1}-{b + 1}", "class": cls, "first": a + 1, "last": b + 1,
                        "n_frames": b - a + 1, "seconds": round((b - a + 1) * spf, 1),
                        "level_pct_max": round(float(lvl[np.argmax(np.abs(lvl))]), 1),
                        "shape_pct_median": round(shape_med, 1), "jump_pct_max": round(float(jmp.max()), 1),
                        "default": default, "default_why": why})
    out.sort(key=lambda d: d["first"])
    return out


# ----------------------------------------------------------------------------- decisions

DECISION_MEANING = ("keep = nothing of this stretch leaves the statistics; exclude = only its FLAGGED frames "
                    "(the frames that differ from their neighbours by the core rule, plus frames no longer on "
                    "disk) leave the statistics (events starting in them +/- a margin, and their time); "
                    "exclude-all = every frame of the stretch; frames:68,76-77 = exactly these recorded frames "
                    "(1-based, within the stretch). Frames stay in the movie either way unless you drop them "
                    "explicitly.")


def parse_frame_ranges(txt: str) -> list[int]:
    """'68, 76-77' -> [68, 76, 77] (1-based, sorted, unique)."""
    out = set()
    for part in str(txt).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a, b = int(a.strip()), int(b.strip())
            if b < a:
                a, b = b, a
            out.update(range(a, b + 1))
        else:
            out.add(int(part))
    return sorted(out)


def frames_str(frames, max_runs: int = 12) -> str:
    """[67, 68, 76, 77, 99] -> '67-68, 76-77, 99'."""
    frames = sorted({int(f) for f in (frames if frames is not None else [])})
    if not frames:
        return ""
    runs = []
    s = p = frames[0]
    for f in frames[1:]:
        if f == p + 1:
            p = f
            continue
        runs.append((s, p)); s = p = f
    runs.append((s, p))
    txt = ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b in runs[:max_runs])
    return txt + (" ..." if len(runs) > max_runs else "")


def normalise_decision(d) -> str:
    """-> 'keep' | 'exclude' | 'exclude-all' | 'frames:<ranges>'."""
    s = str(d or "").strip().lower()
    if s.startswith("frames:") or s.startswith("frames="):
        try:
            return "frames:" + frames_str(parse_frame_ranges(s[7:])).replace(" ", "")
        except ValueError:
            return "keep"
    if s.replace("_", "-").replace(" ", "-") in ("exclude-all", "all", "exclude-stretch", "exclude-whole-stretch"):
        return "exclude-all"
    if s.startswith("ex"):
        return "exclude"
    return "keep"


def decision_frames(st: dict, decision) -> list[int]:
    """The 1-based recorded frames a decision excludes in this stretch."""
    d = normalise_decision(decision)
    lo, hi = int(st["first"]), int(st["last"])
    if d == "keep":
        return []
    if d == "exclude-all":
        return list(range(lo, hi + 1))
    if d.startswith("frames:"):
        return [f for f in parse_frame_ranges(d[7:]) if lo <= f <= hi]
    core = st.get("core_frames")
    if core is None:                    # results of an older version without the per-frame check
        return list(range(lo, hi + 1))
    return sorted({int(f) for f in core if lo <= int(f) <= hi})


def describe_decision(st: dict, decision) -> str:
    """'exclude (flagged frames 67-68, 76-77)' — for logs and the GUI."""
    d = normalise_decision(decision)
    fr = decision_frames(st, d)
    if d == "keep":
        return "keep"
    if d == "exclude-all":
        return f"exclude-all ({len(fr)} frames)"
    if d.startswith("frames:"):
        return f"frames {frames_str(fr)}" if fr else "frames: none inside the stretch (= keep)"
    if st.get("core_frames") is None:
        return "exclude (whole stretch: no per-frame check yet — re-run with force)"
    return f"exclude (flagged frames {frames_str(fr)})" if fr else "exclude (no flagged frames = keep)"


def decisions_path(tif_path: Path) -> Path:
    return Path(str(Path(tif_path).with_suffix("")) + "_frames_decisions.json")


def read_decisions(tif_path: Path) -> dict:
    p = decisions_path(tif_path)
    if not p.exists():
        return {}
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        return {str(k): str(v) for k, v in d.get("decisions", {}).items()}
    except Exception:
        return {}


def write_decisions(tif_path: Path, decisions: dict, notes: dict | None = None):
    from datetime import datetime
    p = decisions_path(tif_path)
    old = {}
    if p.exists():
        try:
            old = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            old = {}
    rec = {"decisions": {str(k): normalise_decision(v) for k, v in decisions.items()},
           "notes": {**old.get("notes", {}), **(notes or {})},
           "meaning": DECISION_MEANING,
           "updated": datetime.now().isoformat(timespec="seconds")}
    p.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    return p


def excluded_frames(stretches: list[dict], decisions: dict, T: int) -> tuple[np.ndarray, list[str]]:
    """Bool over recorded frames + the (normalised) decision applied to each
    stretch. Sets st['excluded_frames'] (1-based) on every stretch."""
    ex = np.zeros(T, bool)
    applied = []
    for st in stretches:
        d = normalise_decision(decisions.get(st["id"], st["default"]))
        applied.append(d)
        fr = [f for f in decision_frames(st, d) if 1 <= f <= T]
        st["excluded_frames"] = fr
        for f in fr:
            ex[f - 1] = True
    return ex, applied


# ----------------------------------------------------------------------------- per-frame csv

def write_frames_csv(path: Path, spf: float, tr: dict, artefact, state, stretch_of, decision_of, excluded,
                     reasons, gain, tiff_frame, info: dict | None = None, replaced=None, core=None,
                     check: dict | None = None):
    """One row per RECORDED frame. flagged = artefact class, state = state
    class, stretch = the stretch id the frame belongs to, decision = what the
    statistics do with that stretch, flagged_frame = 1 when the frame is one
    of the stretch's FLAGGED frames (differs from its neighbours by the core
    rule), excluded = 1 when its events leave the statistics, tiff_frame =
    1-based index in the movie (0 = dropped), replaced = 1 when the movie
    holds a filler instead of the frame, check_diff_pct / check_ratio = the
    per-frame check (RMS % difference from the clean neighbours, and that
    relative to the clean floor)."""
    import csv
    T = tr["T"]
    info = info or {}
    base = info.get("base"); shape = info.get("shape"); sbase = info.get("shape_base"); jump = info.get("jump_rel")
    tiff_frame = np.asarray(tiff_frame, int)
    artefact = np.asarray(artefact, bool); state = np.asarray(state, bool); excluded = np.asarray(excluded, bool)
    replaced = np.zeros(T, bool) if replaced is None else np.asarray(replaced, bool)
    core = np.zeros(T, bool) if core is None else np.asarray(core, bool)
    check = check or {}
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "time_s", "level_median", "level_mean", "p10", "p90", "corr_ref",
                    "gain_applied", "flagged", "state", "flag_class", "stretch", "decision", "flagged_frame",
                    "excluded", "tiff_frame", "replaced", "reason", "level_baseline", "p10_over_median",
                    "shape_baseline", "jump_pct", "check_field_pct", "check_diff_pct", "check_ratio"])
        for t in range(T):
            cls = "artefact" if artefact[t] else "state" if state[t] else ""
            c = check.get(t)
            w.writerow([t + 1, round(t * spf, 3), round(tr["median"][t], 1), round(tr["mean"][t], 1),
                        round(tr["p10"][t], 1), round(tr["p90"][t], 1), round(tr["corr_ref"][t], 4),
                        round(float(gain[t]), 4), int(artefact[t]), int(state[t]), cls,
                        stretch_of[t] if stretch_of is not None else "",
                        decision_of[t] if decision_of is not None else "",
                        int(core[t]), int(excluded[t]), int(tiff_frame[t]), int(replaced[t]), reasons[t],
                        round(float(base[t]), 1) if base is not None else "",
                        round(float(shape[t]), 4) if shape is not None else "",
                        round(float(sbase[t]), 4) if sbase is not None and np.isfinite(sbase[t]) else "",
                        round(float(jump[t]) * 100, 2) if jump is not None else "",
                        round(c[1] * 100, 1) if c else "", round(c[0] * 100, 2) if c else "",
                        round(c[2], 2) if c else ""])


def read_frames_csv(path: Path):
    """-> dict over RECORDED frames: frame (1-based), flagged (artefact),
    state, excluded (events leave the statistics; older files: = flagged),
    core (the stretch's flagged frames; None when the file predates the
    per-frame check), stretch, decision, tiff_frame (1-based index in the
    movie, 0 = dropped), level, p10, base, shape, sbase, check_ratio, ... or
    None."""
    path = Path(path)
    if not path.exists():
        return None
    import csv
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if not rows:
        return None

    def col(name, cast=float, default=np.nan):
        return np.array([cast(r[name]) if r.get(name, "") not in ("", None) else default for r in rows])

    frame = col("frame", int, 0)
    flagged = np.array([r.get("flagged", "0") == "1" for r in rows])
    tiff_frame = col("tiff_frame", int, 0) if "tiff_frame" in rows[0] else frame.copy()
    excluded = (np.array([r.get("excluded", "0") == "1" for r in rows]) if "excluded" in rows[0]
                else flagged.copy())
    core = np.array([r.get("flagged_frame", "0") == "1" for r in rows]) if "flagged_frame" in rows[0] else None
    return {"frame": frame, "flagged": flagged,
            "state": np.array([r.get("state", "0") == "1" for r in rows]),
            "excluded": excluded, "core": core,
            "stretch": [r.get("stretch", "") for r in rows],
            "decision": [r.get("decision", "") for r in rows],
            "replaced": np.array([r.get("replaced", "0") == "1" for r in rows]),
            "tiff_frame": tiff_frame, "dropped": tiff_frame <= 0,
            "level": col("level_median"), "mean": col("level_mean"),
            "p10": col("p10"), "p90": col("p90"), "corr": col("corr_ref"),
            "gain": col("gain_applied", float, 1.0),
            "base": col("level_baseline") if "level_baseline" in rows[0] else None,
            "shape": col("p10_over_median") if "p10_over_median" in rows[0] else None,
            "sbase": col("shape_baseline") if "shape_baseline" in rows[0] else None,
            "check_field": col("check_field_pct") if "check_field_pct" in rows[0] else None,
            "check_diff": col("check_diff_pct") if "check_diff_pct" in rows[0] else None,
            "check_ratio": col("check_ratio") if "check_ratio" in rows[0] else None,
            "reason": [r.get("reason", "") for r in rows]}


def frame_map(frames: dict | None, n_tiff: int | None = None) -> np.ndarray | None:
    """Original 1-based frame number for each movie frame (index 0 = movie
    frame 1), from read_frames_csv output. None when nothing was dropped."""
    if frames is None or not frames["dropped"].any():
        return None
    present = np.flatnonzero(~frames["dropped"])
    orig = frames["frame"][present]
    if n_tiff is not None and orig.size != n_tiff:
        return None
    return orig


def core_from_check_csv(path: Path, stretches: list[dict], core_ratio: float = 3.0, core_min_pct: float = 10.0,
                        gone=None) -> bool:
    """Fill st['core_frames'] from an existing <base>_flagged_frames.csv (its
    per-frame ratio / difference columns) — for results made before the
    per-frame check stored them. Returns True when the file was usable."""
    import csv
    path = Path(path)
    if not path.exists():
        return False
    try:
        rows = list(csv.DictReader(path.open(encoding="utf-8")))
    except Exception:
        return False
    if not rows or "ratio" not in rows[0]:
        return False
    per = {}
    for r in rows:
        try:
            per[int(r["frame"])] = (float(r["diff_from_neighbours_rms_pct"]), float(r["ratio"]))
        except (ValueError, KeyError):
            continue
    gone = set(int(g) for g in (gone or []))
    for st in stretches:
        if st.get("class") != "artefact":
            st["core_frames"] = []
            continue
        core = [f for f in range(int(st["first"]), int(st["last"]) + 1)
                if f in gone or (f in per and per[f][1] >= core_ratio and per[f][0] >= core_min_pct)]
        st["core_frames"] = core
    return True


def apply_decisions(tif_path: Path, decisions: dict | None = None, log=print) -> dict | None:
    """Change what the statistics do with flagged stretches without touching
    the movie: merge `decisions` ({stretch id: 'keep' | 'exclude' |
    'exclude-all' | 'frames:68,76-77'}) into the sidecar, recompute the
    decision / excluded columns of <base>_frames.csv and the json summary.
    Returns the stretch list with the decisions applied."""
    tif_path = Path(tif_path)
    base = tif_path.with_suffix("")
    csv_path = Path(str(base) + "_frames.csv")
    json_path = Path(str(base) + ".json")
    fr = read_frames_csv(csv_path)
    if fr is None:
        log(f"  {tif_path.name}: no _frames.csv — run the brightness step first")
        return None
    meta = {}
    if json_path.exists():
        try:
            meta = json.loads(json_path.read_text(encoding="utf-8"))
        except Exception:
            meta = {}
    gc = meta.get("global_correction") or {}
    stretches = gc.get("stretches")
    T = fr["frame"].size
    if not stretches:                               # rebuild from the csv
        stretches = []
        for cls, flag in (("artefact", fr["flagged"]), ("state", fr["state"])):
            for (a, b) in _runs(flag):
                stretches.append({"id": f"{a + 1}-{b + 1}", "class": cls, "first": a + 1, "last": b + 1,
                                  "n_frames": b - a + 1, "default": "exclude" if cls == "artefact" else "keep"})
        stretches.sort(key=lambda d: d["first"])
    if any("core_frames" not in st for st in stretches):
        gone = [int(f) for f in fr["frame"][fr["dropped"] | fr["replaced"]]]
        if fr.get("core") is not None:
            for st in stretches:
                st["core_frames"] = [int(f) for f in fr["frame"][fr["core"]]
                                     if int(st["first"]) <= int(f) <= int(st["last"])]
        elif core_from_check_csv(Path(str(base) + "_flagged_frames.csv"), stretches,
                                 gone=gone, **{k: v for k, v in (gc.get("params") or {}).items()
                                               if k in ("core_ratio", "core_min_pct")}):
            log(f"  {tif_path.name}: flagged frames taken from the existing check csv (older result — "
                f"run the brightness step with force to get the check stacks and the current defaults)")
        else:
            log(f"  {tif_path.name}: no per-frame check on file — 'exclude' means the whole stretch here; run the "
                f"brightness step with force to check frame by frame")
        for st in stretches:
            if st.get("class") == "artefact" and st.get("core_frames") is not None and st.get("default") == "exclude" \
                    and not st["core_frames"]:
                st["default"] = "keep"
                st["default_why"] = "no frame of it differs clearly from its neighbours"
            if st.get("class") == "state":
                st["core_frames"] = []
    dec = read_decisions(tif_path)
    if decisions:
        dec.update({str(k): normalise_decision(v) for k, v in decisions.items()})
        write_decisions(tif_path, dec)
    excluded, applied = excluded_frames(stretches, dec, T)
    core = np.zeros(T, bool)
    for st, d in zip(stretches, applied):
        st["decision"] = d
        for f in st.get("core_frames") or []:
            if 1 <= int(f) <= T:
                core[int(f) - 1] = True
    # rewrite the csv columns decision / flagged_frame / excluded
    import csv
    rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    fields = list(rows[0].keys())
    for extra in ("decision", "flagged_frame", "excluded"):
        if extra not in fields:
            fields.append(extra)
    dec_of = [""] * T
    for st, d in zip(stretches, applied):
        for t in range(st["first"] - 1, st["last"]):
            dec_of[t] = d
    for t, r in enumerate(rows[:T]):
        r["decision"] = dec_of[t]
        r["flagged_frame"] = int(core[t])
        r["excluded"] = int(excluded[t])
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    # the check csv's decision / excluded columns, when it exists
    chk_csv = Path(str(base) + "_flagged_frames.csv")
    if chk_csv.exists():
        try:
            crows = list(csv.DictReader(chk_csv.open(encoding="utf-8")))
            if crows:
                cfields = list(crows[0].keys())
                for extra in ("flagged_frame", "decision", "excluded"):
                    if extra not in cfields:
                        cfields.append(extra)
                for r in crows:
                    try:
                        f = int(r["frame"])
                    except (ValueError, KeyError):
                        continue
                    if 1 <= f <= T:
                        r["decision"] = dec_of[f - 1]; r["flagged_frame"] = int(core[f - 1]); r["excluded"] = int(excluded[f - 1])
                with open(chk_csv, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=cfields)
                    w.writeheader(); w.writerows(crows)
        except Exception as e:
            log(f"  warning: could not update {chk_csv.name}: {e}")
    gc["stretches"] = stretches
    gc["n_excluded"] = int(excluded.sum())
    gc["excluded_frames_1based"] = [int(i) + 1 for i in np.flatnonzero(excluded)]
    meta["global_correction"] = gc
    try:
        json_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    except Exception as e:
        log(f"  warning: could not update {json_path.name}: {e}")
    n_st = sum(1 for st in stretches if st.get("excluded_frames"))
    log(f"  {tif_path.name}: {int(excluded.sum())} frame(s) excluded from the statistics "
        f"({frames_str(np.flatnonzero(excluded) + 1) or 'none'}; {n_st} of {len(stretches)} stretches)")
    for st in stretches:
        if st["id"] in (decisions or {}):
            log(f"    {st['class']:8s} {st['id']:>9s}: {describe_decision(st, st['decision'])}")
    return stretches


# ----------------------------------------------------------------------------- QC plot

def qc_plot(png: Path, spf: float, tr: dict, artefact, state, excluded, info: dict,
            behavior_csv: Path | None, title: str, mode: str = "keep"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    T = tr["T"]
    t = np.arange(T) * spf / 60.0
    artefact = np.asarray(artefact, bool); state = np.asarray(state, bool); excluded = np.asarray(excluded, bool)
    fig, axes = plt.subplots(3, 1, figsize=(13, 8.5), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1.3, 1]})
    ax = axes[0]
    ax.plot(t, tr["median"], color="k", lw=0.9, label="frame median (global level)")
    ax.plot(t, info["base"], color="tab:blue", lw=0.9, ls="--", label="1-min moving median")
    for s in np.flatnonzero(artefact):
        ax.axvspan(t[s] - spf / 120, t[s] + spf / 120, color="#d62728", alpha=0.35, lw=0)
    for s in np.flatnonzero(state):
        ax.axvspan(t[s] - spf / 120, t[s] + spf / 120, color="#ff9f1c", alpha=0.30, lw=0)
    for s in np.flatnonzero(excluded):
        ax.axvspan(t[s] - spf / 120, t[s] + spf / 120, ymin=0.0, ymax=0.04, color="#333333", alpha=0.9, lw=0)
    ax.set_ylabel("intensity")
    sthr = info.get("shape_threshold", float("nan"))
    what = {"keep": "frames kept in the movie", "drop": "excluded frames left out of the movie",
            "fill": "excluded frames replaced by a filler"}.get(mode, mode)
    ax.set_title(f"{title}\n"
                 f"red = artefact ({int(artefact.sum())} frames: flashes, jumps > {info['threshold'] * 100:.1f}%, content), "
                 f"amber = state change ({int(state.sum())} frames: plateaus, p10/median change > {sthr * 100:.1f}%); "
                 f"black bar = excluded from the statistics ({int(excluded.sum())}); {what}", fontsize=9)
    ax.legend(fontsize=8, loc="upper right")
    ax = axes[1]
    ax.plot(t, info["rel"] * 100, color="k", lw=0.8, label="level vs moving median")
    if info.get("shape") is not None:
        ax.plot(t, info["shape_rel"] * 100, color="tab:orange", lw=0.8, alpha=0.9,
                label="p10/median vs clean-frame moving median")
        ax.axhline(sthr * 100, color="tab:orange", ls=":", lw=0.8)
        ax.axhline(-sthr * 100, color="tab:orange", ls=":", lw=0.8)
    ax.axhline(info["threshold"] * 100, color="tab:red", ls="--", lw=0.8)
    ax.axhline(-info["threshold"] * 100, color="tab:red", ls="--", lw=0.8)
    lim = max(10.0, float(np.nanpercentile(np.abs(info["rel"]) * 100, 99.5)) * 1.2)
    ax.set_ylim(-lim, lim)
    ax.set_ylabel("deviation (%)")
    ax.legend(fontsize=7, loc="upper right")
    ax2 = ax.twinx()
    ax2.plot(t, tr["corr_ref"], color="tab:green", lw=0.7, alpha=0.6)
    ax2.set_ylabel("corr with\nreference", color="tab:green")
    ax = axes[2]
    shown = False
    if behavior_csv is not None and Path(behavior_csv).exists():
        try:
            import csv
            rows = list(csv.DictReader(Path(behavior_csv).open()))
            if rows and "velocity_mean" in rows[0]:
                v = np.array([float(r["velocity_mean"]) for r in rows])
                ax.plot(np.arange(v.size) * spf / 60.0, np.abs(v), color="tab:purple", lw=0.8)
                ax.set_ylabel("running speed\n(counts/s)")
                shown = True
        except Exception:
            pass
    if not shown:
        ax.text(0.5, 0.5, "no behaviour csv (unpack_bin writes one when a quadrature file exists)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_yticks([])
    ax.set_xlabel("time (min)")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- the per-frame check

def _read_speed(behavior_csv):
    if behavior_csv is None or not Path(behavior_csv).exists():
        return None
    try:
        import csv
        rows = list(csv.DictReader(Path(behavior_csv).open()))
        if rows and "velocity_mean" in rows[0]:
            return np.array([abs(float(r["velocity_mean"])) for r in rows])
    except Exception:
        pass
    return None


def _ds_page(a: np.ndarray, ds: int) -> np.ndarray:
    """Block mean by ds x ds (float32)."""
    if ds <= 1:
        return a
    H, W = a.shape[:2]
    Hd, Wd = H // ds, W // ds
    return a[:Hd * ds, :Wd * ds].reshape(Hd, ds, Wd, ds).mean(axis=(1, 3))


class _FrameSource:
    """Recorded-brightness frames of a (corrected) movie — page / gain — read
    on demand; the downsampled copies are cached."""

    def __init__(self, tif_path: Path, present: np.ndarray, file_gain: np.ndarray, ds: int):
        import tifffile
        self.tf = tifffile.TiffFile(str(tif_path))
        self.pages = self.tf.pages
        self.H, self.W = self.pages[0].shape[:2]
        self.present = np.asarray(present, bool)
        self.page_of = np.cumsum(self.present) - 1
        self.gain = np.asarray(file_gain, float)
        self.ds = int(ds)
        self.cache: dict[int, np.ndarray] = {}

    def raw(self, t: int) -> np.ndarray:
        return self.pages[int(self.page_of[t])].asarray().astype(np.float32) / np.float32(self.gain[t])

    def get(self, t: int) -> np.ndarray:
        if t not in self.cache:
            self.cache[t] = _ds_page(self.raw(t), self.ds)
        return self.cache[t]

    def close(self):
        try:
            self.tf.close()
        except Exception:
            pass


def _ref_at(ref: dict, t: int) -> np.ndarray:
    """Clean reference image for frame t: the median of the kept frames before
    the stretch blended into the median of those after it."""
    A, B = ref["A"], ref["B"]
    if A is None:
        return B
    if B is None:
        return A
    a, b = int(ref["lo"][-1]), int(ref["hi"][0])
    w = min(1.0, max(0.0, (t - a) / float(b - a))) if b > a else 0.5
    return (1 - w) * A + w * B


def check_frames(tif_path: Path, stretches: list[dict], present: np.ndarray, file_gain: np.ndarray, spf: float,
                 info: dict, reasons, behavior_csv: Path | None = None, raw_lost: np.ndarray | None = None,
                 ds: int = 4, n_ref: int = 10, smooth_px: float = 4.0, core_ratio: float = 3.0,
                 core_min_pct: float = 10.0, log=print) -> dict:
    """The frame-by-frame check of every flagged stretch, made while its frames
    are on disk. Every frame of a stretch is compared with the kept frames
    around it (up to n_ref on each side — for an artefact stretch its
    immediate surroundings, so that a flash inside a plateau is judged against
    the plateau; for a state stretch the frames outside every stretch):
    whole-field change, the RMS of the smoothed % difference from those
    neighbours, and that RMS relative to how much two clean frames differ
    (the ratio, 'x1.0' = like any clean frame).

    FLAGGED FRAMES of an artefact stretch (st['core_frames'], 1-based): the
    frames whose difference is at least core_min_pct % of the frame AND
    core_ratio times the clean floor — a flash, a shutter frame, a shifted
    frame — plus frames no longer on disk. They are what the decision
    'exclude' excludes; the other frames of the stretch (the tails, margins
    and mild frames around them) are kept. A state stretch never has flagged
    frames. Sets verdict / max_ratio / core_frames and the final default
    decision on every stretch (artefact with flagged frames -> exclude, else
    keep; state -> keep). Returns the check data used by write_check_csv,
    check_sheet and check_stacks."""
    from scipy.ndimage import gaussian_filter
    T = int(info["rel"].size)
    present = np.asarray(present, bool)
    raw_lost = np.zeros(T, bool) if raw_lost is None else np.asarray(raw_lost, bool)
    avail = present & ~raw_lost
    is_art = np.zeros(T, bool); is_state = np.zeros(T, bool)
    for st in stretches:
        (is_art if st["class"] == "artefact" else is_state)[st["first"] - 1:st["last"]] = True
    clean_for = {"artefact": np.flatnonzero(~is_art & avail),
                 "state": np.flatnonzero(~is_art & ~is_state & avail)}
    chk = {"per": {}, "refs": {}, "speed": _read_speed(behavior_csv), "ds": int(ds), "smooth_px": float(smooth_px),
           "frames": {}, "T": T, "core_ratio": float(core_ratio), "core_min_pct": float(core_min_pct),
           "display": None, "avail": avail, "shape_ds": None}
    for st in stretches:
        st["core_frames"] = []
        st["verdict"] = ""
        st["stands_out"] = None
        st["max_ratio"] = None
    if not stretches:
        return chk
    src = _FrameSource(tif_path, present, file_gain, ds)
    chk["frames"] = src.cache
    chk["shape_ds"] = (src.H // ds, src.W // ds)

    def smooth_dev(img, ref):
        rel = (img - ref) / np.maximum(ref, 1e-6)
        sm = gaussian_filter(rel, smooth_px)
        return float(np.sqrt(np.mean(sm ** 2))), float(np.mean(rel))

    disp = []
    try:
        for st in stretches:
            s, e = st["first"] - 1, st["last"] - 1
            cls = st["class"]
            inside = [t for t in range(s, e + 1) if avail[t]]
            lost = [t + 1 for t in range(s, e + 1) if not avail[t]]
            if cls == "artefact":
                st["core_frames"] = list(lost)          # content gone: nothing to check, stays out
            pool = clean_for[cls]
            lo = pool[pool < s][-n_ref:]
            hi = pool[pool > e][:n_ref]
            refs = np.concatenate([lo, hi])
            if refs.size < 2 or not inside:
                st["verdict"] = ("no recorded frames left to check (removed by an earlier run)" if not inside
                                 else "no clean neighbours to compare with")
                continue
            A = np.median(np.stack([src.get(int(i)) for i in lo]), axis=0) if lo.size else None
            B = np.median(np.stack([src.get(int(i)) for i in hi]), axis=0) if hi.size else None
            M_all = np.median(np.stack([src.get(int(i)) for i in refs]), axis=0)
            floor = float(np.median([smooth_dev(src.get(int(i)), M_all)[0] for i in refs])) + 1e-9
            ref = {"lo": lo, "hi": hi, "A": A, "B": B, "M_all": M_all, "floor": floor}
            chk["refs"][st["id"]] = ref
            for t in inside:
                dev, mean_rel = smooth_dev(src.get(t), _ref_at(ref, t))
                chk["per"][t] = (dev, mean_rel, dev / floor)
            for i in refs[:4]:
                disp.append(np.percentile(src.get(int(i)), [0.5, 99.8]))
            ratios = {t: chk["per"][t][2] for t in inside}
            t_max = max(ratios, key=ratios.get)
            mx = float(ratios[t_max])
            st["max_ratio"] = round(mx, 2)
            st["max_diff_pct"] = round(float(chk["per"][t_max][0]) * 100, 1)
            st["stands_out"] = bool(mx >= 1.3)
            st["compared_with"] = f"frames {int(lo[-1]) + 1 if lo.size else '-'} / {int(hi[0]) + 1 if hi.size else '-'}"
            if cls == "artefact":
                core = [t + 1 for t in inside
                        if chk["per"][t][2] >= core_ratio and chk["per"][t][0] * 100 >= core_min_pct]
                st["core_frames"] = sorted(set(core) | set(lost))
                if core:
                    st["verdict"] = (f"{len(core)} frame(s) clearly different from the neighbours: {frames_str(core)} "
                                     f"(up to {st['max_diff_pct']:.0f}%, x{mx:.1f})")
                elif mx >= 1.3:
                    st["verdict"] = (f"mildly different from its neighbours (at most {st['max_diff_pct']:.0f}%, x{mx:.1f} "
                                     f"at f{t_max + 1}) — below the rule, kept")
                else:
                    st["verdict"] = "looks like its neighbours — kept"
                if lost:
                    st["verdict"] += f"; {len(lost)} frame(s) not on disk stay out"
            else:
                how = "clearly" if mx >= 2.0 else "mildly" if mx >= 1.3 else "hardly"
                st["verdict"] = (f"state change: {how} different from the frames outside the stretches "
                                 f"(at most {st['max_diff_pct']:.0f}%, x{mx:.1f}) — kept")
    finally:
        src.close()
    if disp:
        d = np.median(np.stack(disp), axis=0)
        chk["display"] = (float(d[0]), float(d[1]))
    # final defaults
    for st in stretches:
        if st["class"] == "artefact":
            if st["core_frames"]:
                st["default"] = "exclude"
                st["default_why"] = (f"flagged frames {frames_str(st['core_frames'])} differ from their neighbours "
                                     f"by >= {core_min_pct:.0f}% and >= x{core_ratio:.0f} the clean floor; the rest kept")
            else:
                st["default"] = "keep"
                st["default_why"] = "no frame of it differs clearly from its neighbours"
        else:
            st["default"] = "keep"
    return chk


def write_check_csv(out_csv: Path, stretches: list[dict], chk: dict, reasons, info: dict, excluded: np.ndarray):
    """<base>_flagged_frames.csv: one row per flagged frame with the rule that
    fired, the per-frame check, whether it is a flagged frame of its stretch,
    the stretch verdict, the decision and whether it is excluded."""
    import csv
    rel_level = info.get("rel"); rel_shape = info.get("shape_rel")
    speed = chk.get("speed")
    excluded = np.asarray(excluded, bool)
    rows = []
    for st in stretches:
        core = set(int(f) for f in (st.get("core_frames") or []))
        decision = st.get("decision", st["default"])
        for t in range(st["first"] - 1, st["last"]):
            c = chk["per"].get(t)
            rows.append([t + 1, st["id"], st["class"], reasons[t],
                         round(float(rel_level[t]) * 100, 1) if rel_level is not None else "",
                         round(float(rel_shape[t]) * 100, 1) if rel_shape is not None else "",
                         round(c[1] * 100, 1) if c else "", round(c[0] * 100, 2) if c else "",
                         round(chk["refs"][st["id"]]["floor"] * 100, 2) if c and st["id"] in chk["refs"] else "",
                         round(c[2], 2) if c else "",
                         round(float(speed[t]), 1) if speed is not None and t < speed.size else "",
                         int((t + 1) in core), st.get("verdict", "") if c else "recorded frame not on disk",
                         decision, int(excluded[t]) if t < excluded.size else ""])
    rows.sort(key=lambda r: r[0])
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "stretch", "class", "reason", "level_dev_pct", "shape_dev_pct", "whole_field_change_pct",
                    "diff_from_neighbours_rms_pct", "clean_neighbour_floor_rms_pct", "ratio", "speed",
                    "flagged_frame", "stretch_verdict", "decision", "excluded"])
        w.writerows(rows)


def check_sheet(out_png: Path, title: str, stretches: list[dict], chk: dict, info: dict, spf: float,
                excluded: np.ndarray, n_show: int = 8, context: int = 25, log=print):
    """The overview sheet (<base>_flagged_frames.png): per stretch its traces
    (stretch shaded, flagged frames ticked), a few of its frames between kept
    neighbours and their smoothed % difference from the clean neighbours.
    The frame stacks (check_stacks) are the frame-by-frame view; this is the
    one-glance summary."""
    from scipy.ndimage import gaussian_filter
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
        from matplotlib.gridspec import GridSpec
    except Exception as e:
        log(f"    warning: matplotlib not usable, no sheet: {e}")
        return
    T = chk["T"]
    frames = chk["frames"]; per = chk["per"]; smooth_px = chk["smooth_px"]
    rel_level = info.get("rel"); rel_shape = info.get("shape_rel")
    excluded = np.asarray(excluded, bool)
    sheet = []
    for st in stretches:
        ref = chk["refs"].get(st["id"])
        if ref is None:
            continue
        s, e = st["first"] - 1, st["last"] - 1
        inside = [t for t in range(s, e + 1) if t in per]
        core = [int(f) - 1 for f in (st.get("core_frames") or []) if (int(f) - 1) in per]
        pick = list(core[:n_show])
        rest = [t for t in inside if t not in pick]
        n_more = max(0, min(n_show, len(inside)) - len(pick))
        if rest and n_more:
            pick += [rest[int(round(x))] for x in np.linspace(0, len(rest) - 1, min(n_more, len(rest)))]
        if inside:
            t_max = max(inside, key=lambda t: per[t][2])
            if t_max not in pick:
                pick.append(t_max)
        pick = sorted(set(pick))[:n_show + 1]
        panels = [("kept", int(i), frames[int(i)], ref["M_all"]) for i in ref["lo"][-2:]]
        panels += [(st["class"], t, frames[t], _ref_at(ref, t)) for t in pick]
        panels += [("kept", int(i), frames[int(i)], ref["M_all"]) for i in ref["hi"][:2]]
        sheet.append((st, panels, ref))
    if not sheet:
        return
    try:
        colours = {"artefact": "#d62728", "state": "#e08a00", "kept": "#1a5e1a"}
        Hd, Wd = chk["shape_ds"]
        ncol = max(len(p) for (_, p, _) in sheet)
        panel_w = 1.45
        panel_h = panel_w * (Hd / max(Wd, 1)) + 0.05
        trace_h = 1.15
        fig_w = panel_w * ncol + 1.7
        fig_h = len(sheet) * (trace_h + 2 * panel_h + 1.0) + 0.7
        fig = plt.figure(figsize=(fig_w, fig_h))
        gs = GridSpec(3 * len(sheet), ncol, figure=fig, height_ratios=[trace_h, panel_h, panel_h] * len(sheet),
                      left=0.05, right=0.9, top=1 - 0.55 / fig_h, bottom=0.25 / fig_h, hspace=0.55, wspace=0.06)
        t_min = np.arange(T) * spf / 60.0
        for j, (st, panels, ref) in enumerate(sheet):
            s, e = st["first"] - 1, st["last"] - 1
            cls = st["class"]; col = colours[cls]
            decision = st.get("decision", st["default"])
            core = set(int(f) - 1 for f in (st.get("core_frames") or []))
            tax = fig.add_subplot(gs[3 * j, :])
            a0, b0 = max(0, s - context), min(T - 1, e + context)
            seg = slice(a0, b0 + 1)
            tax.plot(t_min[seg], rel_level[seg] * 100, color="k", lw=1.0, label="level vs moving median")
            if rel_shape is not None:
                tax.plot(t_min[seg], rel_shape[seg] * 100, color="tab:orange", lw=1.0, label="p10/median vs baseline")
            tax.axvspan(t_min[s] - spf / 120, t_min[e] + spf / 120, color=col, alpha=0.18, lw=0)
            for t in range(a0, b0 + 1):
                if excluded[t]:
                    tax.axvspan(t_min[t] - spf / 120, t_min[t] + spf / 120, ymin=0, ymax=0.06, color="#222222", lw=0)
                if t in core:
                    tax.axvspan(t_min[t] - spf / 120, t_min[t] + spf / 120, ymin=0.94, ymax=1.0, color="#d62728", lw=0)
            for (kind, t, _, _) in panels:
                tax.axvline(t_min[t], color=colours[kind], lw=0.6, alpha=0.5)
            tax.axhline(0, color="#999999", lw=0.5)
            tax.set_xlim(t_min[a0], t_min[b0])
            tax.set_ylabel("%", fontsize=7); tax.tick_params(labelsize=7)
            tax.legend(fontsize=6, loc="upper right", ncol=2, frameon=False)
            head = (f"{cls.upper()} frames {st['id']} ({st['n_frames']} frames, {st['seconds']:.0f} s)  —  "
                    f"{st.get('verdict', '')}  —  decision: {describe_decision(st, decision).upper()}"
                    + ("" if normalise_decision(decision) == normalise_decision(st["default"]) else " (yours)"))
            tax.set_title(head, fontsize=8.5, fontweight="bold", loc="left", color=col)
            vmin, vmax = np.percentile(ref["M_all"], [1, 99.5])
            sms, lim = [], 20.0
            for kind, t, img, rimg in panels:
                sm = gaussian_filter((img - rimg) / np.maximum(rimg, 1e-6), smooth_px) * 100
                sms.append(sm)
                if kind != "kept":
                    lim = max(lim, min(100.0, float(np.percentile(np.abs(sm), 99)) * 1.1))
            norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
            im = None
            for i, (kind, t, img, rimg) in enumerate(panels):
                ax = fig.add_subplot(gs[3 * j + 1, i])
                ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax, interpolation="nearest")
                ax.set_xticks([]); ax.set_yticks([])
                c = colours[kind]
                if kind == "kept":
                    ax.set_title(f"neighbour f{t + 1}", fontsize=7, color=c, pad=2)
                else:
                    dev, mean_rel, r = per.get(t, (np.nan, np.nan, np.nan))
                    mark = " FLAGGED" if t in core else ""
                    ax.set_title(f"f{t + 1}{mark}: field {mean_rel * 100:+.0f}%, x{r:.1f}", fontsize=7, color=c, pad=2,
                                 fontweight="bold" if t in core else "normal")
                for sp in ax.spines.values():
                    sp.set_edgecolor(c); sp.set_linewidth(2.2 if t in core and kind != "kept" else 1.3)
                ax2 = fig.add_subplot(gs[3 * j + 2, i])
                im = ax2.imshow(sms[i], cmap="RdBu_r", norm=norm, interpolation="nearest")
                ax2.set_xticks([]); ax2.set_yticks([])
                for sp in ax2.spines.values():
                    sp.set_edgecolor("#bbbbbb")
            if im is not None:
                bb = fig.axes[-1].get_position()
                cax = fig.add_axes([bb.x1 + 0.012, bb.y0, 0.012, bb.height])
                cb = fig.colorbar(im, cax=cax)
                cb.set_label("% vs clean neighbours", fontsize=7)
                cb.ax.tick_params(labelsize=7)
        fig.suptitle(f"{title} — flagged stretches (overview; the _check_frames.tif stack has every frame)\n"
                     f"trace: stretch shaded, red ticks = flagged frames, black bar = excluded; middle: recorded "
                     f"brightness (one grey scale per stretch); bottom: smoothed % difference from the clean neighbours "
                     f"(red brighter, blue darker); x = that difference relative to how much two clean frames differ",
                     fontsize=8.5, y=0.998)
        fig.savefig(out_png, dpi=90)
        plt.close(fig)
    except Exception as e:
        log(f"    warning: could not draw {Path(out_png).name}: {e}")


# ----------------------------------------------------------------------------- the ImageJ check stacks

def _font(size: int):
    from PIL import ImageFont
    try:
        from matplotlib import font_manager
        return ImageFont.truetype(font_manager.findfont("DejaVu Sans"), size)
    except Exception:
        try:
            return ImageFont.load_default(size=size)
        except TypeError:
            return ImageFont.load_default()


def _text_size(draw, txt, font):
    try:
        l, t, r, b = draw.textbbox((0, 0), txt, font=font)
        return r - l, b - t
    except Exception:
        return draw.textsize(txt, font=font)


def _wrap(draw, txt: str, font, width: int, sep: str = "  |  ") -> list[str]:
    """Split a ' | '-separated line into lines that fit `width` pixels."""
    parts = txt.split(sep)
    lines, cur = [], ""
    for p in parts:
        cand = p if not cur else cur + sep + p
        if cur and _text_size(draw, cand, font)[0] > width:
            lines.append(cur); cur = p
        else:
            cur = cand
    if cur:
        lines.append(cur)
    return lines or [""]


def _merge_segments(stretches: list[dict], T: int, context: int) -> list[tuple[int, int]]:
    segs = sorted((max(0, st["first"] - 1 - context), min(T - 1, st["last"] - 1 + context)) for st in stretches)
    out: list[list[int]] = []
    for a, b in segs:
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def _trace_panel(draw, x0, y0, w, h, a, b, t_cur, rel, srel, stretches, core, excluded, font_small):
    """Level / shape traces of frames a..b (0-based) with the stretches shaded,
    flagged frames ticked, excluded frames barred and the cursor on t_cur."""
    pal = {"artefact": (247, 212, 212), "state": (251, 230, 194)}
    n = b - a + 1
    px = w / max(n, 1)

    def xf(t):
        return x0 + (t - a + 0.5) * px

    vals = np.concatenate([rel[a:b + 1] * 100, (srel[a:b + 1] * 100) if srel is not None else np.zeros(0)])
    vals = vals[np.isfinite(vals)]
    lo_v = min(-5.0, float(vals.min()) - 2) if vals.size else -5.0
    hi_v = max(5.0, float(vals.max()) + 2) if vals.size else 5.0

    def yf(v):
        return y0 + (hi_v - v) / (hi_v - lo_v) * h

    draw.rectangle([x0, y0, x0 + w, y0 + h], fill=(255, 255, 255), outline=(160, 160, 160))
    for st in stretches:
        s, e = st["first"] - 1, st["last"] - 1
        if e < a or s > b:
            continue
        s, e = max(s, a), min(e, b)
        draw.rectangle([xf(s) - px / 2, y0 + 1, xf(e) + px / 2, y0 + h - 1], fill=pal[st["class"]])
    for t in range(a, b + 1):
        if excluded[t]:
            draw.rectangle([xf(t) - px / 2, y0 + h - 6, xf(t) + px / 2, y0 + h - 1], fill=(30, 30, 30))
        if core[t]:
            draw.polygon([(xf(t), y0 + 8), (xf(t) - 4, y0 + 1), (xf(t) + 4, y0 + 1)], fill=(214, 39, 40))
    draw.line([x0, yf(0), x0 + w, yf(0)], fill=(170, 170, 170), width=1)
    pts = [(xf(t), yf(float(rel[t]) * 100)) for t in range(a, b + 1) if np.isfinite(rel[t])]
    if srel is not None:
        pts2 = [(xf(t), yf(float(srel[t]) * 100)) for t in range(a, b + 1) if np.isfinite(srel[t])]
        if len(pts2) > 1:
            draw.line(pts2, fill=(224, 138, 0), width=2)
    if len(pts) > 1:
        draw.line(pts, fill=(0, 0, 0), width=2)
    draw.line([xf(t_cur), y0, xf(t_cur), y0 + h], fill=(42, 120, 214), width=2)
    step = 5 if px >= 12 else 10 if px >= 6 else 20
    for t in range(a, b + 1):
        if t == t_cur or ((t + 1) % step == 0 and abs(t - t_cur) * px >= 22):
            col = (42, 120, 214) if t == t_cur else (90, 90, 90)
            draw.line([xf(t), y0 + h, xf(t), y0 + h + 3], fill=col, width=1)
            txt = str(t + 1)
            tw, th = _text_size(draw, txt, font_small)
            draw.text((xf(t) - tw / 2, y0 + h + 4), txt, fill=col, font=font_small)
    for v in (lo_v, 0.0, hi_v):
        txt = f"{v:+.0f}%" if v else "0"
        tw, th = _text_size(draw, txt, font_small)
        draw.text((x0 - tw - 3, yf(v) - th / 2), txt, fill=(90, 90, 90), font=font_small)


def check_stacks(tif_path: Path, stretches: list[dict], present: np.ndarray, file_gain: np.ndarray, spf: float,
                 rel: np.ndarray, srel: np.ndarray | None, chk: dict, excluded: np.ndarray, out_frames: Path,
                 out_movie: Path, context: int = 6, ds: int = 1, raw_lost: np.ndarray | None = None,
                 title: str = "", movie_ds: int = 3, log=print) -> dict:
    """Two ImageJ stacks of the flagged stretches, each with `context` frames on
    both sides (overlapping stretches merged); slice i of both is the same
    recorded frame:
      out_frames  16-bit: the recorded frames at their recorded brightness
                  (full resolution, or block-averaged by ds), each with an
                  ImageJ slice label and a burned-in text line on top:
                  frame number, stretch, class, whole-field change, difference
                  from the neighbours (x the clean floor), flagged or not,
                  EXCLUDED or kept;
      out_movie   RGB composite (frame block-averaged by movie_ds): the frame
                  (one fixed grey scale for the whole stack, so a flash shows
                  as brighter), its smoothed % difference from the clean
                  neighbours (red brighter, blue darker, +/-20 %), and the
                  level / shape traces of the segment with the cursor on this
                  frame, stretches shaded, flagged frames ticked and excluded
                  frames barred.
    Scroll through them in ImageJ; the label tells which recorded frame the
    slice is."""
    import tifffile
    from PIL import Image, ImageDraw
    from scipy.ndimage import gaussian_filter
    T = int(len(rel))
    present = np.asarray(present, bool)
    raw_lost = np.zeros(T, bool) if raw_lost is None else np.asarray(raw_lost, bool)
    excluded = np.asarray(excluded, bool)
    core = np.zeros(T, bool)
    of_stretch: dict[int, dict] = {}
    for st in stretches:
        for f in st.get("core_frames") or []:
            if 1 <= int(f) <= T:
                core[int(f) - 1] = True
        for t in range(st["first"] - 1, st["last"]):
            of_stretch[t] = st
    segs = _merge_segments(stretches, T, context)
    order = [(seg, t) for seg in segs for t in range(seg[0], seg[1] + 1) if present[t]]
    if not order:
        return {"n_slices": 0}
    ds_chk = int(chk.get("ds", 4)); smooth_px = float(chk.get("smooth_px", 4.0))
    per = chk.get("per", {}); refs = chk.get("refs", {})

    def nearest_ref(t):
        best, dist = None, None
        for st in stretches:
            if st["id"] not in refs:
                continue
            d = 0 if st["first"] - 1 <= t <= st["last"] - 1 else min(abs(t - (st["first"] - 1)), abs(t - (st["last"] - 1)))
            if dist is None or d < dist:
                best, dist = st["id"], d
        return best

    font = _font(13); font_small = _font(10); font_cap = _font(11)
    src = _FrameSource(tif_path, present, file_gain, ds_chk)
    src.cache.update(chk.get("frames", {}))
    try:
        H, W = src.H, src.W
        Hs, Ws = (H // ds, W // ds) if ds > 1 else (H, W)
        strip = 20
        movie_ds = max(1, int(movie_ds))
        Hm, Wm = H // movie_ds, W // movie_ds
        disp = chk.get("display") or (0.0, 1.0)
        vmin, vmax = float(disp[0]), max(float(disp[1]), float(disp[0]) + 1.0)
        lim = 20.0
        pad = 8; th = 16; ch = 14; Htr = 96; ax_h = 18; n_text = 4
        Wc = max(2 * Wm + 3 * pad + 30, 600)              # room for the text lines on small movies
        Hc = pad + n_text * th + 4 + ch + Hm + pad + Htr + ax_h + pad
        frames16 = np.zeros((len(order), Hs + strip, Ws), np.uint16)
        movie = np.zeros((len(order), Hc, Wc, 3), np.uint8)
        labels = []
        for i, (seg, t) in enumerate(order):
            raw = src.raw(t)
            st = of_stretch.get(t)
            c = per.get(t)
            if c is None:
                rid = nearest_ref(t)
                if rid is not None:
                    img_ds = src.get(t)
                    r_img = _ref_at(refs[rid], t)
                    relimg = (img_ds - r_img) / np.maximum(r_img, 1e-6)
                    sm = gaussian_filter(relimg, smooth_px)
                    c = (float(np.sqrt(np.mean(sm ** 2))), float(np.mean(relimg)), float(np.sqrt(np.mean(sm ** 2)) / refs[rid]["floor"]))
                else:
                    img_ds = src.get(t); sm = np.zeros_like(img_ds)
            else:
                img_ds = src.get(t)
                rid = st["id"] if st is not None and st["id"] in refs else nearest_ref(t)
                r_img = _ref_at(refs[rid], t) if rid is not None else img_ds
                sm = gaussian_filter((img_ds - r_img) / np.maximum(r_img, 1e-6), smooth_px)
            numbers = f"field {c[1] * 100:+.1f}%, diff {c[0] * 100:.1f}% = x{c[2]:.1f} clean" if c is not None else "not compared"
            if st is not None:
                what = f"{st['id']} {st['class'].upper()} ({st['n_frames']} frames)"
                flag = "FLAGGED FRAME" if core[t] else "not flagged"
                fate = "EXCLUDED" if excluded[t] else "kept"
                line1 = f"f{t + 1}  |  {what}  |  {numbers}  |  {flag}  ->  {fate}"
                dec = st.get("decision", st["default"])
                line2 = (f"{st.get('verdict', '')}  |  decision: {describe_decision(st, dec)}"
                         + ("" if normalise_decision(dec) == normalise_decision(st["default"]) else " (yours)"))
                label = (f"f{t + 1} | {st['id']} {st['class']} | {c[1] * 100:+.0f}% x{c[2]:.1f} | " if c is not None
                         else f"f{t + 1} | {st['id']} {st['class']} | ") + ("flagged " if core[t] else "") + fate
            else:
                line1 = f"f{t + 1}  |  context  |  {numbers}"
                line2 = "kept frame between / around the flagged stretches"
                label = f"f{t + 1} | context"
            if raw_lost[t]:
                line1 += "  (movie holds a filler here — content gone)"
                label += " | filler"
            labels.append(label)
            # --- 16-bit frame with a burned-in text line
            fr = _ds_page(raw, ds) if ds > 1 else raw
            frames16[i, strip:, :] = np.clip(np.round(fr), 0, 65535).astype(np.uint16)
            bright = float(np.percentile(fr, 99.9))
            tim = Image.new("L", (Ws, strip), 0)
            d = ImageDraw.Draw(tim)
            d.text((4, 2), line1, fill=255, font=font)
            frames16[i, :strip, :] = (np.asarray(tim, np.float32) / 255.0 * bright).astype(np.uint16)
            # --- composite
            im = Image.new("RGB", (Wc, Hc), (255, 255, 255))
            d = ImageDraw.Draw(im)
            l1 = _wrap(d, line1, font, Wc - 2 * pad)
            lines = l1 + _wrap(d, line2, font_cap, Wc - 2 * pad)
            y = pad
            for k, ln in enumerate(lines[:n_text]):
                first = k < len(l1)
                d.text((pad, y), ln, fill=(180, 0, 0) if (first and excluded[t]) else (0, 0, 0) if first else (60, 60, 60),
                       font=font if first else font_cap)
                y += th
            y_cap = pad + n_text * th + 4
            y_img = y_cap + ch
            x1 = pad + 30; x2 = x1 + Wm + pad
            if Wm >= 200:
                d.text((x1, y_cap), "recorded frame (fixed grey scale)", fill=(60, 60, 60), font=font_cap)
                d.text((x2, y_cap), f"% vs clean neighbours (red +, blue -, ±{lim:.0f}%)", fill=(60, 60, 60), font=font_cap)
            else:
                d.text((x1, y_cap), f"frame | % vs neighbours (±{lim:.0f}%)", fill=(60, 60, 60), font=font_cap)
            img_m = _ds_page(raw, movie_ds) if movie_ds > 1 else raw
            g = np.clip((img_m - vmin) / (vmax - vmin), 0, 1)
            g8 = (g * 255).astype(np.uint8)
            im.paste(Image.fromarray(np.stack([g8] * 3, axis=-1), "RGB"), (x1, y_img))
            v = np.clip(sm * 100 / lim, -1, 1)
            rgb = np.empty(v.shape + (3,), np.uint8)
            pos = np.clip(v, 0, 1); neg = np.clip(-v, 0, 1)
            rgb[..., 0] = (255 * (1 - neg)).astype(np.uint8)
            rgb[..., 1] = (255 * (1 - np.maximum(pos, neg))).astype(np.uint8)
            rgb[..., 2] = (255 * (1 - pos)).astype(np.uint8)
            dimg = Image.fromarray(rgb, "RGB")
            if dimg.size != (Wm, Hm):
                dimg = dimg.resize((Wm, Hm), Image.BILINEAR)
            im.paste(dimg, (x2, y_img))
            if st is not None:
                colr = (214, 39, 40) if st["class"] == "artefact" else (224, 138, 0)
                d.rectangle([x1 - 2, y_img - 2, x1 + Wm + 1, y_img + Hm + 1], outline=colr, width=3 if core[t] else 2)
            y_tr = y_img + Hm + pad
            _trace_panel(d, x1, y_tr, Wc - x1 - pad, Htr, seg[0], seg[1], t, rel, srel, stretches, core, excluded, font_small)
            movie[i] = np.asarray(im, np.uint8)
    finally:
        src.close()
    gb16 = frames16.nbytes / 1e9
    kw16 = dict(imagej=True, metadata={"axes": "TYX", "Labels": labels}) if gb16 < 3.9 else dict(bigtiff=True)
    tifffile.imwrite(str(out_frames), frames16, **kw16)
    kwm = (dict(imagej=True, photometric="rgb", metadata={"axes": "TYXS", "Labels": labels})
           if movie.nbytes / 1e9 < 3.9 else dict(bigtiff=True, photometric="rgb"))
    tifffile.imwrite(str(out_movie), movie, **kwm)
    return {"n_slices": len(order), "segments": [(a + 1, b + 1) for a, b in segs],
            "mb": round((frames16.nbytes + movie.nbytes) / 1e6)}


def make_check_stacks(tif_path: Path, context: int = 6, ds: int = 1, log=print) -> dict | None:
    """Rebuild the two check stacks of a corrected movie from its <base>.json /
    _frames.csv (no re-flagging; the current decisions are shown). Reads the
    flagged frames from the movie again."""
    tif_path = Path(tif_path)
    base = tif_path.with_suffix("")
    fr = read_frames_csv(Path(str(base) + "_frames.csv"))
    json_path = Path(str(base) + ".json")
    if fr is None or not json_path.exists():
        log(f"  {tif_path.name}: no _frames.csv / .json — run the brightness step first")
        return None
    try:
        meta = json.loads(json_path.read_text(encoding="utf-8")) or {}
    except Exception:
        meta = {}
    gc = meta.get("global_correction") or {}
    stretches = gc.get("stretches") or []
    if not stretches:
        log(f"  {tif_path.name}: nothing flagged — no stacks to make")
        return {"n_slices": 0}
    T = fr["frame"].size
    spf = float(meta.get("seconds_per_frame") or 2.0)
    base_lvl = fr.get("base"); shape = fr.get("shape"); sbase = fr.get("sbase")
    if base_lvl is None or not np.isfinite(base_lvl).any():
        base_lvl = _moving_median(fr["level"], 31)
    rel = fr["level"] / np.maximum(base_lvl, 1e-9) - 1.0
    srel = None
    if shape is not None and sbase is not None and np.isfinite(sbase).any():
        srel = shape / np.maximum(np.where(np.isfinite(sbase), sbase, np.nanmedian(sbase)), 1e-9) - 1.0
    present = ~fr["dropped"]
    raw_lost = fr["dropped"] | fr["replaced"]
    import tifffile
    with tifffile.TiffFile(str(tif_path)) as tf:
        n_pages = len(tf.pages)
    if n_pages != int(present.sum()):
        log(f"  {tif_path.name}: {n_pages} pages but {int(present.sum())} frames in _frames.csv — re-run the brightness step")
        return None
    params = gc.get("params") or {}
    info = {"rel": rel, "shape_rel": srel}
    reasons = fr["reason"]
    chk = check_frames(tif_path, stretches, present, fr["gain"], spf, info, reasons,
                       behavior_csv=Path(str(base) + "_behavior.csv"), raw_lost=raw_lost,
                       core_ratio=float(params.get("core_ratio", 3.0)), core_min_pct=float(params.get("core_min_pct", 10.0)),
                       log=log)
    dec = read_decisions(tif_path)
    excluded, applied = excluded_frames(stretches, dec, T)
    for st, d in zip(stretches, applied):
        st["decision"] = d
    out_frames = Path(str(base) + "_check_frames.tif"); out_movie = Path(str(base) + "_check_movie.tif")
    r = check_stacks(tif_path, stretches, present, fr["gain"], spf, rel, srel, chk, excluded, out_frames, out_movie,
                     context=context, ds=ds, raw_lost=raw_lost, title=tif_path.name, log=log)
    log(f"  {tif_path.name}: {out_frames.name} + {out_movie.name}: {r.get('n_slices', 0)} slices "
        f"({r.get('mb', 0)} MB), segments {', '.join(f'{a}-{b}' for a, b in r.get('segments', []))}")
    return r


# ----------------------------------------------------------------------------- the whole step

FLAG_PARAMS = ("win_s", "k", "min_rel", "corr_k", "min_corr_drop", "grow", "shape_rel", "shape_win_s",
               "min_clean_s", "max_spike_s", "core_ratio", "core_min_pct")


def apply(tif_path: Path, spf: float | None = None, win_s: float = 60.0, k: float = 4.0,
          min_rel: float = 0.03, corr_k: float = 5.0, min_corr_drop: float = 0.05, grow: int = 1,
          shape_rel: float = 0.05, shape_win_s: float = 300.0, min_clean_s: float = 10.0,
          max_spike_s: float = 10.0, n_fill: int = 10, mode: str = "keep", normalise: bool = True,
          force: bool = False, decisions: dict | None = None, report_only: bool = False,
          core_ratio: float = 3.0, core_min_pct: float = 10.0, stacks: bool = True, stack_context: int = 6,
          stack_ds: int = 1, log=print) -> dict | None:
    """The brightness step on one binned tiff, in place. Returns the info
    written to the json sidecar, or None when skipped.

    mode "keep" (default): every frame stays in the movie; the level is
    flattened (normalise) and the flags, the per-frame check, the decisions
    and the check files (sheet + the two ImageJ stacks) are written.
    mode "drop": the excluded frames are left out of the movie (tiff_frame
    maps them); "fill": replaced by a frozen filler.
    report_only: flags, csv, QC plot and the check files only; the movie and
    its json are not touched.
    core_ratio / core_min_pct: the rule for the FLAGGED frames of an artefact
    stretch (the frames 'exclude' excludes): differ from the clean neighbours
    by >= core_min_pct % of the frame and >= core_ratio times the clean floor.
    stacks: write <base>_check_frames.tif / _check_movie.tif (stretches +/-
    stack_context frames; the 16-bit stack block-averaged by stack_ds).
    On a tiff corrected before, `force` re-flags using the raw per-frame
    statistics recorded in <base>_frames.csv; frames dropped or replaced
    earlier stay out (their content is gone). A companion <base>_ch2.tif loses
    the same frames in drop mode."""
    if mode not in ("drop", "fill", "keep"):
        raise ValueError(f"mode must be keep, drop or fill (got {mode!r})")
    tif_path = Path(tif_path)
    base = tif_path.with_suffix("")
    json_path = Path(str(base) + ".json")
    frames_csv = Path(str(base) + "_frames.csv")
    meta = {}
    if json_path.exists():
        try:
            meta = json.loads(json_path.read_text())
        except Exception:
            meta = {}
    already = meta.get("global_correction")
    if already and not force and not report_only:
        log(f"  {tif_path.name}: already corrected on {already.get('written', '?')} — skipping (force re-does it)")
        return None
    spf = float(spf or meta.get("seconds_per_frame") or 2.0)
    prev = read_frames_csv(frames_csv) if already else None
    import tifffile
    with tifffile.TiffFile(str(tif_path)) as tf:
        n_pages = len(tf.pages)
    ok = bool(already) and prev is not None and np.isfinite(prev["level"]).all() \
        and int((~prev["dropped"]).sum()) == n_pages
    tr = None if ok else global_trace(tif_path, log=log)        # the raw statistics are on file: no full read
    T = n_pages
    old_gain = np.ones(T)
    present = np.ones(T, bool)
    gone = np.zeros(T, bool)
    if already:
        if ok:
            T = prev["frame"].size
            present = ~prev["dropped"]
            gone = prev["dropped"] | prev["replaced"]
            old_gain = prev["gain"]
            tr = dict(T=T, median=prev["level"], mean=prev["mean"], p10=prev["p10"], p90=prev["p90"],
                      corr_ref=prev["corr"] if np.isfinite(prev["corr"]).all() else np.full(T, np.nan))
            if not np.isfinite(tr["corr_ref"]).all():
                c = np.full(T, np.nan); c[present] = global_trace(tif_path, log=log)["corr_ref"]
                tr["corr_ref"] = np.where(np.isfinite(c), c, np.nanmedian(c))
            log(f"  {tif_path.name}: corrected on {already.get('written', '?')} — re-flagging from the raw levels "
                f"recorded in {frames_csv.name}")
        else:
            log(f"  {tif_path.name}: corrected on {already.get('written', '?')} but {frames_csv.name} is missing or "
                f"does not match the movie — re-flagging from the movie itself (frames listed in the json stay out)")
            if already.get("bad_frames", "filled" if already.get("replaced_flagged") else "kept") != "kept":
                for i in already.get("excluded_frames_1based", already.get("flagged_frames_1based", [])):
                    if 1 <= int(i) <= T:
                        gone[int(i) - 1] = True
    artefact, reasons, info = flag_frames(tr["median"], tr["corr_ref"], spf, win_s=win_s, k=k, min_rel=min_rel,
                                          corr_k=corr_k, min_corr_drop=min_corr_drop, grow=grow, p10=tr["p10"],
                                          shape_rel=shape_rel, shape_win_s=shape_win_s, min_clean_s=min_clean_s,
                                          max_spike_s=max_spike_s)
    state = np.asarray(info["state"], bool)
    if gone.any():
        for i in np.flatnonzero(gone & ~artefact):
            reasons[i] = "removed in an earlier run (content no longer on disk)"
        artefact = artefact | gone
        state = state & ~artefact
        log(f"  {tif_path.name}: {int(gone.sum())} frame(s) are no longer on disk (dropped or replaced by an earlier "
            f"run) and stay excluded; to get them back restore {Path(str(base) + '_raw.tif').name} "
            f"(--restore-raw) or re-unpack this run")
    stretches = build_stretches(artefact, state, info, spf)
    # the frame-by-frame check, while every frame is still here: which frames of
    # each artefact stretch really differ from their neighbours (-> defaults)
    beh = Path(str(base) + "_behavior.csv")
    chk = None
    if stretches:
        try:
            chk = check_frames(tif_path, stretches, present, old_gain, spf, info, reasons,
                               behavior_csv=beh if beh.exists() else None, raw_lost=gone,
                               core_ratio=core_ratio, core_min_pct=core_min_pct, log=log)
        except Exception as e:
            log(f"    warning: the frame-by-frame check failed ({e}); artefact stretches are excluded whole")
            for st in stretches:
                st.pop("core_frames", None)
            chk = None
    dec = read_decisions(tif_path)
    if decisions:
        dec.update({str(kk): normalise_decision(v) for kk, v in decisions.items()})
    if dec and not report_only:
        write_decisions(tif_path, dec)
    excluded, applied = excluded_frames(stretches, dec, T)
    stretch_of = [""] * T
    decision_of = [""] * T
    core = np.zeros(T, bool)
    for st, d in zip(stretches, applied):
        st["decision"] = d
        for t in range(st["first"] - 1, st["last"]):
            stretch_of[t] = st["id"]; decision_of[t] = d
        for f in st.get("core_frames") or []:
            if 1 <= int(f) <= T:
                core[int(f) - 1] = True
    G = float(np.median(tr["median"][~artefact])) if (~artefact).any() else float(np.median(tr["median"]))
    gain = (G / np.maximum(tr["median"], 1e-9)) if normalise else np.ones(T)
    n_art, n_state, n_ex = int(artefact.sum()), int(state.sum()), int(excluded.sum())
    log(f"  {tif_path.name}: level {G:.0f}, robust SD {info['robust_sd'] * 100:.2f}%; "
        f"artefact {n_art} frame(s)" + (f" [{_runs_str(artefact)}]" if n_art else "")
        + f"; state change {n_state} frame(s)" + (f" [{_runs_str(state)}]" if n_state else "")
        + f"; excluded from the statistics: {n_ex} frame(s)" + (f" [{frames_str(np.flatnonzero(excluded) + 1)}]" if n_ex else ""))
    for st in stretches:
        own = normalise_decision(st["decision"]) != normalise_decision(st["default"])
        log(f"    {st['class']:8s} {st['id']:>9s} ({st['n_frames']} frames): {st.get('verdict', '')}"
            f" -> {describe_decision(st, st['decision'])}" + (" (your decision)" if own else ""))
    # the check files
    check_files = []
    if stretches and chk is not None:
        try:
            write_check_csv(Path(str(base) + "_flagged_frames.csv"), stretches, chk, reasons, info, excluded)
            check_sheet(Path(str(base) + "_flagged_frames.png"), tif_path.name, stretches, chk, info, spf, excluded, log=log)
            check_files += [Path(str(base) + "_flagged_frames.png").name, Path(str(base) + "_flagged_frames.csv").name]
        except Exception as e:
            log(f"    warning: could not write the check sheet / csv: {e}")
        if stacks:
            try:
                r = check_stacks(tif_path, stretches, present, old_gain, spf, info["rel"], info.get("shape_rel"), chk,
                                 excluded, Path(str(base) + "_check_frames.tif"), Path(str(base) + "_check_movie.tif"),
                                 context=stack_context, ds=stack_ds, raw_lost=gone, title=tif_path.name, log=log)
                check_files += [Path(str(base) + "_check_frames.tif").name, Path(str(base) + "_check_movie.tif").name]
                log(f"    check stacks: {r.get('n_slices', 0)} slices ({r.get('mb', 0)} MB) — "
                    f"{Path(str(base) + '_check_frames.tif').name} (16-bit, labelled) and _check_movie.tif (composite)")
            except Exception as e:
                log(f"    warning: could not write the check stacks: {e}")
    check_files = [f for f in check_files if (tif_path.parent / f).exists()]
    if report_only:
        tiff_frame = np.zeros(T, int); tiff_frame[present] = np.arange(1, int(present.sum()) + 1)
        write_frames_csv(frames_csv, spf, tr, artefact, state, stretch_of, decision_of, excluded, reasons,
                         old_gain, tiff_frame, info, replaced=prev["replaced"] if prev is not None else None,
                         core=core, check=chk["per"] if chk else None)
        qc_plot(Path(str(base) + "_globalfix.png"), spf, tr, artefact, state, excluded, info,
                beh if beh.exists() else None, tif_path.name + "  (report only — movie not modified)",
                mode=already.get("bad_frames", "keep") if already else "keep")
        return {"report_only": True, "n_artefact": n_art, "n_state": n_state, "n_excluded": n_ex, "stretches": stretches}
    factor = gain / np.maximum(old_gain, 1e-9)
    if already:
        factor[np.abs(factor - 1.0) < 5e-3] = 1.0
        gain = old_gain * factor
    remove = excluded if mode in ("drop", "fill") else np.zeros(T, bool)
    need_write = (not np.allclose(factor[present & ~remove], 1.0)
                  or (mode == "drop" and bool((remove & present).any()))
                  or (mode == "fill" and bool((remove | ~present).any()))
                  or (mode == "keep" and bool((~present).any())))
    if need_write:
        cinfo = correct_tiff(tif_path, factor, remove, mode=mode, present=present, n_fill=n_fill, log=log)
        tiff_frame = cinfo["tiff_frame"]
        what = {"drop": f"{cinfo['n_dropped']} frame(s) dropped from the movie, {cinfo['n_written']} written",
                "fill": f"{cinfo['n_filled']} frame(s) replaced by a filler",
                "keep": "all frames kept in the movie (level flattened)"}[mode]
        log(f"    {what}")
    else:
        tiff_frame = np.zeros(T, int); tiff_frame[present] = np.arange(1, int(present.sum()) + 1)
        cinfo = {"n_written": int(present.sum()), "n_dropped": int((~present).sum()), "n_filled": 0}
    ch2 = Path(str(base) + "_ch2.tif")
    ch2_note = None
    if mode == "drop" and ch2.exists() and cinfo["n_dropped"]:
        try:
            import tifffile
            with tifffile.TiffFile(str(ch2)) as tf:
                n2 = len(tf.pages)
            keep = None
            if n2 == T:
                keep = tiff_frame > 0
            elif n2 == int(present.sum()):
                keep = (tiff_frame > 0)[present]
            if keep is not None and not keep.all():
                drop_frames(ch2, keep, log=log)
                ch2_note = f"{int((~keep).sum())} frame(s) dropped from {ch2.name} to match"
                log(f"    {ch2_note}")
            elif keep is None:
                log(f"    warning: {ch2.name} has {n2} frames, expected {T}: not touched")
        except Exception as e:
            log(f"    warning: could not shorten {ch2.name}: {e}")
    # frames whose content on disk is a filler: earlier fillers (gone but present) + new ones in fill mode
    replaced = (gone & present & (tiff_frame > 0))
    if mode == "fill":
        replaced |= (remove | ~present) & (tiff_frame > 0)
    write_frames_csv(frames_csv, spf, tr, artefact, state, stretch_of, decision_of, excluded, reasons, gain,
                     tiff_frame, info, replaced=replaced, core=core, check=chk["per"] if chk else None)
    from datetime import datetime
    qc_plot(Path(str(base) + "_globalfix.png"), spf, tr, artefact, state, excluded, info,
            beh if beh.exists() else None, tif_path.name, mode=mode)
    rec = {"method": ("per-frame division by the spatial median (level flattened to the median of the clean frames)"
                      if normalise else "flagging only"),
           "rules": {"artefact": f"level jump between frames or spike (<= {max_spike_s:.0f} s, returning) > "
                                 f"max({k} robust SD, {min_rel * 100:.0f}%) vs a {win_s:.0f}-s moving median, or "
                                 f"corr_ref < median - max({corr_k} robust SD, {min_corr_drop}); tails while still "
                                 f"moving, grown by {grow}, clean islands < {min_clean_s:.0f} s",
                     "state": f"level plateau beyond the floor for > {max_spike_s:.0f} s or not returning, or "
                              f"|p10/median vs clean-frame moving median ({shape_win_s:.0f} s) - 1| > "
                              f"max({k} robust SD, {shape_rel * 100:.0f}%)",
                     "flagged_frames": f"a frame of an artefact stretch whose smoothed difference from the clean "
                                       f"neighbours is >= {core_min_pct:.0f}% of the frame and >= {core_ratio:.1f} x "
                                       f"the clean floor (plus frames no longer on disk)",
                     "default_decision": "artefact with flagged frames -> exclude (those frames only); artefact "
                                         "without -> keep; state -> keep",
                     "decisions": DECISION_MEANING},
           "bad_frames": {"drop": "dropped", "fill": "filled", "keep": "kept"}[mode],
           "n_frames": int(T), "frames_in_movie": int(cinfo["n_written"]),
           "n_artefact": n_art, "n_state": n_state, "n_excluded": n_ex,
           "artefact_frames_1based": [int(i) + 1 for i in np.flatnonzero(artefact)],
           "state_frames_1based": [int(i) + 1 for i in np.flatnonzero(state)],
           "core_frames_1based": [int(i) + 1 for i in np.flatnonzero(core)],
           "excluded_frames_1based": [int(i) + 1 for i in np.flatnonzero(excluded)],
           "stretches": stretches,
           "decisions_file": decisions_path(tif_path).name,
           "check_files": check_files or None,
           "movie_frame_of_original_frame": "column tiff_frame in " + frames_csv.name,
           "reference_level": G, "gain_range": [float(gain.min()), float(gain.max())],
           "params": {"win_s": win_s, "k": k, "min_rel": min_rel, "corr_k": corr_k, "min_corr_drop": min_corr_drop,
                      "grow": grow, "shape_rel": shape_rel, "shape_win_s": shape_win_s, "min_clean_s": min_clean_s,
                      "max_spike_s": max_spike_s, "core_ratio": core_ratio, "core_min_pct": core_min_pct,
                      "stack_context": stack_context},
           "written": datetime.now().isoformat(timespec="seconds")}
    # backwards-compatible summary fields
    rec["n_flagged"] = n_art
    rec["flagged_frames_1based"] = rec["artefact_frames_1based"]
    if ch2_note:
        rec["ch2"] = ch2_note
    if already:
        rec["previous"] = {"written": already.get("written"), "n_flagged": already.get("n_flagged"),
                           "bad_frames": already.get("bad_frames", "filled" if already.get("replaced_flagged") else "kept")}
    meta["global_correction"] = rec
    if cinfo["n_dropped"]:
        meta["frames_out"] = int(cinfo["n_written"])
        meta["frames_recorded"] = int(T)
    try:
        json_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    except Exception as e:
        log(f"  warning: could not update {json_path.name}: {e}")
    return rec


def restore_raw(tif_path: Path, log=print) -> bool:
    """Put the uncorrected movie (<base>_raw.tif, written by unpack --keep-raw)
    back in place of the corrected one and forget the correction, so the
    brightness step can be run again from scratch. Returns True when done."""
    tif_path = Path(tif_path)
    base = tif_path.with_suffix("")
    raw = Path(str(base) + "_raw.tif")
    if not raw.exists():
        log(f"  {tif_path.name}: no {raw.name} — re-unpack that run instead (unpack_bin.py --overwrite)")
        return False
    import shutil
    shutil.copy2(raw, tif_path)
    json_path = Path(str(base) + ".json")
    if json_path.exists():
        try:
            meta = json.loads(json_path.read_text())
            meta.pop("global_correction", None)
            if "frames_recorded" in meta:
                meta["frames_out"] = meta.pop("frames_recorded")
            json_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        except Exception as e:
            log(f"  warning: could not update {json_path.name}: {e}")
    for suffix in ("_frames.csv", "_globalfix.png", "_flagged_frames.png", "_flagged_frames.csv",
                   "_removed_frames.png", "_removed_frames.csv", "_check_frames.tif", "_check_movie.tif"):
        p = Path(str(base) + suffix)
        if p.exists():
            p.unlink()
    log(f"  {tif_path.name}: restored from {raw.name} (decisions file kept)")
    return True
