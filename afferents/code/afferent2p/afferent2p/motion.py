"""motion — how much the field moves in a run, what AQuA2's registration leaves
behind, and whether the detected events coincide with motion.

Why: AQuA2 registers the binned movie itself (reg.regCrossCorrelation: one
rigid, whole-pixel translation per frame, found by cross-correlation against
the mean of the first 10 frames, then the border is cropped), so the detection
runs on a steadied movie, not on the raw tiff. What that cannot undo: motion
WITHIN a 2-s bin (the raw frames are averaged before AQuA2 sees them — a
moving bin is a blurred frame), deformation that is not one translation of
the whole field (non-rigid), and z-motion. Any of these makes the brightness
at fibre edges flicker in step with the movement, and a flicker that is
coherent over a fibre-shaped patch can be called an event. This module
measures all of that per run from the files that exist:

  * the binned tiff (what AQuA2 was given): a rigid shift per frame against
    AQuA2's own reference (mean of the first 10 frames) by phase correlation
    (frames smoothed, background removed, 2x downsampled for speed, parabolic
    sub-pixel peak) -> dy, dx, |shift|, the |jump| between consecutive frames,
    the correlation to the reference and a sharpness index (mean gradient
    over mean brightness: a blurred frame scores lower);
  * the registered movie (res.datOrg1 of the _AQuA2.mat, when readable): the
    rigid shift that is still there after AQuA2 (residual) and a non-rigid
    index: the field is cut into 3 x 4 blocks, each block's own shift is
    measured, and the spread of the block shifts around the whole-field shift
    is the deformation AQuA2's single translation could not remove (its
    level), while the median movement of the blocks since the previous frame
    is the deformation that happened at that frame (the flicker-maker);
  * the behaviour csv: all of the above split into still and running frames,
    and the correlation of |shift| with running speed;
  * the events csv: how many event onsets fall within +-1 frame of a jump
    (>= 2 px between consecutive raw frames, or more when the run is restless:
    median + 3 MAD) or of a deformation jump (the blocks of the registered
    movie moved >= 1 px, or median + 3 MAD, since the previous frame), against
    the number expected if events and motion were unrelated (the share of
    such frames). Noise floor of the indices on these movies: ~0.2
    px (white noise at the real level gives 0.2; a 0.5 % stretch of the field
    gives 1.1). A ratio well above 1
    with several events says the detections follow the movement.

Files, next to the tiff:  <base>_motion.csv (one row per movie frame) and
<base>_motion.json (the summary) — rebuilt when the tiff or the AQuA2 result
is newer. aqua_analyze.py writes them for every run and gathers the summaries
in aqua_analysis/motion_summary.csv; the QC report shows the figure and the
numbers per recording and flags runs whose events follow the motion.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

JUMP_PX = 2.0            # |shift| change between consecutive frames that counts as a jump
NONRIGID_PX = 1.5        # non-rigid residual (level) worth a note
NONRIGID_JUMP_PX = 1.0   # change of the block shifts between consecutive registered frames that counts as motion
NEAR_FRAMES = 1          # +-frames around an onset that count as coincident
REF_FRAMES = 10          # AQuA2's reference: mean of the first 10 frames
MIN_STATE_FRAMES = 10    # still / running statistics need at least this many frames


# ----------------------------------------------------------------------------- shifts

def _down(img: np.ndarray, k: int) -> np.ndarray:
    if k <= 1:
        return np.asarray(img, np.float32)
    H, W = img.shape
    H2, W2 = (H // k) * k, (W // k) * k
    return np.asarray(img[:H2, :W2], np.float32).reshape(H2 // k, k, W2 // k, k).mean(axis=(1, 3))


def _prep(img: np.ndarray, sigma: float, hp: float, down: int = 1) -> np.ndarray:
    """Downsample by `down` (mean), smooth the shot noise, remove the slow
    background (vessels' shadow, vignetting); sigma/hp are in original pixels."""
    from scipy.ndimage import gaussian_filter
    f = _down(img, down)
    f = gaussian_filter(f, max(sigma / down, 0.5))
    f = f - gaussian_filter(f, max(hp / down, 2.0))
    return f


class _Matcher:
    """Phase correlation of frames against one reference (prepared once).
    `prepared=True` means the reference (and the frames given to shift()) are
    already prepared (_prep) and downsampled — for block matching, where one
    prepared frame serves every block."""

    def __init__(self, ref: np.ndarray, max_shift: float = 40.0, sigma: float = 1.5, hp: float = 20.0, down: int = 2,
                 prepared: bool = False):
        self.sigma, self.hp, self.down = sigma, hp, max(1, int(down))
        self.prepared = prepared
        r = np.asarray(ref, np.float32) if prepared else _prep(ref, sigma, hp, self.down)
        H, W = r.shape
        self.H, self.W = H, W
        wy = np.hanning(H).astype(np.float32); wx = np.hanning(W).astype(np.float32)
        self.win = wy[:, None] * wx[None, :]
        r = r * self.win
        r = r - r.mean()
        self.rn = float(np.linalg.norm(r)) + 1e-9
        self.FR = np.conj(np.fft.rfft2(r))
        self.ms = int(max(2, round(max_shift / self.down)))

    def prep(self, frame: np.ndarray) -> np.ndarray:
        return _prep(frame, self.sigma, self.hp, self.down)

    def shift(self, frame: np.ndarray) -> tuple[float, float, float]:
        """(dy, dx, corr) in pixels of the original frame: the frame's content sits
        dy rows and dx columns away from where it is in the reference."""
        m = (np.asarray(frame, np.float32) if self.prepared else self.prep(frame)) * self.win
        m = m - m.mean()
        mn = float(np.linalg.norm(m)) + 1e-9
        cc = np.fft.irfft2(np.fft.rfft2(m) * self.FR, s=(self.H, self.W))
        cc = np.fft.fftshift(cc) / (self.rn * mn)
        cy, cx = self.H // 2, self.W // 2
        ms = min(self.ms, cy - 1, cx - 1)
        sub = cc[cy - ms:cy + ms + 1, cx - ms:cx + ms + 1]
        p = np.unravel_index(int(np.argmax(sub)), sub.shape)
        dy, dx = float(p[0] - ms), float(p[1] - ms)
        peak = float(sub[p])

        def refine(a, b, c):
            d = a - 2 * b + c
            return 0.0 if d >= 0 else float(np.clip(0.5 * (a - c) / d, -1, 1))

        y0, x0 = p
        if 0 < y0 < sub.shape[0] - 1:
            dy += refine(sub[y0 - 1, x0], sub[y0, x0], sub[y0 + 1, x0])
        if 0 < x0 < sub.shape[1] - 1:
            dx += refine(sub[y0, x0 - 1], sub[y0, x0], sub[y0, x0 + 1])
        return dy * self.down, dx * self.down, peak


def whole_shift(matcher: "_Matcher", prepared_frame: np.ndarray) -> tuple[float, float, float]:
    """shift() for a frame that was prepared with matcher.prep() already."""
    saved = matcher.prepared
    matcher.prepared = True
    try:
        return matcher.shift(prepared_frame)
    finally:
        matcher.prepared = saved


def sharpness(img: np.ndarray, sigma: float = 1.0) -> float:
    """Mean gradient magnitude of the lightly smoothed frame over its mean
    brightness — a blurred (moving) frame scores lower than a steady one."""
    from scipy.ndimage import gaussian_filter
    f = gaussian_filter(np.asarray(img, np.float32), sigma)
    gy, gx = np.gradient(f)
    return float(np.hypot(gy, gx).mean() / max(float(f.mean()), 1e-6))


def frame_shifts(movie, ref: np.ndarray | None = None, max_shift: float = 40.0, down: int = 2,
                 with_sharpness: bool = True, progress=None) -> dict:
    """Rigid shift of every frame of `movie` ((T, H, W) array-like) against
    `ref` (default: mean of its first REF_FRAMES frames, as AQuA2 does).
    Returns {dy, dx, corr, sharp} arrays of length T."""
    T = int(movie.shape[0])
    if ref is None:
        n = min(REF_FRAMES, T)
        ref = np.mean(np.stack([np.asarray(movie[i], np.float32) for i in range(n)]), axis=0)
    M = _Matcher(ref, max_shift=max_shift, down=down)
    dy = np.zeros(T); dx = np.zeros(T); corr = np.zeros(T); sharp = np.full(T, np.nan)
    for t in range(T):
        fr = np.asarray(movie[t], np.float32)
        dy[t], dx[t], corr[t] = M.shift(fr)
        if with_sharpness:
            sharp[t] = sharpness(fr)
        if progress and t % 100 == 99:
            progress(t + 1, T)
    return {"dy": dy, "dx": dx, "corr": corr, "sharp": sharp}


def block_residual(movie, ref: np.ndarray | None = None, rows: int = 3, cols: int = 4, max_shift: float = 12.0,
                   down: int = 2) -> dict:
    """Whole-field shift and the non-rigid residual of every frame of a movie
    that is already registered: each of rows x cols blocks is matched on its
    own; nonrigid = median over blocks of |block shift - whole-field shift|,
    nonrigid_max = the 90th percentile over blocks. Blocks without structure
    (correlation < 0.2) are ignored."""
    T, H, W = (int(v) for v in movie.shape)
    if ref is None:
        n = min(REF_FRAMES, T)
        ref = np.mean(np.stack([np.asarray(movie[i], np.float32) for i in range(n)]), axis=0)
    whole = _Matcher(ref, max_shift=max_shift, down=down)
    ref_p = whole.prep(ref)                                   # prepared + downsampled once; blocks are cut from it
    Hp, Wp = ref_p.shape
    edges_y = np.linspace(0, Hp, rows + 1).astype(int); edges_x = np.linspace(0, Wp, cols + 1).astype(int)
    blocks = [(edges_y[i], edges_y[i + 1], edges_x[j], edges_x[j + 1]) for i in range(rows) for j in range(cols)]
    matchers = [_Matcher(ref_p[y0:y1, x0:x1], max_shift=max_shift, down=down, prepared=True) for (y0, y1, x0, x1) in blocks]
    rdy = np.zeros(T); rdx = np.zeros(T); nr = np.full(T, np.nan); nr_max = np.full(T, np.nan); bjump = np.full(T, np.nan)
    nb = len(blocks)
    prev = np.full((nb, 2), np.nan)
    for t in range(T):
        fr = np.asarray(movie[t], np.float32)
        fp = whole.prep(fr)
        rdy[t], rdx[t], _ = whole_shift(whole, fp)
        dev, cur = [], np.full((nb, 2), np.nan)
        for j, ((y0, y1, x0, x1), m) in enumerate(zip(blocks, matchers)):
            by, bx, bc = m.shift(fp[y0:y1, x0:x1])
            if bc >= 0.2:
                dev.append(float(np.hypot(by - rdy[t], bx - rdx[t]))); cur[j] = (by, bx)
        if len(dev) >= 3:
            nr[t] = float(np.median(dev)); nr_max[t] = float(np.percentile(dev, 90))
        ok = np.isfinite(cur[:, 0]) & np.isfinite(prev[:, 0])
        if ok.sum() >= 3:                      # how much the blocks moved since the previous frame (what makes edges flicker)
            bjump[t] = float(np.median(np.hypot(cur[ok, 0] - prev[ok, 0], cur[ok, 1] - prev[ok, 1])))
        prev = cur
    return {"res_dy": rdy, "res_dx": rdx, "nonrigid": nr, "nonrigid_max": nr_max, "block_jump": bjump}


def aqua2_style_register(movie, gap: int = 5, max_shift: float = 40.0, log=None) -> tuple[np.ndarray, dict]:
    """What AQuA2's registrateCorrect=2 does, in numpy (for tests and the fake
    MATLAB): crop `gap` pixels off every border (regMaskGap), shift every frame
    by a whole-pixel translation that best matches the mean of the first 10
    frames (cross-correlation), then crop the motion margin (max / min shift).
    Returns (registered uint16 movie (T, H', W'), {dy, dx, crop_offset (row,
    col of the crop in the raw frame)})."""
    T, H, W = (int(v) for v in movie.shape)
    def frame(t):
        return np.asarray(movie[t], np.float32)[gap:H - gap, gap:W - gap]
    ref = np.mean(np.stack([frame(t) for t in range(min(REF_FRAMES, T))]), axis=0)
    M = _Matcher(ref, max_shift=max_shift, down=1, sigma=1.0, hp=30.0)
    dy = np.zeros(T, int); dx = np.zeros(T, int)
    for t in range(T):
        a, b, _ = M.shift(frame(t))
        dy[t], dx[t] = int(round(a)), int(round(b))
    y0, y1 = int(max(dy.max(), 0)), int(min(dy.min(), 0)); x0, x1 = int(max(dx.max(), 0)), int(min(dx.min(), 0))
    Hc, Wc = H - 2 * gap, W - 2 * gap
    out = np.zeros((T, Hc - y0 + y1, Wc - x0 + x1), np.uint16)
    for t in range(T):
        f = frame(t)
        sh = np.zeros_like(f)
        # content sits (dy, dx) away from the reference: move it back
        ys, ye = max(0, dy[t]), min(Hc, Hc + dy[t]); xs, xe = max(0, dx[t]), min(Wc, Wc + dx[t])
        sh[ys - dy[t]:ye - dy[t], xs - dx[t]:xe - dx[t]] = f[ys:ye, xs:xe]
        out[t] = np.clip(np.round(sh[y0:Hc + y1, x0:Wc + x1]), 0, 65535).astype(np.uint16)
    if log:
        log(f"  registered {T} frames: shifts dy {dy.min()}..{dy.max()}, dx {dx.min()}..{dx.max()} -> {out.shape[1]}x{out.shape[2]}")
    return out, {"dy": dy, "dx": dx, "crop_offset": (gap + y0, gap + x0)}


def write_fake_aqua_mat(mat_path: Path, movie_reg: np.ndarray, opts: dict | None = None) -> Path:
    """An HDF5 file shaped like a MATLAB -v7.3 AQuA2 result: res/datOrg1 (H x W x 1 x T
    uint16, stored as h5py sees a MATLAB array: (T, 1, W, H)) and a few opts."""
    import h5py
    T, H, W = movie_reg.shape
    with h5py.File(mat_path, "w") as f:
        g = f.create_group("res")
        ds = g.create_dataset("datOrg1", data=np.ascontiguousarray(np.swapaxes(movie_reg, 1, 2))[:, None, :, :], dtype=np.uint16)
        ds.attrs["MATLAB_class"] = np.bytes_("uint16")
        o = g.create_group("opts")
        for k, v in {"maxValueDat1": float(movie_reg.max()), "minValueDat1": float(movie_reg.min()), "regMaskGap": 5.0,
                     "frameRate": 1.984, "spatialRes": 0.792, **(opts or {})}.items():
            d = o.create_dataset(k, data=np.array([[float(v)]]))
            d.attrs["MATLAB_class"] = np.bytes_("double")
    return Path(mat_path)


# ----------------------------------------------------------------------------- per run

def _running(beh, T: int, still_thresh: float = 0.5):
    """(running bool per frame, speed per frame) from the behaviour csv, padded/cut to T."""
    run = np.zeros(T, bool); speed = np.full(T, np.nan)
    if beh is None:
        return run, speed
    try:
        sf = beh["still_fraction"].to_numpy(float) if "still_fraction" in beh.columns else None
        v = beh["velocity_mean"].to_numpy(float) if "velocity_mean" in beh.columns else None
    except Exception:
        return run, speed
    n = min(T, len(beh))
    if sf is not None:
        run[:n] = sf[:n] < still_thresh
    elif v is not None:
        run[:n] = np.abs(v[:n]) >= 5
    if v is not None:
        speed[:n] = np.abs(v[:n])
    return run, speed


def _onsets_from_events(events) -> np.ndarray:
    """0-based movie frames of the event onsets (events csv: t0_frame 1-based)."""
    if events is None or len(events) == 0 or "t0_frame" not in getattr(events, "columns", []):
        return np.zeros(0, int)
    t0 = events["t0_frame"].to_numpy(float)
    return (np.round(t0[np.isfinite(t0)]).astype(int) - 1)


def coincidence(onsets: np.ndarray, flagged: np.ndarray, near: int = NEAR_FRAMES) -> dict:
    """How many onsets fall within +-near frames of a flagged frame, and how many
    would by chance (share of frames within +-near of a flagged one)."""
    T = flagged.size
    if T == 0:
        return {"n": 0, "near": 0, "expected": np.nan, "ratio": np.nan, "share_frames": np.nan}
    window = flagged.copy()
    for k in range(1, near + 1):
        window[k:] |= flagged[:-k]; window[:-k] |= flagged[k:]
    share = float(window.mean())
    on = onsets[(onsets >= 0) & (onsets < T)]
    n = int(on.size)
    nearc = int(window[on].sum()) if n else 0
    exp = share * n
    ratio = (nearc / exp) if exp > 0 else (np.nan if nearc == 0 else np.inf)
    return {"n": n, "near": nearc, "expected": float(exp), "ratio": float(ratio), "share_frames": share}


def motion_paths(tif: Path) -> tuple[Path, Path]:
    base = str(Path(tif).with_suffix(""))
    return Path(base + "_motion.csv"), Path(base + "_motion.json")


def _mtime(p: Path) -> float:
    try:
        return p.stat().st_mtime
    except OSError:
        return 0.0


def run_motion(tif: Path, beh=None, events=None, force: bool = False, registered: bool = True, log=print) -> dict | None:
    """Measure one run; write <base>_motion.csv / .json; return the summary dict
    (None when the movie cannot be read). Cached: nothing is recomputed when
    the files are newer than the tiff and the AQuA2 result."""
    import pandas as pd
    from .aqua_io import aqua_files, open_registered
    tif = Path(tif)
    csv_p, json_p = motion_paths(tif)
    f = aqua_files(tif)
    newest_src = max(_mtime(tif), _mtime(f["res"]), _mtime(Path(str(tif.with_suffix("")) + "_AQuA2_events.csv")))
    if not force and json_p.exists() and csv_p.exists() and _mtime(json_p) >= newest_src:
        try:
            summ = json.loads(json_p.read_text(encoding="utf-8"))
            if summ.get("version") == 1:
                summ["verdict"] = verdict(summ)          # wording follows the current rules
                return summ
        except Exception:
            pass
    if not tif.exists() or tif.stat().st_size == 0:
        return None
    import tifffile as tfl
    t_start = time.time()
    try:
        mov = tfl.memmap(str(tif), mode="r")
    except Exception:
        mov = tfl.imread(str(tif))
    if mov.ndim == 2:
        mov = mov[None]
    T = int(mov.shape[0])
    spf = 2.0
    if f["json"].exists():
        try:
            spf = float(json.loads(f["json"].read_text())["seconds_per_frame"])
        except Exception:
            pass
    if beh is None and f["behavior"].exists():
        try:
            beh = pd.read_csv(f["behavior"])
        except Exception:
            beh = None
    if events is None:
        ev_csv = Path(str(tif.with_suffix("")) + "_AQuA2_events.csv")
        if ev_csv.exists():
            try:
                events = pd.read_csv(ev_csv)
            except Exception:
                events = None
    log(f"  motion: {tif.name} ({T} frames) …")
    sh = frame_shifts(mov)
    amp = np.hypot(sh["dy"], sh["dx"])
    jump = np.concatenate([[0.0], np.hypot(np.diff(sh["dy"]), np.diff(sh["dx"]))])
    run, speed = _running(beh, T)
    df = pd.DataFrame({"frame": np.arange(T), "time_s": (np.arange(T) + 0.5) * spf, "dy": sh["dy"], "dx": sh["dx"],
                       "shift": amp, "jump": jump, "corr": sh["corr"], "sharp": sh["sharp"],
                       "running": run.astype(int), "speed": speed})
    reg_info = {"available": False}
    if registered:
        reg = open_registered(tif, log=log)
        if reg is not None:
            try:
                Tr = min(int(reg.shape[0]), T)
                br = block_residual(reg)
                for k in ("res_dy", "res_dx", "nonrigid", "nonrigid_max", "block_jump"):
                    col = np.full(T, np.nan); col[:min(Tr, br[k].size)] = br[k][:min(Tr, br[k].size)]
                    df[k] = col
                reg_info = {"available": True, "frames": int(reg.shape[0]), "shape": [int(reg.shape[1]), int(reg.shape[2])]}
            except Exception as e:
                log(f"    registered movie: residual not measured ({e})")
            finally:
                reg.close()
    # summary
    still = ~run
    fin = np.isfinite
    summ = {"version": 1, "run": tif.stem, "frames": T, "seconds_per_frame": spf, "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "shift_median": float(np.median(amp)), "shift_p95": float(np.percentile(amp, 95)), "shift_max": float(amp.max()),
            "dy_range": [float(sh["dy"].min()), float(sh["dy"].max())], "dx_range": [float(sh["dx"].min()), float(sh["dx"].max())],
            "corr_median": float(np.median(sh["corr"])), "corr_min": float(sh["corr"].min()),
            "jump_p95": float(np.percentile(jump, 95)), "jump_max": float(jump.max()),
            "jumps_ge2_per_min": float((jump >= JUMP_PX).sum() / (T * spf / 60.0)),
            "running_pct": float(100.0 * run.mean()),
            "shift_sd_still": float(np.hypot(sh["dy"][still].std(), sh["dx"][still].std())) if still.sum() >= MIN_STATE_FRAMES else np.nan,
            "shift_sd_run": float(np.hypot(sh["dy"][run].std(), sh["dx"][run].std())) if run.sum() >= MIN_STATE_FRAMES else np.nan,
            "jumps_ge2_per_min_still": float((jump[still] >= JUMP_PX).sum() / (still.sum() * spf / 60.0)) if still.sum() >= MIN_STATE_FRAMES else np.nan,
            "jumps_ge2_per_min_run": float((jump[run] >= JUMP_PX).sum() / (run.sum() * spf / 60.0)) if run.sum() >= MIN_STATE_FRAMES else np.nan,
            "corr_shift_speed": float(np.corrcoef(amp[fin(speed)], speed[fin(speed)])[0, 1]) if fin(speed).sum() > 10 and np.nanstd(speed) > 0 else np.nan,
            "sharp_still": float(np.nanmedian(sh["sharp"][still])) if still.sum() >= MIN_STATE_FRAMES else np.nan,
            "sharp_run": float(np.nanmedian(sh["sharp"][run])) if run.sum() >= MIN_STATE_FRAMES else np.nan,
            "registered": reg_info}
    s_st, s_ru = summ["sharp_still"], summ["sharp_run"]
    summ["sharp_drop_run_pct"] = float(100.0 * (1 - s_ru / s_st)) if (fin(s_st) and fin(s_ru) and s_st > 0) else np.nan
    if "nonrigid" in df.columns:
        nr = df["nonrigid"].to_numpy(float); rres = np.hypot(df["res_dy"].to_numpy(float), df["res_dx"].to_numpy(float))
        summ.update({"residual_rigid_p95": float(np.nanpercentile(rres, 95)), "residual_rigid_max": float(np.nanmax(rres)),
                     "nonrigid_median": float(np.nanmedian(nr)), "nonrigid_p95": float(np.nanpercentile(nr, 95)),
                     "nonrigid_median_still": float(np.nanmedian(nr[still])) if still.sum() >= MIN_STATE_FRAMES else np.nan,
                     "nonrigid_median_run": float(np.nanmedian(nr[run])) if run.sum() >= MIN_STATE_FRAMES else np.nan,
                     "frames_nonrigid_ge1p5_pct": float(100.0 * np.nanmean(nr >= NONRIGID_PX))})
    # events vs motion
    onsets = _onsets_from_events(events)
    medj = float(np.median(jump)); madj = float(np.median(np.abs(jump - medj))) * 1.4826
    thr_jump = max(JUMP_PX, medj + 3.0 * madj)                   # a jump unusual for this run (at least JUMP_PX)
    summ["jump_threshold"] = thr_jump
    flag_jump = jump >= thr_jump
    c_jump = coincidence(onsets, flag_jump)
    summ["events"] = {"n_onsets": c_jump["n"], "near_jump": c_jump["near"], "expected_near_jump": c_jump["expected"],
                      "ratio_jump": c_jump["ratio"], "share_frames_near_jump": c_jump["share_frames"]}
    if "block_jump" in df.columns:
        bj = np.nan_to_num(df["block_jump"].to_numpy(float))
        med = float(np.median(bj)); mad = float(np.median(np.abs(bj - med))) * 1.4826
        thr_nr = max(NONRIGID_JUMP_PX, med + 3.0 * mad)          # a deformation step unusual for this run
        summ["nonrigid_jump_threshold"] = thr_nr
        flag_nr = bj >= thr_nr
        c_nr = coincidence(onsets, flag_nr)
        summ["events"].update({"near_nonrigid": c_nr["near"], "expected_near_nonrigid": c_nr["expected"], "ratio_nonrigid": c_nr["ratio"]})
        summ["nonrigid_jumps_per_min"] = float(flag_nr.sum() / (T * spf / 60.0))
        flag_any = flag_jump | flag_nr
    else:
        flag_any = flag_jump
    c_any = coincidence(onsets, flag_any)
    summ["events"].update({"near_motion": c_any["near"], "expected_near_motion": c_any["expected"], "ratio_motion": c_any["ratio"],
                           "share_frames_near_motion": c_any["share_frames"]})
    df["motion_frame"] = flag_any.astype(int)
    summ["verdict"] = verdict(summ)
    summ["seconds"] = round(time.time() - t_start, 1)
    try:
        df.to_csv(csv_p, index=False, float_format="%.4g")
        json_p.write_text(json.dumps(_jsonable(summ), indent=1), encoding="utf-8")
    except OSError as e:
        log(f"    motion files not written ({e})")
    log(f"    shift median {summ['shift_median']:.1f} px (max {summ['shift_max']:.1f}), jumps >= {JUMP_PX:g} px: {summ['jumps_ge2_per_min']:.1f}/min"
        + (f", deformation after AQuA2 median {summ['nonrigid_median']:.1f} px, {summ.get('nonrigid_jumps_per_min', 0):.1f} deformation jumps/min" if "nonrigid_median" in summ else "")
        + f"; events near motion {summ['events']['near_motion']} of {summ['events']['n_onsets']} (expected {summ['events']['expected_near_motion']:.1f})"
        + f" — {summ['verdict'][0]} ({summ['seconds']} s)")
    return summ


def _fin(v) -> bool:
    try:
        return v is not None and bool(np.isfinite(float(v)))
    except (TypeError, ValueError):
        return False


def verdict(summ: dict) -> list:
    """[level, explanation]: 'ok', 'check' or 'warn' about the run's motion."""
    e = summ.get("events", {})
    n, near, ratio = e.get("n_onsets", 0) or 0, e.get("near_motion", 0) or 0, e.get("ratio_motion", np.nan)
    level, reasons, notes = "ok", [], []
    if n >= 5 and near >= 3 and _fin(ratio) and ratio >= 2.0:
        level = "warn"; reasons.append(f"{near} of {n} event onsets fall on moving frames ({ratio:.1f}x the chance level) — look at those events")
    elif n >= 5 and _fin(ratio) and ratio >= 1.5 and near >= 2:
        level = "check"; reasons.append(f"{near} of {n} event onsets on moving frames ({ratio:.1f}x chance)")
    elif n >= 5:
        notes.append(f"the event onsets do not follow the movement ({near} of {n} on moving frames, {(e.get('expected_near_motion') or 0):.1f} expected)")
    elif n:
        notes.append(f"too few events ({n}) to test whether they follow the movement")
    else:
        notes.append("no events to test against the movement")

    def raise_to(lv):
        nonlocal level
        if lv == "warn" or level == "ok":
            level = lv
    if (summ.get("shift_max") or 0) >= 15:
        raise_to("check"); reasons.append(f"the field moved up to {summ['shift_max']:.0f} px from the reference (AQuA2 cropped that margin away)")
    if (summ.get("jumps_ge2_per_min") or 0) >= 10:
        raise_to("check"); reasons.append(f"{summ['jumps_ge2_per_min']:.0f} jumps of >= {JUMP_PX:g} px per minute between consecutive frames")
    nj = summ.get("nonrigid_jumps_per_min", np.nan)
    if _fin(nj) and nj >= 10:
        raise_to("check"); reasons.append(f"{nj:.0f} frames per minute where the field deformed by >= {summ.get('nonrigid_jump_threshold', NONRIGID_JUMP_PX):.1f} px since the previous frame")
    sdr = summ.get("sharp_drop_run_pct", np.nan)
    if _fin(sdr) and sdr >= 3 and (summ.get("running_pct") or 0) >= 5:
        raise_to("check"); reasons.append(f"frames are {sdr:.0f} % less sharp while running (motion within the 2-s bins)")
    nr = summ.get("nonrigid_median", np.nan)
    if _fin(nr) and nr >= NONRIGID_PX:
        notes.append(f"after AQuA2's registration the field is still deformed by {nr:.1f} px (median) relative to its reference frames — "
                     "a slow warp, not a flicker; it blurs the CFU footprints rather than inventing events")
    text = "; ".join(reasons + notes) if reasons else "; ".join(notes) or "movement within the usual range"
    return [level, text]


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    return o


def read_motion(tif: Path) -> tuple:
    """(summary dict or None, per-frame DataFrame or None) from the cached files."""
    import pandas as pd
    csv_p, json_p = motion_paths(tif)
    summ = df = None
    if json_p.exists():
        try:
            summ = json.loads(json_p.read_text(encoding="utf-8"))
        except Exception:
            summ = None
    if csv_p.exists():
        try:
            df = pd.read_csv(csv_p)
        except Exception:
            df = None
    return summ, df


# ----------------------------------------------------------------------------- figure

def motion_figure(df, summ: dict | None, events=None, title: str = "") -> str | None:
    """Shift traces (dy, dx), the jump between frames, the non-rigid residual
    (when measured) and the running speed on one clock, with the event onsets
    marked; returns a PNG data URI (afferent2p.report._fig_uri)."""
    if df is None or len(df) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .report import _fig_uri, _style, INK, INK2, SERIES, GRID
    t = df["time_s"].to_numpy(float) / 60.0
    have_nr = "nonrigid" in df.columns and np.isfinite(df["nonrigid"].to_numpy(float)).any()
    n_ax = 3 + int(have_nr)
    fig, axes = plt.subplots(n_ax, 1, figsize=(11, 1.25 * n_ax + 0.9), sharex=True,
                             gridspec_kw=dict(hspace=0.12, left=0.07, right=0.99, top=0.9, bottom=0.11))
    axes = np.atleast_1d(axes)
    onsets = _onsets_from_events(events)
    on_t = t[onsets[(onsets >= 0) & (onsets < len(t))]] if onsets.size else np.zeros(0)
    ax = axes[0]; _style(ax)
    ax.plot(t, df["dy"], color=SERIES[0], lw=0.9, label="dy (rows)")
    ax.plot(t, df["dx"], color=SERIES[1], lw=0.9, label="dx (columns)")
    ax.axhline(0, color=GRID, lw=0.7)
    ax.set_ylabel("shift of the raw\nframe (px)", fontsize=8, color=INK2)
    ax.legend(fontsize=7, loc="upper right", frameon=False, ncol=2)
    ax = axes[1]; _style(ax)
    ax.fill_between(t, 0, df["jump"], color=INK, alpha=0.6, lw=0, step="mid")
    ax.axhline((summ or {}).get("jump_threshold", JUMP_PX), color="#d64545", lw=0.7, ls="--")
    ax.set_ylabel("jump between\nframes (px)", fontsize=8, color=INK2)
    k = 2
    if have_nr:
        ax = axes[2]; _style(ax)
        if "block_jump" in df.columns:
            ax.fill_between(t, 0, df["block_jump"].fillna(0), color="#7a5195", alpha=0.6, lw=0, step="mid", label="blocks moved since the previous frame")
            ax.axhline((summ or {}).get("nonrigid_jump_threshold", NONRIGID_JUMP_PX), color="#d64545", lw=0.7, ls="--")
        ax.plot(t, df["nonrigid"], color="#ffa600", lw=0.8, label="deformation vs the reference")
        ax.set_ylabel("after AQuA2:\nnon-rigid (px)", fontsize=8, color=INK2)
        ax.legend(fontsize=7, loc="upper right", frameon=False)
        k = 3
    ax = axes[k]; _style(ax)
    sp = df["speed"].to_numpy(float)
    if np.isfinite(sp).any():
        ax.fill_between(t, 0, np.nan_to_num(sp), color=SERIES[0], alpha=0.5, lw=0, step="mid")
        ax.set_ylabel("running\nspeed", fontsize=8, color=INK2)
    else:
        ax.text(0.005, 0.5, "no behaviour csv", transform=ax.transAxes, fontsize=7, color=INK2, va="center")
    ax.set_xlabel("minutes", fontsize=8, color=INK2)
    for a in axes:
        for x in on_t:
            a.axvline(x, color="#2a9d8f", lw=0.6, alpha=0.55, zorder=0)
    mf = df["motion_frame"].to_numpy(int) if "motion_frame" in df.columns else None
    if mf is not None and mf.any():
        for a in axes:
            a.fill_between(t, 0, 1, where=mf > 0, transform=a.get_xaxis_transform(), color="#f3c6c6", alpha=0.45, lw=0, step="mid", zorder=0)
    head = title or (summ or {}).get("run", "")
    if summ:
        v = summ.get("verdict", ["", ""])
        head += f" — motion {v[0]}: {v[1]}"
    fig.suptitle(head, fontsize=9, x=0.01, ha="left", color=INK)
    axes[0].set_xlim(t.min() if t.size else 0, t.max() if t.size else 1)
    return _fig_uri(fig, dpi=100)
