"""Choosing AQuA2 detection settings by measurement instead of by eye.

There is no ground truth for afferent events, but for any setting three
things can be measured without anyone watching:

  false alarms   the setting run on a NULL copy of the movie — every pixel's
                 time series circularly shifted by its own random offset, so
                 each pixel keeps its brightness, noise and slow drift while
                 nothing is left that is spatially coherent in time. Every
                 event AQuA2 finds there is a false positive.
  recall         the setting run on an INJECTED copy — the real movie plus
                 synthetic transients of known amplitude, size and duration
                 placed on the fibres. The fraction recovered, per amplitude,
                 is what the setting can see (corrected for the chance that a
                 real event happens to coincide).
  specificity    the fraction of the real movie's events that sit on the
                 fibres of the mean image, and how the count changes with the
                 threshold (a plateau is a robust setting, a cliff is not).

aqua_tune.py writes the three movies per chosen run, a grid of settings as
AQuA2 presets (the FOV's own Opt file with the swept values replaced), one
one-row batch file per job and a MATLAB script that runs the jobs one by one,
exports each result's events and mean image and deletes the 1.2 GB .mat again.
Afterwards --report reads the exports, scores every setting, recommends one
per movie (best recall among the settings whose false-alarm share and
on-fibre fraction pass), and --apply writes it as an Opt file the batch uses.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

DEFAULT_GRID = {"thrARScl": [1.5, 2.0, 2.5], "minSize": [20, 50, 100]}        # v37: minSize 20 in, thrARScl 3 and 4 out
AMPLITUDES = (0.3, 0.6, 1.0)          # nominal dF/F of the injected transients (of each pixel's own baseline). AQuA2 grows
                                      # the footprint 2-3x and reports the footprint-average dF/F, so a nominal 0.3 / 0.6 /
                                      # 1.0 is reported at ~0.15 / 0.3 / 0.5 — below, at and above the real events, which
                                      # are reported at 0.2-0.45 dF/F (2608-1 thr 2; 2609-1 FOV 2 by hand: p10-p90
                                      # 0.18-0.41). v37: the earlier default 0.1 / 0.2 / 0.4 probed BELOW the real signal
                                      # (reported ~0.05-0.2): every setting looked blind (recall <= 0.3), and the choice
                                      # among blind settings was meaningless (the Oct 21 sweeps of all three animals)
SIZES_PX = (20, 60, 200)              # footprint sizes in pixels on the fibre mask (detected 2-3x larger): the real events
                                      # of 2608-1 (thick bright fibres) are reported at ~150 (p10) to ~1200 px (p90), those
                                      # of 2609-1 FOV 2 (a fine network) at 29 to 194 px, median 42 — the small class is
                                      # what tells whether minSize lets such events through. v37: 20 in, 400 out
OLD_AMPLITUDES = ((0.1, 0.2, 0.4),)   # class lists of earlier versions, replaced by the defaults when a plan is re-made
OLD_SIZES = ((80, 200, 400),)
MASK_VERSION = 2                      # v37: the fibre network is kept whole (see _elongated_only); a saved mask of an
                                      # older version is recomputed by the plan
# which way is "more permissive" for the swept variables (lower = more events, unless listed in the second set)
PERMISSIVE_LOWER = {"thrARScl", "minSize", "minDur", "sigThr", "sourceSzRatio", "seedSzRatio", "circularityThr", "maxSize"}
PERMISSIVE_HIGHER = {"sourceSensitivity", "maxDelay", "gapExt"}
KINDS = ("short", "long")             # ~4 and ~12 frames above half amplitude (8 s / 24 s at 2 s per frame;
                                      # the real events' interquartile range is 4-14 frames)


# ----------------------------------------------------------------------------- movies

def _open(tif: Path):
    import tifffile as tf
    try:
        arr = tf.memmap(str(tif), mode="r")
    except Exception:
        arr = tf.imread(str(tif))
    if arr.ndim == 2:
        arr = arr[None]
    return arr


def mean_image(tif: Path, n: int = 150) -> np.ndarray:
    """Mean image of a movie: unpack_bin's <stem>_meanproj.tif (the mean of every
    frame) when it exists, else the mean of up to n evenly spaced frames."""
    import tifffile as tf
    tif = Path(tif)
    for cand in (tif.with_name(tif.stem + "_meanproj.tif"),):
        if cand.exists():
            try:
                m = np.asarray(tf.imread(str(cand)), np.float32)
                if m.ndim == 2:
                    return m
            except Exception:
                pass
    arr = _open(tif)
    idx = np.unique(np.linspace(0, arr.shape[0] - 1, min(n, arr.shape[0])).astype(int))
    return np.mean(np.stack([np.asarray(arr[int(i)], np.float32) for i in idx]), axis=0)


def manual_mask(tif: Path, mouse_dir: Path | None, fov=None, shape=None):
    """A hand-drawn fibre mask, if the user made one: <stem>_fibre_mask.tif next
    to the movie, or <mouse>\\FOV<k>_fibre_mask.tif (any non-zero pixel = fibre;
    same size as the movie frame, e.g. drawn in ImageJ on <stem>_meanproj.tif).
    Returns (mask, description) or (None, '')."""
    import tifffile as tf
    tif = Path(tif)
    cands = [tif.with_name(tif.stem + "_fibre_mask.tif")]
    if mouse_dir is not None and fov is not None:
        cands.append(Path(mouse_dir) / f"FOV{fov}_fibre_mask.tif")
    for c in cands:
        if not c.exists():
            continue
        try:
            m = np.asarray(tf.imread(str(c)))
            if m.ndim == 3:
                m = m.max(axis=0)
            m = m > 0
            if shape is not None and m.shape != tuple(shape):
                return None, f"ignored {c.name}: {m.shape[0]}x{m.shape[1]} px but the movie frame is {shape[0]}x{shape[1]}"
            if m.sum() == 0:
                return None, f"ignored {c.name}: empty"
            return m, f"hand-drawn mask {c.name}"
        except Exception as e:
            return None, f"ignored {c.name}: {e}"
    return None, ""


def _elongated_only(mask: np.ndarray, um_per_px: float, min_len_um: float = 10.0, min_ratio: float = 2.5):
    """Drop the round, compact pieces of a mask (autofluorescent specks, cell
    bodies) and keep the fibres. A piece is kept when its skeleton is at least
    min_len_um long (a fibre segment or a whole branched network — v37: before,
    a dense network of crossing fibres formed ONE big connected piece that is
    about as wide as it is long, failed the elongation test below and was
    dropped, which left only isolated specks for 2609-1 FOV 2) or when it is
    min_len_um long and min_ratio times longer than wide. Returns (mask,
    number dropped)."""
    from scipy.ndimage import label
    from skimage.measure import regionprops
    from skimage.morphology import skeletonize
    lab, n = label(mask, structure=np.ones((3, 3)))
    if n == 0:
        return mask.astype(bool), 0
    keep = np.zeros(n + 1, bool)
    min_len = min_len_um / um_per_px
    skel = skeletonize(mask.astype(bool))
    skel_len = np.bincount(lab[skel], minlength=n + 1)
    for rp in regionprops(lab):
        major = rp.axis_major_length; minor = max(rp.axis_minor_length, 1.0)
        if skel_len[rp.label] >= min_len or (major >= min_len and major / minor >= min_ratio):
            keep[rp.label] = True
    out = keep[lab]
    return out, int(n - keep[1:].sum())


def long_pieces(mask: np.ndarray, min_skel_px: int = 25) -> np.ndarray:
    """The pieces of the mask whose skeleton is at least min_skel_px long — the
    fibres proper, without the specks and short bits (used to place the
    injection sites, v37)."""
    from scipy.ndimage import label
    from skimage.morphology import skeletonize
    m = np.asarray(mask, bool)
    lab, n = label(m, structure=np.ones((3, 3)))
    if n == 0:
        return m
    sk = skeletonize(m)
    skel_len = np.bincount(lab[sk], minlength=n + 1)
    keep = skel_len >= min_skel_px
    keep[0] = False
    return keep[lab]


def injection_mask(mean_img: np.ndarray, um_per_px: float, min_px: int = 5000, min_px_relaxed: int = 1500, log=print):
    """Where the fibres are, for placing injected transients and for the on-fibre
    fraction: the ridge-filter fibre mask of spatial.py; if that is small (dim or
    noisy movie) the same with a lower ridge threshold; if even that is small,
    the brightest 5 % of the smoothed, flattened mean image. Returns (mask, description)."""
    from afferent2p.spatial import fibre_mask
    fm = fibre_mask(mean_img, um_per_px)
    mask, n_blob = _elongated_only(fm["mask"], um_per_px)
    how = f"fibre mask, {fm['n_segments']} segments" + (f", {n_blob} round specks dropped" if n_blob else "")
    b = 6                                                   # the frame's own border is never a fibre
    for m_ in (mask,):
        m_[:b, :] = False; m_[-b:, :] = False; m_[:, :b] = False; m_[:, -b:] = False
    if mask.sum() < min_px:
        fm2 = fibre_mask(mean_img, um_per_px, z_thr=2.0)
        m2, n_blob2 = _elongated_only(fm2["mask"], um_per_px)
        if m2.sum() > mask.sum():
            mask, how = m2, f"relaxed fibre mask (ridge z > 2; {fm2['n_segments']} segments" + (f", {n_blob2} round specks dropped)" if n_blob2 else ")")
    if mask.sum() < min_px_relaxed:
        from scipy.ndimage import gaussian_filter
        from skimage.morphology import remove_small_objects
        im = np.asarray(mean_img, np.float32)
        flat = gaussian_filter(im, 1.5) - gaussian_filter(im, 25)      # pixel noise smoothed, vignetting removed
        bright = flat > np.percentile(flat, 95)
        try:
            bright = remove_small_objects(bright, max_size=19)
        except TypeError:
            bright = remove_small_objects(bright, min_size=20)
        if bright.sum() > mask.sum():
            mask, how = bright, "brightest 5 % of the flattened mean image (fibres too faint for the ridge filter)"
    return mask.astype(bool), how


def _new_movie(dst: Path, shape, sec_per_frame: float | None):
    """An empty movie file to fill frame by frame, written the way unpack_bin.py
    writes the binned movies (ImageJ-style TIFF below 4 GB, BigTIFF above).
    Returns (memmap, partial path); rename the partial path to dst when done."""
    import tifffile as tf
    T, H, W = shape
    gb = T * H * W * 2 / 1e9
    md = {"axes": "TYX"}
    if sec_per_frame:
        md["finterval"] = float(sec_per_frame)
    kw = dict(imagej=True, metadata=md) if gb < 3.9 else dict(bigtiff=True)
    part = dst.with_name(dst.name + ".partial")
    return tf.memmap(str(part), shape=(T, H, W), dtype=np.uint16, **kw), part


def make_null_movie(src: Path, dst: Path, seed: int = 0, chunk_rows: int = 32, sec_per_frame: float | None = None,
                    log=print) -> Path:
    """Per-pixel circular time shift (offsets drawn from T/4..3T/4): same pixels,
    same noise and drift per pixel, no spatially coherent transients."""
    arr = _open(src)
    T, H, W = arr.shape
    rng = np.random.default_rng(seed)
    out, part = _new_movie(dst, (T, H, W), sec_per_frame)
    tt = np.arange(T)[:, None, None]
    for r0 in range(0, H, chunk_rows):
        r1 = min(H, r0 + chunk_rows)
        blk = np.asarray(arr[:, r0:r1, :])
        shift = rng.integers(T // 4, 3 * T // 4, size=(1, r1 - r0, W))
        idx = (tt + shift) % T
        out[:, r0:r1, :] = np.take_along_axis(blk, idx, axis=0)
    out.flush(); del out
    part.replace(dst)
    log(f"  null movie: {dst.name}")
    return dst


def _kernel(kind: str) -> np.ndarray:
    """Transient shape in frames (binned, ~2 s/frame): rise, plateau, exponential decay.
    short: 4 frames above half amplitude; long: 12 frames."""
    if kind == "short":
        rise, plateau, tau, n = 1, 1, 2.0, 10
    else:
        rise, plateau, tau, n = 2, 6, 6.0, 30
    k = np.zeros(n)
    for i in range(n):
        if i < rise:
            k[i] = (i + 1) / (rise + 1)
        elif i < rise + plateau:
            k[i] = 1.0
        else:
            k[i] = np.exp(-(i - rise - plateau) / tau)
    return k


def make_injected_movie(src: Path, dst: Path, truth_csv: Path, fmask: np.ndarray, n_sites: int = 72,
                        amplitudes=AMPLITUDES, sizes=SIZES_PX, kinds=KINDS, seed: int = 1,
                        min_sep_px: int = 25, sec_per_frame: float | None = None, log=print) -> int:
    """Real movie + synthetic transients on the fibre mask. Each site gets one
    event; amplitude (dF/F of the pixel's own baseline), footprint size and
    kind cycle through the given classes so every class gets the same number.
    Writes the movie and a truth table (1-based frames; 0-based raw-frame
    pixel coordinates, i.e. numpy indices of the input movie)."""
    import pandas as pd
    arr = _open(src)
    T, H, W = arr.shape
    rng = np.random.default_rng(seed)
    # per-pixel baseline (temporal median of a sample of frames)
    idx = np.unique(np.linspace(0, T - 1, min(60, T)).astype(int))
    f0 = np.median(np.stack([np.asarray(arr[int(i)], np.float32) for i in idx]), axis=0)
    # sites: fibre pixels at least min_sep apart, away from the border — on the long pieces of the mask (the
    # fibres proper) when those offer enough room, so that specks hanging on the mask get no transients (v37)
    fm_sites = long_pieces(fmask)
    if fm_sites.sum() < 2000:
        fm_sites = np.asarray(fmask, bool)
    cand = np.argwhere(fm_sites)
    cand = cand[(cand[:, 0] > 20) & (cand[:, 0] < H - 20) & (cand[:, 1] > 20) & (cand[:, 1] < W - 20)]
    if cand.size == 0:
        raise ValueError("no fibre pixels to inject on (the fibre mask is empty)")
    rng.shuffle(cand)
    sites = []
    for p in cand:
        if all(np.hypot(*(p - q)) >= min_sep_px for q in sites):
            sites.append(p)
        if len(sites) >= n_sites:
            break
    if len(sites) < n_sites:
        log(f"  note: only {len(sites)} injection sites fit on the fibre mask at {min_sep_px} px separation")
    # footprints are built locally around each site
    dst_arr, part = _new_movie(dst, (T, H, W), sec_per_frame)
    rows = []
    classes = [(a, s, k) for a in amplitudes for s in sizes for k in kinds]
    rng.shuffle(classes)
    adds = []                      # (t0, kernel, pixel index array, amplitude)
    for j, (y, x) in enumerate(sites):
        amp, size, kind = classes[j % len(classes)]
        # grow a disc on the fibre mask until it holds ~size pixels (cap radius 40)
        r = 2
        while True:
            y0, y1 = max(0, y - r), min(H, y + r + 1); x0, x1 = max(0, x - r), min(W, x + r + 1)
            yy, xx = np.mgrid[y0:y1, x0:x1]
            disc = ((yy - y) ** 2 + (xx - x) ** 2) <= r * r
            foot = disc & fmask[y0:y1, x0:x1]
            if foot.sum() >= size or r >= 40:
                break
            r += 1
        fy, fx = np.nonzero(foot); fy = fy + y0; fx = fx + x0
        k = _kernel(kind)
        t0 = int(rng.integers(10, T - k.size - 10))
        adds.append((t0, k, fy, fx, amp))
        rows.append({"event": j + 1, "site_y": int(y), "site_x": int(x), "y_px": float(fy.mean()), "x_px": float(fx.mean()),
                     "n_pixels": int(fy.size), "radius_px": r, "amplitude": amp, "size_class": size, "kind": kind,
                     "t0_frame": t0 + 1, "t_peak_frame": t0 + 1 + int(np.argmax(k)),
                     "duration_frames": int((k >= 0.5).sum()), "n_frames_kernel": int(k.size)})
    # write frame by frame (memory-light)
    by_frame: dict[int, list] = {}
    for t0, k, fy, fx, amp in adds:
        for i in range(k.size):
            by_frame.setdefault(t0 + i, []).append((fy, fx, amp * k[i]))
    for t in range(T):
        fr = np.asarray(arr[t], np.float32)
        for fy, fx, g in by_frame.get(t, []):
            fr[fy, fx] += g * f0[fy, fx]
        dst_arr[t] = np.clip(np.round(fr), 0, 65535).astype(np.uint16)
    dst_arr.flush(); del dst_arr
    part.replace(dst)
    pd.DataFrame(rows).to_csv(truth_csv, index=False)
    log(f"  injected movie: {dst.name} — {len(rows)} transients on {len(sites)} sites "
        f"({', '.join(f'{a:g}' for a in amplitudes)} dF/F; {', '.join(str(s) for s in sizes)} px; {', '.join(kinds)})")
    return len(rows)


# ----------------------------------------------------------------------------- grid

def parse_grid(items: list[str] | None) -> dict:
    """['thrARScl=1.5,2,3', 'minSize=50,100'] -> {var: [values]}."""
    if not items:
        return {k: list(v) for k, v in DEFAULT_GRID.items()}
    grid = {}
    for it in items:
        if "=" not in it:
            raise ValueError(f"grid entries look like thrARScl=1.5,2,3 (got {it!r})")
        k, v = it.split("=", 1)
        grid[k.strip()] = [float(x) for x in v.replace(";", ",").split(",") if x.strip()]
    return grid


def _fmt(v) -> str:
    return str(int(v)) if float(v).is_integer() else f"{float(v):g}"


def settings_from_grid(base: dict, grid: dict) -> list[dict]:
    """Setting 0 = the FOV's current values; then every combination of the grid
    (the current combination is not repeated)."""
    import itertools
    keys = list(grid)
    out = [{"setting": 0, "label": "current", **{k: _fmt(base.get(k, "nan")) for k in keys}}]
    for combo in itertools.product(*[grid[k] for k in keys]):
        vals = {k: _fmt(v) for k, v in zip(keys, combo)}
        if all(vals[k] == out[0][k] for k in keys):
            out[0]["label"] = "current (in grid)"
            continue
        out.append({"setting": len(out), "label": " ".join(f"{k}={vals[k]}" for k in keys), **vals})
    return out


# ----------------------------------------------------------------------------- MATLAB

TUNE_SCRIPT = r"""%% run_aqua_tune.m — generated by aqua_tune.py
% Runs every job of tune_jobs.csv (one movie x one setting) through AQuA2's
% aqua_cmd_batch, exports the events and the mean image of each result and
% deletes the big _AQuA2.mat again (tune_one_job.m does one job). Jobs whose
% events file exists are skipped, so the script can be stopped and restarted.
% nWorkers > 1 runs that many jobs at once in a parallel pool (needs the
% Parallel Computing Toolbox; each job takes several GB of memory, 0 = choose
% from the machine). Without the toolbox, run_aqua_tune_parallel.bat starts
% several MATLAB processes, each with its own share of the jobs.
nWorkers = {n_workers};
if ~exist('tuneSplit', 'var'), tuneSplit = [1 1]; end     % set by run_aqua_tune_<k>of<n>.m
aquaDir = '{aqua_dir}';
if ~isfolder(aquaDir)                      % the folder set in the GUI is not visible here: find AQuA2 on the MATLAB path
    w = which('aqua_cmd_batch');
    if ~isempty(w), aquaDir = fileparts(w); else, error('AQuA2 not found: set the AQuA2 folder in the GUI (top row) or addpath it in MATLAB'); end
end
tuneDir = '{tune_dir}';
cd(aquaDir);
startup;
addpath(tuneDir);
jobs = readtable('{jobs_csv}', 'TextType', 'string', 'Delimiter', ',');
paramsCsv = '{params_csv}';
logFile = '{log_file}';
if tuneSplit(2) > 1, logFile = strrep(logFile, '.txt', sprintf('_%dof%d.txt', tuneSplit(1), tuneSplit(2))); end
try
    diary(logFile);      % aqua_cmd_batch clears the command window; the log keeps everything
catch                    % the old log cannot be appended to (written by another account, e.g. an
    logFile = strrep(logFile, '.txt', ['_' datestr(now, 'yyyymmdd_HHMM') '.txt']);   % Administrator window): a fresh file
    diary(logFile);
end
nJobs = height(jobs);
jobList = tuneSplit(1):tuneSplit(2):nJobs;
nDone = 0;
for j = jobList
    if exist(char(jobs.events_csv(j)), 'file'), nDone = nDone + 1; end
end
havePCT = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
if tuneSplit(2) > 1
    nWorkers = 1;                                   % one process per share
elseif nWorkers == 0
    if havePCT
        nWorkers = max(1, min(6, floor(feature('numcores') / 2)));
        try
            m = memory; nWorkers = max(1, min(nWorkers, floor(m.MemAvailableAllArrays / 8e9)));
        catch
        end
    else
        nWorkers = 1;
    end
elseif nWorkers > 1 && ~havePCT
    fprintf('no Parallel Computing Toolbox here - one job at a time (run_aqua_tune_parallel.bat starts several MATLAB processes instead)\n');
    nWorkers = 1;
end
fprintf('\n%s  AQuA2 tuning: %d job(s) in this share, %d done already (about %.1f h left at 1.4 min each, %d at a time)\n', ...
    datestr(now), numel(jobList), nDone, (numel(jobList) - nDone) * 1.4 / 60 / max(1, nWorkers), nWorkers);
tAll = tic;
if nWorkers > 1
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= nWorkers
        delete(gcp('nocreate'));
        pool = parpool(nWorkers);
    end
    addAttachedFiles(pool, {{fullfile(tuneDir, 'tune_one_job.m'), fullfile(tuneDir, 'tune_worker_init.m')}});
    f = parfevalOnAll(pool, @addpath, 0, tuneDir); wait(f);
    f = parfevalOnAll(pool, @tune_worker_init, 0, aquaDir); wait(f);
    parfor jj = 1:numel(jobList)
        tune_one_job(jobList(jj), jobs, paramsCsv, nJobs, aquaDir);
    end
else
    for jj = 1:numel(jobList)
        tune_one_job(jobList(jj), jobs, paramsCsv, nJobs, aquaDir);
    end
end
fprintf('%s  tuning done: %.1f h - now: python aqua_tune.py "<mouse folder>" --report\n', datestr(now), toc(tAll) / 3600);
diary off;
"""

TUNE_JOB_FUNCTION = r"""function tune_one_job(j, jobs, paramsCsv, nJobs, aquaDir)
%TUNE_ONE_JOB one job of run_aqua_tune.m: AQuA2 on one movie with one setting,
% export of the events and the mean image, removal of the big result.
if exist('aqua_cmd_batch', 'file') ~= 2
    tune_worker_init(aquaDir);
end
outDir = char(jobs.out_dir(j)); evCsv = char(jobs.events_csv(j)); matPath = char(jobs.mat_path(j));
if exist(evCsv, 'file')
    return;
end
if ~isfolder(outDir), mkdir(outDir); end
fprintf('[%d/%d] %s  %s  setting %d (Preset%d)  %s\n', j, nJobs, char(jobs.run(j)), char(jobs.variant(j)), ...
    jobs.setting(j), jobs.preset(j), datestr(now, 'HH:MM'));
tJob = tic;
try
    aqua_cmd_batch(char(jobs.batch_csv(j)), paramsCsv);
catch err
    fprintf('   [%d] AQuA2 FAILED: %s\n', j, err.message);
    fid = fopen(fullfile(outDir, 'FAILED.txt'), 'w'); fprintf(fid, '%s', err.message); fclose(fid);
    return;
end
try
    S = load(matPath, 'res'); res = S.res; clear S;
    H = size(res.datOrg1, 1); W = size(res.datOrg1, 2);
    [~, mname, ~] = fileparts(matPath);
    imwrite(uint16(round(mean(single(res.datOrg1), 4))), fullfile(outDir, [mname, '_meanimg.tif']));
    fts = res.fts1; nEvt = 0;
    if isstruct(fts) && isfield(fts, 'basic') && isfield(fts.basic, 'area'), nEvt = numel(fts.basic.area); end
    ctr = zeros(nEvt, 2); npx = zeros(nEvt, 1); bb = zeros(nEvt, 4);
    t0 = zeros(nEvt, 1); t1 = zeros(nEvt, 1); dffMax = zeros(nEvt, 1); area = zeros(nEvt, 1); pv = zeros(nEvt, 1);
    for k = 1:nEvt
        c = fts.basic.center{k}; ctr(k, :) = c(1:2);
        if isfield(fts.loc, 'xSpa')
            pix = fts.loc.xSpa{k};
        else
            [ih0, iw0, ~, ~] = ind2sub(size(res.datOrg1), res.evt1{k});
            pix = unique(sub2ind([H, W], ih0, iw0));
        end
        npx(k) = numel(pix);
        [ih, iw] = ind2sub([H, W], pix);
        bb(k, :) = [min(ih) max(ih) min(iw) max(iw)];
        t0(k) = fts.loc.t0(k); t1(k) = fts.loc.t1(k);
        dffMax(k) = fts.curve.dffMax(k); area(k) = fts.basic.area(k);
        if isfield(fts.curve, 'dffMaxPval'), pv(k) = fts.curve.dffMaxPval(k); else, pv(k) = NaN; end
    end
    evT = table((1:nEvt)', t0, t1, ctr(:, 1), ctr(:, 2), npx, bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4), dffMax, area, pv, ...
        'VariableNames', {'event', 't0_frame', 't1_frame', 'y_px', 'x_px', 'n_pixels', 'y0', 'y1', 'x0', 'x1', 'dff_max', 'area_um2', 'p_value'});
    writetable(evT, evCsv);
    fprintf('   [%d] %d events, %.1f min\n', j, nEvt, toc(tJob) / 60);
    clear res;
catch err
    fprintf('   [%d] export FAILED: %s\n', j, err.message);
    fid = fopen(fullfile(outDir, 'FAILED.txt'), 'w'); fprintf(fid, '%s', err.message); fclose(fid);
end
% the big files go; the small feature table (_Ch1.csv) and the events stay
if exist(matPath, 'file'), delete(matPath); end
d = dir(fullfile(outDir, '*_risingMaps_*'));
for q = 1:numel(d)
    if d(q).isdir, rmdir(fullfile(outDir, d(q).name), 's'); end
end
d = dir(fullfile(outDir, '*.xlsx'));
for q = 1:numel(d), delete(fullfile(outDir, d(q).name)); end
end
"""

TUNE_INIT_FUNCTION = r"""function tune_worker_init(aquaDir)
%TUNE_WORKER_INIT AQuA2 on the path of this (worker) process.
if ~isfolder(aquaDir)
    w = which('aqua_cmd_batch');
    if ~isempty(w), aquaDir = fileparts(w); else, error('AQuA2 not found: set the AQuA2 folder in the GUI (top row) or addpath it in MATLAB'); end
end
cd(aquaDir);
run('startup');
end
"""

TUNE_SPLIT_SCRIPT = r"""%% run_aqua_tune_{k}of{n}.m — share {k} of {n} of the tuning jobs (generated by aqua_tune.py)
% Start every share in its own MATLAB (run_aqua_tune_parallel.bat does), each
% takes every {n}th job; a share can be restarted at any time.
tuneSplit = [{k} {n}];
run('{script}');
"""

TUNE_BAT = r"""@echo off
rem run_aqua_tune_parallel.bat — {n} MATLAB processes, one share of the tuning jobs each (generated by aqua_tune.py).
rem Needs "matlab" on the PATH (else open {n} MATLAB windows and run run_aqua_tune_<k>of{n}.m in each).
{lines}
"""


def _mat(s) -> str:
    return str(s).replace("'", "''")


def write_tune_batch(tune_dir: Path, jobs: list[dict], presets: list[dict], template_rows: list, aqua_dir: str,
                     seen=lambda p: str(p), n_workers: int = 0) -> dict:
    """parameters_for_tune.csv (one Preset column per entry of `presets`), one
    one-row batch csv per job, tune_jobs.csv and run_aqua_tune.m. Each job has
    run, variant, setting (the movie's own setting index), preset (1-based
    column of parameters_for_tune.csv), input, out_dir, mat_path, events_csv.
    `seen` maps a path to how the MATLAB machine sees it."""
    tune_dir.mkdir(exist_ok=True)
    (tune_dir / "batch").mkdir(exist_ok=True)
    lines = ["Name,Variable,Type," + ",".join(f"Preset{i + 1}" for i in range(len(presets)))]
    for name, var, typ, default in template_rows:
        if not var:
            lines.append("," * (2 + len(presets)))
            continue
        vals = [str(p.get(var, default)) for p in presets]
        lines.append(",".join([f'"{name}"' if "," in name else name, var, typ] + vals))
    from .aqua_io import write_text_free
    # (a csv left open in Excel or by a stuck reader is locked on Windows: written under the next free
    # name then, and the MATLAB script points at that name — v40)
    params_csv = write_text_free(tune_dir / "parameters_for_tune.csv", "\n".join(lines) + "\n")
    import pandas as pd
    rows = []
    for j, job in enumerate(jobs):
        bpath = tune_dir / "batch" / f"job_{j + 1:04d}.csv"
        bpath.write_text("inputPath,preset,outputPath,cellMaskPath,landmarkMaskPath\n"
                         f"{seen(job['input'])},{int(job['preset'])},{seen(job['out_dir'])},,\n", encoding="utf-8")
        rows.append({"job": j + 1, "run": job["run"], "variant": job["variant"], "setting": job["setting"],
                     "preset": int(job["preset"]),
                     "batch_csv": seen(bpath), "mat_path": seen(job["mat_path"]), "events_csv": seen(job["events_csv"]),
                     "out_dir": seen(job["out_dir"]), "input": seen(job["input"])})
    jobs_csv = write_text_free(tune_dir / "tune_jobs.csv", pd.DataFrame(rows).to_csv(index=False).replace("\r\n", "\n"))
    script = write_text_free(tune_dir / "run_aqua_tune.m", TUNE_SCRIPT.format(
        n_workers=int(n_workers), aqua_dir=_mat(aqua_dir), tune_dir=_mat(seen(tune_dir)),
        jobs_csv=_mat(seen(jobs_csv)), params_csv=_mat(seen(params_csv)),
        log_file=_mat(seen(tune_dir / "tune_log.txt"))))
    (tune_dir / "tune_one_job.m").write_text(TUNE_JOB_FUNCTION, encoding="utf-8")
    (tune_dir / "tune_worker_init.m").write_text(TUNE_INIT_FUNCTION, encoding="utf-8")
    # several MATLAB processes instead of a parallel pool: one script per share + a .bat that starts them
    for old in list(tune_dir.glob("run_aqua_tune_*of*.m")) + [tune_dir / "run_aqua_tune_parallel.bat"]:
        if old.exists():
            old.unlink()
    n_split = max(2, int(n_workers)) if int(n_workers) != 1 else 0
    if n_split:
        lines = []
        for k in range(1, n_split + 1):
            sp = tune_dir / f"run_aqua_tune_{k}of{n_split}.m"
            sp.write_text(TUNE_SPLIT_SCRIPT.format(k=k, n=n_split, script=_mat(seen(script))), encoding="utf-8")
            lines.append(f'start "aqua_tune {k} of {n_split}" matlab -batch "run(\'{_mat(seen(sp))}\')"')
        (tune_dir / "run_aqua_tune_parallel.bat").write_text(
            TUNE_BAT.format(n=n_split, lines="\n".join(lines)) + "\n", encoding="utf-8")
    return {"params_csv": params_csv, "jobs_csv": jobs_csv, "script": script}


# ----------------------------------------------------------------------------- scoring

def crop_offset(raw_mean: np.ndarray, exported_mean: np.ndarray, min_match: float = 0.5):
    """Where AQuA2's registered, cropped mean image sits in the raw frame: (dy, dx, match) or None."""
    try:
        from skimage.feature import match_template
        tmpl = np.asarray(exported_mean, np.float32); ref = np.asarray(raw_mean, np.float32)
        if tmpl.shape[0] > ref.shape[0] or tmpl.shape[1] > ref.shape[1]:
            return None
        cc = match_template(ref, tmpl)
        dy, dx = np.unravel_index(int(np.argmax(cc)), cc.shape)
        if float(cc[dy, dx]) < min_match:
            return None
        return int(dy), int(dx), float(cc[dy, dx])
    except Exception:
        return None


def match_truth(truth, det, t_before: int = 2, t_after: int = 6, pad_px: int = 4, max_dist_px: float = 15.0,
                return_index: bool = False):
    """Greedy one-to-one match of injected transients to detected events:
    detection onset within [t0 - t_before, t0 + t_after] (or the detection
    spans the onset) and the truth centre inside the detection's box (+pad) or
    within max_dist of its centre. Returns a bool array per truth row (with
    return_index: the matched detection's row index per truth row, -1 = none)."""
    hit = np.zeros(len(truth), bool)
    idx = np.full(len(truth), -1, int)
    if det is None or len(det) == 0 or len(truth) == 0:
        return idx if return_index else hit
    used = np.zeros(len(det), bool)
    d_t0 = det["t0_frame"].to_numpy(float); d_t1 = det["t1_frame"].to_numpy(float)
    d_y = det["y_px"].to_numpy(float); d_x = det["x_px"].to_numpy(float)
    have_box = all(c in det.columns for c in ("y0", "y1", "x0", "x1"))
    for i, tr in enumerate(truth.itertuples(index=False)):
        t0 = float(tr.t0_frame); y = float(tr.y_px); x = float(tr.x_px)
        ok_t = ((d_t0 >= t0 - t_before) & (d_t0 <= t0 + t_after)) | ((d_t0 <= t0) & (d_t1 >= t0))
        if have_box:
            inside = ((det["y0"].to_numpy(float) - pad_px <= y) & (det["y1"].to_numpy(float) + pad_px >= y) &
                      (det["x0"].to_numpy(float) - pad_px <= x) & (det["x1"].to_numpy(float) + pad_px >= x))
        else:
            inside = np.zeros(len(det), bool)
        near = np.hypot(d_y - y, d_x - x) <= max_dist_px
        cand = np.flatnonzero(ok_t & (inside | near) & ~used)
        if cand.size:
            k = cand[np.argmin(np.abs(d_t0[cand] - t0))]
            used[k] = True; hit[i] = True; idx[i] = k
    return idx if return_index else hit


def _on_mask(det, fmask_raw, offset) -> np.ndarray:
    """Which detected events (AQuA2 1-based coordinates in its crop at `offset`
    in the raw frame) have their centre on the raw-frame mask."""
    if det is None or len(det) == 0:
        return np.zeros(0, bool)
    yy = np.clip(np.round(det["y_px"].to_numpy(float)).astype(int) - 1 + int(offset[0]), 0, fmask_raw.shape[0] - 1)
    xx = np.clip(np.round(det["x_px"].to_numpy(float)).astype(int) - 1 + int(offset[1]), 0, fmask_raw.shape[1] - 1)
    return fmask_raw[yy, xx]


def _truth_in_crop(truth, offset, crop_shape):
    """Truth table (0-based raw-frame coordinates) in a variant's cropped, 1-based frame; rows outside dropped."""
    tr = truth.copy()
    tr["y_px"] = tr["y_px"] - offset[0] + 1; tr["x_px"] = tr["x_px"] - offset[1] + 1
    if crop_shape is not None:
        inside = ((tr["y_px"] >= 1) & (tr["y_px"] <= crop_shape[0]) & (tr["x_px"] >= 1) & (tr["x_px"] <= crop_shape[1])).to_numpy()
        tr = tr[inside]
    return tr.reset_index(drop=True)


def shift_events(det, from_off, to_off):
    """An exported events table moved from one variant's cropped frame to another's (both offsets = (dy, dx, ...)
    in the raw frame)."""
    if det is None:
        return None
    d = det.copy()
    ddy, ddx = from_off[0] - to_off[0], from_off[1] - to_off[1]
    for c in ("y_px", "y0", "y1"):
        if c in d:
            d[c] = d[c] + ddy
    for c in ("x_px", "x0", "x1"):
        if c in d:
            d[c] = d[c] + ddx
    return d


def score_reference(real, ref, real_off, ref_off, tol_frames: int = 3) -> dict:
    """A setting's real-movie events against a hand-made AQuA2 result of the same movie (v37): ref_recall = share
    of the hand events a setting recovers (onset within tol_frames, centre within the detection's box or 15 px),
    ref_extra = share of the setting's events that match no hand event, ref_n = number of hand events."""
    out = {"ref_n": 0 if ref is None else int(len(ref))}
    if ref is None or not len(ref):
        return out
    ref_r = shift_events(ref, ref_off, real_off)                         # the hand events in this setting's frame
    tr = ref_r.rename(columns={"t0_frame": "t0_frame"})[["t0_frame", "y_px", "x_px"]].reset_index(drop=True)
    if real is None or not len(real):
        out.update(ref_recall=0.0, ref_extra=np.nan, ref_hit=np.zeros(len(tr), bool))
        return out
    hit = match_truth(tr, real, t_before=tol_frames, t_after=tol_frames)
    back = match_truth(real[["t0_frame", "y_px", "x_px"]].reset_index(drop=True), ref_r, t_before=tol_frames, t_after=tol_frames)
    out.update(ref_recall=float(hit.mean()), ref_extra=float(1.0 - back.mean()), ref_hit=hit)
    return out


def score_setting(real, null, inj, truth, minutes: float, fmask_raw=None, offsets=None) -> dict:
    """Metrics of one (movie, setting). `truth` has 0-based raw-frame pixel
    coordinates; the exported events of each variant have AQuA2's 1-based
    coordinates in that variant's cropped, registered frame, which sits at
    offsets[variant] = (dy, dx, crop_shape) in the raw frame. `fmask_raw` is
    the (dilated) fibre mask in the raw frame: everything that matters for the
    afferents is counted on it — on-fibre event rate, on-fibre false alarms —
    and what AQuA2 finds off the fibres is reported separately."""
    offsets = offsets or {}
    out = {}

    def off(variant):
        o = offsets.get(variant) or offsets.get("real") or (0, 0, None)
        return o

    n_real = 0 if real is None else len(real)
    n_null = 0 if null is None else len(null)
    out["n_real"] = n_real; out["rate_real"] = n_real / minutes
    out["n_null"] = n_null if null is not None else np.nan
    out["rate_null"] = (n_null / minutes) if null is not None else np.nan
    out["false_share"] = (min(1.0, n_null / n_real) if (null is not None and n_real > 0) else
                          (0.0 if (null is not None and n_null == 0) else np.nan))
    if real is not None and n_real:
        out["dff_median"] = float(np.nanmedian(real["dff_max"])) if "dff_max" in real else np.nan
        out["duration_median_frames"] = float(np.nanmedian(real["t1_frame"] - real["t0_frame"] + 1))
        out["size_median_px"] = float(np.nanmedian(real["n_pixels"])) if "n_pixels" in real else np.nan
    if fmask_raw is not None:
        on_r = _on_mask(real, fmask_raw, off("real")[:2]) if n_real else np.zeros(0, bool)
        n_on = int(on_r.sum())
        out["n_real_on"] = n_on; out["rate_real_on"] = n_on / minutes; out["rate_off"] = (n_real - n_on) / minutes
        out["on_fibre"] = float(on_r.mean()) if n_real else np.nan
        if null is not None:
            on_n = _on_mask(null, fmask_raw, off("null")[:2]) if n_null else np.zeros(0, bool)
            n_null_on = int(on_n.sum())
            out["n_null_on"] = n_null_on; out["rate_null_on"] = n_null_on / minutes
            out["false_share_on"] = (min(1.0, n_null_on / n_on) if n_on > 0 else (0.0 if n_null_on == 0 else np.nan))
    if truth is not None and len(truth):
        o_inj = off("inject"); o_real = off("real")
        tr = _truth_in_crop(truth, o_inj[:2], o_inj[2])
        out["n_injected_outside_crop"] = int(len(truth) - len(tr))
        out["n_injected"] = int(len(tr))
        if not len(tr):
            return out
        hit = match_truth(tr, inj) if inj is not None else np.zeros(len(tr), bool)
        # chance coincidence with the real movie's own detections, in the real variant's frame
        tr_real = _truth_in_crop(truth, o_real[:2], None)
        keep = tr_real["event"].isin(tr["event"]).to_numpy()
        tr_real = tr_real[keep].reset_index(drop=True)
        coin = match_truth(tr_real, real) if real is not None else np.zeros(len(tr), bool)
        if coin.size != hit.size:
            coin = np.zeros(hit.size, bool)
        out["recall_raw"] = float(hit.mean())
        c = float(coin.mean())
        out["coincidence"] = c
        out["recall"] = float(max(0.0, (hit.mean() - c) / max(1 - c, 1e-9)))
        # what AQuA2 reports for the injected transients it found (footprint-average dF/F, v37): tells whether the
        # nominal classes sit below, within or above the real events of this movie
        if inj is not None and len(inj) and "dff_max" in inj and hit.any():
            try:
                mi = match_truth(tr, inj, return_index=True)
                got = mi[hit & (mi >= 0)]
                out["inj_dff_reported"] = float(np.nanmedian(inj["dff_max"].to_numpy(float)[got])) if got.size else np.nan
                for a in sorted(tr["amplitude"].unique()):
                    m = ((tr["amplitude"] == a).to_numpy()) & hit & (mi >= 0)
                    out[f"inj_dff_reported_amp{a:g}"] = float(np.nanmedian(inj["dff_max"].to_numpy(float)[mi[m]])) if m.any() else np.nan
            except Exception:
                pass
        for a in sorted(tr["amplitude"].unique()):
            m = (tr["amplitude"] == a).to_numpy()
            ca = float(coin[m].mean()) if m.any() else 0.0
            out[f"recall_amp{a:g}"] = float(max(0.0, (hit[m].mean() - ca) / max(1 - ca, 1e-9))) if m.any() else np.nan
        for k in sorted(tr["kind"].unique()):
            m = (tr["kind"] == k).to_numpy()
            ck = float(coin[m].mean()) if m.any() else 0.0
            out[f"recall_{k}"] = float(max(0.0, (hit[m].mean() - ck) / max(1 - ck, 1e-9))) if m.any() else np.nan
        for s in sorted(tr["size_class"].unique()):
            m = (tr["size_class"] == s).to_numpy()
            cs = float(coin[m].mean()) if m.any() else 0.0
            out[f"recall_size{int(s)}"] = float(max(0.0, (hit[m].mean() - cs) / max(1 - cs, 1e-9))) if m.any() else np.nan
    return out


def false_key(rows: list[dict]) -> str:
    """The false-alarm measure to judge by: on the fibres when it was measured, else overall."""
    return "false_share_on" if any(np.isfinite(r.get("false_share_on", np.nan)) for r in rows) else "false_share"


GRACE_NULL_EVENTS = 1


def is_clean(r: dict, fk: str, max_false: float, grace: int = GRACE_NULL_EVENTS) -> bool:
    """Does a setting keep the false-alarm share at or below max_false? With one null event of grace (v38): the
    share is a ratio of small counts (4 null events among 36 real = 11 %, 1 among 4 = 25 %), and a single null
    event should never decide — the setting also passes when the share without `grace` null events would."""
    v = r.get(fk, np.nan)
    if not np.isfinite(v):
        return False
    if v <= max_false:
        return True
    on = fk == "false_share_on"
    n_null = r.get("n_null_on" if on else "n_null", np.nan); n_real = r.get("n_real_on" if on else "n_real", np.nan)
    try:
        n_null = float(n_null); n_real = float(n_real)
    except (TypeError, ValueError):
        return False
    # the grace needs a few real events to lean on: with 1-3 real events one null event is 33-100 % of them, and a
    # setting that finds almost nothing must not be called clean on the strength of a forgiven event
    if np.isfinite(n_null) and np.isfinite(n_real) and n_real >= 4 and n_null >= 1:
        return max(0.0, n_null - grace) / n_real <= max_false
    return False



def recommend(rows: list[dict], max_false: float = 0.1, min_on_fibre: float = 0.8, target_amp: float | None = None):
    """Pick the setting: highest recall (at the target amplitude, else overall)
    among those whose false-alarm share — on the fibres when measured — is
    <= max_false and whose on_fibre fraction is >= min_on_fibre (that rule is
    skipped when no setting reaches it). When no setting keeps the false-alarm
    share low enough, the best trade-off — highest recall x (1 - false share)
    — is returned instead, flagged as not passing.
    Returns (chosen row or None, reason, passed)."""
    if not rows:
        return None, "no scored settings", False
    key = f"recall_amp{target_amp:g}" if target_amp is not None and any(f"recall_amp{target_amp:g}" in r for r in rows) else "recall"
    fk = false_key(rows)
    fname = "false-alarm share on the fibres" if fk == "false_share_on" else "false-alarm share"
    have_recall = [r for r in rows if np.isfinite(r.get(key, np.nan))]
    if not have_recall:
        return None, "no injected-movie results yet", False
    scored = [r for r in have_recall if np.isfinite(r.get(fk, np.nan))]
    if not scored:
        return None, "no null-movie results yet", False
    ok = [r for r in scored if is_clean(r, fk, max_false)]
    if not ok:
        best = max(scored, key=lambda r: (r.get(key, 0.0) * (1 - r[fk]), -r[fk]))
        return best, (f"no setting keeps the {fname} at or below {max_false:.0%} (lowest: "
                      f"{min(r[fk] for r in scored):.0%}); shown: the best trade-off, highest {key} x (1 - false share) "
                      f"— accept it with a higher --max-false, or treat this field as not reliably detectable at these settings"), False
    note = ""
    ok2 = [r for r in ok if np.isfinite(r.get("on_fibre", np.nan)) and r["on_fibre"] >= min_on_fibre]
    if ok2:
        ok = ok2
    else:
        note = f" (no setting has >= {min_on_fibre:.0%} of its events on the fibres — rule skipped)"
    best = max(ok, key=lambda r: (r.get(key, 0.0), -r.get(fk, 1.0), r.get("n_real_on", r.get("n_real", 0))))
    return best, f"highest {key} with {fname} <= {max_false:.0%}{note}", True


def _step_beyond(key: str, v: float, permissive: bool) -> float | None:
    """One grid step beyond value v of variable `key`, towards more (permissive=True) or fewer events."""
    lower = (key in PERMISSIVE_LOWER) or (key not in PERMISSIVE_HIGHER)       # unknown variables: lower = more
    towards_lower = lower == permissive
    if key == "thrARScl":
        nv = (v - 0.5 if v > 2.0 else v - 0.25) if towards_lower else (v + 0.5)
        nv = max(1.0, nv)
    elif key == "minSize":
        nv = max(10, int(round(v / 2.5 / 5.0)) * 5) if towards_lower else int(v * 2)
    elif key == "minDur":
        nv = max(1, int(v) - 1) if towards_lower else int(v) + 1
    else:
        nv = v * 0.8 if towards_lower else v * 1.25
    if towards_lower and nv >= v or (not towards_lower and nv <= v):
        return None
    return round(float(nv), 4)


def grid_verdict(rows: list[dict], best: dict | None, passed: bool, grid_keys: list[str], max_false: float = 0.1,
                 min_recall: float = 0.5, edge_recall: float = 0.7) -> dict:
    """Does the grid bracket the optimum of this movie, or does it only touch it (v37)?
    'too strict': the best clean setting is the most permissive value swept of some variable and its recall is
    below min_recall — the real optimum lies beyond the grid (2609-1 FOV 2: recall 0.22 at thrARScl 1.5 /
    minSize 50, both the lowest swept, while the hand setting minSize 20 finds 9 events/min); 'edge': the same
    with a recall below edge_recall — more sensitivity is probably available one step further; 'too permissive':
    no setting keeps the false alarms low; 'weak': low recall although the grid brackets it (a hard field);
    'ok' otherwise. `suggest` = {variable: [value]} one step beyond the edge for the variables concerned."""
    out = {"status": "unscored", "suggest": {}, "text": "", "recall": float("nan")}
    if not rows or best is None:
        return out
    key = "recall"
    rec = float(best.get(key, np.nan))
    out["recall"] = rec
    swept = {k: sorted({float(r[k]) for r in rows if np.isfinite(_num(r.get(k)))}) for k in grid_keys}
    def edge_val(k, permissive):
        vals = swept.get(k) or []
        if not vals:
            return None
        lower = (k in PERMISSIVE_LOWER) or (k not in PERMISSIVE_HIGHER)
        return vals[0] if lower == permissive else vals[-1]
    if not passed:
        # nothing clean: stricter values are needed beyond the strict edge
        sug = {}
        for k in grid_keys:
            e = edge_val(k, permissive=False)
            if e is not None and float(best.get(k, np.nan)) == e:
                nv = _step_beyond(k, e, permissive=False)
                if nv is not None:
                    sug[k] = [nv]
        if not sug:
            for k in grid_keys:
                e = edge_val(k, permissive=False)
                nv = _step_beyond(k, e, permissive=False) if e is not None else None
                if nv is not None:
                    sug[k] = [nv]
        out.update(status="too permissive", suggest=sug,
                   text=f"GRID TOO PERMISSIVE: no setting keeps the false alarms <= {max_false:.0%}"
                        + (" — sweep stricter values: " + ", ".join(f"{k}={_fmt(v[0])}" for k, v in sug.items()) if sug else ""))
        return out
    at_edge = [k for k in grid_keys if edge_val(k, True) is not None and float(best.get(k, np.nan)) == edge_val(k, True)]
    sug = {}
    for k in at_edge:
        nv = _step_beyond(k, float(best[k]), permissive=True)
        if nv is not None:
            sug[k] = [nv]
    if at_edge and np.isfinite(rec) and rec < min_recall:
        out.update(status="too strict", suggest=sug,
                   text=f"GRID TOO STRICT: the best clean setting is the most permissive one swept ({', '.join(f'{k}={_fmt(best[k])}' for k in at_edge)}) "
                        f"and still finds only {rec:.0%} of the injected transients — the optimum lies beyond the grid"
                        + (": sweep " + ", ".join(f"{k}={_fmt(v[0])}" for k, v in sug.items()) if sug else ""))
    elif at_edge and np.isfinite(rec) and rec < edge_recall:
        out.update(status="edge", suggest=sug,
                   text=f"at the edge of the grid ({', '.join(f'{k}={_fmt(best[k])}' for k in at_edge)}, recall {rec:.2f}): more sensitivity "
                        f"may be available one step further" + (" — " + ", ".join(f"{k}={_fmt(v[0])}" for k, v in sug.items()) if sug else ""))
    elif np.isfinite(rec) and rec < min_recall:
        out.update(status="weak", text=f"WEAK: only {rec:.0%} of the injected transients found although more permissive settings were swept "
                                        f"(they fail the false-alarm or on-fibre rule) — a hard field; check the fibre mask and the injected classes")
    else:
        out.update(status="ok", text="")
    return out


def _num(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def recommend_common(per_movie: dict, grid_keys: list[str], max_false: float = 0.1, min_on_fibre: float = 0.8,
                     target_amp: float | None = None):
    """One setting for several movies (the days of one FOV): among the grid
    settings scored on every movie, the highest mean recall whose false-alarm
    share stays <= max_false on EVERY movie (on-fibre rule as in recommend,
    skipped when no setting meets it on all movies). Falls back to the best
    trade-off, flagged as not passing. Returns (setting tuple, per-movie rows,
    reason, passed) or (None, {}, reason, False)."""
    runs = [r for r in per_movie if per_movie[r]]
    if not runs:
        return None, {}, "no scored settings", False
    key = f"recall_amp{target_amp:g}" if target_amp is not None and any(f"recall_amp{target_amp:g}" in x for r in runs for x in per_movie[r]) else "recall"
    fk = false_key([x for r in runs for x in per_movie[r]])
    fname = "false-alarm share on the fibres" if fk == "false_share_on" else "false-alarm share"
    table: dict = {}
    for r in runs:
        for x in per_movie[r]:
            if np.isfinite(x.get(key, np.nan)) and np.isfinite(x.get(fk, np.nan)):
                table.setdefault(tuple(x[k] for k in grid_keys), {})[r] = x
    cands = {c: rows for c, rows in table.items() if len(rows) == len(runs)}
    if not cands:
        return None, {}, "no setting was scored on every movie yet", False
    ok = {c: rows for c, rows in cands.items() if all(is_clean(x, fk, max_false) for x in rows.values())}
    if not ok:
        best = max(cands, key=lambda c: (np.mean([x.get(key, 0.0) * (1 - x[fk]) for x in cands[c].values()]),
                                         -max(x[fk] for x in cands[c].values())))
        return best, cands[best], (f"no setting keeps the {fname} at or below {max_false:.0%} on every day; shown: the best "
                                   f"trade-off (mean recall x (1 - false share))"), False
    note = ""
    ok2 = {c: rows for c, rows in ok.items()
           if all(np.isfinite(x.get("on_fibre", np.nan)) and x["on_fibre"] >= min_on_fibre for x in rows.values())}
    if ok2:
        ok = ok2
    else:
        note = f" (no setting has >= {min_on_fibre:.0%} of its events on the fibres on every day — rule skipped)"
    best = max(ok, key=lambda c: (np.mean([x.get(key, 0.0) for x in ok[c].values()]), -max(x[fk] for x in ok[c].values())))
    return best, ok[best], f"highest mean {key} over the days with {fname} <= {max_false:.0%} on every day{note}", True


def _fp_rate(row: dict) -> float:
    """The false-alarm rate (null-movie events per minute, on the fibres when measured)."""
    v = row.get("rate_null_on", np.nan)
    if not np.isfinite(v):
        v = row.get("rate_null", np.nan)
    return float(v)


def recommend_matched(per_movie: dict, grid_keys: list[str], max_false: float = 0.1, min_on_fibre: float = 0.8,
                      tol_recall: float = 0.05, tol_fp_frac: float = 0.10, fp_floor: float = 0.05,
                      target_amp: float | None = None, min_level: float = 0.3):
    """Matched-sensitivity settings for the movies of one day (pre and post, v36):
    one setting PER MOVIE, chosen so that the movies detect the same injected
    transients with the same recall (within tol_recall) and invent false events
    at the same rate (null-movie events/min on the fibres within tol_fp =
    max(fp_floor, tol_fp_frac x the lower real on-fibre rate)), every setting
    clean on its own movie (false share <= max_false, with one null event of
    grace; on-fibre >= min_on_fibre when reachable). Among the matched
    combinations the highest common sensitivity (max of the lowest recall)
    wins; ties: the smaller recall gap, then the fewer false alarms. Movies
    that pass the rule with the same best setting get that one setting. v38:
    when the on-fibre rule leaves only a match below min_level, it is relaxed
    for the movies it limits (noted in info["on_fibre_relaxed"]) — a setting
    that finds the fibre events and some off-fibre ones besides is better
    than one that finds nothing; the off-fibre events are for the curation.
    Returns (choice, rows, reason, passed, info) with choice = {run: setting
    tuple} (None when nothing can be matched), rows = {run: scored row},
    info = {level, recall_gap, fp_gap, fp_tol, same_setting, ...}."""
    import itertools
    runs = [r for r in per_movie if per_movie[r]]
    if len(runs) < 2:
        return None, {}, "fewer than two scored movies", False, {}
    if len(runs) > 3:
        return None, {}, "matched sensitivity is defined for the pre and post movies of one day (2-3 movies)", False, {}
    key = f"recall_amp{target_amp:g}" if target_amp is not None and any(f"recall_amp{target_amp:g}" in x for r in runs for x in per_movie[r]) else "recall"
    fk = false_key([x for r in runs for x in per_movie[r]])
    clean: dict = {}; strict: dict = {}
    for r in runs:
        ok = [x for x in per_movie[r] if np.isfinite(x.get(key, np.nan)) and is_clean(x, fk, max_false)]
        ok2 = [x for x in ok if np.isfinite(x.get("on_fibre", np.nan)) and x["on_fibre"] >= min_on_fibre]
        clean[r] = ok; strict[r] = ok2 or ok
        if not ok:
            return None, {}, f"{r}: no setting keeps the false-alarm share at or below {max_false:.0%}", False, {}

    def search(cands):
        real_rates = [max(x.get("rate_real_on", x.get("rate_real", 0.0)) for x in cands[r]) for r in runs]
        fp_tol = max(fp_floor, tol_fp_frac * min(real_rates)) if real_rates else fp_floor
        best = None; best_same = None
        for combo in itertools.product(*[cands[r] for r in runs]):
            rec = [float(x.get(key, 0.0)) for x in combo]
            fps = [_fp_rate(x) for x in combo]
            gap_r = max(rec) - min(rec)
            gap_f = (max(fps) - min(fps)) if all(np.isfinite(fps)) else 0.0
            if gap_r > tol_recall or gap_f > fp_tol:
                continue
            score = (min(rec), -gap_r, -sum(f for f in fps if np.isfinite(f)))
            if best is None or score > best[0]:
                best = (score, combo)
            if len({tuple(x[k] for k in grid_keys) for x in combo}) == 1 and (best_same is None or score > best_same[0]):
                best_same = (score, combo)
        # the same setting for every movie whenever it matches within the tolerance (v38): equal recall of the injected
        # transients makes two settings equivalent by our measure, and one setting cannot favour either condition on
        # the real events (a sham day with pre thr 1.5/50 and post thr 2/20 showed a 2x "increase" that one setting
        # for both did not)
        if best_same is not None and best is not None and best_same[0][0] >= best[0][0] - tol_recall:
            best = best_same
        return best, fp_tol

    best, fp_tol = search(strict)
    relaxed = []
    if (best is None or best[0][0] < min_level) and any(len(clean[r]) > len(strict[r]) for r in runs):
        best2, fp_tol2 = search(clean)
        if best2 is not None and best2[0][0] >= min_level and (best is None or best2[0][0] > best[0][0]):
            best, fp_tol = best2, fp_tol2
            relaxed = [r for r, x in zip(runs, best[1]) if x not in strict[r]]
    if best is None:
        tops = {r: max(strict[r], key=lambda x: x.get(key, 0.0)) for r in runs}
        reach = "; ".join(f"{r.split('_')[-3] + '_' + r.split('_')[-2] if r.count('_') >= 3 else r}: at most {tops[r].get(key, 0.0):.2f}" for r in runs)
        return None, tops, (f"no combination of clean settings detects the movies with the same sensitivity (within {tol_recall:.2f}) "
                            f"and the same false-alarm rate (within {fp_tol:.2f}/min) — {reach}"), False, {"fp_tol": fp_tol}
    combo = best[1]
    choice = {r: tuple(x[k] for k in grid_keys) for r, x in zip(runs, combo)}
    rows = {r: x for r, x in zip(runs, combo)}
    rec = [float(x.get(key, 0.0)) for x in combo]
    fps = [_fp_rate(x) for x in combo]
    info = {"level": min(rec), "recall_gap": max(rec) - min(rec), "fp_gap": (max(fps) - min(fps)) if all(np.isfinite(fps)) else np.nan,
            "fp_tol": fp_tol, "same_setting": len(set(choice.values())) == 1, "key": key, "on_fibre_relaxed": relaxed}
    if min(rec) < min_level:
        # matched, but at a sensitivity too low to mean anything (v37): both movies miss most transients — the grid
        # is too strict for this field (or the field is not detectable); the day is not comparable as it stands
        info["weak"] = True
        tops = {r: max(strict[r], key=lambda x: x.get(key, 0.0)) for r in runs}
        reach = "; ".join(f"{r.split('_')[-3] + '_' + r.split('_')[-2] if r.count('_') >= 3 else r}: at most {tops[r].get(key, 0.0):.2f}" for r in runs)
        return choice, rows, (f"TOO WEAK — the best common sensitivity is {min(rec):.2f} (< {min_level:.2f}): the movies can be matched only "
                              f"at a level where most injected transients are missed ({reach}) — extend the grid towards more "
                              f"permissive values (the report's grid verdict says which) before trusting this day"), False, info
    why = (f"highest common sensitivity ({min(rec):.2f}) with recalls within {tol_recall:.2f} and false-alarm rates within "
           f"{fp_tol:.2f}/min, every setting with the {('false-alarm share on the fibres' if fk == 'false_share_on' else 'false-alarm share')} "
           f"<= {max_false:.0%}")
    if relaxed:
        why += ("; the on-fibre rule (>= " + f"{min_on_fibre:.0%}" + ") was relaxed for "
                + ", ".join(f"{r.split('_')[-3] + '_' + r.split('_')[-2] if r.count('_') >= 3 else r} ({rows[r].get('on_fibre', np.nan):.0%} on fibres)" for r in relaxed)
                + " — without it the day could only be matched below the useful level; the off-fibre events are for the curation")
    return choice, rows, why, True, info


# ----------------------------------------------------------------------------- plain words, candidates, summary figure (v38)

def setting_short(r: dict, grid_keys: list[str]) -> str:
    """'2/20' for thrARScl=2 minSize=20 (the grid values joined with '/')."""
    return "/".join(_fmt(r[k]) if np.isfinite(_num(r.get(k))) else str(r.get(k)) for k in grid_keys)


def setting_words(r: dict, grid_keys: list[str]) -> str:
    """'thrARScl 2 / minSize 20'."""
    return " / ".join(f"{k} {_fmt(r[k]) if np.isfinite(_num(r.get(k))) else r.get(k)}" for k in grid_keys)


def noise_pct(r: dict, fk: str) -> float:
    v = r.get(fk, np.nan)
    return float(min(1.0, v) * 100) if np.isfinite(v) else np.nan


def runner_up(rows: list[dict], best: dict | None, fk: str, max_false: float, min_on_fibre: float, min_gain: float = 0.05):
    """The most sensitive setting the rule rejected although it is not absurd: recall at least min_gain above the
    chosen one, noise at most 50 % — what a person tuning by eye would be tempted by. None when there is none."""
    if best is None:
        return None
    cands = [r for r in rows if np.isfinite(r.get("recall", np.nan)) and r.get("recall", 0) >= best.get("recall", 0) + min_gain
             and np.isfinite(r.get(fk, np.nan)) and r[fk] <= 0.5 and r is not best]
    if not cands:
        return None
    return max(cands, key=lambda r: (r.get("recall", 0.0), -r.get(fk, 1.0)))


def candidates(rows: list[dict], best: dict | None, grid_keys: list[str], fk: str, max_false: float, min_on_fibre: float,
               n_max: int = 4) -> list[tuple]:
    """The few settings worth a look by eye (v38): the rule's choice, the more sensitive one it rejected, the current
    setting in the batch, and the next cleaner clean setting. Returns [(row, tag, colour)]."""
    out = []
    seen = set()

    def add(r, tag, col):
        if r is None:
            return
        key = tuple(r.get(k) for k in grid_keys)
        if key in seen:
            return
        seen.add(key); out.append((r, tag, col))
    add(best, "chosen by the rule", "#1a7f2e")
    ru = runner_up(rows, best, fk, max_false, min_on_fibre)
    if ru is not None:
        why = (f"{noise_pct(ru, fk):.0f} % noise" if not is_clean(ru, fk, max_false) else f"only {ru.get('on_fibre', np.nan):.0%} on fibres")
        add(ru, f"more sensitive, rejected: {why}", "#d35400")
    cur = next((r for r in rows if r.get("setting") == 0), None)
    add(cur, "in the batch now", "#7f7f7f")
    if best is not None:
        cleaner = [r for r in rows if is_clean(r, fk, max_false) and np.isfinite(r.get(fk, np.nan)) and r[fk] < best.get(fk, 0)
                   and r.get("recall", 0) >= 0.6 * best.get("recall", 0) and r is not best]
        if cleaner:
            add(max(cleaner, key=lambda r: r.get("recall", 0.0)), "cleaner, less sensitive", "#2b6cb0")
    return out[:n_max]


def plain_movie_words(rows: list[dict], best: dict | None, passed: bool, grid_keys: list[str], fk: str, max_false: float,
                      min_on_fibre: float) -> list[str]:
    """Two or three plain sentences about one movie's result."""
    if best is None:
        return ["no setting scored yet"]
    minutes = float(best.get("minutes", np.nan))
    n_on = best.get("n_real_on", best.get("n_real", 0)) or 0
    rate = best.get("rate_real_on", best.get("rate_real", np.nan))
    fp = _fp_rate(best)
    out = [(f"{'Chosen' if passed else 'Best trade-off (rule not met)'}: {setting_words(best, grid_keys)} — "
            f"{rate:.2f} events/min on the fibres ({int(n_on)} in {minutes:.0f} min), of which {fp:.2f}/min ({noise_pct(best, fk):.0f} %) "
            f"would also appear in a movie without signal; catches {100 * best.get('recall', np.nan):.0f} % of the planted events"
            + (f"; {best.get('on_fibre', np.nan):.0%} of its events on the fibres" if np.isfinite(best.get("on_fibre", np.nan)) else "") + ".")]
    if int(best.get("n_real", 0) or 0) == 0:
        out.append("It finds no event in this movie: the field is silent at any clean setting, or the detection is too strict for it.")
    ru = runner_up(rows, best, fk, max_false, min_on_fibre)
    if ru is not None:
        r2 = ru.get("rate_real_on", ru.get("rate_real", np.nan)); n2 = ru.get("n_real_on", ru.get("n_real", 0)) or 0
        nn = ru.get("n_null_on", ru.get("n_null", 0)) or 0
        if not is_clean(ru, fk, max_false):
            out.append(f"More sensitive but rejected: {setting_words(ru, grid_keys)} — catches {100 * ru.get('recall', np.nan):.0f} % and finds "
                       f"{r2:.2f}/min ({int(n2)} events), but {noise_pct(ru, fk):.0f} % of them would be noise ({int(nn)} events in the scrambled "
                       f"movie; the limit is {max_false:.0%}). To take it anyway: max false-alarm share {min(1.0, np.ceil(ru[fk] * 20) / 20):.2f} "
                       f"(the analysis subtracts the expected noise from the rates, not from the CFU maps).")
        else:
            out.append(f"More sensitive but rejected: {setting_words(ru, grid_keys)} — catches {100 * ru.get('recall', np.nan):.0f} %, {r2:.2f}/min, "
                       f"noise {noise_pct(ru, fk):.0f} %, but only {ru.get('on_fibre', np.nan):.0%} of its events lie on the fibres (limit {min_on_fibre:.0%}).")
    return out


def plain_day_words(conds: dict, chosen: dict, mchoice, mrows, mpass: bool, minfo: dict, mwhy: str, grid_keys: list[str], fk: str,
                    common=None, crow=None) -> list[str]:
    """Plain sentences about a day's pre/post choice."""
    out = []
    runs = list(conds)
    if mpass and mchoice:
        same = minfo.get("same_setting")
        parts = []
        for r in runs:
            x = mrows.get(r)
            if x is None:
                continue
            parts.append(f"{conds[r]} {x.get('rate_real_on', x.get('rate_real', np.nan)):.2f}/min")
        if same:
            st = setting_words(mrows[runs[0]], grid_keys)
            out.append(f"Pre and post both at {st}: catch " + " / ".join(f"{100 * mrows[r].get('recall', np.nan):.0f} %" for r in runs if r in mrows)
                       + ", noise " + " / ".join(f"{noise_pct(mrows[r], fk):.0f} %" for r in runs if r in mrows)
                       + "; events on the fibres " + " → ".join(f"{mrows[r].get('rate_real_on', mrows[r].get('rate_real', np.nan)):.2f}" for r in runs if r in mrows) + " per min.")
        else:
            out.append("Matched at different settings — " + "; ".join(f"{conds[r]} at {setting_words(mrows[r], grid_keys)} (catch {100 * mrows[r].get('recall', np.nan):.0f} %, "
                                                                   f"noise {noise_pct(mrows[r], fk):.0f} %)" for r in runs if r in mrows)
                       + "; events on the fibres " + " → ".join(f"{mrows[r].get('rate_real_on', mrows[r].get('rate_real', np.nan)):.2f}" for r in runs if r in mrows) + " per min.")
        if minfo.get("on_fibre_relaxed"):
            out.append("The on-fibre rule was relaxed for " + ", ".join(f"{conds.get(r, r)} ({mrows[r].get('on_fibre', np.nan):.0%} on fibres)" for r in minfo["on_fibre_relaxed"] if r in mrows)
                       + " — its off-fibre events are for the curation.")
    else:
        out.append("NOT COMPARABLE: " + (mwhy.split(" — ")[0] if mwhy else "pre and post cannot be detected with the same sensitivity") + ".")
        if common is not None and crow:
            out.append(f"The day gets the common setting {' / '.join(f'{k} {_fmt(v)}' for k, v in zip(grid_keys, common))}: catch "
                       + " / ".join(f"{100 * crow[r].get('recall', np.nan):.0f} %" for r in runs if r in crow)
                       + "; events on the fibres " + " → ".join(f"{crow[r].get('rate_real_on', crow[r].get('rate_real', np.nan)):.2f}" for r in runs if r in crow) + " per min.")
    return out


def tune_summary_figure(groups: list, per_movie: dict, chosen: dict, passed: dict, grid_keys: list[str], max_false: float,
                        min_on_fibre: float, day_matched: dict, conds: dict, out_png: Path, dpi: int = 150, n_planted: int = 72,
                        classes: str = "", day_common: dict | None = None) -> Path | None:
    """The tuning at a glance (v38): one panel per experiment (pre and post together) or per movie, every setting a
    dot at (noise, catch) — right = invents more, up = finds more of the planted events — sized by the events/min it
    reports on the fibres, labelled with its values, the rule's pick starred, the pair the day uses ringed in green,
    the noise limit as a red line; plain sentences under each panel. `groups` = [(title, exp, [runs])]."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if not groups:
        return None
    ncol = 2 if len(groups) > 1 else 1
    nrow = int(np.ceil(len(groups) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(8.6 * ncol, 6.3 * nrow + 1.6), squeeze=False)
    cond_cols = {"pre": "#3b7dd8", "post": "#e8702a"}
    fallback_cols = ["#3b7dd8", "#e8702a", "#2ca02c"]
    for gi, (title, exp, runs) in enumerate(groups):
        ax = axes[gi // ncol][gi % ncol]
        rows_all = [r for run in runs for r in per_movie.get(run, [])]
        if not rows_all:
            ax.set_axis_off(); ax.set_title(title + " — not scored yet", fontsize=10, loc="left"); continue
        fk = false_key(rows_all)
        max_rate = max([r.get("rate_real_on", r.get("rate_real", 0.0)) or 0.0 for r in rows_all] + [1e-6])
        dm = day_matched.get(exp) if day_matched else None
        mchoice = dm[0] if (dm and dm[2] and dm[0]) else None
        ring_pts = []
        for ri, run in enumerate(runs):
            rows = per_movie.get(run, [])
            if not rows:
                continue
            cond = conds.get(run, "")
            col = cond_cols.get(cond, fallback_cols[ri % 3])
            best = chosen.get(run)
            for r in rows:
                x = noise_pct(r, fk); y = 100 * r.get("recall", np.nan)
                if not (np.isfinite(x) and np.isfinite(y)):
                    continue
                rate = r.get("rate_real_on", r.get("rate_real", 0.0)) or 0.0
                size = 30 + 420 * np.sqrt(rate / max_rate)
                hollow = np.isfinite(r.get("on_fibre", np.nan)) and r["on_fibre"] < min_on_fibre and (r.get("n_real", 0) or 0) >= 5
                ax.scatter([x], [y], s=size, facecolor="none" if hollow else col, edgecolor=col, linewidths=1.4, alpha=0.85, zorder=3)
                ax.annotate(setting_short(r, grid_keys), (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7.5, color="0.25", zorder=4)
                if r.get("setting") == 0:
                    ax.scatter([x], [y], s=size + 160, facecolor="none", edgecolor="#d62728", linestyle=":", linewidths=1.2, zorder=4)
                if best is not None and r is best:
                    ax.scatter([x], [y], marker="*", s=170, color="black", zorder=6)
                if mchoice is not None and run in mchoice and tuple(r[k] for k in grid_keys) == tuple(mchoice[run]):
                    ax.scatter([x], [y], s=size + 260, facecolor="none", edgecolor="#2ecc40", linewidths=2.4, zorder=5)
                    ring_pts.append((x, y))
        if len(ring_pts) == 2:
            ax.plot([ring_pts[0][0], ring_pts[1][0]], [ring_pts[0][1], ring_pts[1][1]], "--", color="#2ecc40", lw=1.2, zorder=2)
        ax.axvline(100 * max_false, color="#d62728", ls="--", lw=1.2)
        ax.text(100 * max_false + 1.5, 97, f"noise limit {max_false:.0%}", color="#d62728", fontsize=8, va="top")
        ax.set_xlim(-3, 103); ax.set_ylim(-3, 103)
        ax.set_xlabel("NOISE — events it invents in the scrambled movie, as % of the events it reports in the real one", fontsize=8.5)
        ax.set_ylabel(f"CATCH — % of the {n_planted} planted events it finds", fontsize=8.5)
        ax.grid(alpha=0.3); ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=10.5, loc="left", fontweight="bold")
        handles = []
        for run in runs:
            cond = conds.get(run, "")
            handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=cond_cols.get(cond, fallback_cols[runs.index(run) % 3]), markersize=9,
                                  label=f"{cond or run} ({run.split('_')[2] if run.count('_') >= 2 else run})"))
        handles += [Line2D([0], [0], marker="*", color="w", markerfacecolor="black", markersize=12, label="the rule's pick"),
                    Line2D([0], [0], marker="o", color="#2ecc40", markerfacecolor="none", markersize=11, markeredgewidth=2, linestyle="none", label="what the day uses (matched pair)"),
                    Line2D([0], [0], marker="o", color="#d62728", markerfacecolor="none", markersize=11, linestyle="none", label="in the batch now")]
        ax.legend(handles=handles, fontsize=7.5, loc="lower right", frameon=True, framealpha=0.85)
        # plain words under the panel
        words = []
        if len(runs) >= 2 and dm is not None:
            mchoice_, mrows_, mpass_, minfo_, mwhy_, conds_ = dm
            dc = (day_common or {}).get(exp)
            common_, crow_ = (dc[0], dc[1]) if dc else (None, None)
            words += plain_day_words(conds_, chosen, mchoice_, mrows_, mpass_, minfo_, mwhy_, grid_keys, fk, common=common_, crow=crow_)
        for run in runs:
            rows = per_movie.get(run, [])
            if rows and chosen.get(run) is not None:
                w = plain_movie_words(rows, chosen[run], passed.get(run, False), grid_keys, fk, max_false, min_on_fibre)
                words.append(f"{conds.get(run, '') or run}: " + w[0])
                words += ["    " + x for x in w[1:]]
        import textwrap
        wrapped = []
        for w in words:
            wrapped += textwrap.wrap(w, 118, subsequent_indent="      ") or [""]
        ax.text(0, -0.17, "\n".join(wrapped[:9]), transform=ax.transAxes, fontsize=7.6, va="top", ha="left", family="DejaVu Sans")
    for gi in range(len(groups), nrow * ncol):
        axes[gi // ncol][gi % ncol].set_axis_off()
    fig.suptitle("HOW A SETTING IS MEASURED — on each tuned movie, two numbers per setting. CATCH: 72 fake transients are planted on the fibres "
                 f"({classes}); catch = the share the setting finds.\nNOISE: in a copy of the movie with every pixel's time course scrambled "
                 "(no real events left) everything found is invented; noise = those events as a share of what the setting reports in the real movie.\n"
                 f"THE PICK (star): the highest dot left of the red line — most catch with noise <= {max_false:.0%} (and >= {min_on_fibre:.0%} of its events on "
                 "the fibres; hollow dots fail that). Dot size = events/min it reports on the fibres. A day uses one setting for pre and post whenever "
                 "one matches both (green rings).", fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.02, 1, 0.955), h_pad=6.0)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return out_png


# ----------------------------------------------------------------------------- figure

def tune_figure(per_movie: dict, out_png: Path, grid_keys: list[str], settings: list[dict], chosen: dict, titles: dict | None = None,
                dpi: int = 150, notes: dict | None = None):
    """One row per movie, four panels: recall vs the first grid variable (one
    line per value of the second), recall per injected amplitude (at the
    chosen setting's value of the second variable), false-alarm share, and the
    real movie's event rate with the on-fibre fraction. `notes` = {run: text}
    printed under the row title (the grid verdict, v38); `dpi` sets the file's
    resolution (150 default; the GUI's 'Redraw at 2x' uses 260)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    movies = list(per_movie)
    if not movies:
        return None
    k1 = grid_keys[0]; k2 = grid_keys[1] if len(grid_keys) > 1 else None
    fig, axes = plt.subplots(len(movies), 4, figsize=(18, 3.9 * len(movies)), squeeze=False)
    cols = ["#3b7dd8", "#e8702a", "#2ca02c", "#9467bd", "#8c564b", "#17becf"]
    # amplitude classes: light to dark by rank, whatever the class values (v38 — before, only 0.1/0.2/0.4 had colours)
    all_amps = sorted({float(c[len("recall_amp"):]) for rows_ in per_movie.values() for r in rows_ for c in r if c.startswith("recall_amp")})
    shades = ["#b0bec5", "#5b8def", "#1a237e", "#000000"]
    amp_cols = {a: shades[min(i, len(shades) - 1)] if len(all_amps) <= 3 else plt.cm.viridis(0.85 - 0.8 * i / max(len(all_amps) - 1, 1))
                for i, a in enumerate(all_amps)}
    for i, mv in enumerate(movies):
        rows = list(per_movie[mv])                     # the current setting joins its grid line
        cur = [r for r in per_movie[mv] if r["setting"] == 0]
        ch = chosen.get(mv)
        groups = sorted({r[k2] for r in rows}, key=float) if k2 else [None]
        ax_r, ax_a, ax_f, ax_n = axes[i]
        fk = false_key(rows)
        on_rates = any(np.isfinite(r.get("rate_real_on", np.nan)) for r in rows)
        rk = "rate_real_on" if on_rates else "rate_real"
        for gi, g in enumerate(groups):
            rr = sorted([r for r in rows if (k2 is None or r[k2] == g)], key=lambda r: float(r[k1]))
            x = [float(r[k1]) for r in rr]
            lab = f"{k2}={g}" if k2 else "grid"
            ax_r.plot(x, [r.get("recall", np.nan) for r in rr], "o-", color=cols[gi % 6], label=lab)
            ax_f.plot(x, [r.get(fk, np.nan) for r in rr], "o-", color=cols[gi % 6], label=lab)
            ax_n.plot(x, [r.get(rk, np.nan) for r in rr], "o-", color=cols[gi % 6], label=lab)
            if on_rates:
                ax_n.plot(x, [r.get("rate_off", np.nan) for r in rr], "v:", color=cols[gi % 6], ms=4, lw=0.9, alpha=0.8)
        # recall per amplitude, at the chosen setting's value of the second variable (else the first group)
        g_amp = ch[k2] if (ch and k2) else (groups[0] if k2 else None)
        rr = sorted([r for r in rows if (k2 is None or r[k2] == g_amp)], key=lambda r: float(r[k1]))
        x = [float(r[k1]) for r in rr]
        amps = sorted({float(c[len("recall_amp"):]) for r in rr for c in r if c.startswith("recall_amp")})
        amp_markers = ["o", "s", "^", "D"]
        for ai, a in enumerate(amps):
            tag = " (dimmest)" if ai == 0 and len(amps) > 1 else (" (brightest)" if ai == len(amps) - 1 and len(amps) > 1 else "")
            ax_a.plot(x, [r.get(f"recall_amp{a:g}", np.nan) for r in rr], amp_markers[ai % 4] + "-", color=amp_cols.get(a, "k"),
                      ms=7 if ai else 8, label=f"{a:g} dF/F planted{tag}")
        for kind, ls in (("short", ":"), ("long", "--")):
            if any(f"recall_{kind}" in r for r in rr):
                ax_a.plot(x, [r.get(f"recall_{kind}", np.nan) for r in rr], ls, color="0.3", lw=1, label=f"{kind} (all amplitudes)")
        ax_a.set_title(f"by amplitude and duration, {k2}={g_amp}" if k2 else "by amplitude and duration", fontsize=9, loc="right")
        ax2 = ax_n.twinx()
        for gi, g in enumerate(groups):
            rr = sorted([r for r in rows if (k2 is None or r[k2] == g)], key=lambda r: float(r[k1]))
            ax2.plot([float(r[k1]) for r in rr], [r.get("on_fibre", np.nan) for r in rr], "s--", color=cols[gi % 6], ms=3, lw=0.8, alpha=0.6)
        ax2.set_ylim(0, 1.05); ax2.set_ylabel("on-fibre fraction (squares)", fontsize=7, color="0.4"); ax2.tick_params(labelsize=7, colors="0.4")
        for ax in axes[i]:
            if ch:
                ax.axvline(float(ch[k1]), color="k", lw=1, alpha=0.5)
            for c in cur:
                ax.axvline(float(c[k1]), color="#d62728", lw=1, ls=":", alpha=0.7)
            ax.tick_params(labelsize=9); ax.spines[["top"]].set_visible(False); ax.grid(alpha=0.3)
            ax.set_xlabel(k1, fontsize=9)
        if ch:
            ax_r.plot([float(ch[k1])], [ch.get("recall", np.nan)], "*", color="k", ms=13, zorder=5)
            ax_f.plot([float(ch[k1])], [ch.get(fk, np.nan)], "*", color="k", ms=13, zorder=5)
        ax_r.set_ylim(0, 1.02); ax_a.set_ylim(0, 1.02); ax_f.set_ylim(0, 1.02)
        ax_r.set_ylabel("recall of injected transients (all)", fontsize=9)
        ax_a.set_ylabel("recall", fontsize=9)
        ax_f.set_ylabel("false-alarm share on the fibres\n(null events / real events, both on the fibre mask)" if fk == "false_share_on"
                        else "false-alarm share (null events / real events)", fontsize=9)
        ax_n.set_ylabel("real movie: events / min on the fibres (circles)\nand off the fibres (triangles, dotted)" if on_rates
                        else "real movie: events / min", fontsize=9)
        title = (titles.get(mv, mv) if titles else mv) + ("  — chosen: " + ch["label"] if ch else "  — no setting chosen")
        note = (notes or {}).get(mv, "")
        if note:
            title += "\n" + note
        ax_r.set_title(title, fontsize=9, loc="left", fontweight="bold")
        ax_r.legend(frameon=False, fontsize=8); ax_a.legend(frameon=False, fontsize=8)
    fig.suptitle("AQuA2 settings measured, one row per tuned movie.  Left: the share of the injected transients each setting recovers (recall), "
                 "one line per second grid value.  2nd: the same by injected amplitude class (light = dimmest) and duration.\n"
                 "3rd: false alarms — null-movie events on the fibres as a share of the real events (0 = clean; 1 = as many invented as real).  "
                 "Right: events/min the setting reports in the real movie on the fibres (circles), off them (triangles), on-fibre fraction (squares, right axis).\n"
                 "Red dotted line = the setting in the batch now; black line and star = the chosen one (highest recall with false alarms <= 10 % and >= 80 % on fibres).",
                 fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return out_png


def _why_panel(ax, bg: np.ndarray, fmask_raw: np.ndarray, vname: str, det, off, truth, reference=None, ref_off=None) -> None:
    """One picture of the 'why' figures: the mean image with what one setting found in one variant movie —
    real: events as circles of their size (green on the fibres, orange off); null: every detection red (a false
    alarm); inject: the injected transients found (filled cyan) or missed (hollow magenta). With a hand-made
    result (reference, v37) its events are drawn on the real panel as thin white dashed circles."""
    from matplotlib.patches import Circle
    H, W = bg.shape
    ax.imshow(bg, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    try:
        ax.contour(fmask_raw.astype(float), levels=[0.5], colors=["#00d0d0"], linewidths=0.3, alpha=0.3)
    except Exception:
        pass
    if vname == "real" and reference is not None and len(reference):
        ro = ref_off or off
        ry = reference["y_px"].to_numpy(float) - 1 + ro[0]; rx = reference["x_px"].to_numpy(float) - 1 + ro[1]
        rr = np.sqrt(np.clip(reference["n_pixels"].to_numpy(float), 1, None) / np.pi) if "n_pixels" in reference else np.full(len(reference), 6.0)
        for y, x, r_ in zip(ry, rx, rr):
            ax.add_patch(Circle((x, y), max(r_, 3.0) + 1.5, fill=False, lw=0.6, color="white", alpha=0.8, linestyle=(0, (2, 2))))
    if vname in ("real", "null") and det is not None and len(det):
        yy = det["y_px"].to_numpy(float) - 1 + off[0]; xx = det["x_px"].to_numpy(float) - 1 + off[1]
        rad = np.sqrt(np.clip(det["n_pixels"].to_numpy(float), 1, None) / np.pi) if "n_pixels" in det else np.full(len(det), 6.0)
        on = _on_mask(det, fmask_raw, off[:2]) if vname == "real" else np.zeros(len(det), bool)
        for y, x, rr_, o in zip(yy, xx, rad, on):
            col = "#d62728" if vname == "null" else ("#2ecc40" if o else "#ff9f1a")
            ax.add_patch(Circle((x, y), max(rr_, 3.0), fill=False, lw=0.9, color=col, alpha=0.9))
    if vname == "inject" and truth is not None and len(truth):
        tr = _truth_in_crop(truth, off[:2], off[2])
        hit = match_truth(tr, det) if det is not None else np.zeros(len(tr), bool)
        ty = tr["y_px"].to_numpy(float) - 1 + off[0]; tx = tr["x_px"].to_numpy(float) - 1 + off[1]
        for y, x, h in zip(ty, tx, hit):
            if h:
                ax.add_patch(Circle((x, y), 7.0, fill=True, lw=0, color="#00e5ff", alpha=0.8))
            else:
                ax.add_patch(Circle((x, y), 7.0, fill=False, lw=1.4, color="#ff3fd0", alpha=1.0))
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_xticks([]); ax.set_yticks([])


def _bg_image(raw_mean: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(raw_mean, [1, 99.7])
    return np.clip((raw_mean - lo) / max(hi - lo, 1e-6), 0, 1) ** 0.6


def why_day_figure(day_title: str, movies: list[dict], grid_keys: list[str], matched: dict | None, common: tuple | None,
                   max_false: float, out_png: Path, max_cols: int = 6):
    """The pre and post movies of one day side by side (v36): one column per setting (the first grid variable's
    values at the matched / common second value), rows = pre real, pre null, post real, post null (and the injected
    rows when there are at most 4 columns). The matched-sensitivity choice is framed green in its own movie's rows
    (pre may get one column, post another); the common setting is framed with a dashed grey line; a red tag marks a
    setting that fails the false-alarm rule on that movie. `movies` = [{cond, rows, raw_mean, fmask_raw, offsets,
    truth, read_events}] in pre/post order; `matched` = {cond: setting tuple} or None; `common` = setting tuple."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if len(movies) < 2 or not all(m.get("rows") for m in movies):
        return None
    k1 = grid_keys[0]; k2 = grid_keys[1] if len(grid_keys) > 1 else None
    fk = false_key([x for m in movies for x in m["rows"]])
    # the second grid value to show: the matched/common choice's (first movie), else the most common value
    if k2:
        if matched:
            g = next(iter(matched.values()))[1]
        elif common:
            g = common[1]
        else:
            g = sorted({r[k2] for m in movies for r in m["rows"]}, key=float)[0]
    # the columns: every value of the first grid variable scored at that second value (the current setting included
    # when it lies on the grid — a grid point and the current setting with the same values are the same column)
    cols = sorted({r[k1] for m in movies for r in m["rows"] if (k2 is None or r[k2] == g)}, key=float)
    must = {tuple(st)[0] for st in ((matched or {}).values())} | ({common[0]} if common else set())
    if len(cols) > max_cols:
        keep = {i for i, c in enumerate(cols) if c in must}
        for i in np.linspace(0, len(cols) - 1, max_cols).round().astype(int):
            if len(keep) < max_cols:
                keep.add(int(i))
        cols = [c for i, c in enumerate(cols) if i in keep]
    if not cols:
        return None
    show_inject = len(cols) <= 4
    var_rows = [("real", "real"), ("null", "null")] + ([("inject", "injected")] if show_inject else [])
    nrow = len(movies) * len(var_rows)
    n = len(cols)
    fig_w = min(3.4 * n + 0.8, 24)
    head_in = 0.85                                   # inches reserved for the four title lines
    fig_h = nrow * fig_w / n * 0.66 + head_in + 0.4
    fig, axes = plt.subplots(nrow, n, figsize=(fig_w, fig_h), squeeze=False)
    for mi, m in enumerate(movies):
        bg = _bg_image(m["raw_mean"])
        by_val: dict = {}
        for r in sorted(m["rows"], key=lambda r: (r["setting"] == 0, r["setting"])):   # a grid point before its 'current' twin
            if (k2 is None or r[k2] == g) and r[k1] not in by_val:
                by_val[r[k1]] = r
        m_key = matched.get(m["cond"]) if matched else None
        for ci, val in enumerate(cols):
            r = by_val.get(val)
            for vi, (vname, vlabel) in enumerate(var_rows):
                ax = axes[mi * len(var_rows) + vi, ci]
                if r is None:
                    ax.imshow(bg, cmap="gray", vmin=0, vmax=1); ax.set_xticks([]); ax.set_yticks([])
                    ax.text(0.5, 0.5, "not scored", transform=ax.transAxes, ha="center", color="white", fontsize=8)
                    continue
                det = m["read_events"](r["setting"], vname)
                off = (m["offsets"].get(vname) or m["offsets"].get("real") or (0, 0, None))
                _why_panel(ax, bg, m["fmask_raw"], vname, det, off, m["truth"], reference=m.get("reference"), ref_off=m.get("ref_off"))
                key = tuple(r[k] for k in grid_keys)
                is_m = m_key is not None and key == tuple(m_key)
                is_c = common is not None and key == tuple(common)
                fails = bool(np.isfinite(r.get(fk, np.nan)) and not is_clean(r, fk, max_false))
                for sp in ax.spines.values():
                    if is_m:
                        sp.set_color("#2ecc40"); sp.set_linewidth(3.2)
                    elif is_c:
                        sp.set_linestyle((0, (4, 3))); sp.set_color("0.25"); sp.set_linewidth(2.0)
                    else:
                        sp.set_color("#d62728" if fails else "0.6"); sp.set_linewidth(1.2 if fails else 0.8)
                if vname == "real":
                    txt = f"{int(r.get('n_real', 0))} events, {r.get('on_fibre', np.nan):.0%} on fibres"
                    # the verdict of this setting on this movie, top left
                    if is_m:
                        tag, tcol = f"{m['cond'].upper()}'S SETTING (matched)", "#2ecc40"
                    elif fails:
                        tag, tcol = f"fails: false alarms {r[fk]:.0%} > {max_false:.0%}", "#ff5252"
                    elif is_c:
                        tag, tcol = "common setting of the day", "#dddddd"
                    else:
                        tag, tcol = "passes the false-alarm rule", "#bbbbbb"
                    ax.text(0.01, 0.98, tag, transform=ax.transAxes, fontsize=7, color=tcol, va="top", ha="left",
                            fontweight="bold" if is_m else "normal",
                            bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="none", alpha=0.6))
                elif vname == "null":
                    txt = (f"{int(r.get('n_null_on', r.get('n_null', 0)) or 0)} false alarms = {r[fk]:.0%} of the real ones, {_fp_rate(r):.2f}/min"
                           if np.isfinite(r.get(fk, np.nan)) else "null not scored")
                else:
                    txt = f"recall {r.get('recall', np.nan):.2f}"
                ax.text(0.01, 0.02, txt, transform=ax.transAxes, fontsize=7, color="white", va="bottom", ha="left",
                        bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="none", alpha=0.6))
            if mi == 0:
                axes[0, ci].set_title(f"{k1} {val}" + (f"  {k2} {g}" if k2 else ""), fontsize=9)
        for vi, (_, vlabel) in enumerate(var_rows):
            axes[mi * len(var_rows) + vi, 0].set_ylabel(f"{m['cond'].upper()}\n{vlabel} movie", fontsize=8, labelpad=6)
    if matched:
        head = ("matched sensitivity: " + "; ".join(f"{c} {' '.join(f'{k}={v}' for k, v in zip(grid_keys, st))}" for c, st in matched.items())
                + "  (green frames)")
    else:
        head = "no matched pair of clean settings — the day is not comparable" + (
            f"; common setting {' '.join(f'{k}={v}' for k, v in zip(grid_keys, common))} (dashed)" if common else "")
    fig.suptitle(f"{day_title} — pre (top) and post (bottom) at the same settings, side by side.\n{head}\n"
                 "Rows per movie: real (events: green on fibres / orange off), null (every detection a false alarm, red)"
                 + (", injected (found cyan / missed magenta)" if show_inject else "") + ".\n"
                 "Green frame = that movie's own matched setting; dashed = the one common setting of the day; red = fails the false-alarm rule on that movie.",
                 fontsize=8.5, x=0.01, y=0.995, ha="left", va="top")
    fig.tight_layout(rect=(0, 0, 1, 1 - head_in / fig_h))
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return out_png


def why_figure(run_title: str, rows: list[dict], chosen: dict | None, passed: bool, grid_keys: list[str], raw_mean: np.ndarray,
               fmask_raw: np.ndarray, offsets: dict, truth, read_events, max_false: float, min_on_fibre: float, out_png: Path,
               max_cols: int = 6, reference=None, ref_off=None):
    """The sanity check behind the recommendation, as pictures instead of
    curves (v35): one column per setting (the thresholds at the chosen
    setting's second grid value, plus the current setting when it is not
    among them), three rows — the REAL movie with every detected event drawn
    as a circle of its size on the mean image (green on the fibres, orange off
    them), the NULL movie (time-scrambled copy: every detection is a false
    alarm, red) and the INJECTED movie (the 72 synthetic transients: found =
    filled cyan, missed = hollow magenta). The numbers under each column are
    the ones the rule uses; the chosen column is framed in green, columns
    that fail the false-alarm rule are tagged red. `read_events(setting,
    variant)` returns the exported events table (or None)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not rows:
        return None
    k1 = grid_keys[0]; k2 = grid_keys[1] if len(grid_keys) > 1 else None
    fk = false_key(rows)
    # the settings to show: the sweep of the first grid variable at the chosen (else the most common) value of the second
    if k2:
        g = chosen[k2] if chosen else sorted({r[k2] for r in rows}, key=float)[0]
        sel = [r for r in rows if r[k2] == g and r["setting"] != 0]
    else:
        sel = [r for r in rows if r["setting"] != 0]
    sel = sorted(sel, key=lambda r: float(r[k1]))
    cur = next((r for r in rows if r["setting"] == 0), None)
    if cur is not None and not any(tuple(r[k] for k in grid_keys) == tuple(cur[k] for k in grid_keys) for r in sel):
        sel.append(cur)
    if len(sel) > max_cols:                         # keep the extremes, the chosen and an even spread
        keep = {0, len(sel) - 1}
        if chosen:
            keep |= {i for i, r in enumerate(sel) if r is chosen or tuple(r[k] for k in grid_keys) == tuple(chosen[k] for k in grid_keys)}
        for i in np.linspace(0, len(sel) - 1, max_cols).round().astype(int):
            keep.add(int(i))
        sel = [r for i, r in enumerate(sel) if i in keep][:max_cols]
    if not sel:
        return None
    bg = _bg_image(raw_mean)
    chosen_key = tuple(chosen[k] for k in grid_keys) if chosen else None
    n = len(sel)
    fig_w = min(3.4 * n + 0.6, 24)
    fig, axes = plt.subplots(3, n, figsize=(fig_w, 3 * fig_w / n * 0.68 + 1.2), squeeze=False)
    variants = (("real", "REAL movie\nevents found\n(green on fibres, orange off" + (";\nwhite dashed = hand result)" if reference is not None else ")")),
                ("null", "NULL movie (scrambled)\nevery detection\n= a false alarm (red)"),
                ("inject", "INJECTED movie\nsynthetic transients\nfound (cyan) / missed (magenta)"))
    for ci, r in enumerate(sel):
        key = tuple(r[k] for k in grid_keys)
        is_chosen = chosen_key is not None and key == chosen_key
        is_cur = r["setting"] == 0 or (cur is not None and key == tuple(cur[k] for k in grid_keys))
        fails_false = np.isfinite(r.get(fk, np.nan)) and not is_clean(r, fk, max_false)
        for ri, (vname, _) in enumerate(variants):
            ax = axes[ri, ci]
            det = read_events(r["setting"], vname)
            off = offsets.get(vname) or offsets.get("real") or (0, 0, None)
            _why_panel(ax, bg, fmask_raw, vname, det, off, truth, reference=reference, ref_off=ref_off)
            frame_col = "#2ecc40" if is_chosen else ("#d62728" if fails_false else "0.5")
            for sp in ax.spines.values():
                sp.set_color(frame_col); sp.set_linewidth(3.0 if is_chosen else (1.6 if fails_false else 0.8))
            # the numbers of this row
            if vname == "real":
                n_all = r.get("n_real", np.nan)
                txt = f"{int(n_all) if np.isfinite(n_all) else '-'} events, {r.get('on_fibre', np.nan):.0%} on fibres"
                if np.isfinite(r.get("ref_recall", np.nan)):
                    txt += f"\nhand result: {r['ref_recall']:.0%} of {int(r.get('ref_n', 0))} found" + (
                        f", {r['ref_extra']:.0%} extra" if np.isfinite(r.get("ref_extra", np.nan)) else "")
            elif vname == "null":
                n_null = r.get("n_null_on", r.get("n_null", np.nan))
                txt = (f"{int(n_null) if np.isfinite(n_null) else '-'} false alarms on fibres = {r[fk]:.0%} of the real ones"
                       if np.isfinite(r.get(fk, np.nan)) else "null movie not scored")
            else:
                nt = r.get("n_injected", np.nan)
                txt = (f"recall {r.get('recall', np.nan):.2f} of {int(nt) if np.isfinite(nt) else '-'} injected"
                       + (f"  (0.1 / 0.2 / 0.4 dF/F: {r.get('recall_amp0.1', np.nan):.2f} / {r.get('recall_amp0.2', np.nan):.2f} / "
                          f"{r.get('recall_amp0.4', np.nan):.2f})" if np.isfinite(r.get('recall_amp0.1', np.nan)) else ""))
            ax.text(0.01, 0.02, txt, transform=ax.transAxes, fontsize=7, color="white", va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="none", alpha=0.6))
        # column header
        head = "  ".join(f"{k} {r[k]}" for k in grid_keys)
        if is_chosen:
            tag = "CHOSEN" if passed else "best trade-off (rule not met)"; tcol = "#1a7f2e" if passed else "#b8860b"
        elif fails_false:
            tag = f"too permissive: false alarms {r[fk]:.0%} > {max_false:.0%}"; tcol = "#b30000"
        elif np.isfinite(r.get("on_fibre", np.nan)) and r["on_fibre"] < min_on_fibre and chosen is not None:
            tag = f"off the fibres: only {r['on_fibre']:.0%} on fibres"; tcol = "#b36b00"
        elif chosen is not None and np.isfinite(r.get("recall", np.nan)) and np.isfinite(chosen.get("recall", np.nan)) and r["recall"] < chosen["recall"] - 1e-9:
            tag = f"too strict: recall {r.get('recall', np.nan):.2f} < {chosen.get('recall', np.nan):.2f}"; tcol = "0.35"
        else:
            tag = "passes the rule" if chosen is not None else ""; tcol = "0.35"
        if is_cur:
            tag = (tag + "  — the current setting") if tag else "the current setting"
        axes[0, ci].set_title(head + "\n" + tag, fontsize=8.5, color=tcol, fontweight="bold" if is_chosen else "normal")
    for ri, (_, label) in enumerate(variants):
        axes[ri, 0].set_ylabel(label, fontsize=7.5, labelpad=6)
    rule = (f"rule: the most sensitive setting (highest recall of the injected transients) whose false alarms on the fibres stay <= {max_false:.0%} "
            f"of the real events and whose events are >= {min_on_fibre:.0%} on the fibres.\n"
            "Lower thresholds add red (false alarms) faster than green (real events); higher ones turn cyan into magenta (missed transients).")
    fig.suptitle(f"{run_title} — why this setting?  Each column = one setting on the same movie.\n{rule}", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return out_png


# ----------------------------------------------------------------------------- opt file writer

def write_opt_csv(path: Path, values: dict, template_rows: list, note: str = "", skip_unchanged: bool = True) -> tuple[Path, bool]:
    """An AQuA2 parameter file (Name,Variable,Type,Value,Notes) the batch reads.
    Returns (path, written). When the file exists with the same values it is
    left untouched and written is False: its mtime decides which runs the batch
    queues again, and a recomputed common setting that did not change must not
    send every run through AQuA2 once more (adding a day to an animal)."""
    lines = ["Name,Variable,Type,Value,Notes"]
    new_vals = {}
    for name, var, typ, default in template_rows:
        if not var:
            lines.append(",,,,")
            continue
        v = str(values.get(var, default))
        new_vals[var] = v
        n = note.replace('"', "'") if var == "thrARScl" else ""
        lines.append(",".join([f'"{name}"' if "," in name else name, var, typ, v, f'"{n}"' if "," in n else n]))
    path = Path(path)
    if skip_unchanged and path.exists():
        try:
            from .aqua_io import read_opt_csv
            old = {k: str(v) for k, v in read_opt_csv(path).items()}
            if all(old.get(k) == v for k, v in new_vals.items()):
                return path, False
        except Exception:
            pass
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path, True
