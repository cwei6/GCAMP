#!/usr/bin/env python3
"""cfu_curate — look through the CFUs AQuA2 found, mark each one real or noise,
and add CFUs it missed. Decisions are read back by aqua_analyze.py.

    python cfu_curate.py <mouse folder>            # GUI
    python cfu_curate.py <mouse folder> --list      # what is there, per run
    python cfu_curate.py <run tiff> --reject 3,7    # no GUI: mark CFU 3 and 7 as noise
    python cfu_curate.py <run tiff> --keep-all
    python cfu_curate.py <mouse folder> --propagate # copy decisions to matched CFUs
                                                    # in the other runs of each FOV

Nothing in AQuA2's own files is modified. Each run gets a sidecar

    <run>_bin31f_1p98s_AQuA2_res_cfu_curated.json
        {"rejected": [3, 7], "added": [ {...} ], "notes": "..."}

so re-running the CFU batch never loses the curation, and `aqua_analyze.py`
drops rejected CFUs (and, by default, their events) and includes added ones.

GUI
  * left: the runs ("kept / total" per run) and, below them, the CFUs of the
    current run (number, events, kept / noise) — click one to select it; this
    is the easy way to pick one of several overlapping CFUs
  * top left: the mean image with the CFU outlines. Every CFU has its own
    colour, chosen so that overlapping or touching CFUs never share one; a CFU
    marked noise is thin grey, a hand-drawn region is cyan, the selected one is
    drawn thick in white over black. Click an outline to select it (clicking
    again where several overlap cycles through them).
  * top right: the MOVIE with the same outlines — by default AQuA2's own
    REGISTERED movie, read out of the run's _AQuA2.mat (res.datOrg1: the
    frames the events were detected on, after AQuA2's rigid registration,
    border crop and bleach correction, so the outlines sit exactly on it and
    the field does not wander); "as recorded" switches to the binned tiff as
    it was recorded (the crop offset applied, the motion not) — playable:
    space = play / pause, left / right = one frame back / forward (shift: 10),
    home / end, the slider, a speed box, raw or dF/F view (v), a contrast
    slider, and [ / ] = jump to the previous / next event of the selected CFU
    (3 frames before its onset, so one keypress on space shows the event).
    Click in the trace or in the raster to jump to that time.
  * zoom: scroll wheel over either image (around the cursor; both images move
    together), + / - keys, 0 = whole field; move the zoomed view: drag with
    the right (or middle) mouse button, or tick "move view" (h) and drag with the
    left button, or ctrl + arrow keys (a quarter of the view per press). The
    view is kept while you step through CFUs and through the runs of the same
    field.
  * middle: the selected CFU's dF/F curve with its event onsets marked and the
    movie's time cursor; right: its numbers
  * bottom: event onsets of every CFU (raster) with the time cursor
  * several CFUs at once (v31): drag a box with the left button over either
    image — every CFU whose centre (or half of whose footprint) lies in the box
    is selected (shift-drag adds to the selection, ctrl-click toggles one, Esc
    clears; ctrl/shift-click also works in the list); k / x then act on all of
    them. The title says how many are selected.
  * "looks" column and "Select suggested" (u): each CFU is checked against the
    fibre mask of the registered mean image (afferent2p/cfu_shape.py) —
    "fibre 85%" = that share of its footprint lies on fibres; "off-fibre 5%"
    = it sits beside the fibres; "blob 12%" = a compact bright spot off the
    fibres (debris, a cell body, autofluorescence). Rows that look non-fibre
    are tinted; u selects them all so one x marks them noise — a suggestion,
    not a decision: a bright varicosity ON a fibre stays "fibre", and a faint
    fibre the mask misses can be flagged wrongly, so look before you press x.
    f shows the fibre mask (orange) on the image. `--list --suggest` prints
    the suggestions per run without the GUI.
  * keys: k = keep, x = noise, a = add region (click vertices, Enter to close,
    Esc to cancel), tab / shift-tab = next / previous CFU, n / p = next /
    previous run, s = save
  * "Propagate to this FOV" copies the decisions of the current run to the CFUs
    of the other runs of the same FOV that overlap them (same matching as
    aqua_analyze.py: mean images aligned, then intersection-over-union).

AQuA2 registers the binned movie itself (one whole-pixel translation per
frame against the mean of its first 10 frames) and crops regMaskGap plus the
motion margin; the registered movie is what it detects on and what the panel
shows by default (needs h5py and a -v7.3 _AQuA2.mat; older results without
res.datOrg1 fall back to the tiff as recorded, with the crop offset). The dF/F
view is (F - F0) / F0 per pixel with F0 the median of ~60 frames spread over
the run.

Added CFUs are measured from the movie shown (the registered one when it is
available — same frame as the maps — else the binned tiff), not by AQuA2: the
trace is the mean over the drawn region, dF/F against a moving 10th-percentile
baseline, and events are threshold crossings (z > add_z, default 3, on the
MAD-scaled trace).
They are flagged as `source: manual` everywhere, so they never silently mix
with AQuA2's own event statistics.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

try:
    import tifffile  # noqa: F401
    import scipy  # noqa: F401
except ImportError as e:
    raise SystemExit(f"\nMissing Python package: {e.name}\n"
                     f'Install it: "{sys.executable}" -m pip install tifffile scipy h5py pandas\n')

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.aqua_io import (aqua_files, pick_binned_tif, read_cfu, relocate,      # noqa: E402
                                read_curation, write_curation, curation_path)
from afferent2p.plan import plan_mouse                                      # noqa: E402


# ----------------------------------------------------------------------------- loading

def discover_runs(target: Path, sec: float = 2.0, tol: float = 0.25, log=print):
    """[{experiment, fov, condition, run, tif}] for a mouse folder (via the
    aqua_batch manifest, else by planning) or for a single tiff."""
    target = Path(target)
    if target.is_file():
        return [{"experiment": "(single file)", "fov": None, "condition": "?",
                 "run": target.stem, "tif": target}]
    man = target / "aqua_batch" / "manifest.json"
    out = []
    if man.exists():
        m = json.loads(man.read_text())
        for e in m["experiments"]:
            for r in e.get("runs", []):
                out.append({"experiment": e["name"], "fov": e.get("fov"),
                            "condition": r["condition"], "run": r["run"],
                            "tif": relocate(r["tif"], target)})
        if out:
            return out
    for e in plan_mouse(target):
        for sbx_name, cond in e.runs:
            hits = list(target.rglob(sbx_name))
            if not hits:
                continue
            tif = pick_binned_tif(hits[0].parent, sec, tol=tol, log=log)
            if tif:
                out.append({"experiment": e.name, "fov": e.fov, "condition": cond,
                            "run": tif.stem, "tif": tif})
    return out


def load_run(tif: Path, log=print) -> dict | None:
    """CFUs + curation + display image for one run."""
    f = aqua_files(tif)
    if not f["cfu"].exists():
        return None
    c = read_cfu(f["cfu"])
    img = c["datPro"]
    if img is None and f["meanproj"].exists():
        import tifffile as tf
        img = np.asarray(tf.imread(str(f["meanproj"])), np.float32)
        img = (img - img.min()) / max(float(np.ptp(img)), 1e-9)
    spf = 2.0
    if f["json"].exists():
        try:
            spf = float(json.loads(f["json"].read_text())["seconds_per_frame"])
        except Exception:
            pass
    cur = read_curation(tif) or {"rejected": [], "added": []}
    # v43: decisions made for an earlier CFU result (e.g. another overlap threshold) come back emptied
    # with stale=True — the run starts clean; the old file is kept as *_superseded_* on the first save
    stale = bool(cur.get("stale"))
    n_ev = [int(np.asarray(e).size) for e in c["events"]]
    # AQuA2 crops regMaskGap pixels off every border (burst.prep1), so CFU maps
    # and datPro are smaller than the input tiff; a region drawn on the map must
    # be shifted by that gap to read the right pixels out of the tiff
    gap = 0
    ref = c["maps"].shape[1:] if c["maps"].size else (img.shape if img is not None else None)
    if ref is not None:
        try:
            import tifffile as tf
            with tf.TiffFile(str(tif)) as t:
                th, tw = t.pages[0].shape[:2]
            gap = int(round((th - ref[0]) / 2))
            if gap != int(round((tw - ref[1]) / 2)) or gap < 0:
                log(f"  note: {tif.name}: tiff {th}x{tw} vs CFU maps {ref[0]}x{ref[1]} — AQuA2's crop is not "
                    "symmetric; the curation window measures the offset from the movie")
                gap = max(0, gap)
        except Exception:
            gap = 0
    mean_raw = None
    mi = Path(str(tif.with_suffix("")) + "_AQuA2_meanimg.tif")
    if mi.exists():
        try:
            import tifffile as tf
            mean_raw = np.asarray(tf.imread(str(mi)), np.float32)
        except Exception:
            mean_raw = None
    return {"tif": tif, "cfu": c, "image": img, "mean_raw": mean_raw, "seconds_per_frame": spf, "crop_gap": gap,
            "rejected": set(int(i) for i in cur["rejected"]), "added": list(cur["added"]),
            "n_events": n_ev, "notes": "" if stale else cur.get("notes", ""),
            "stale": stale, "stale_reason": cur.get("stale_reason", "") if stale else "",
            "curated": curation_path(tif).exists() and not stale}


# ----------------------------------------------------------------------------- tuning-sweep events (overlay)

TUNE_COLOUR = "#ff4fd8"


def tune_events_dir(mouse_dir: Path, run: str, setting: int, variant: str = "real") -> Path | None:
    d = Path(mouse_dir) / "aqua_tune" / run / f"s{int(setting):02d}" / variant
    return d if d.is_dir() else None


def tune_events_records(mouse_dir: Path, run: str, setting: int, variant: str = "real", movie=None,
                        spf: float = 2.0, log=print) -> tuple[list, np.ndarray | None, str]:
    """The events AQuA2 found for one setting of the tuning sweep (tab 3), as
    added-style records the curation window can show: one box per event,
    its onset = AQuA2's t0, its trace read from the movie in a single pass.
    Returns (records, exported mean image or None, label of the setting)."""
    import pandas as pd
    d = tune_events_dir(mouse_dir, run, setting, variant)
    if d is None:
        return [], None, ""
    hits = sorted(d.glob("*_AQuA2_events.csv"))
    if not hits:
        return [], None, ""
    ev = pd.read_csv(hits[0])
    label = f"s{int(setting):02d}"
    try:
        man = json.loads((Path(mouse_dir) / "aqua_tune" / "tune_manifest.json").read_text())
        for mv in man.get("movies", []):
            if mv["run"] == run:
                for st in mv["settings"]:
                    if int(st["setting"]) == int(setting):
                        label = f"s{int(setting):02d} ({st['label']})"
    except Exception:
        pass
    mean = None
    mi = sorted(d.glob("*_meanimg.tif"))
    if mi:
        try:
            import tifffile as tf
            mean = np.asarray(tf.imread(str(mi[0])), np.float32)
        except Exception:
            mean = None
    recs = []
    boxes = []
    for r in ev.itertuples(index=False):
        x0, x1, y0, y1 = float(r.x0) - 1, float(r.x1) - 1, float(r.y0) - 1, float(r.y1) - 1
        poly = [[x0 - 0.5, y0 - 0.5], [x1 + 0.5, y0 - 0.5], [x1 + 0.5, y1 + 0.5], [x0 - 0.5, y1 + 0.5]]
        boxes.append((int(y0), int(y1) + 1, int(x0), int(x1) + 1))
        recs.append({"source": "tune", "label": f"T{int(r.event)} {label}",
                     "note": (f"\nAQuA2 event {int(r.event)} of setting {label}: frames {int(r.t0_frame)}-{int(r.t1_frame)}, "
                              f"dF/F max {float(r.dff_max):.2f}, {int(r.n_pixels)} px (box shown)"),
                     "polygon_xy": poly, "npix": int(r.n_pixels), "crop_gap": [0, 0],
                     "n_events": 1, "onset_frames": [int(r.t0_frame)], "t1_frame": int(r.t1_frame),
                     "dff": [], "z_thresh": float("nan"), "peak_dff": float(r.dff_max), "seconds_per_frame": float(spf)})
    # traces: one pass over the movie, every box at once
    if movie is not None and recs:
        try:
            T = movie.T
            acc = np.zeros((len(boxes), T), np.float32)
            for t in range(T):
                fr = movie.frame(t, "raw")
                for k, (ya, yb, xa, xb) in enumerate(boxes):
                    acc[k, t] = float(fr[max(0, ya):yb, max(0, xa):xb].mean())
            for k, rec in enumerate(recs):
                dff, _, _ = trace_to_events(acc[k], spf)
                rec["dff"] = [round(float(v), 4) for v in dff]
        except Exception as e:
            log(f"  note: traces of the tuning events not computed ({e})")
    log(f"  {run}: {len(recs)} events of tuning setting {label} overlaid (magenta boxes T1..; not saved with the curation)")
    return recs, mean, label


# ----------------------------------------------------------------------------- hand-drawn CFUs

def polygon_mask(points, shape) -> np.ndarray:
    """Boolean mask of a polygon given as [(x, y), ...] (matplotlib order)."""
    from matplotlib.path import Path as MplPath
    H, W = shape
    yy, xx = np.mgrid[0:H, 0:W]
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    inside = MplPath(np.asarray(points, float)).contains_points(pts)
    return inside.reshape(H, W)


def split_gap(gap) -> tuple[int, int]:
    """AQuA2's border crop as (rows, columns) — an int means the same on both axes."""
    if gap is None:
        return 0, 0
    if np.isscalar(gap):
        return int(gap), int(gap)
    gy, gx = gap
    return int(gy), int(gx)


def region_trace(tif: Path, mask: np.ndarray, gap=0, movie=None) -> np.ndarray:
    """Mean intensity inside a mask, frame by frame, read page by page (the
    binned movie is small but this keeps memory flat for any length).

    `gap` shifts the mask from AQuA2's cropped coordinates into the tiff's.
    With `movie` (a frame source with .shape (T, H, W) and movie[i], e.g.
    AQuA2's registered movie) the pixels are read from it instead and `gap`
    is ignored (its frame IS the CFU maps' frame)."""
    import tifffile as tf
    if movie is not None:
        T, th, tw = (int(v) for v in movie.shape)
        full = np.zeros((th, tw), bool)
        h = min(mask.shape[0], th); w = min(mask.shape[1], tw)
        full[:h, :w] = mask[:h, :w]
        idx = np.flatnonzero(full.ravel())
        return np.asarray([float(np.asarray(movie[i]).ravel()[idx].mean()) for i in range(T)], np.float32)
    with tf.TiffFile(str(tif)) as t:
        th, tw = t.pages[0].shape[:2]
    gy, gx = split_gap(gap)
    if gy or gx or mask.shape != (th, tw):
        full = np.zeros((th, tw), bool)
        h = min(mask.shape[0], th - gy)
        w = min(mask.shape[1], tw - gx)
        full[gy:gy + h, gx:gx + w] = mask[:h, :w]
        mask = full
    idx = np.flatnonzero(mask.ravel())
    out = []
    with tf.TiffFile(str(tif)) as t:
        for page in t.pages:
            out.append(float(np.asarray(page.asarray()).ravel()[idx].mean()))
    return np.asarray(out, np.float32)


def trace_to_events(curve: np.ndarray, spf: float, z_thresh: float = 3.0,
                    baseline_frames: int = 25, baseline_pct: float = 10.0):
    """dF/F + event onsets for a hand-drawn region.

    Baseline: moving percentile (a plain moving mean, as AQuA2 uses, is dragged
    up by the events themselves). Events: z (MAD-scaled) crossing above
    z_thresh; onset = the first frame of each crossing.
    """
    from scipy.ndimage import percentile_filter
    win = max(3, int(baseline_frames) | 1)
    f0 = percentile_filter(curve.astype(np.float32), baseline_pct, size=win, mode="nearest")
    dff = (curve - f0) / np.maximum(f0, 1e-6)
    med = float(np.median(dff))
    mad = float(np.median(np.abs(dff - med))) * 1.4826 + 1e-9
    z = (dff - med) / mad
    above = z > z_thresh
    onsets = np.flatnonzero(np.diff(np.r_[0, above.astype(int)]) == 1)
    return dff, z, onsets


def add_region(tif: Path, points, shape, spf: float, z_thresh: float = 3.0,
               label: str | None = None, gap=0, movie=None) -> dict:
    """Measure a hand-drawn region and return the record stored in the sidecar.

    `points` are in the coordinates of the displayed image (AQuA2's cropped
    frame); `gap` is the border AQuA2 removed, so the pixels are read from the
    right place in the tiff — unless `movie` is AQuA2's registered movie, which
    is read directly (same frame as the maps, and steadied)."""
    mask = polygon_mask(points, shape)
    if mask.sum() == 0:
        raise ValueError("the drawn region contains no pixels")
    curve = region_trace(tif, mask, gap=gap, movie=movie)
    dff, z, onsets = trace_to_events(curve, spf, z_thresh)
    return {
        "source": "manual", "label": label or "added",
        "polygon_xy": [[float(x), float(y)] for x, y in points],
        "npix": int(mask.sum()), "crop_gap": list(split_gap(gap)), "movie": "registered" if movie is not None else "raw",
        "n_events": int(onsets.size),
        "onset_frames": [int(i) + 1 for i in onsets],          # 1-based, like AQuA2
        "dff": [round(float(v), 4) for v in dff],
        "z_thresh": float(z_thresh),
        "peak_dff": round(float(np.nanmax(dff)), 4),
        "seconds_per_frame": float(spf),
    }


def added_mask(rec: dict, shape) -> np.ndarray:
    return polygon_mask(rec["polygon_xy"], shape)


# ----------------------------------------------------------------------------- propagation

def propagate(runs: list[dict], iou: float = 0.3, log=print) -> dict:
    """Copy each run's decisions onto the overlapping CFUs of the other runs of
    the same FOV. A CFU identity counts as noise everywhere as soon as it was
    rejected in one run. Returns {tif: number of CFUs newly rejected}."""
    from aqua_analyze import match_cfus                                    # noqa: E402
    by_exp: dict[str, list[dict]] = {}
    for r in runs:
        by_exp.setdefault(r["experiment"], []).append(r)
    changed = {}
    for exp, group in by_exp.items():
        loaded = [(g, load_run(g["tif"], log=log)) for g in group]
        loaded = [(g, d) for g, d in loaded if d is not None]
        if len(loaded) < 2:
            continue
        cfu_list = [d["cfu"] for _, d in loaded]
        ids, n_matched, transforms, n_id = match_cfus(cfu_list, iou_thr=iou, log=log)
        bad_ids = set()
        for (_, d), idv in zip(loaded, ids):
            for k, ident in enumerate(idv):
                if (k + 1) in d["rejected"]:
                    bad_ids.add(int(ident))
        total_new = 0
        for (g, d), idv in zip(loaded, ids):
            new = {k + 1 for k, ident in enumerate(idv) if int(ident) in bad_ids}
            n_new = len(new - d["rejected"])
            total_new += n_new
            if n_new:
                write_curation(g["tif"], sorted(d["rejected"] | new), d["added"],
                               notes=(d["notes"] + " | propagated within " + exp).strip(" |"))
            changed[str(g["tif"])] = n_new
        log(f"  {exp}: {len(bad_ids)} CFU identity(ies) marked noise, "
            f"{total_new} new rejection(s) written across {len(loaded)} runs")
    return changed


# ----------------------------------------------------------------------------- outlines, colours, movie

# 19 colours that stay apart on a grey image (no grey, no white); a CFU gets the
# first one that none of the CFUs it overlaps or touches already has
PALETTE = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4",
           "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff", "#9a6324", "#fffac8",
           "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075"]
NOISE_COLOUR = "0.55"
ADDED_COLOUR = "cyan"


def cfu_masks(maps: np.ndarray, thr: float = 0.1) -> list[np.ndarray]:
    return [maps[k] > thr for k in range(maps.shape[0] if maps.size else 0)]


def mask_outline(mask: np.ndarray) -> list[np.ndarray]:
    """Closed outline(s) of a boolean mask as (N, 2) arrays of (x, y) in pixel
    coordinates (every CFU, however ragged, in one pass — far faster than one
    ax.contour per CFU)."""
    from skimage import measure
    if not mask.any():
        return []
    padded = np.pad(mask.astype(np.float32), 1)
    return [c[:, ::-1] - 1.0 for c in measure.find_contours(padded, 0.5)]


def polygon_outline(points) -> list[np.ndarray]:
    p = np.asarray(points, float)
    return [np.vstack([p, p[:1]])] if len(p) else []


def assign_colours(masks: list[np.ndarray], n_colours: int = len(PALETTE), touch_px: int = 3) -> list[int]:
    """Palette index per CFU such that CFUs that overlap or come within
    `touch_px` of each other get different colours (greedy graph colouring;
    only CFUs whose bounding boxes meet are compared)."""
    from scipy.ndimage import binary_dilation
    n = len(masks)
    boxes, grown = [], []
    for m in masks:
        if m.any():
            ys, xs = np.nonzero(m)
            boxes.append((ys.min() - touch_px, ys.max() + touch_px, xs.min() - touch_px, xs.max() + touch_px))
            grown.append(binary_dilation(m, iterations=touch_px))
        else:
            boxes.append(None); grown.append(m)
    colour = [-1] * n
    for k in range(n):
        used = set()
        bk = boxes[k]
        if bk is not None:
            for j in range(k):
                bj = boxes[j]
                if bj is None or bk[0] > bj[1] or bj[0] > bk[1] or bk[2] > bj[3] or bj[2] > bk[3]:
                    continue
                if (grown[k] & masks[j]).any():
                    used.add(colour[j])
        colour[k] = next((c for c in range(n_colours) if c not in used), k % n_colours)
    return colour


class Movie:
    """Frame access for the movie panel, raw or dF/F, in the CFU maps' frame:
    source "registered" = AQuA2's own registered movie read out of the
    _AQuA2.mat (res.datOrg1: the frames the events were detected on, already
    cropped like the maps — no offset, no motion beyond what AQuA2 left);
    source "raw" = the binned tiff as recorded, lazily (memmap or page by
    page), cropped to the maps' frame at the crop offset."""

    def __init__(self, tif: Path, gap: int = 0, shape=None, source: str = "registered", log=print):
        import tifffile as tf
        self.tif = Path(tif)
        self.pages = None
        self.arr = None
        self.source = "raw"
        self._reg = None
        if source == "registered":
            from afferent2p.aqua_io import open_registered
            reg = open_registered(self.tif, log=log)
            if reg is not None:
                want = tuple(int(v) for v in shape) if shape is not None else None
                if want is not None and tuple(reg.shape[1:]) != want:
                    log(f"  note: {self.tif.name}: the registered movie is {reg.shape[1]}x{reg.shape[2]} but the CFU maps are "
                        f"{want[0]}x{want[1]} (result files from different AQuA2 runs?) — showing the binned tiff instead")
                    reg.close()
                else:
                    self._reg = reg; self.arr = reg; self.source = "registered"
        if self.arr is None:
            try:
                self.arr = tf.memmap(str(self.tif), mode="r")
                if self.arr.ndim == 2:
                    self.arr = self.arr[None]
            except Exception:
                self.arr = None
                self._tf = tf.TiffFile(str(self.tif))
                self.pages = self._tf.pages
        self.T = int(self.arr.shape[0] if self.arr is not None else len(self.pages))
        th, tw = (self.arr.shape[1:] if self.arr is not None else self.pages[0].shape[:2])
        self.full_shape = (int(th), int(tw))
        self.gap = split_gap(gap) if self.source == "raw" else (0, 0)
        self.shape = tuple(int(v) for v in shape) if shape is not None else self.full_shape
        self.sample_mean = None
        self.f0 = None
        self.vmax_dff = 0.5
        sample = self._sample(10)
        self.vmin_raw = float(np.percentile(sample, 1))
        self.vmax_raw = float(np.percentile(sample, 99.7))
        if self.vmax_raw <= self.vmin_raw:
            self.vmax_raw = self.vmin_raw + 1.0

    def _read(self, i: int) -> np.ndarray:
        if self.arr is not None:
            return np.asarray(self.arr[i])
        return np.asarray(self.pages[i].asarray())

    def crop(self, frame: np.ndarray) -> np.ndarray:
        """The tiff frame in the CFU maps' coordinates (AQuA2's border crop)."""
        H, W = self.shape
        gy, gx = self.gap
        out = np.zeros((H, W), np.float32)
        h = min(H, frame.shape[0] - gy); w = min(W, frame.shape[1] - gx)
        if h > 0 and w > 0:
            out[:h, :w] = frame[gy:gy + h, gx:gx + w]
        return out

    def _sample(self, n: int, full: bool = False) -> np.ndarray:
        idx = np.unique(np.linspace(0, self.T - 1, min(n, self.T)).astype(int))
        frames = [np.asarray(self._read(int(i)), np.float32) for i in idx]
        if full:
            return np.stack(frames)
        return np.stack([self.crop(f) for f in frames])

    def align_to(self, mean_img, min_match: float = 0.5):
        """Find where AQuA2's (registered, cropped) mean image sits in the raw
        frame by normalised cross-correlation and set the crop offset to it;
        returns (dy, dx, match) or None when the match is poor. AQuA2 crops
        regMaskGap plus the motion margin, which need not be symmetric.
        Nothing to do for the registered movie (offset 0)."""
        if mean_img is None or self.source == "registered":
            return None
        try:
            from skimage.feature import match_template
            ref = np.mean(self._sample(10, full=True), axis=0)
            tmpl = np.asarray(mean_img, np.float32)
            if tmpl.shape[0] > ref.shape[0] or tmpl.shape[1] > ref.shape[1]:
                return None
            cc = match_template(ref, tmpl)
            dy, dx = np.unravel_index(int(np.argmax(cc)), cc.shape)
            score = float(cc[dy, dx])
            if score < min_match:
                return None
            self.gap = (int(dy), int(dx))
            return int(dy), int(dx), score
        except Exception:
            return None

    def ensure_f0(self, n: int = 60):
        if self.f0 is None:
            s = self._sample(n)
            self.f0 = np.maximum(np.median(s, axis=0), 1.0).astype(np.float32)
            d = (s - self.f0) / self.f0
            self.vmax_dff = max(float(np.percentile(d, 99.9)), 0.2)

    def frame(self, i: int, mode: str = "raw") -> np.ndarray:
        i = int(min(max(i, 0), self.T - 1))
        f = self.crop(self._read(i))
        if mode == "dff":
            self.ensure_f0()
            return (f - self.f0) / self.f0
        return f

    def clim(self, mode: str, contrast: float = 1.0):
        if mode == "dff":
            self.ensure_f0()
            return (-0.05, self.vmax_dff * contrast)
        return (self.vmin_raw, self.vmin_raw + (self.vmax_raw - self.vmin_raw) * contrast)

    def close(self):
        if self.pages is not None:
            try:
                self._tf.close()
            except Exception:
                pass
        if self._reg is not None:
            try:
                self._reg.close()
            except Exception:
                pass
        self.arr = None; self.pages = None; self._reg = None


# ----------------------------------------------------------------------------- GUI

def launch_gui(runs: list[dict], log=print, fps: float = 8.0, tune: int | None = None, tune_variant: str = "real",
               mouse_dir: Path | None = None, um_per_px: float = 0.792):
    import tkinter as tk
    from tkinter import ttk, messagebox
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.collections import LineCollection
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    SPEEDS = ("2", "4", "8", "15", "30")
    TYPING = ("Entry", "TEntry", "TCombobox", "TSpinbox", "Spinbox", "Text", "TScale", "Scale")

    def fixed_map(style, option):
        # Tk 8.6.9 ignores Treeview tag colours unless the style map is rewritten
        return [e for e in style.map("Treeview", query_opt=option) if e[:2] != ("!disabled", "!selected")]

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("cfu_curate — mark AQuA2 CFUs real or noise")
            self.geometry("1560x940")
            self.runs = runs
            self.i = 0
            self.data = None
            self.tune = tune           # tuning setting to overlay (tab 3), or None
            self.tune_variant = tune_variant
            self.mouse_dir = Path(mouse_dir) if mouse_dir else None
            self.sel = None            # (kind, index): ("cfu", k) or ("added", j) — the CFU whose trace is shown
            self.multi = set()         # several CFUs selected at once (drag a box, ctrl-click, "Select suggested"): k, k, ...
            self.adding = None         # list of clicked points while adding
            self.um_per_px = um_per_px
            self.feats = []            # per CFU: shape / position features and a suggestion (afferent2p.cfu_shape)
            self.fmask = None          # fibre mask of this run's registered mean image
            self.show_fibre = False
            self._drag = None          # left-button press: (x_px, y_px, xdata, ydata, shift, ctrl)
            self._band = None          # the rubber band (two Rectangle patches) while dragging
            self.movie = None
            self.frame = 0
            self.playing = False
            self.mode = "raw"
            self.view = None           # (xlim, ylim) kept across redraws / runs of the same shape
            self.view_shape = None
            self.masks, self.outlines, self.colours = [], [], []
            self.centroids, self.onsets = [], []
            self.event_note = ""
            self._bg = None
            self._anim = []
            self._pan = None
            self.hand = False              # h: left-drag moves the view instead of drawing a selection box
            self._slider_lock = False
            self._ctree_lock = False
            self.mov_im = self.mov_txt = self.cur_line = self.all_line = None
            self.mov_lines = self.mov_sel = None
            self._build()
            self.load(0)
            self.after(300, self.report_stale)

        def report_stale(self):
            """v43: say once which runs have decisions from an earlier CFU result (they are ignored)."""
            stale, lost, why = [], [], ""
            for r in self.runs:
                cu = read_curation(r["tif"])
                if cu and cu.get("stale"):
                    stale.append(r["run"].replace("_bin31f_1p98s", ""))
                    why = why or cu.get("stale_reason", "")
                    try:
                        raw = json.loads(curation_path(r["tif"]).read_text(encoding="utf-8"))
                        if raw.get("rejected") or raw.get("added"):
                            lost.append(stale[-1])
                    except Exception:
                        pass
            if not stale:
                return
            log(f"  {len(stale)} run(s) have curation files from an earlier CFU result ({why}); they are ignored"
                + (f" — {len(lost)} of them hold decisions: " + ", ".join(lost[:8]) if lost else ""))
            if lost:
                messagebox.showinfo("earlier decisions ignored",
                                    f"{len(lost)} run(s) have CFU decisions made for an earlier CFU result\n({why}).\n\n"
                                    "The CFU numbers changed with the result, so those decisions are NOT applied: these runs "
                                    "start clean. The old file is kept beside the new one as\n"
                                    "*_res_cfu_curated_superseded_<date>.json the first time you save the run.\n\n"
                                    + "\n".join(lost[:12]) + ("\n…" if len(lost) > 12 else ""))

        # ---- layout
        def _build(self):
            style = ttk.Style(self)
            try:
                style.map("Treeview", foreground=fixed_map(style, "foreground"),
                          background=fixed_map(style, "background"))
            except Exception:
                pass
            top = ttk.Frame(self); top.pack(fill="x", padx=6, pady=4)
            ttk.Button(top, text="< prev run (p)", command=lambda: self.step(-1)).pack(side="left")
            ttk.Button(top, text="next run (n) >", command=lambda: self.step(1)).pack(side="left", padx=4)
            ttk.Button(top, text="Keep (k)", command=lambda: self.mark(True)).pack(side="left", padx=(16, 2))
            ttk.Button(top, text="Noise (x)", command=lambda: self.mark(False)).pack(side="left", padx=2)
            ttk.Button(top, text="< CFU", width=6, command=lambda: self.step_cfu(-1)).pack(side="left", padx=(16, 1))
            ttk.Button(top, text="CFU >", width=6, command=lambda: self.step_cfu(1)).pack(side="left", padx=1)
            ttk.Button(top, text="Add region (a)", command=self.start_add).pack(side="left", padx=(16, 2))
            ttk.Button(top, text="Delete added", command=self.delete_added).pack(side="left", padx=2)
            ttk.Button(top, text="Keep all", command=lambda: self.mark_all(True)).pack(side="left", padx=(16, 2))
            ttk.Button(top, text="Reject all", command=lambda: self.mark_all(False)).pack(side="left", padx=2)
            ttk.Button(top, text="Select suggested (u)", command=self.select_suggested).pack(side="left", padx=(16, 2))
            self.fibre_var = tk.BooleanVar(value=False)
            ttk.Checkbutton(top, text="fibre mask (f)", variable=self.fibre_var, command=self.toggle_fibre).pack(side="left", padx=2)
            ttk.Button(top, text="Save (s)", command=self.save).pack(side="right")
            ttk.Button(top, text="Propagate to this FOV", command=self.propagate_fov).pack(side="right", padx=6)

            mv = ttk.Frame(self); mv.pack(fill="x", padx=6, pady=(0, 4))
            ttk.Label(mv, text="movie:").pack(side="left")
            ttk.Button(mv, text="|<", width=3, command=lambda: self.set_frame(0)).pack(side="left")
            ttk.Button(mv, text="< 1", width=4, command=lambda: self.set_frame(self.frame - 1)).pack(side="left")
            self.play_btn = ttk.Button(mv, text="play (space)", width=13, command=self.toggle_play)
            self.play_btn.pack(side="left", padx=2)
            ttk.Button(mv, text="1 >", width=4, command=lambda: self.set_frame(self.frame + 1)).pack(side="left")
            ttk.Button(mv, text=">|", width=3, command=lambda: self.set_frame(10 ** 9)).pack(side="left")
            ttk.Label(mv, text="  fps").pack(side="left")
            self.fps_var = tk.StringVar(value=f"{fps:g}")
            cb = ttk.Combobox(mv, textvariable=self.fps_var, values=SPEEDS, width=4, state="readonly")
            cb.pack(side="left"); cb.bind("<<ComboboxSelected>>", lambda _e: self.focus_canvas())
            self.mode_var = tk.StringVar(value="raw")
            ttk.Radiobutton(mv, text="raw", variable=self.mode_var, value="raw",
                            command=self.set_mode).pack(side="left", padx=(10, 0))
            ttk.Radiobutton(mv, text="dF/F (v)", variable=self.mode_var, value="dff",
                            command=self.set_mode).pack(side="left")
            ttk.Label(mv, text="  contrast").pack(side="left")
            self.contrast = tk.DoubleVar(value=1.0)
            sc = ttk.Scale(mv, from_=0.3, to=3.0, variable=self.contrast, length=90,
                           command=lambda _v: self.apply_clim())
            sc.pack(side="left"); sc.bind("<ButtonRelease-1>", lambda _e: self.focus_canvas())
            self.outl_var = tk.BooleanVar(value=True)
            ttk.Checkbutton(mv, text="outlines on movie", variable=self.outl_var,
                            command=self.show_frame).pack(side="left", padx=(10, 0))
            ttk.Button(mv, text="[ < event", width=9, command=lambda: self.jump_event(-1)).pack(side="left", padx=(10, 1))
            ttk.Button(mv, text="event > ]", width=9, command=lambda: self.jump_event(1)).pack(side="left", padx=1)
            ttk.Button(mv, text="whole field (0)", command=self.reset_view).pack(side="left", padx=(10, 0))
            sl = ttk.Frame(self); sl.pack(fill="x", padx=6, pady=(0, 4))
            ttk.Label(sl, text="show:").pack(side="left")
            self.src_var = tk.StringVar(value="registered")
            ttk.Radiobutton(sl, text="AQuA2-registered", variable=self.src_var, value="registered",
                            command=self.set_source).pack(side="left")
            ttk.Radiobutton(sl, text="as recorded", variable=self.src_var, value="raw",
                            command=self.set_source).pack(side="left")
            self.hand_var = tk.BooleanVar(value=False)
            ttk.Checkbutton(sl, text="move view (h)", variable=self.hand_var, command=self.toggle_hand).pack(side="left", padx=(10, 0))
            ttk.Label(sl, text="  frame:").pack(side="left")
            self.frame_lbl = ttk.Label(sl, text="", width=30, anchor="e"); self.frame_lbl.pack(side="right")
            self.slider = ttk.Scale(sl, from_=0, to=1, orient="horizontal", command=self.on_slider)
            self.slider.pack(side="left", fill="x", expand=True, padx=8)
            self.slider.bind("<ButtonRelease-1>", lambda _e: self.focus_canvas())

            body = ttk.Frame(self); body.pack(fill="both", expand=True, padx=6)
            left = ttk.Frame(body); left.pack(side="left", fill="y")
            lf = ttk.LabelFrame(left, text="runs"); lf.pack(fill="x")
            self.tree = ttk.Treeview(lf, columns=("kept",), height=9, show="tree headings")
            self.tree.heading("#0", text="run"); self.tree.heading("kept", text="kept/total")
            self.tree.column("#0", width=250); self.tree.column("kept", width=96, anchor="center")
            self.tree.pack(side="left", fill="x", expand=True)
            sb = ttk.Scrollbar(lf, orient="vertical", command=self.tree.yview); sb.pack(side="right", fill="y")
            self.tree.configure(yscrollcommand=sb.set)
            self.tree.bind("<<TreeviewSelect>>", self.on_tree)
            for k, r in enumerate(self.runs):
                name = r["run"].replace("_bin31f_1p98s", "")
                self.tree.insert("", "end", iid=str(k),
                                 text=f"[{r['condition']}] {name[-30:]}", values=("?",))
            cf = ttk.LabelFrame(left, text="CFUs of this run (click to select; ctrl/shift-click several)"); cf.pack(fill="both", expand=True, pady=(6, 0))
            self.ctree = ttk.Treeview(cf, columns=("ev", "looks", "st"), show="tree headings", selectmode="extended")
            self.ctree.heading("#0", text="CFU"); self.ctree.heading("ev", text="events"); self.ctree.heading("looks", text="looks")
            self.ctree.heading("st", text="status")
            self.ctree.column("#0", width=96); self.ctree.column("ev", width=48, anchor="center")
            self.ctree.column("looks", width=104, anchor="center"); self.ctree.column("st", width=52, anchor="center")
            self.ctree.pack(side="left", fill="both", expand=True)
            sb2 = ttk.Scrollbar(cf, orient="vertical", command=self.ctree.yview); sb2.pack(side="right", fill="y")
            self.ctree.configure(yscrollcommand=sb2.set)
            self.ctree.bind("<<TreeviewSelect>>", self.on_ctree)
            self.ctree.tag_configure("added", foreground="#0097a7")
            self.ctree.tag_configure("noise", foreground="#8a8a8a")
            self.ctree.tag_configure("sugg", background="#fff3d6")

            self.fig = Figure(figsize=(12, 8))
            self.ax_img = self.fig.add_axes([0.02, 0.44, 0.47, 0.53])
            self.ax_mov = self.fig.add_axes([0.51, 0.44, 0.47, 0.53], sharex=self.ax_img, sharey=self.ax_img)
            self.ax_cur = self.fig.add_axes([0.05, 0.255, 0.66, 0.145])
            self.ax_txt = self.fig.add_axes([0.74, 0.19, 0.24, 0.21]); self.ax_txt.axis("off")
            self.ax_all = self.fig.add_axes([0.05, 0.05, 0.66, 0.125])
            self.canvas = FigureCanvasTkAgg(self.fig, master=body)
            self.canvas.get_tk_widget().pack(side="left", fill="both", expand=True)
            self.canvas.mpl_connect("button_press_event", self.on_click)
            self.canvas.mpl_connect("button_release_event", self.on_release)
            self.canvas.mpl_connect("motion_notify_event", self.on_motion)
            self.canvas.mpl_connect("scroll_event", self.on_scroll)
            self.canvas.mpl_connect("draw_event", self.on_draw)
            self.bind("<Key>", self.on_key)
            self.bind("<Tab>", lambda e: self.on_key(e) or "break")
            self.bind("<Shift-Tab>", lambda e: self.on_key(e) or "break")
            try:
                self.bind("<ISO_Left_Tab>", lambda e: self.on_key(e) or "break")   # shift-tab on X11
            except tk.TclError:
                pass
            for bar in (top, mv, sl):                 # buttons never keep the keyboard focus
                for w in bar.winfo_children():
                    if isinstance(w, (ttk.Button, ttk.Checkbutton, ttk.Radiobutton)):
                        w.configure(takefocus=0)
            self.protocol("WM_DELETE_WINDOW", self.on_close)
            self.focus_canvas()

        def focus_canvas(self):
            try:
                self.canvas.get_tk_widget().focus_set()
            except Exception:
                pass

        def on_close(self):
            self.playing = False
            self.save(quiet=True)
            if self.movie is not None:
                self.movie.close()
            self.destroy()

        # ---- data
        def shape(self):
            d = self.data
            if d is None:
                return (1, 1)
            if d["cfu"]["maps"].size:
                return tuple(d["cfu"]["maps"].shape[1:])
            return tuple(d["image"].shape) if d["image"] is not None else (1, 1)

        def load(self, i):
            self.playing = False
            self.i = max(0, min(i, len(self.runs) - 1))
            r = self.runs[self.i]
            self.data = load_run(r["tif"], log=log)
            self.sel = None; self.adding = None; self.event_note = ""
            if self.movie is not None:
                self.movie.close(); self.movie = None
            self.masks, self.outlines, self.colours = [], [], []
            self.centroids, self.onsets = [], []
            tune_recs, tune_mean, tune_label = [], None, ""
            if self.tune is not None and self.mouse_dir is not None:
                tune_recs, tune_mean, tune_label = tune_events_records(self.mouse_dir, r["run"], self.tune, self.tune_variant,
                                                                       movie=None, spf=2.0, log=log)
                if not tune_recs and tune_events_dir(self.mouse_dir, r["run"], self.tune, self.tune_variant) is None:
                    log(f"  {r['run']}: no tuning results for setting {self.tune} (only the runs chosen by --plan have them)")
            if self.data is None and tune_mean is not None:
                # no CFU result, but a tuning result: show its mean image and its events alone
                img = (tune_mean - tune_mean.min()) / max(float(np.ptp(tune_mean)), 1e-9)
                spf = 2.0
                try:
                    spf = float(json.loads(aqua_files(r["tif"])["json"].read_text())["seconds_per_frame"])
                except Exception:
                    pass
                self.data = {"tif": r["tif"], "image": img, "seconds_per_frame": spf, "crop_gap": 0,
                             "cfu": {"maps": np.zeros((0,) + img.shape, np.float32), "events": [], "occurrence": np.zeros((0, 0), bool),
                                     "dff": np.zeros((0, 0), np.float32), "datPro": img, "freq": [], "n": 0,
                                     "note": "no CFU result — tuning events only"},
                             "rejected": set(), "added": [], "n_events": [], "notes": "", "tune_only": True}
            if self.data is None:
                messagebox.showwarning("no CFU file",
                                       f"{r['run']}\n\nhas no *_res_cfu.mat yet — run run_cfu_batch.m"
                                       + ("" if self.tune is None else f"\n(and no tuning result for setting {self.tune} either)"))
            else:
                maps = self.data["cfu"]["maps"]
                self.masks = cfu_masks(maps)
                self.outlines = [mask_outline(m) for m in self.masks]
                self.colours = [PALETTE[c] for c in assign_colours(self.masks)]
                self.centroids = [tuple(np.argwhere(m).mean(axis=0)) if m.any() else None for m in self.masks]
                self.multi = set(); self.feats = []; self.fmask = None
                try:
                    from afferent2p.cfu_shape import cfu_features, fibre_mask_for
                    ref = self.data.get("mean_raw")
                    if ref is None or ref.shape != self.shape():
                        ref = self.data["image"]
                    if ref is not None and ref.shape == self.shape() and self.masks:
                        self.fmask, _ = fibre_mask_for(ref, self.um_per_px)
                        self.feats = cfu_features(self.masks, ref, self.um_per_px, fmask=self.fmask)
                        n_sugg = sum(1 for f in self.feats if f["suggest"] != "fibre")
                        if n_sugg:
                            log(f"  {r['run']}: {n_sugg} of {len(self.feats)} CFUs look like non-fibre structures (column 'looks'; "
                                f"'Select suggested' picks them, x marks them noise)")
                except Exception as e:
                    log(f"  note: CFU shape suggestions not available ({e})")
                occ = self.data["cfu"]["occurrence"]
                self.onsets = [np.flatnonzero(occ[k]) for k in range(occ.shape[0])] if occ.size else []
                self.open_movie(r, log_it=True)
            if self.movie is not None:
                self.frame = min(self.frame, self.movie.T - 1)
                self._slider_lock = True
                self.slider.configure(to=max(self.movie.T - 1, 1)); self.slider.set(self.frame)
                self._slider_lock = False
            if self.data is not None and tune_recs:
                # the tuning events ride along as read-only "added" records (never written to the curation file)
                self.data["added"] = [x for x in self.data["added"] if x.get("source") != "tune"]
                if self.movie is not None:
                    tune_recs, _, _ = tune_events_records(self.mouse_dir, r["run"], self.tune, self.tune_variant,
                                                          movie=self.movie, spf=self.data["seconds_per_frame"], log=lambda *_: None)
                self.data["added"].extend(tune_recs)
                self.data["tune_label"] = tune_label
            if self.view_shape != self.shape():
                self.view = None
            self.refresh_tree(); self.fill_ctree(); self.draw()

        def open_movie(self, r, log_it=False):
            """The movie panel's frames: AQuA2's registered movie (default, when the
            _AQuA2.mat can be read) or the binned tiff as recorded, per the radio."""
            if self.movie is not None:
                self.movie.close(); self.movie = None
            want = self.src_var.get() if hasattr(self, "src_var") else "registered"
            try:
                self.movie = Movie(r["tif"], gap=self.data.get("crop_gap", 0), shape=self.shape(), source=want, log=log)
                if self.movie.source == "raw":
                    fit = self.movie.align_to(self.data["image"])
                    if fit is not None:
                        self.data["crop_gap"] = (fit[0], fit[1])      # hand-drawn regions use the same offset
                        if log_it:
                            log(f"  {r['run']}: AQuA2's crop sits at row {fit[0]}, column {fit[1]} of the tiff (match {fit[2]:.2f})")
                    if want == "registered" and log_it:
                        log(f"  {r['run']}: no readable registered movie in the _AQuA2.mat — showing the binned tiff as recorded")
                elif log_it:
                    log(f"  {r['run']}: movie = AQuA2's registered movie (res.datOrg1 of the _AQuA2.mat)")
            except Exception as e:
                log(f"  note: {r['tif'].name}: movie not available ({e})")
                self.movie = None
            self.data["movie_source"] = self.movie.source if self.movie is not None else None

        def set_source(self):
            """Radio: registered (AQuA2) / as recorded — reopen the movie, keep the frame."""
            if self.data is None:
                return
            self.playing = False
            self.open_movie(self.runs[self.i], log_it=True)
            if self.movie is not None:
                self.frame = min(self.frame, self.movie.T - 1)
            self.focus_canvas(); self.draw()

        def toggle_hand(self):
            self.hand = bool(self.hand_var.get())
            try:
                self.canvas.get_tk_widget().configure(cursor="fleur" if self.hand else "")
            except Exception:
                pass
            self.focus_canvas()

        def pan_by(self, fx: float, fy: float):
            """Move the view by a fraction of its width / height (keyboard panning)."""
            if self.data is None:
                return
            x0, x1 = self.ax_img.get_xlim(); y0, y1 = self.ax_img.get_ylim()
            dx = (x1 - x0) * fx; dy = (y0 - y1) * fy
            self.set_view(x0 + dx, x1 + dx, y0 + dy, y1 + dy)

        def refresh_tree(self):
            for k, r in enumerate(self.runs):
                if k == self.i and self.data is not None:
                    n = self.data["cfu"]["n"]
                    kept = n - len(self.data["rejected"])
                    extra = f"+{len(self.data['added'])}" if self.data["added"] else ""
                    self.tree.item(str(k), values=(f"{kept}/{n}{extra}" + (" (old decisions ignored)" if self.data.get("stale") else ""),))
            self.tree.selection_set(str(self.i))

        def on_tree(self, _):
            s = self.tree.selection()
            if s and int(s[0]) != self.i:
                self.load(int(s[0]))

        def fill_ctree(self):
            self._ctree_lock = True
            self.ctree.delete(*self.ctree.get_children())
            d = self.data
            if d is not None:
                for k in range(d["cfu"]["n"]):
                    self.ctree.tag_configure(f"c{k}", foreground=self.colours[k])
                    rej = (k + 1) in d["rejected"]
                    self.ctree.insert("", "end", iid=f"cfu{k}", text=f"CFU {k + 1}",
                                      values=(d["n_events"][k], self.looks(k), "noise" if rej else "kept"),
                                      tags=self.ctags(k))
                self.ctree.tag_configure("tune", foreground=TUNE_COLOUR)
                for j, rec in enumerate(d["added"]):
                    if rec.get("source") == "tune":
                        self.ctree.insert("", "end", iid=f"add{j}", text=rec["label"].split(" ")[0] + " tuning",
                                          values=(f"fr {rec['onset_frames'][0]}", f"dF/F {rec.get('peak_dff', 0):.2f}"), tags=("tune",))
                    else:
                        self.ctree.insert("", "end", iid=f"add{j}", text=f"+{j + 1} (added)",
                                          values=(rec["n_events"], "added"), tags=("added",))
            self._ctree_lock = False
            self.sync_ctree()

        def refresh_ctree_status(self):
            d = self.data
            if d is None:
                return
            for k in range(d["cfu"]["n"]):
                rej = (k + 1) in d["rejected"]
                if self.ctree.exists(f"cfu{k}"):
                    self.ctree.item(f"cfu{k}", values=(d["n_events"][k], self.looks(k), "noise" if rej else "kept"),
                                    tags=self.ctags(k))

        def looks(self, k) -> str:
            if k < len(self.feats):
                from afferent2p.cfu_shape import looks_label
                return looks_label(self.feats[k])
            return ""

        def ctags(self, k):
            rej = (k + 1) in self.data["rejected"]
            sugg = k < len(self.feats) and self.feats[k]["suggest"] != "fibre" and not rej
            return (("noise",) if rej else (f"c{k}",)) + (("sugg",) if sugg else ())

        def sync_ctree(self):
            self._ctree_lock = True
            try:
                iids = [f"cfu{k}" for k in sorted(self.multi) if self.ctree.exists(f"cfu{k}")]
                if self.sel is not None:
                    iid = f"cfu{self.sel[1]}" if self.sel[0] == "cfu" else f"add{self.sel[1]}"
                    if self.ctree.exists(iid) and iid not in iids:
                        iids.append(iid)
                if iids:
                    self.ctree.selection_set(iids); self.ctree.see(iids[-1])
                else:
                    self.ctree.selection_remove(*self.ctree.selection())
            finally:
                self._ctree_lock = False

        def on_ctree(self, _):
            if self._ctree_lock or self.data is None:
                return
            s = self.ctree.selection()
            if not s:
                return
            cfus = {int(i[3:]) for i in s if i.startswith("cfu")}
            iid = s[-1]
            new = ("cfu", int(iid[3:])) if iid.startswith("cfu") else ("added", int(iid[3:]))
            multi = cfus if len(cfus) > 1 else set()
            if new != self.sel or multi != self.multi:
                self.sel = new; self.multi = multi; self.draw()

        # ---- drawing
        def sel_paths(self):
            if self.data is None:
                return []
            paths = []
            for k in sorted(self.multi):
                if k < len(self.outlines):
                    paths += self.outlines[k]
            if self.sel is not None:
                kind, k = self.sel
                if kind == "cfu":
                    if k not in self.multi and k < len(self.outlines):
                        paths += self.outlines[k]
                else:
                    paths += polygon_outline(self.data["added"][k]["polygon_xy"])
            return paths

        def outline_collections(self, animated=False):
            """(all outlines, selected-black, selected-white) LineCollections."""
            d = self.data
            segs, cols, lws = [], [], []
            for k, paths in enumerate(self.outlines):
                rejected = (k + 1) in d["rejected"]
                col = NOISE_COLOUR if rejected else self.colours[k]
                for p in paths:
                    segs.append(p); cols.append(col); lws.append(0.7 if rejected else 1.3)
            for rec in d["added"]:
                tune_rec = rec.get("source") == "tune"
                for p in polygon_outline(rec["polygon_xy"]):
                    segs.append(p); cols.append(TUNE_COLOUR if tune_rec else ADDED_COLOUR); lws.append(0.9 if tune_rec else 1.4)
            lc = LineCollection(segs, colors=cols, linewidths=lws, animated=animated)
            sp = self.sel_paths()
            lb = LineCollection(sp, colors=["black"], linewidths=4.0, animated=animated)
            lw = LineCollection(sp, colors=["white"], linewidths=1.8, animated=animated)
            return lc, lb, lw

        def apply_view(self):
            H, W = self.shape()
            if self.view is not None and self.view_shape == (H, W):
                self.ax_img.set_xlim(*self.view[0]); self.ax_img.set_ylim(*self.view[1])
            else:
                self.ax_img.set_xlim(-0.5, W - 0.5); self.ax_img.set_ylim(H - 0.5, -0.5)

        def draw(self):
            d = self.data
            for ax in (self.ax_img, self.ax_mov, self.ax_cur, self.ax_all):
                ax.clear()
            self.ax_txt.clear(); self.ax_txt.axis("off")
            self._anim = []; self._bg = None
            self.mov_im = self.mov_txt = self.cur_line = self.all_line = self.mov_lines = self.mov_sel = None
            r = self.runs[self.i]
            if d is None:
                self.ax_img.text(0.5, 0.5, "no CFU result for this run", ha="center")
                self.ax_img.axis("off"); self.ax_mov.axis("off"); self.canvas.draw(); return
            img = d["image"]
            H, W = self.shape()
            if img is not None:
                self.ax_img.imshow(img, cmap="gray", interpolation="nearest",
                                   vmin=float(np.percentile(img, 1)), vmax=float(np.percentile(img, 99.5)))
            if self.show_fibre and self.fmask is not None and self.fmask.shape == (H, W) and self.fmask.any():
                try:
                    self.ax_img.contour(self.fmask.astype(float), levels=[0.5], colors=["#eda100"], linewidths=0.5, alpha=0.8)
                except Exception:
                    pass
            lc, lb, lw = self.outline_collections()
            self.ax_img.add_collection(lc); self.ax_img.add_collection(lb); self.ax_img.add_collection(lw)
            for k, cen in enumerate(self.centroids):
                if cen is None:
                    continue
                rejected = (k + 1) in d["rejected"]
                selected = self.sel == ("cfu", k) or k in self.multi
                cy, cx = cen
                self.ax_img.text(cx, cy, str(k + 1), color="white" if selected else (NOISE_COLOUR if rejected else self.colours[k]),
                                 fontsize=9 if selected else 7, fontweight="bold" if selected else "normal",
                                 ha="center", va="center", clip_on=True,
                                 bbox=dict(boxstyle="round,pad=0.15", fc="black", ec="none", alpha=0.6) if selected else None)
            for j, rec in enumerate(d["added"]):
                poly = np.asarray(rec["polygon_xy"], float)
                tune_rec = rec.get("source") == "tune"
                self.ax_img.text(poly[:, 0].mean(), poly[:, 1].min(), rec["label"].split(" ")[0] if tune_rec else f"+{j + 1}",
                                 color=TUNE_COLOUR if tune_rec else ADDED_COLOUR, fontsize=6 if tune_rec else 7,
                                 ha="center", va="bottom", clip_on=True)
            if self.adding:
                p = np.asarray(self.adding, float)
                self.ax_img.plot(p[:, 0], p[:, 1], "o-", color="lime", ms=3, lw=1)
            n = d["cfu"]["n"]
            n_sugg = sum(1 for f in self.feats if f["suggest"] != "fibre")
            import re as _re
            short = (_re.search(r"(\d{6}_\d{3})", r["run"]) or [r["run"]])[0]
            rules = d["cfu"].get("rules", "")                             # v41: which AQuA2 CFU rules made this file
            rules_txt = ""
            if rules and not self.multi:
                rules_txt = ("   ·   AQuA2 newer rules (basin splitting + post-merge): " if "post-merge" in rules
                             else "   ·   AQuA2 rule: ") + rules
            self.ax_img.set_title(
                f"[{r['condition']}] {short}   —   {n - len(d['rejected'])} of {n} CFUs kept"
                + rules_txt
                + (f"   ·   {len(self.multi)} selected (k / x act on all)" if self.multi else "")
                + (f"   ·   {n_sugg} look non-fibre (u selects them)" if n_sugg and not self.multi else "")
                + ("   ·   earlier decisions (other CFU result) ignored — curate afresh" if d.get("stale") and not self.multi else "")
                + ((lambda na, nt: (f", {na} added" if na else "") + (f", {nt} tuning events ({d.get('tune_label', '')}) in magenta" if nt else ""))(
                    sum(1 for x in d["added"] if x.get("source") != "tune"), sum(1 for x in d["added"] if x.get("source") == "tune")))
                + ("\nNO CFUs in this result: no place had enough repeated events (CFU min events in tab 4); "
                   "'a' draws a region by hand" if n == 0 and not d["added"] else "")
                + ("\nADD MODE: click vertices, Enter to close, Esc to cancel" if self.adding is not None else ""),
                fontsize=9)
            self.ax_img.set_xticks([]); self.ax_img.set_yticks([])
            self.ax_mov.set_xticks([]); self.ax_mov.set_yticks([])
            self.draw_movie_panel()
            self.apply_view()
            self.draw_curve(); self.draw_raster()
            self.canvas.draw()

        def draw_movie_panel(self):
            ax = self.ax_mov
            if self.movie is None:
                ax.text(0.5, 0.5, "movie not available\n(the binned tiff could not be opened)", ha="center",
                        va="center", transform=ax.transAxes, fontsize=9)
                ax.set_title("movie", fontsize=9); return
            f = self.movie.frame(self.frame, self.mode)
            vmin, vmax = self.movie.clim(self.mode, float(self.contrast.get()))
            self.mov_im = ax.imshow(f, cmap="gray", vmin=vmin, vmax=vmax, interpolation="nearest", animated=True)
            lc, lb, lw = self.outline_collections(animated=True)
            lb.set_linewidth(2.6); lw.set_linewidth(1.1)
            ax.add_collection(lc); ax.add_collection(lb); ax.add_collection(lw)
            self.mov_lines, self.mov_sel = lc, (lb, lw)
            self.mov_txt = ax.text(0.01, 0.99, "", transform=ax.transAxes, color="white", fontsize=8,
                                   va="top", ha="left", animated=True,
                                   bbox=dict(boxstyle="round,pad=0.2", fc="black", ec="none", alpha=0.55))
            self._anim += [(ax, self.mov_im), (ax, lc), (ax, lb), (ax, lw), (ax, self.mov_txt)]
            src = ("movie: AQuA2's REGISTERED frames (what it detected on; same frame as the maps)"
                   if self.movie is not None and self.movie.source == "registered"
                   else "movie: binned tiff AS RECORDED (AQuA2's registration not applied, its crop offset is)")
            ax.set_title(src + "\nspace play · ←/→ frame · [ ] events of this CFU · scroll zoom · "
                         "move: right-drag / h+drag / ctrl+arrows", fontsize=7.5)
            self.update_frame_artists()

        def draw_curve(self):
            d = self.data
            spf = d["seconds_per_frame"]
            if self.sel is None:
                self.ax_cur.text(0.5, 0.5, "click a CFU (or pick one in the list on the left)", ha="center",
                                 transform=self.ax_cur.transAxes)
                self.ax_cur.set_xticks([]); self.ax_cur.set_yticks([])
                self.ax_txt.text(0, 1, self.help_text(), va="top", fontsize=8, family="monospace")
                self.add_cursor(self.ax_cur, "cur_line"); return
            kind, k = self.sel
            if kind == "cfu":
                dff = d["cfu"]["dff"][k]
                onsets = self.onsets[k] if k < len(self.onsets) else []
                nev, npix = d["n_events"][k], int(self.masks[k].sum())
                title = f"CFU {k + 1}" + ("  (noise)" if (k + 1) in d["rejected"] else "")
                extra = ""
                if k < len(self.feats):
                    f_ = self.feats[k]
                    extra = f"\nlooks: {f_['suggest']} — {f_['reason']}"
                col = self.colours[k]
            else:
                rec = d["added"][k]
                dff = np.asarray(rec["dff"], float)
                onsets = np.asarray(rec["onset_frames"], int) - 1
                nev, npix = rec["n_events"], rec["npix"]
                if rec.get("source") == "tune":
                    title = rec["label"]
                    extra = rec.get("note", "") + "\n(trace = mean of the box in the movie; not part of the curation)"
                    col = TUNE_COLOUR
                else:
                    title = f"added region {k + 1}"
                    extra = f"\nthreshold z > {rec['z_thresh']:g} (not AQuA2)"
                    col = ADDED_COLOUR
            t = np.arange(len(dff)) * spf / 60.0
            if len(dff) == 0:
                self.ax_cur.text(0.5, 0.5, "no trace (movie not available)", ha="center", transform=self.ax_cur.transAxes)
            self.ax_cur.plot(t, dff, lw=0.8, color="k")
            for o in np.atleast_1d(onsets):
                if 0 <= o < len(dff):
                    self.ax_cur.plot(t[o], dff[o], "v", color=col, mec="black", mew=0.4, ms=6)
            self.ax_cur.set_xlim(0, max(t[-1] if len(t) else 1, 1e-3))
            self.ax_cur.set_xlabel("min"); self.ax_cur.set_ylabel("dF/F")
            self.ax_cur.set_title(f"{title}  —  click in the trace to move the movie there", fontsize=9, color=col if kind == "cfu" else "k")
            minutes = len(dff) * spf / 60.0
            peak = np.nanmax(dff) if len(dff) else float("nan")
            self.ax_txt.text(0, 1, f"{title}\n{nev} events, {nev / max(minutes, 1e-9):.2f} events/min\n"
                                   f"{npix} px, peak dF/F {peak:.2f}{extra}\n\n" + self.help_text(),
                             va="top", fontsize=8, family="monospace")
            self.add_cursor(self.ax_cur, "cur_line")

        @staticmethod
        def help_text():
            return ("k keep   x noise   a add region\n"
                    "tab / shift-tab  next / previous CFU\n"
                    "click: select (again = next overlapping)\n"
                    "drag a box: select several (shift adds,\n"
                    "ctrl-click toggles, Esc clears); k/x act on all\n"
                    "u select suggested non-fibre   f fibre mask\n"
                    "space play   <- -> frame (shift: 10)\n"
                    "[ ] previous / next event of this CFU\n"
                    "v raw / dF/F   scroll zoom   0 whole field\n"
                    "move the view: right-drag, h + left-drag,\n"
                    "ctrl+arrows        n p runs   s save")

        def add_cursor(self, ax, attr):
            if self.movie is None:
                return
            t = self.frame * self.data["seconds_per_frame"] / 60.0
            line = ax.axvline(t, color="#1f77d8", lw=1.2, alpha=0.9, animated=True)
            setattr(self, attr, line)
            self._anim.append((ax, line))

        def draw_raster(self):
            d = self.data
            occ = d["cfu"]["occurrence"]
            if occ.size == 0:
                self.ax_all.axis("off"); return
            spf = d["seconds_per_frame"]
            xs, ys, cs = [], [], []
            for k, fr in enumerate(self.onsets):
                if fr.size == 0:
                    continue
                col = "0.8" if (k + 1) in d["rejected"] else self.colours[k]
                xs.append(fr * spf / 60.0); ys.append(np.full(fr.size, k + 1)); cs += [col] * fr.size
            if xs:
                self.ax_all.scatter(np.concatenate(xs), np.concatenate(ys), marker="|", c=cs, s=30, linewidths=1.0)
            if self.sel is not None and self.sel[0] == "cfu":
                self.ax_all.axhspan(self.sel[1] + 0.5, self.sel[1] + 1.5, color="#ffd24d", alpha=0.75, lw=0)
            self.ax_all.set_xlim(0, occ.shape[1] * spf / 60.0)
            self.ax_all.set_ylim(0.5, occ.shape[0] + 0.5)
            self.ax_all.set_xlabel("time (min)"); self.ax_all.set_ylabel("CFU")
            self.ax_all.set_title("event onsets per CFU (grey = marked noise); click to move the movie there", fontsize=8)
            self.add_cursor(self.ax_all, "all_line")

        # ---- movie
        def on_draw(self, _event):
            """after every full draw: remember the background, paint the animated artists"""
            self._bg = self.canvas.copy_from_bbox(self.fig.bbox)
            for ax, a in self._anim:
                ax.draw_artist(a)

        def update_frame_artists(self):
            if self.movie is None or self.mov_im is None:
                return
            self.mov_im.set_data(self.movie.frame(self.frame, self.mode))
            spf = self.data["seconds_per_frame"]
            t_min = self.frame * spf / 60.0
            self.mov_txt.set_text(f"frame {self.frame + 1} / {self.movie.T}   {t_min:.2f} min   "
                                  f"{'dF/F' if self.mode == 'dff' else 'raw'}"
                                  + (f"\n{self.event_note}" if self.event_note else ""))
            show = bool(self.outl_var.get())
            self.mov_lines.set_visible(show)
            for a in self.mov_sel:
                a.set_visible(show)
            for line in (self.cur_line, self.all_line):
                if line is not None:
                    line.set_xdata([t_min, t_min])
            self.frame_lbl.config(text=f"frame {self.frame + 1} / {self.movie.T}   {t_min:.2f} min")

        def show_frame(self):
            if self.movie is None or self.mov_im is None:
                return
            self.update_frame_artists()
            if self._bg is None:
                self.canvas.draw(); return
            self.canvas.restore_region(self._bg)
            for ax, a in self._anim:
                ax.draw_artist(a)
            self.canvas.blit(self.fig.bbox)
            self.canvas.flush_events()

        def set_frame(self, i, note_keep=False):
            if self.movie is None:
                return
            i = int(min(max(i, 0), self.movie.T - 1))
            if not note_keep:
                self.event_note = ""
            self.frame = i
            self._slider_lock = True
            try:
                self.slider.set(i)
            finally:
                self._slider_lock = False
            self.show_frame()

        def on_slider(self, value):
            if self._slider_lock:
                return
            self.set_frame(int(round(float(value))))

        def toggle_play(self):
            if self.movie is None:
                return
            self.playing = not self.playing
            self.play_btn.config(text="pause (space)" if self.playing else "play (space)")
            if self.playing:
                if self.frame >= self.movie.T - 1:
                    self.frame = -1
                self._tick()

        def _tick(self):
            if not self.playing or self.movie is None:
                return
            if self.frame + 1 >= self.movie.T:
                self.playing = False; self.play_btn.config(text="play (space)"); return
            self.set_frame(self.frame + 1, note_keep=True)
            try:
                fps = max(0.5, float(self.fps_var.get()))
            except ValueError:
                fps = 8.0
            self.after(max(5, int(1000 / fps)), self._tick)

        def set_mode(self):
            self.mode = self.mode_var.get()
            if self.movie is None or self.mov_im is None:
                return
            if self.mode == "dff" and self.movie.f0 is None:
                self.frame_lbl.config(text="computing F0 ..."); self.update_idletasks()
            self.apply_clim()

        def apply_clim(self):
            if self.movie is None or self.mov_im is None:
                return
            self.mov_im.set_clim(*self.movie.clim(self.mode, float(self.contrast.get())))
            self.show_frame()

        def sel_onsets(self):
            d = self.data
            if d is None or self.sel is None:
                return np.zeros(0, int)
            kind, k = self.sel
            if kind == "cfu":
                return self.onsets[k] if k < len(self.onsets) else np.zeros(0, int)
            return np.asarray(d["added"][k]["onset_frames"], int) - 1

        def jump_event(self, direction):
            if self.movie is None or self.sel is None:
                return
            onsets = np.sort(self.sel_onsets())
            if onsets.size == 0:
                self.event_note = "this CFU has no events"; self.show_frame(); return
            cur = self.frame + 3
            if direction > 0:
                cand = onsets[onsets > cur]
                target = cand[0] if cand.size else onsets[0]
            else:
                cand = onsets[onsets < cur]
                target = cand[-1] if cand.size else onsets[-1]
            idx = int(np.searchsorted(onsets, target))
            if self.sel[0] == "cfu":
                name = f"CFU {self.sel[1] + 1}"
            else:
                rec = self.data["added"][self.sel[1]]
                name = rec["label"] if rec.get("source") == "tune" else f"added region {self.sel[1] + 1}"
            self.event_note = f"event {idx + 1} / {onsets.size} of {name}: onset at frame {int(target) + 1}"
            self.playing = False; self.play_btn.config(text="play (space)")
            self.set_frame(int(target) - 3, note_keep=True)

        # ---- zoom and pan (both image panels share their axes)
        def set_view(self, x0, x1, y0, y1):
            H, W = self.shape()
            w = min(max(x1 - x0, 6.0), W); h = min(max(y0 - y1, 6.0), H)      # y runs downwards
            x0 = min(max(x0, -0.5), W - 0.5 - w); y1 = min(max(y1, -0.5), H - 0.5 - h)
            self.ax_img.set_xlim(x0, x0 + w); self.ax_img.set_ylim(y1 + h, y1)
            self.view = (self.ax_img.get_xlim(), self.ax_img.get_ylim()); self.view_shape = (H, W)
            self._bg = None; self.canvas.draw_idle()

        def zoom_at(self, x, y, factor):
            if self.data is None:
                return
            x0, x1 = self.ax_img.get_xlim(); y0, y1 = self.ax_img.get_ylim()
            H, W = self.shape()
            if factor >= 1 and (x1 - x0) >= W and (y0 - y1) >= H:
                return
            if x is None or y is None:
                x, y = (x0 + x1) / 2, (y0 + y1) / 2
            nx0 = x - (x - x0) * factor; nx1 = x + (x1 - x) * factor
            ny1 = y - (y - y1) * factor; ny0 = y + (y0 - y) * factor
            self.set_view(nx0, nx1, ny0, ny1)

        def reset_view(self):
            self.view = None
            self.apply_view()
            self._bg = None; self.canvas.draw_idle()

        def on_scroll(self, event):
            if event.inaxes not in (self.ax_img, self.ax_mov):
                return
            self.zoom_at(event.xdata, event.ydata, 1 / 1.3 if event.button == "up" else 1.3)

        # ---- interaction
        def on_click(self, event):
            if self.data is None or event.inaxes is None:
                return
            self.focus_canvas()
            if event.inaxes in (self.ax_cur, self.ax_all):
                if event.button == 1 and self.movie is not None and event.xdata is not None:
                    self.set_frame(int(round(event.xdata * 60.0 / self.data["seconds_per_frame"])))
                return
            if event.inaxes not in (self.ax_img, self.ax_mov):
                return
            if event.button in (2, 3) or (event.button == 1 and self.hand and self.adding is None):
                self._pan = (event.x, event.y, self.ax_img.get_xlim(), self.ax_img.get_ylim()); return
            if event.button != 1 or event.xdata is None:
                return
            x, y = event.xdata, event.ydata
            if self.adding is not None:
                self.adding.append((x, y)); self.draw(); return
            key = (event.key or "").lower()
            self._drag = (event.x, event.y, x, y, "shift" in key, "control" in key or "ctrl" in key)

        def _click_select(self, x, y, ctrl: bool):
            """A plain click: the CFU under the cursor (again = the next overlapping one); ctrl toggles it in the multi-selection."""
            xi, yi = int(round(x)), int(round(y))
            H, W = self.shape()
            inside = 0 <= yi < H and 0 <= xi < W
            hits = [k for k, m in enumerate(self.masks) if inside and m[yi, xi]]
            if ctrl:
                if hits:
                    cur = self.sel[1] if (self.sel and self.sel[0] == "cfu") else -1
                    k = hits[(hits.index(cur) + 1) % len(hits)] if cur in hits and len(hits) > 1 else hits[0]
                    if self.sel == ("cfu", k) or k in self.multi:
                        self.multi.discard(k)
                        if self.sel == ("cfu", k):
                            self.sel = ("cfu", next(iter(sorted(self.multi)))) if self.multi else None
                    else:
                        if self.sel is not None and self.sel[0] == "cfu":
                            self.multi.add(self.sel[1])
                        self.multi.add(k); self.sel = ("cfu", k)
                return
            self.multi = set()
            if hits:
                cur = self.sel[1] if (self.sel and self.sel[0] == "cfu") else -1
                nxt = hits[(hits.index(cur) + 1) % len(hits)] if cur in hits else hits[0]
                self.sel = ("cfu", nxt)
            else:
                for j, rec in enumerate(self.data["added"]):
                    if inside and added_mask(rec, (H, W))[yi, xi]:
                        self.sel = ("added", j); break
                else:
                    self.sel = None

        def _box_select(self, x0, y0, x1, y1, add: bool):
            """Every CFU whose centroid lies in the box, or half of whose footprint does."""
            xa, xb = sorted((x0, x1)); ya, yb = sorted((y0, y1))
            H, W = self.shape()
            r0, r1 = int(np.clip(np.floor(ya), 0, H)), int(np.clip(np.ceil(yb) + 1, 0, H))
            c0, c1 = int(np.clip(np.floor(xa), 0, W)), int(np.clip(np.ceil(xb) + 1, 0, W))
            hits = set()
            for k, m in enumerate(self.masks):
                cen = self.centroids[k]
                if cen is not None and ya <= cen[0] <= yb and xa <= cen[1] <= xb:
                    hits.add(k); continue
                tot = int(m.sum())
                if tot and r1 > r0 and c1 > c0 and int(m[r0:r1, c0:c1].sum()) >= 0.5 * tot:
                    hits.add(k)
            if add:
                if self.sel is not None and self.sel[0] == "cfu":
                    hits.add(self.sel[1])
                self.multi |= hits
            else:
                self.multi = hits
            if self.multi:
                self.sel = ("cfu", min(self.multi)) if not (self.sel and self.sel[0] == "cfu" and self.sel[1] in self.multi) else self.sel
            elif not add:
                self.sel = None

        def _band_update(self, x1, y1):
            x0, y0 = self._drag[2], self._drag[3]
            if self._band is None:
                self._band = []
                for ax in (self.ax_img, self.ax_mov):
                    rect = Rectangle((x0, y0), 0, 0, fill=False, ec="#00e5ff", lw=1.2, ls="--", animated=True)
                    ax.add_patch(rect); self._band.append(rect)
                    self._anim.append((ax, rect))
            for rect in self._band:
                rect.set_xy((min(x0, x1), min(y0, y1))); rect.set_width(abs(x1 - x0)); rect.set_height(abs(y1 - y0))
            if self._bg is None:
                self.canvas.draw(); return
            self.canvas.restore_region(self._bg)
            for ax, a in self._anim:
                ax.draw_artist(a)
            self.canvas.blit(self.fig.bbox)
            self.canvas.flush_events()

        def on_release(self, event):
            if self._pan is not None and (event.button in (2, 3) or (event.button == 1 and self.hand)):
                self._pan = None
                return
            if self._drag is not None and event.button == 1:
                x0px, y0px, x0, y0, shift, ctrl = self._drag
                self._drag = None
                moved = abs(event.x - x0px) > 4 or abs(event.y - y0px) > 4
                if self._band is not None:
                    for rect in self._band:
                        try:
                            rect.remove()
                        except Exception:
                            pass
                    self._anim = [(ax, a) for ax, a in self._anim if a not in self._band]
                    self._band = None
                if moved and event.inaxes in (self.ax_img, self.ax_mov) and event.xdata is not None:
                    self._box_select(x0, y0, event.xdata, event.ydata, add=shift or ctrl)
                elif not moved:
                    self._click_select(x0, y0, ctrl)
                self.sync_ctree(); self.draw()

        def on_motion(self, event):
            if self._drag is not None and event.inaxes in (self.ax_img, self.ax_mov) and event.xdata is not None:
                if abs(event.x - self._drag[0]) > 4 or abs(event.y - self._drag[1]) > 4:
                    self._band_update(event.xdata, event.ydata)
                return
            if self._pan is None:
                return
            x0px, y0px, (xl0, xl1), (yl0, yl1) = self._pan
            inv = self.ax_img.transData.inverted()
            (xa, ya), (xb, yb) = inv.transform([(x0px, y0px), (event.x, event.y)])
            dx, dy = xb - xa, yb - ya
            self.set_view(xl0 - dx, xl1 - dx, yl0 - dy, yl1 - dy)

        def on_key(self, event):
            w = self.focus_get()
            try:
                if w is not None and w.winfo_class() in TYPING:
                    return
            except Exception:
                pass
            k = (event.keysym or "").lower()
            k = {"enter": "return", "kp_enter": "return", "iso_left_tab": "tab", "kp_add": "plus", "kp_subtract": "minus"}.get(k, k)
            shift = bool(event.state & 0x1)
            ctrl = bool(event.state & 0x4)
            if ctrl and k in ("left", "right", "up", "down"):
                f = 0.25
                self.pan_by({"left": -f, "right": f}.get(k, 0.0), {"up": -f, "down": f}.get(k, 0.0)); return
            if self.adding is not None:
                if k == "return":
                    self.finish_add()
                elif k == "escape":
                    self.adding = None; self.draw()
                return
            if k == "k":
                self.mark(True)
            elif k in ("x", "r"):
                self.mark(False)
            elif k == "escape":
                if self.multi:
                    self.multi = set(); self.sync_ctree(); self.draw()
            elif k == "u":
                self.select_suggested()
            elif k == "f":
                self.fibre_var.set(not self.fibre_var.get()); self.toggle_fibre()
            elif k == "h":
                self.hand_var.set(not self.hand_var.get()); self.toggle_hand()
            elif k == "n":
                self.step(1)
            elif k == "p":
                self.step(-1)
            elif k == "a":
                self.start_add()
            elif k == "s":
                self.save()
            elif k == "space":
                self.toggle_play()
            elif k == "left":
                self.set_frame(self.frame - (10 if shift else 1))
            elif k == "right":
                self.set_frame(self.frame + (10 if shift else 1))
            elif k == "home":
                self.set_frame(0)
            elif k == "end":
                self.set_frame(10 ** 9)
            elif k == "bracketleft":
                self.jump_event(-1)
            elif k == "bracketright":
                self.jump_event(1)
            elif k == "v":
                self.mode_var.set("dff" if self.mode == "raw" else "raw"); self.set_mode()
            elif k == "0":
                self.reset_view()
            elif k in ("plus", "equal"):
                self.zoom_at(None, None, 1 / 1.3)
            elif k == "minus":
                self.zoom_at(None, None, 1.3)
            elif k == "tab":
                self.step_cfu(-1 if shift else 1)
            elif k == "period":
                self.step_cfu(1)
            elif k == "comma":
                self.step_cfu(-1)

        def step(self, d):
            self.save(quiet=True); self.load(self.i + d)

        def step_cfu(self, d):
            if self.data is None:
                return
            n = self.data["cfu"]["n"]; na = len(self.data["added"])
            items = [("cfu", k) for k in range(n)] + [("added", j) for j in range(na)]
            if not items:
                return
            if self.sel in items:
                self.sel = items[(items.index(self.sel) + d) % len(items)]
            else:
                self.sel = items[0 if d > 0 else -1]
            self.sync_ctree(); self.draw()

        def mark(self, keep):
            if self.data is None:
                return
            targets = set(self.multi)
            if self.sel is not None and self.sel[0] == "cfu":
                targets.add(self.sel[1])
            if not targets:
                return
            for k in targets:
                if keep:
                    self.data["rejected"].discard(k + 1)
                else:
                    self.data["rejected"].add(k + 1)
            self.refresh_tree(); self.refresh_ctree_status(); self.draw()

        def select_suggested(self):
            """Select every CFU the shape rules call non-fibre (off the fibre mask, or a compact bright blob)."""
            if self.data is None or not self.feats:
                return
            hits = {k for k, f in enumerate(self.feats) if f["suggest"] != "fibre" and (k + 1) not in self.data["rejected"]}
            if not hits:
                log("  no CFU looks like a non-fibre structure (or all of them are marked noise already)")
                self.multi = set()
            else:
                self.multi = hits; self.sel = ("cfu", min(hits))
                log(f"  {len(hits)} CFU(s) selected: " + ", ".join(f"CFU {k + 1} ({self.feats[k]['reason']})" for k in sorted(hits)))
            self.sync_ctree(); self.draw()

        def toggle_fibre(self):
            self.show_fibre = bool(self.fibre_var.get())
            self.draw()

        def mark_all(self, keep):
            if self.data is None:
                return
            n = self.data["cfu"]["n"]
            self.data["rejected"] = set() if keep else set(range(1, n + 1))
            self.refresh_tree(); self.refresh_ctree_status(); self.draw()

        def start_add(self):
            if self.data is None:
                return
            self.adding = []; self.sel = None; self.sync_ctree(); self.draw()

        def finish_add(self):
            pts = self.adding or []
            self.adding = None
            if len(pts) < 3:
                self.draw(); return
            d = self.data
            shape = self.shape()
            if shape == (1, 1):
                messagebox.showerror("cannot add a region",
                                     "this run has neither CFU maps nor a mean image")
                return
            try:
                reg = self.movie.arr if (self.movie is not None and self.movie.source == "registered") else None
                rec = add_region(d["tif"], pts, shape, d["seconds_per_frame"],
                                 gap=d.get("crop_gap", 0), movie=reg)
            except Exception as e:
                messagebox.showerror("could not measure the region", str(e)); return
            d["added"].append(rec)
            self.sel = ("added", len(d["added"]) - 1)
            self.refresh_tree(); self.fill_ctree(); self.draw()

        def delete_added(self):
            if self.data is None or not self.sel or self.sel[0] != "added":
                return
            self.data["added"].pop(self.sel[1]); self.sel = None
            self.refresh_tree(); self.fill_ctree(); self.draw()

        def save(self, quiet=False):
            if self.data is None:
                return
            if self.data.get("tune_only"):
                log("nothing to save: this run has no CFU result (tuning events only)"); return
            p = write_curation(self.data["tif"], self.data["rejected"],
                               [x for x in self.data["added"] if x.get("source") != "tune"],
                               notes=self.data.get("notes", ""))
            if not quiet:
                log(f"saved {p}")
                self.ax_img.set_xlabel(f"saved {p.name}", fontsize=8); self._bg = None; self.canvas.draw()

        def propagate_fov(self):
            self.save(quiet=True)
            exp = self.runs[self.i]["experiment"]
            group = [r for r in self.runs if r["experiment"] == exp]
            changed = propagate(group, log=log)
            messagebox.showinfo("propagated",
                                "\n".join(f"{Path(k).stem}: +{v} rejected" for k, v in changed.items())
                                or "nothing to copy")
            self.load(self.i)

    app = App()
    app.mainloop()
    return app


# ----------------------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", type=Path, help="mouse folder, or one binned tiff")
    ap.add_argument("--list", action="store_true", help="print what is there and exit")
    ap.add_argument("--reject", default=None, help="CFU numbers to mark as noise, e.g. 3,7,12")
    ap.add_argument("--keep", default=None, help="CFU numbers to mark as real again")
    ap.add_argument("--keep-all", action="store_true")
    ap.add_argument("--propagate", action="store_true",
                    help="copy each run's decisions to the matched CFUs of the other runs of its FOV")
    ap.add_argument("--iou", type=float, default=0.3)
    ap.add_argument("--sec", type=float, default=2.0)
    ap.add_argument("--fps", type=float, default=8.0, help="movie playback speed in the GUI (frames per second)")
    ap.add_argument("--tune", type=int, default=None, metavar="N",
                    help="overlay the events AQuA2 found with setting N of the tuning sweep (tab 3; aqua_tune\\<run>\\sNN\\real) "
                         "as magenta boxes T1.. with their onsets and traces; read-only")
    ap.add_argument("--tune-variant", default="real", choices=["real", "null", "inject"], help="which movie of the sweep (default real)")
    ap.add_argument("--um-per-px", type=float, default=0.792, help="pixel size, for the fibre mask behind the 'looks' column")
    ap.add_argument("--suggest", action="store_true", help="with --list: print per CFU what it looks like (fibre / off-fibre / blob) and exit")
    a = ap.parse_args()
    runs = discover_runs(a.target, sec=a.sec)
    if not runs:
        raise SystemExit(f"{a.target}: no runs found (run aqua_batch.py first, or pass a tiff)")

    if a.list:
        n_all = n_cur = 0
        for r in runs:
            d = load_run(r["tif"])
            if d is None:
                print(f"  [{r['condition']:>7}] {r['run']}: no _res_cfu.mat")
                continue
            n = d["cfu"]["n"]
            n_all += 1; n_cur += int(bool(d["curated"]))
            print(f"  [{r['condition']:>7}] {r['run']}: {n} CFUs, "
                  f"{n - len(d['rejected'])} kept, {len(d['added'])} added, "
                  f"events {sum(d['n_events'])}"
                  + (f"   (rules: {d['cfu']['rules']})" if d["cfu"].get("rules") else "")
                  + ("   (NOT curated under these rules: earlier decisions ignored — " + d["stale_reason"] + ")" if d["stale"]
                     else "" if d["curated"] else "   (not curated yet)"))
            if a.suggest and n:
                from afferent2p.cfu_shape import cfu_features, looks_label
                ref = d.get("mean_raw") if d.get("mean_raw") is not None else d["image"]
                feats = cfu_features(cfu_masks(d["cfu"]["maps"]), ref, a.um_per_px)
                for k, f in enumerate(feats):
                    if f["suggest"] != "fibre":
                        print(f"        CFU {k + 1}: {looks_label(f)} — {f['reason']}")
        print(f"  {n_cur} of {n_all} runs with a CFU result are curated under the current CFU rules")
        return
    if a.propagate:
        propagate(runs, iou=a.iou)
        return
    if a.reject or a.keep or a.keep_all:
        for r in runs:
            d = load_run(r["tif"])
            if d is None:
                continue
            rej = set() if a.keep_all else set(d["rejected"])
            if a.reject:
                rej |= {int(x) for x in a.reject.replace(" ", "").split(",") if x}
            if a.keep:
                rej -= {int(x) for x in a.keep.replace(" ", "").split(",") if x}
            p = write_curation(r["tif"], rej, d["added"], notes="set from the command line")
            print(f"  {r['run']}: {d['cfu']['n'] - len(rej)}/{d['cfu']['n']} kept -> {p.name}")
        return
    launch_gui(runs, fps=a.fps, tune=a.tune, tune_variant=a.tune_variant, um_per_px=a.um_per_px,
               mouse_dir=a.target if Path(a.target).is_dir() else Path(a.target).parent.parent)


if __name__ == "__main__":
    main()
