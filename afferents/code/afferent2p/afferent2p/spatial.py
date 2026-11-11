"""Where the activity is, and whether it moved.

An afferent does not respond at the same spot in every recording (a different
functional unit of the same fibre can be the active one), and the field of
view is never exactly the same twice. So the comparison between conditions is
made on three things that do not depend on matching hotspots one to one:

  * the COMMON REGION — every run of a FOV is aligned to the first run
    (rotation + shift of its mean image, afferent2p.prealign) and only the
    part of the field present in all runs is compared;
  * the FIBRE ATLAS — the mean images of all runs are averaged after
    alignment and the fibres (bright ridges) are segmented from that average;
    every event is assigned to a tile of a coarse grid restricted to fibre
    pixels, and to the fibre SEGMENT it lies on. Tiles are large enough
    (25 um) that a hotspot moving along a fibre by a few um stays in place;
    segments are the afferent-level unit: the same segment in every run
    whatever part of it was active;
  * MAP COMPARISONS — per condition an events/min map on the tiles; between
    conditions the tile-wise ratio, the correlation of the two maps (how much
    of the spatial pattern is kept), and the fraction of active tiles that
    switched (active in one condition only) — 'the afferents respond
    elsewhere' made into numbers.

Inputs per run (written by run_cfu_batch.m next to the AQuA2 result):
  <base>_AQuA2_events.csv    one row per event: t0_frame, t1_frame, y_px, x_px
                             (1-based centre), n_pixels, box, dff_max, area_um2
  <base>_AQuA2_meanimg.tif   mean of the registered movie (datPro of the CFU
                             file is used when the tif is missing)
"""
from __future__ import annotations
from pathlib import Path

import numpy as np

from .aqua_io import AQUA_SUFFIX


# ----------------------------------------------------------------------------- files

def event_files(tif_path: Path) -> dict:
    b = str(Path(tif_path).with_suffix("")) + AQUA_SUFFIX
    return {"events": Path(b + "_events.csv"), "eventmap": Path(b + "_eventmap.tif"),
            "meanimg": Path(b + "_meanimg.tif")}


def read_events(tif_path: Path):
    """The per-event location table as a DataFrame with 0-based float y, x
    (columns y, x), or None when the export is missing."""
    import pandas as pd
    f = event_files(tif_path)["events"]
    if not f.exists():
        return None
    ev = pd.read_csv(f)
    if len(ev) == 0:                                   # a result without events
        return empty_events()
    if "y_px" not in ev.columns:
        return ev.assign(y=np.nan, x=np.nan)
    ev["y"] = ev["y_px"].astype(float) - 1.0
    ev["x"] = ev["x_px"].astype(float) - 1.0
    return ev


def empty_events():
    """The per-event location table of a run without events."""
    import pandas as pd
    cols = ["event", "t0_frame", "t1_frame", "y_px", "x_px", "n_pixels", "y0", "y1", "x0", "x1",
            "dff_max", "area_um2", "y", "x"]
    return pd.DataFrame({c: pd.Series(dtype=float) for c in cols})


def read_meanimg(tif_path: Path, fallback=None):
    import tifffile
    f = event_files(tif_path)["meanimg"]
    if f.exists():
        return tifffile.imread(str(f)).astype(np.float32)
    return None if fallback is None else np.asarray(fallback, np.float32)


# ----------------------------------------------------------------------------- alignment

def transform_points(pts: np.ndarray, t: dict, shape) -> np.ndarray:
    """Event centres (N, 2) as (y, x) of a run -> the reference frame, with
    the transform prealign.estimate_transform found for that run's mean
    image (the inverse of what apply_transform does to an image)."""
    pts = np.asarray(pts, float)
    if pts.size == 0:
        return pts.reshape(0, 2)
    H, W = shape[:2]
    th = np.deg2rad(t.get("angle_deg", 0.0))
    c, s = np.cos(th), np.sin(th)
    M = np.array([[c, -s], [s, c]])                  # apply_transform: input = M(out - centre - shift) + centre
    centre = np.array([H / 2.0, W / 2.0])
    out = (M.T @ (pts - centre).T).T + centre + np.array([t.get("dy", 0), t.get("dx", 0)])
    return out


def warp_image(img: np.ndarray, t: dict, cval: float = 0.0) -> np.ndarray:
    from .prealign import apply_transform
    return apply_transform(img[None].astype(np.float32), t.get("angle_deg", 0.0), int(t.get("dy", 0)),
                           int(t.get("dx", 0)), cval=cval)[0]


def align_runs(images: list, max_angle_deg: float = 4.0, min_gain: float = 0.05, good_match: float = 0.3,
               log=print):
    """Transforms of every run's mean image onto the first run's. A transform
    is used only when it clearly improves the image match (else the identity);
    a run counts as aligned when its match to the first run ends up >=
    good_match, whatever the transform."""
    from .prealign import estimate_transform
    ref = images[0]
    out = [{"angle_deg": 0.0, "dy": 0, "dx": 0, "ncc": 1.0, "ncc_before": 1.0, "aligned": True, "moved": False}]
    for i, im in enumerate(images[1:], 1):
        if im is None or im.shape != ref.shape:
            out.append({"angle_deg": 0.0, "dy": 0, "dx": 0, "ncc": np.nan, "ncc_before": np.nan, "aligned": False,
                        "moved": False})
            continue
        t = estimate_transform(ref, im, max_angle_deg=max_angle_deg)
        gain = t["ncc"] - t["ncc_before"]
        if gain <= min_gain:                       # no transform worth applying
            t = {"angle_deg": 0.0, "dy": 0, "dx": 0, "ncc": t["ncc_before"], "ncc_before": t["ncc_before"]}
            t["moved"] = False
        else:
            t["moved"] = True
        t["aligned"] = bool(t["ncc"] >= good_match)  # the field is the same one (identity or transformed)
        out.append(t)
    return out


def common_region(shape, transforms: list, margin_px: int = 8, covers: list | None = None) -> np.ndarray:
    """Pixels of the reference frame that every run covers after alignment
    (eroded by margin_px). covers: per run the pixels that hold real image
    (None = the whole frame)."""
    from scipy.ndimage import binary_erosion
    H, W = shape[:2]
    ones = np.ones((H, W), np.float32)
    reg = np.ones((H, W), bool)
    for i, t in enumerate(transforms):
        src = ones if covers is None or covers[i] is None else np.asarray(covers[i], np.float32)
        w = warp_image(src, t, cval=0.0)
        reg &= w > 0.5
    if margin_px > 0:
        reg = binary_erosion(reg, iterations=margin_px)
    return reg


# ----------------------------------------------------------------------------- fibres

def fibre_mask(img: np.ndarray, um_per_px: float, sigmas_um=(0.8, 1.6, 2.4), z_thr: float = 3.0,
               min_area_um2: float = 20.0) -> dict:
    """Bright ridges of the mean image -> fibre mask, skeleton and segments.
    Returns dict(mask, skeleton, segments (label image, 0 = none), n_segments,
    ridge (the filter response), segment_length_px (1-based list))."""
    from scipy.ndimage import gaussian_filter, label, distance_transform_edt, binary_closing
    from skimage.filters import sato
    from skimage.morphology import skeletonize, remove_small_objects
    im = np.asarray(img, np.float32)
    im = im - gaussian_filter(im, 25)                                # flatten vignetting
    im = im / (np.percentile(np.abs(im), 99) + 1e-9)
    sig = [max(0.6, s / um_per_px) for s in sigmas_um]
    ridge = sato(im, sigmas=sig, black_ridges=False)
    med = np.median(ridge); mad = np.median(np.abs(ridge - med)) * 1.4826 + 1e-9
    z = (ridge - med) / mad
    mask = z > z_thr
    mask &= im > np.percentile(im, 50)                               # and brighter than the background
    mask = binary_closing(mask, iterations=1)
    min_px = max(8, int(min_area_um2 / um_per_px ** 2))
    try:                                           # scikit-image >= 0.24 (min_size deprecated)
        mask = remove_small_objects(mask, max_size=min_px - 1)
    except TypeError:                              # older scikit-image
        mask = remove_small_objects(mask, min_size=min_px)
    skel = skeletonize(mask)
    # branches: drop junction pixels (more than 2 neighbours), label the rest
    from scipy.ndimage import convolve
    nb = convolve(skel.astype(int), np.ones((3, 3), int), mode="constant") - skel.astype(int)
    branches = skel & (nb <= 2)
    lab, n = label(branches, structure=np.ones((3, 3)))
    lab = _merge_through_junctions(lab, skel & ~branches)
    sizes = np.bincount(lab.ravel())
    min_len = max(3, int(8 / um_per_px))                              # >= 8 um long
    keep = np.flatnonzero(sizes >= min_len); keep = keep[keep > 0]
    remap = np.zeros(sizes.size, int); remap[keep] = np.arange(1, keep.size + 1)
    seg_skel = remap[lab]
    # every fibre pixel belongs to the nearest kept branch
    if keep.size:
        _, (iy, ix) = distance_transform_edt(seg_skel == 0, return_indices=True)
        segments = np.where(mask, seg_skel[iy, ix], 0)
    else:
        segments = np.zeros_like(seg_skel)
    lengths = np.bincount(seg_skel.ravel(), minlength=keep.size + 1)[1:]
    return {"mask": mask, "skeleton": skel, "segments": segments, "n_segments": int(keep.size),
            "ridge": ridge, "segment_length_px": lengths}


def _branch_direction(lab: np.ndarray, b: int, near: tuple, n_px: int = 12):
    """Unit direction of branch b pointing AWAY from the point `near` (y, x),
    from its n_px pixels closest to that point."""
    ys, xs = np.nonzero(lab == b)
    d = (ys - near[0]) ** 2 + (xs - near[1]) ** 2
    order = np.argsort(d)[:n_px]
    if order.size < 2:
        return None
    py, px = ys[order].astype(float), xs[order].astype(float)
    # principal direction of those pixels, oriented away from the junction
    v = np.array([py.mean() - near[0], px.mean() - near[1]])
    c = np.cov(np.vstack([py, px]))
    w, vec = np.linalg.eigh(c)
    u = vec[:, np.argmax(w)]
    if np.dot(u, v) < 0:
        u = -u
    return u / (np.linalg.norm(u) + 1e-9)


def _merge_through_junctions(lab: np.ndarray, junctions: np.ndarray, max_turn_deg: float = 50.0) -> np.ndarray:
    """At every junction of the skeleton, the two branches that continue each
    other in a straight line (turn <= max_turn_deg) get the same label — a
    fibre that crosses another stays ONE segment, side branches stay separate."""
    from scipy.ndimage import label as nd_label, binary_dilation
    # a crossing of two fibres leaves a small cluster of junction pixels, often
    # TWO clusters a few px apart joined by a tiny bridge branch: the bridge and
    # both clusters are one junction, so the crossing has its four arms
    jl, nj = nd_label(binary_dilation(junctions, iterations=2), structure=np.ones((3, 3)))
    if nj == 0:
        return lab
    sizes = np.bincount(lab.ravel())
    bridge_px = max(6, int(round(14 / 0.8)))                        # ~14 um at the usual pixel size
    lab = lab.copy()
    jmask = jl > 0
    for b in np.flatnonzero((sizes > 0) & (sizes <= bridge_px)):
        if b == 0:
            continue
        ys, xs = np.nonzero(lab == b)
        near = np.unique(jl[np.clip(ys[:, None] + np.arange(-2, 3), 0, lab.shape[0] - 1),
                            np.clip(xs[:, None] + np.arange(-2, 3), 0, lab.shape[1] - 1)])
        near = near[near > 0]
        if near.size >= 2:                                           # joins two junctions: part of the crossing
            jmask[ys, xs] = True
            lab[ys, xs] = 0
    jl, nj = nd_label(binary_dilation(jmask, iterations=1), structure=np.ones((3, 3)))
    parent = {}

    def find(a):
        while parent.get(a, a) != a:
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    H, W = lab.shape
    for j in range(1, nj + 1):
        ys, xs = np.nonzero(jl == j)
        cy, cx = float(ys.mean()), float(xs.mean())
        # branches touching this junction (its dilated pixels)
        touch = set(int(v) for v in np.unique(lab[ys, xs]) if v > 0)
        touch = sorted(touch)
        if len(touch) < 2:
            continue
        dirs = {b: _branch_direction(lab, b, (cy, cx)) for b in touch}
        pairs = []
        for i, a in enumerate(touch):
            for b in touch[i + 1:]:
                if dirs[a] is None or dirs[b] is None:
                    continue
                cos = float(np.dot(dirs[a], dirs[b]))            # -1 = straight continuation
                turn = np.degrees(np.arccos(np.clip(-cos, -1, 1)))
                pairs.append((turn, a, b))
        used = set()
        for turn, a, b in sorted(pairs):
            if turn > max_turn_deg or a in used or b in used:
                continue
            union(a, b); used.update((a, b))
    if not parent:
        return lab
    out = lab.copy()
    ids = np.unique(lab); ids = ids[ids > 0]
    for i in ids:
        r = find(int(i))
        if r != i:
            out[lab == i] = r
    return out


def assign_segments(pts: np.ndarray, segments: np.ndarray, max_dist_px: float) -> np.ndarray:
    """Segment id (0 = none) of each event centre: the segment at the pixel,
    else the nearest one within max_dist_px."""
    from scipy.ndimage import distance_transform_edt
    if pts.size == 0:
        return np.zeros(0, int)
    H, W = segments.shape
    d, (iy, ix) = distance_transform_edt(segments == 0, return_indices=True)
    yy = np.clip(np.round(pts[:, 0]).astype(int), 0, H - 1)
    xx = np.clip(np.round(pts[:, 1]).astype(int), 0, W - 1)
    near = segments[iy[yy, xx], ix[yy, xx]]
    return np.where(d[yy, xx] <= max_dist_px, near, 0)


# ----------------------------------------------------------------------------- tiles

def tile_grid(shape, tile_px: int):
    H, W = shape[:2]
    ny, nx = int(np.ceil(H / tile_px)), int(np.ceil(W / tile_px))
    return ny, nx


def tile_of(pts: np.ndarray, tile_px: int, shape):
    ny, nx = tile_grid(shape, tile_px)
    ty = np.clip((pts[:, 0] // tile_px).astype(int), 0, ny - 1)
    tx = np.clip((pts[:, 1] // tile_px).astype(int), 0, nx - 1)
    return ty, tx


def tile_fractions(mask: np.ndarray, tile_px: int) -> np.ndarray:
    """Fraction of each tile covered by a boolean mask."""
    H, W = mask.shape
    ny, nx = tile_grid(mask.shape, tile_px)
    pad = np.zeros((ny * tile_px, nx * tile_px), np.float32)
    pad[:H, :W] = mask
    return pad.reshape(ny, tile_px, nx, tile_px).mean(axis=(1, 3))


def tile_counts(pts: np.ndarray, tile_px: int, shape) -> np.ndarray:
    ny, nx = tile_grid(shape, tile_px)
    out = np.zeros((ny, nx), float)
    if pts.size:
        ty, tx = tile_of(pts, tile_px, shape)
        np.add.at(out, (ty, tx), 1.0)
    return out


# ----------------------------------------------------------------------------- comparisons

def compare_tiles(rate_a: np.ndarray, rate_b: np.ndarray, valid: np.ndarray, n_boot: int = 500,
                  block: int = 3, rng=None) -> dict:
    """Two events/min tile maps (a = baseline, b = condition) on the same
    valid tiles: pooled ratio inside the maps, Spearman correlation of the
    maps, active-tile fractions and the switching fraction, and a spatial
    block bootstrap CI (blocks of block x block tiles) of the mean tile-wise
    log2 ratio (tiles active in either condition)."""
    from scipy.stats import spearmanr
    rng = rng or np.random.default_rng(0)
    v = valid.astype(bool)
    a, b = rate_a[v], rate_b[v]
    out = {"n_tiles": int(v.sum())}
    if v.sum() < 4:
        return out
    out["map_r"] = float(spearmanr(a, b).correlation) if (a.std() > 0 and b.std() > 0) else np.nan
    act_a, act_b = a > 0, b > 0
    out["active_tiles_a"] = int(act_a.sum()); out["active_tiles_b"] = int(act_b.sum())
    either = act_a | act_b; both = act_a & act_b
    out["frac_active_a"] = float(act_a.mean()); out["frac_active_b"] = float(act_b.mean())
    out["frac_switched"] = float((either & ~both).sum() / max(either.sum(), 1))
    out["frac_new"] = float((act_b & ~act_a).sum() / max(either.sum(), 1))
    out["frac_lost"] = float((act_a & ~act_b).sum() / max(either.sum(), 1))
    # tile-wise log2 ratio with a floor of half the smallest non-zero rate
    floor = 0.5 * min(np.min(a[a > 0]) if act_a.any() else 1.0, np.min(b[b > 0]) if act_b.any() else 1.0)
    lr = np.log2((b + floor) / (a + floor))
    out["tile_log2_ratio_mean"] = float(lr[either].mean()) if either.any() else np.nan
    # spatial block bootstrap over the valid tiles' block ids
    ys, xs = np.nonzero(v)
    blk = (ys // block) * 10000 + (xs // block)
    ids = np.unique(blk)
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(ids, ids.size, replace=True)
        sel = np.concatenate([np.flatnonzero(blk == i) for i in pick])
        e = either[sel]
        if e.any():
            boots.append(lr[sel][e].mean())
    if boots:
        out["tile_log2_ratio_ci"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
    return out


def segment_table(seg_ids_by_cond: dict, minutes_by_cond: dict, lengths_px: np.ndarray, um_per_px: float,
                  min_events: int = 3, behaviour_factor: float = 1.0):
    """Per-segment events/min per condition -> DataFrame (segment, length_um,
    n_<cond>, rate_<cond> ...) and a summary dict: how many segments are up /
    down beyond what the behaviour explains (binomial test on the split of a
    segment's events between the conditions, expected split = time x
    behaviour_factor, the FOV-wide rate change the behaviour accounts for),
    and a Wilcoxon over segments."""
    import pandas as pd
    from scipy.stats import poisson, wilcoxon
    conds = list(seg_ids_by_cond)
    n_seg = lengths_px.size
    rows = []
    for s in range(1, n_seg + 1):
        r = {"segment": s, "length_um": float(lengths_px[s - 1] * um_per_px)}
        for c in conds:
            n = int(np.sum(seg_ids_by_cond[c] == s))
            r[f"n_{c}"] = n
            r[f"rate_{c}"] = n / minutes_by_cond[c] if minutes_by_cond[c] > 0 else np.nan
        rows.append(r)
    df = pd.DataFrame(rows)
    summ = {"n_segments": n_seg}
    if len(conds) >= 2 and n_seg:
        base, c = conds[0], conds[1]
        nb, nc = df[f"n_{base}"].to_numpy(), df[f"n_{c}"].to_numpy()
        mb, mc = minutes_by_cond[base], minutes_by_cond[c]
        used = (nb + nc) >= min_events
        # conditional test: given the total, the split is binomial with
        # p = mc*B / (mb + mc*B), B = the rate change the behaviour accounts for
        B = behaviour_factor if np.isfinite(behaviour_factor) and behaviour_factor > 0 else 1.0
        p_c = mc * B / max(mb + mc * B, 1e-9)
        from scipy.stats import binomtest
        up = down = 0
        for i in np.flatnonzero(used):
            pv = binomtest(int(nc[i]), int(nb[i] + nc[i]), p_c).pvalue
            if pv < 0.05:
                if nc[i] / max(mc * B, 1e-9) > nb[i] / max(mb, 1e-9):
                    up += 1
                else:
                    down += 1
        summ.update({"segments_used": int(used.sum()), "segments_up": up, "segments_down": down,
                     "segments_active_base": int((nb > 0).sum()), "segments_active_cond": int((nc > 0).sum())})
        rb, rc = df[f"rate_{base}"].to_numpy(), df[f"rate_{c}"].to_numpy()
        if used.sum() >= 5:
            try:
                summ["segments_wilcoxon_p"] = float(wilcoxon(rc[used], rb[used]).pvalue)
            except ValueError:
                summ["segments_wilcoxon_p"] = np.nan
            floor = 0.5 / max(mb, mc)
            summ["segments_log2_ratio_median"] = float(np.median(np.log2((rc[used] / B + floor) / (rb[used] + floor))))
        summ["behaviour_factor"] = B
    return df, summ


# ----------------------------------------------------------------------------- per experiment

def experiment_maps(exp_name: str, runs: list, conditions: list, um_per_px: float, out_dir: Path,
                    tile_um: float = 25.0, behaviour_factor: dict | None = None, log=print) -> dict | None:
    """runs: [{run, condition, tif, events (DataFrame with y, x, t0_frame),
    minutes_used, meanimg (2D array or None), excluded (bool per recorded
    frame or None)}] in the order of the FOV. Writes tiles_<exp>.csv,
    segments_<exp>.csv, <exp>_maps.png; returns the summary dict."""
    import pandas as pd
    try:
        import skimage  # noqa: F401
    except ImportError:
        log("  maps: scikit-image is not installed in this Python — "
            f"install it (pip install scikit-image) to get the fibre atlas / tile maps; skipped")
        return None
    imgs = [r.get("meanimg") for r in runs]
    if any(im is None for im in imgs):
        log("  maps: no mean image for every run — skipped")
        return None
    # AQuA2's own registration crops every movie by its motion, so the runs'
    # images differ by a few pixels: pad each to one canvas (top-left anchored,
    # which keeps the event coordinates valid) and remember what each covers
    H = max(im.shape[0] for im in imgs); W = max(im.shape[1] for im in imgs)
    shape = (H, W)
    covers = []
    if any(im.shape != shape for im in imgs):
        log(f"  maps: the runs' images differ in size ({', '.join(f'{im.shape[0]}x{im.shape[1]}' for im in imgs)}; "
            f"AQuA2 crops each movie by its motion) — padded to {H}x{W}")
    for i, im in enumerate(imgs):
        cov = np.zeros(shape, np.float32)
        cov[:im.shape[0], :im.shape[1]] = 1.0
        covers.append(cov)
        if im.shape != shape:
            pad = np.full(shape, float(np.median(im)), np.float32)
            pad[:im.shape[0], :im.shape[1]] = im
            imgs[i] = pad
    transforms = align_runs(imgs, log=log)
    for r, t in zip(runs, transforms):
        r["transform"] = t
    n_al = sum(t["aligned"] for t in transforms[1:])
    n_moved = sum(t.get("moved", False) for t in transforms[1:])
    region = common_region(shape, transforms, covers=covers)
    atlas_img = np.mean([warp_image(im, t, cval=float(np.median(im))) for im, t in zip(imgs, transforms)], axis=0)
    fib = fibre_mask(atlas_img, um_per_px)
    fmask = fib["mask"] & region
    segs = np.where(region, fib["segments"], 0)
    tile_px = max(4, int(round(tile_um / um_per_px)))
    ffrac = tile_fractions(fmask, tile_px)
    rfrac = tile_fractions(region, tile_px)
    valid = (ffrac > 0.02) & (rfrac > 0.5)
    # events into the reference frame
    per_cond = {c: {"counts": np.zeros(ffrac.shape), "minutes": 0.0, "n_in": 0, "n_all": 0, "seg": []} for c in conditions}
    for r in runs:
        ev = r.get("events")
        c = r["condition"]
        if c not in per_cond:
            continue
        per_cond[c]["minutes"] += float(r.get("minutes_used", 0.0))
        if ev is None or len(ev) == 0:
            continue
        pts = ev[["y", "x"]].to_numpy(float)
        pts = pts[np.isfinite(pts).all(1)]
        p_ref = transform_points(pts, r["transform"], shape)
        yy = np.clip(np.round(p_ref[:, 0]).astype(int), 0, shape[0] - 1)
        xx = np.clip(np.round(p_ref[:, 1]).astype(int), 0, shape[1] - 1)
        inside = region[yy, xx]
        per_cond[c]["n_all"] += int(pts.shape[0]); per_cond[c]["n_in"] += int(inside.sum())
        per_cond[c]["counts"] += tile_counts(p_ref[inside], tile_px, shape)
        per_cond[c]["seg"].append(assign_segments(p_ref[inside], segs, max_dist_px=max(3, 4 / um_per_px)))
    rates = {c: (per_cond[c]["counts"] / per_cond[c]["minutes"] if per_cond[c]["minutes"] > 0 else per_cond[c]["counts"] * np.nan)
             for c in conditions}
    summary = {"n_runs": len(runs), "runs_aligned": n_al, "runs_not_aligned": len(runs) - 1 - n_al, "runs_moved": n_moved,
               "common_region_frac": float(region.mean()), "fibre_area_um2": float(fmask.sum() * um_per_px ** 2),
               "fibre_frac_of_region": float(fmask.sum() / max(region.sum(), 1)), "n_segments": fib["n_segments"],
               "tile_um": tile_um, "n_valid_tiles": int(valid.sum())}
    base = conditions[0]
    fibre_mm2 = summary["fibre_area_um2"] / 1e6
    for c in conditions:
        pc = per_cond[c]
        summary[f"minutes_{c}"] = pc["minutes"]; summary[f"events_in_region_{c}"] = pc["n_in"]
        summary[f"events_all_{c}"] = pc["n_all"]
        summary[f"rate_region_{c}"] = pc["n_in"] / pc["minutes"] if pc["minutes"] > 0 else np.nan
        summary[f"rate_per_fibre_mm2_{c}"] = summary[f"rate_region_{c}"] / fibre_mm2 if fibre_mm2 > 0 else np.nan
    comp = {}
    bf = behaviour_factor or {}
    for c in conditions[1:]:
        B = bf.get(c, 1.0)
        B = B if np.isfinite(B) and B > 0 else 1.0
        cmp_ = compare_tiles(rates[base], rates[c] / B, valid)      # the condition's map with the behaviour taken out
        rb, rc = summary[f"rate_region_{base}"], summary[f"rate_region_{c}"]
        cmp_["region_ratio"] = rc / rb if rb and np.isfinite(rb) and rb > 0 else np.nan
        cmp_["region_ratio_beyond_behaviour"] = cmp_["region_ratio"] / B if np.isfinite(cmp_["region_ratio"]) else np.nan
        cmp_["behaviour_factor"] = B
        comp[c] = cmp_
    # per-segment rates
    seg_ids = {c: (np.concatenate(per_cond[c]["seg"]) if per_cond[c]["seg"] else np.zeros(0, int)) for c in conditions}
    minutes = {c: per_cond[c]["minutes"] for c in conditions}
    seg_df, seg_summ = segment_table(seg_ids, minutes, fib["segment_length_px"], um_per_px,
                                     behaviour_factor=bf.get(conditions[1], 1.0) if len(conditions) > 1 else 1.0)
    seg_df.to_csv(out_dir / f"segments_{exp_name}.csv", index=False)
    # tiles csv
    ys, xs = np.nonzero(valid)
    tiles = pd.DataFrame({"tile_y": ys, "tile_x": xs, "y_um": (ys + 0.5) * tile_px * um_per_px,
                          "x_um": (xs + 0.5) * tile_px * um_per_px, "fibre_fraction": ffrac[ys, xs]})
    for c in conditions:
        tiles[f"rate_{c}"] = rates[c][ys, xs]
    tiles.to_csv(out_dir / f"tiles_{exp_name}.csv", index=False)
    summary["comparisons"] = comp
    summary["segments"] = seg_summ
    try:
        maps_figure(exp_name, runs, conditions, imgs, transforms, region, fib, fmask, rates, valid, tile_px,
                    um_per_px, seg_df, comp, out_dir / f"{exp_name}_maps.png")
    except Exception as e:                                                   # the numbers are what matters
        log(f"  maps figure failed ({e})")
    return summary


# ----------------------------------------------------------------------------- figure

COND_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def maps_figure(exp_name, runs, conditions, imgs, transforms, region, fib, fmask, rates, valid, tile_px,
                um_per_px, seg_df, comp, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    n_cond = len(conditions)
    colour = {c: COND_COLOURS[i % len(COND_COLOURS)] for i, c in enumerate(conditions)}
    from matplotlib.ticker import FuncFormatter
    H, W = imgs[0].shape
    n_show = min(len(runs), 4)                   # per-run panels in row 1 (all 4 runs of a pre/post day)
    ncols = max(4, n_cond + 2, n_show + 1)
    panel_h = 4.0 * H / W
    fig = plt.figure(figsize=(4.0 * ncols, 1.7 + 2 * panel_h * 1.12))
    gs = fig.add_gridspec(2, ncols, hspace=0.42, wspace=0.2, left=0.04, right=0.99, top=0.87, bottom=0.09)
    um = FuncFormatter(lambda v, _: f"{v * um_per_px:.0f}")

    def um_axes(ax):
        ax.xaxis.set_major_formatter(um); ax.yaxis.set_major_formatter(um)
        ax.tick_params(labelsize=7)
    # row 1: atlas with fibre mask + common region, then per-run aligned images (up to 4)
    ax = fig.add_subplot(gs[0, 0])
    atlas = np.mean([warp_image(im, t, cval=float(np.median(im))) for im, t in zip(imgs, transforms)], axis=0)
    lo, hi = np.percentile(atlas, [1, 99.7])
    ax.imshow(atlas, cmap="gray", vmin=lo, vmax=hi)
    ax.contour(fmask.astype(float), levels=[0.5], colors=["#eda100"], linewidths=0.6)
    ax.contour(region.astype(float), levels=[0.5], colors=["#ffffff"], linewidths=1.0, linestyles="--")
    ax.set_title(f"mean image, all runs aligned — fibres yellow ({fib['n_segments']} segments),\n"
                 f"common region dashed; axes in um", fontsize=7, loc="left")
    um_axes(ax)
    for k, (r, im, t) in enumerate(zip(runs, imgs, transforms)):
        if k >= n_show:
            break
        ax = fig.add_subplot(gs[0, k + 1])
        w = warp_image(im, t, cval=float(np.median(im)))
        lo, hi = np.percentile(w, [1, 99.7])
        ax.imshow(w, cmap="gray", vmin=lo, vmax=hi)
        ev = r.get("events")
        if ev is not None and len(ev):
            pts = ev[["y", "x"]].to_numpy(float); pts = pts[np.isfinite(pts).all(1)]
            p = transform_points(pts, t, (H, W))
            ax.plot(p[:, 1], p[:, 0], ".", ms=2.5, color=colour.get(r["condition"], "#e34948"), alpha=0.7)
        if k == 0:
            al = "reference run"
        elif t.get("moved"):
            al = f"moved {t['angle_deg']:+.1f} deg, ({t['dy']}, {t['dx']}) px; match {t['ncc_before']:.2f} -> {t['ncc']:.2f}"
        elif t["aligned"]:
            al = f"same position (match {t['ncc']:.2f})"
        else:
            al = f"NOT ALIGNED (match {t['ncc']:.2f})"
        ax.set_title(f"{r['run'][-24:]} [{r['condition']}] — {al}\nevents at their centres", fontsize=7, loc="left")
        um_axes(ax)
    # row 2: rate maps per condition, change map, tile scatter, segment scatter
    vmax = max(np.nanmax(rates[c][valid]) if valid.any() else 0 for c in conditions) or 1.0
    for i, c in enumerate(conditions):
        ax = fig.add_subplot(gs[1, i])
        m = np.where(valid, rates[c], np.nan)
        ax.imshow(m, cmap="Blues", vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_title(f"{c}: events/min per {tile_px * um_per_px:.0f}-um tile\n(fibre tiles in the common region; max {vmax:.1f})",
                     fontsize=7, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
    base = conditions[0]
    col = n_cond
    if n_cond >= 2:
        c = conditions[1]
        ax = fig.add_subplot(gs[1, col]); col += 1
        floor = 0.05
        chg = np.where(valid & ((rates[base] > 0) | (rates[c] > 0)), np.log2((rates[c] + floor) / (rates[base] + floor)), np.nan)
        lim = np.nanmax(np.abs(chg)) if np.isfinite(chg).any() else 1
        lim = max(float(lim), 1e-3)                       # identical maps: no change at all
        ax.imshow(chg, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), interpolation="nearest")
        cm = comp.get(c, {})
        ax.set_title(f"log2 {c}/{base} per tile (red = more {c}, +/-{lim:.1f})\nmap r {cm.get('map_r', np.nan):+.2f}, "
                     f"{100 * cm.get('frac_switched', np.nan):.0f} % of active tiles active in one condition only",
                     fontsize=7, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        # segments scatter
        ax = fig.add_subplot(gs[1, col]); col += 1
        B = comp.get(c, {}).get("behaviour_factor", 1.0)
        if len(seg_df):
            rb, rc = seg_df[f"rate_{base}"].to_numpy(), seg_df[f"rate_{c}"].to_numpy() / B
            top = max(np.nanmax(rb), np.nanmax(rc), 0.1)
            ax.plot([0, top], [0, top], color="#999999", lw=0.8)
            ax.plot(rb, rc, "o", ms=4, color="#2a2a2a", alpha=0.6, mec="white", mew=0.5)
            ax.set_xlabel(f"events/min per fibre segment — {base}", fontsize=8)
            ax.set_ylabel(f"— {c}" + (f" (behaviour factor x{B:.2f} taken out)" if abs(B - 1) > 1e-6 else ""), fontsize=8)
        ax.set_title("per fibre segment (the afferent-level unit); line = no change", fontsize=7, loc="left")
        ax.tick_params(labelsize=7)
    fig.suptitle(f"{exp_name} — where the activity is", fontsize=10, x=0.04, y=0.98, ha="left")
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
