"""Manual-ROI workflow: draw ROIs in Fiji on the registered mean image,
then extract traces + run the same event analysis on them.

Phase 1 (no ROI file yet):
    python run_pipeline.py <config> --steps manualroi
  writes  <out>/manual_roi/mean_for_rois.tif  (and struct_for_rois.tif)
  -> open in Fiji, draw ROIs along fiber segments (polygon/freehand/oval),
     add each to the ROI Manager (press 't'), then ROI Manager: More >> Save
     as  RoiSet.zip  into the same manual_roi folder.

Phase 2 (RoiSet.zip present): the same command extracts traces from the
registered movie under each ROI (with ring-neuropil subtraction) and runs the
standard event analysis into <out>/analysis_manual/.

Alternative to RoiSet.zip: a label-mask tiff named roi_mask.tif in the same
folder (integer labels 1..N, 0 = background), same size as the movie.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import binary_dilation

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath


# --------------------------------------------------------------- phase 1

def write_roi_templates(plane0: Path, roi_dir: Path, total_frames: int,
                        Ly: int, Lx: int) -> None:
    roi_dir.mkdir(parents=True, exist_ok=True)
    mov = np.memmap(plane0 / "data.bin", dtype=np.int16, mode="r",
                    shape=(total_frames, Ly, Lx))
    step = max(1, total_frames // 3000)
    sub = np.asarray(mov[::step]).astype(np.float32)
    mean_img = sub.mean(axis=0)
    struct_img = np.percentile(sub, 10, axis=0)
    del mov

    def to16(img):
        lo, hi = np.percentile(img, [0.5, 99.9])
        return (np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1) * 65535).astype(np.uint16)

    tifffile.imwrite(str(roi_dir / "mean_for_rois.tif"), to16(mean_img), imagej=True)
    tifffile.imwrite(str(roi_dir / "struct_for_rois.tif"), to16(struct_img), imagej=True)
    print(f"""
Manual-ROI templates written to: {roi_dir}
  1. Open mean_for_rois.tif in Fiji (struct_for_rois.tif shows static structure).
  2. Draw an ROI along each visible fiber segment (polygon, freehand, or oval;
     avoid line tools) and press 't' to add it to the ROI Manager.
  3. In the ROI Manager: More >> Save...  ->  save as RoiSet.zip in this folder.
  4. Re-run:  python run_pipeline.py <config> --steps manualroi
""")


# --------------------------------------------------------------- ROI loading

def _rasterize_imagej_roi(roi, Ly: int, Lx: int) -> np.ndarray | None:
    import roifile
    t = roi.roitype
    mask = np.zeros((Ly, Lx), bool)
    if t in (roifile.ROI_TYPE.RECT,):
        mask[max(0, roi.top):min(Ly, roi.bottom), max(0, roi.left):min(Lx, roi.right)] = True
        return mask
    if t in (roifile.ROI_TYPE.OVAL,):
        cy, cx = (roi.top + roi.bottom) / 2, (roi.left + roi.right) / 2
        ry, rx = max((roi.bottom - roi.top) / 2, 0.5), max((roi.right - roi.left) / 2, 0.5)
        yy, xx = np.mgrid[0:Ly, 0:Lx]
        mask[((yy - cy) / ry) ** 2 + ((xx - cx) / rx) ** 2 <= 1] = True
        return mask
    if t in (roifile.ROI_TYPE.POLYGON, roifile.ROI_TYPE.FREEHAND, roifile.ROI_TYPE.TRACED):
        coords = roi.coordinates()  # (N, 2) as x, y
        if coords is None or len(coords) < 3:
            return None
        path = MplPath(coords)
        yy, xx = np.mgrid[0:Ly, 0:Lx]
        pts = np.column_stack([xx.ravel(), yy.ravel()])
        mask = path.contains_points(pts).reshape(Ly, Lx)
        # include the outline itself so 1-px-thin ROIs are not empty
        ix = np.clip(coords[:, 0].astype(int), 0, Lx - 1)
        iy = np.clip(coords[:, 1].astype(int), 0, Ly - 1)
        mask[iy, ix] = True
        return mask
    print(f"  manual-roi: WARNING skipping unsupported ROI type {t} ({roi.name}) — "
          "use polygon/freehand/oval")
    return None


def load_rois(roi_dir: Path, Ly: int, Lx: int, explicit: str | None = None):
    """Returns (masks, names). Looks for RoiSet.zip / *.roi / roi_mask*.tif."""
    candidates = []
    if explicit:
        candidates = [Path(explicit)]
    else:
        candidates = (sorted(roi_dir.glob("RoiSet*.zip")) + sorted(roi_dir.glob("*.zip"))
                      + sorted(roi_dir.glob("roi_mask*.tif")) + sorted(roi_dir.glob("*.roi")))
    if not candidates:
        return None, None
    src = candidates[0]
    masks, names = [], []
    if src.suffix in (".zip", ".roi"):
        try:
            import roifile
        except ImportError:
            raise SystemExit("Reading Fiji ROI files needs the 'roifile' package:\n"
                             "  python -m pip install roifile")
        rois = roifile.roiread(str(src))
        if not isinstance(rois, list):
            rois = [rois]
        for r in rois:
            m = _rasterize_imagej_roi(r, Ly, Lx)
            if m is not None and m.sum() > 0:
                masks.append(m)
                names.append(r.name or f"roi{len(masks)}")
    else:  # label mask tiff
        lab = tifffile.imread(str(src))
        if lab.shape != (Ly, Lx):
            raise SystemExit(f"{src.name} is {lab.shape}, expected ({Ly},{Lx}) — "
                             "draw on mean_for_rois.tif")
        for v in np.unique(lab):
            if v == 0:
                continue
            masks.append(lab == v)
            names.append(f"label{int(v)}")
    print(f"  manual-roi: loaded {len(masks)} ROIs from {src.name}")
    return masks, names


# --------------------------------------------------------------- phase 2

def extract_manual_traces(plane0: Path, roi_dir: Path, masks, names,
                          total_frames: int, Ly: int, Lx: int,
                          ring_px: int = 8, chunk: int = 2000) -> None:
    n = len(masks)
    any_roi = np.any(masks, axis=0)
    rings = []
    for m in masks:
        ring = binary_dilation(m, iterations=ring_px) & ~any_roi
        if ring.sum() < 20:  # fallback if ROIs cover the neighborhood
            ring = ~any_roi
        rings.append(ring)

    mov = np.memmap(plane0 / "data.bin", dtype=np.int16, mode="r",
                    shape=(total_frames, Ly, Lx))
    F = np.zeros((n, total_frames), np.float32)
    Fneu = np.zeros((n, total_frames), np.float32)
    flat_masks = [np.flatnonzero(m.ravel()) for m in masks]
    flat_rings = [np.flatnonzero(r.ravel()) for r in rings]
    for s in range(0, total_frames, chunk):
        e = min(s + chunk, total_frames)
        arr = np.asarray(mov[s:e]).reshape(e - s, -1).astype(np.float32)
        for i in range(n):
            F[i, s:e] = arr[:, flat_masks[i]].mean(axis=1)
            Fneu[i, s:e] = arr[:, flat_rings[i]].mean(axis=1)
    del mov

    np.save(roi_dir / "F.npy", F)
    np.save(roi_dir / "Fneu.npy", Fneu)
    iscell = np.ones((n, 2), np.float32)
    np.save(roi_dir / "iscell.npy", iscell)
    (roi_dir / "roi_names.txt").write_text("\n".join(names), encoding="utf-8")
    # suite2p-style stat so the analysis can draw afferent maps for manual ROIs
    stat = []
    for m in masks:
        ys, xs = np.nonzero(m)
        stat.append({"ypix": ys.astype(np.int32), "xpix": xs.astype(np.int32),
                     "lam": np.ones(ys.size, np.float32) / ys.size,
                     "npix": int(ys.size), "med": [float(ys.mean()), float(xs.mean())]})
    np.save(roi_dir / "stat.npy", np.array(stat, dtype=object), allow_pickle=True)

    # overlay for verification
    step = max(1, total_frames // 1500)
    mean_img = np.asarray(np.memmap(plane0 / "data.bin", dtype=np.int16, mode="r",
                                    shape=(total_frames, Ly, Lx))[::step]).astype(np.float32).mean(0)
    lo, hi = np.percentile(mean_img, [0.5, 99.9])
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.imshow(np.clip((mean_img - lo) / max(hi - lo, 1e-9), 0, 1), cmap="gray")
    cmap = plt.get_cmap("tab20")
    for i, m in enumerate(masks):
        ys, xs = np.nonzero(m)
        ax.scatter(xs, ys, s=0.3, color=cmap(i % 20), alpha=0.5)
        ax.text(xs.mean(), ys.mean(), names[i], color=cmap(i % 20), fontsize=8)
    ax.set_title(f"{len(masks)} manual ROIs on registered mean image")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(roi_dir / "roi_overlay.png", dpi=150)
    plt.close(fig)
    print(f"  manual-roi: traces extracted ({n} ROIs x {total_frames} frames); "
          f"check {roi_dir / 'roi_overlay.png'}")
