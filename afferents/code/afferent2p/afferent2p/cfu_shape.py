"""Shape and position of CFUs — which units look like afferent fibre segments and
which look like something else (bright compact blobs, debris, cell bodies, off
the fibres). Used by the curation window (column "looks", "Select suggested")
and available to the analysis.

    feats = cfu_features(masks, mean_img, um_per_px)      # one dict per CFU
    feats[k]["suggest"]  in  {"fibre", "off-fibre", "blob"}

Nothing here decides anything by itself: the suggestions are shown to the
user, who marks the CFUs. Rules (all tunable):
  * on_fibre  = fraction of the CFU's footprint on the ridge-filter fibre mask
                (afferent2p.tune.injection_mask on the registered mean image,
                dilated 2 px); below OFF_FIBRE_MAX -> "off-fibre"
  * elongation = sqrt(ratio of the footprint's principal moments) (1 = round)
  * bright_ratio = mean image inside the footprint / median of a ring around it
  * "blob": compact (elongation < BLOB_ELONG_MAX), bright (bright_ratio >
    BLOB_BRIGHT_MIN) and mostly off the fibres (on_fibre < BLOB_ON_FIBRE_MAX)
    - a bright varicosity ON a fibre stays "fibre"
"""
from __future__ import annotations
import numpy as np

OFF_FIBRE_MAX = 0.25        # on-fibre fraction below this -> "off-fibre"
BLOB_ELONG_MAX = 1.8        # compact
BLOB_BRIGHT_MIN = 1.6       # brighter than its surroundings
BLOB_ON_FIBRE_MAX = 0.5


def fibre_mask_for(mean_img: np.ndarray, um_per_px: float, dilate_px: int = 2) -> tuple[np.ndarray, str]:
    """The fibre mask of a (registered) mean image, dilated a little so that a
    footprint hugging a fibre counts as on it."""
    from scipy.ndimage import binary_dilation
    from .tune import injection_mask
    mask, how = injection_mask(np.asarray(mean_img, np.float32), um_per_px, log=lambda *_: None)
    if dilate_px > 0 and mask.any():
        mask = binary_dilation(mask, iterations=int(dilate_px))
    return mask.astype(bool), how


def _elongation(mask: np.ndarray) -> float:
    pts = np.argwhere(mask)
    if len(pts) < 3:
        return 1.0
    c = np.cov(pts.T.astype(float))
    ev = np.linalg.eigvalsh(c)
    lo, hi = max(float(ev[0]), 1e-6), max(float(ev[-1]), 1e-6)
    return float(np.sqrt(hi / lo))


def cfu_features(masks: list[np.ndarray], mean_img: np.ndarray | None, um_per_px: float = 0.792,
                 fmask: np.ndarray | None = None) -> list[dict]:
    """One dict per CFU: area_px, area_um2, elongation, on_fibre, bright_ratio, suggest, reason."""
    from scipy.ndimage import binary_dilation
    out = []
    if mean_img is not None and fmask is None:
        try:
            fmask, _ = fibre_mask_for(mean_img, um_per_px)
        except Exception:
            fmask = None
    img = np.asarray(mean_img, np.float32) if mean_img is not None else None
    for m in masks:
        m = np.asarray(m, bool)
        area = int(m.sum())
        f = {"area_px": area, "area_um2": area * um_per_px ** 2, "elongation": _elongation(m) if area else 1.0,
             "on_fibre": np.nan, "bright_ratio": np.nan, "suggest": "fibre", "reason": ""}
        if area == 0:
            out.append(f); continue
        if fmask is not None and fmask.shape == m.shape:
            f["on_fibre"] = float((m & fmask).sum() / area)
        if img is not None and img.shape == m.shape:
            inner = float(img[m].mean())
            ring = binary_dilation(m, iterations=7) & ~binary_dilation(m, iterations=2)
            if ring.any():
                base = float(np.median(img[ring]))
                lo = float(np.percentile(img, 1))
                f["bright_ratio"] = (inner - lo) / max(base - lo, 1e-6)
        onf, el, br = f["on_fibre"], f["elongation"], f["bright_ratio"]
        if np.isfinite(onf) and onf < OFF_FIBRE_MAX:
            if np.isfinite(br) and br > BLOB_BRIGHT_MIN and el < BLOB_ELONG_MAX:
                f["suggest"], f["reason"] = "blob", f"compact bright spot off the fibres ({100 * onf:.0f} % on fibre, {br:.1f}x its surroundings)"
            else:
                f["suggest"], f["reason"] = "off-fibre", f"only {100 * onf:.0f} % of the footprint on the fibre mask"
        elif np.isfinite(onf) and onf < BLOB_ON_FIBRE_MAX and np.isfinite(br) and br > BLOB_BRIGHT_MIN and el < BLOB_ELONG_MAX:
            f["suggest"], f["reason"] = "blob", f"compact bright spot, {100 * onf:.0f} % on fibre, {br:.1f}x its surroundings"
        else:
            f["reason"] = (f"{100 * onf:.0f} % on fibre" if np.isfinite(onf) else "no fibre mask") + f", elongation {el:.1f}"
        out.append(f)
    return out


def looks_label(f: dict) -> str:
    """Short text for a list column: 'fibre 85%', 'off-fibre 5%', 'blob 12%'."""
    onf = f.get("on_fibre", np.nan)
    pct = f"{100 * onf:.0f}%" if np.isfinite(onf) else "?"
    return f"{f.get('suggest', 'fibre')} {pct}"
