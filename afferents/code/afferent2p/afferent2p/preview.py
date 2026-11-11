"""Fast raw preview: mean/max projection of the UNregistered sbx data.

This is the first thing to check on a new dataset — if axons are visible here,
the sbx reading is correct and any downstream problem is registration/analysis.
Also writes a small ImageJ-compatible sample tiff (plain TIFF, not BigTIFF).
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import tifffile
from sbxreader import sbx_memmap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .discover import RunSpec


def _norm(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, [1, 99.7])
    return np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1)


def preview_run(
    spec: RunSpec,
    out_dir: Path,
    functional_chan: int = 1,
    block: int = 300,
    sample_frames: int = 200,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    chan = functional_chan - 1
    mmap = sbx_memmap(str(spec.sbx_path))
    n = mmap.shape[0]

    # three blocks: start / middle / end
    starts = sorted({0, max(0, n // 2 - block // 2), max(0, n - block)})
    blocks = []
    for s in starts:
        blocks.append(np.array(mmap[s:min(s + block, n), 0, chan]).astype(np.float32))
    stack = np.concatenate(blocks, axis=0)
    mean_img = stack.mean(axis=0)
    max_img = np.percentile(stack, 99.5, axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    axes[0].imshow(_norm(mean_img), cmap="gray")
    axes[0].set_title(f"RAW mean of {stack.shape[0]} frames (no registration)")
    axes[1].imshow(_norm(max_img), cmap="gray")
    axes[1].set_title("RAW 99.5th percentile (bright events)")
    for ax in axes:
        ax.axis("off")
    fig.suptitle(f"{spec.name} [{spec.condition}] — raw data check: "
                 "axons visible here = sbx reading is correct")
    fig.tight_layout()
    png = out_dir / f"{spec.name}_raw_preview.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)

    # small plain-TIFF sample that vanilla ImageJ can open
    sample = np.array(mmap[:min(sample_frames, n), 0, chan]).astype(np.uint16)
    tif = out_dir / f"{spec.name}_raw_sample.tif"
    tifffile.imwrite(str(tif), sample, imagej=True)
    del mmap
    print(f"  preview: {png.name} and {tif.name} "
          f"(mean of frames {starts} blocks)")


def preview_all(specs: list[RunSpec], out_dir: Path, functional_chan: int = 1) -> None:
    print("Raw previews (pre-registration):")
    for s in specs:
        preview_run(s, out_dir, functional_chan=functional_chan)
