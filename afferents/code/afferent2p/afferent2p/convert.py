"""Convert .sbx runs to int16 tiffs (suite2p input). Replaces UnpackSBX.

We convert to tiff ourselves (rather than using suite2p's sbx reader) so the
pipeline controls run order and channel selection explicitly, and the tiffs
double as a raw, viewable copy of each run.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import tifffile
from sbxreader import sbx_memmap

from .discover import RunSpec


def raw_tiff_name(spec: RunSpec, chan2: int | None = None) -> str:
    return f"{spec.name}_{spec.condition}_raw{'2ch' if chan2 else ''}.tif"


def sbx_run_to_tiff(
    spec: RunSpec,
    out_dir: Path,
    functional_chan: int = 1,
    plane: int = 0,
    chunk: int = 500,
    overwrite: bool = False,
    chan2: int | None = None,
    transform: dict | None = None,
) -> Path:
    """Write one .sbx run to a single int16 tiff: <run>_<condition>_raw.tif.

    With chan2 (1-based, e.g. 2 = red), frames are INTERLEAVED
    (functional, chan2, functional, chan2, ...) — suite2p's 2-channel tiff
    layout (nchannels=2) — as <run>_<condition>_raw2ch.tif.
    Data are divided by 2 (as suite2p does for sbx) so they fit in int16.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / raw_tiff_name(spec, chan2)
    n_per_frame = 2 if chan2 else 1
    tr = transform if (transform and transform.get("apply")) else None
    tr_tag = {"angle_deg": tr["angle_deg"], "dy": tr["dy"], "dx": tr["dx"]} if tr else None
    side = out_path.with_suffix(".json")
    if out_path.exists() and not overwrite:
        # verify the existing tiff is COMPLETE (an interrupted earlier run can
        # leave a truncated file, which would silently shorten the whole analysis)
        try:
            with tifffile.TiffFile(str(out_path)) as tf:
                npages = len(tf.pages)
        except Exception:
            npages = -1
        old_tag = None
        if side.exists():
            try:
                old_tag = json.loads(side.read_text()).get("transform")
            except Exception:
                old_tag = None
        if npages == spec.nframes * n_per_frame and old_tag == tr_tag:
            print(f"  complete ({npages} frames), skipping: {out_path.name}")
            return out_path
        if npages == spec.nframes * n_per_frame and old_tag != tr_tag:
            print(f"  pre-alignment transform changed — re-converting: {out_path.name}")
        else:
            print(f"  INCOMPLETE existing tiff ({npages}/{spec.nframes * n_per_frame} frames) — "
                  f"re-converting: {out_path.name}")

    mmap = sbx_memmap(str(spec.sbx_path))
    nframes = mmap.shape[0]
    chan = functional_chan - 1
    if chan2 and chan2 > mmap.shape[2]:
        raise ValueError(f"{spec.sbx_path.name} has {mmap.shape[2]} channel(s); channel {chan2} not present")
    print(f"  {spec.sbx_path.name} -> {out_path.name} ({nframes} frames"
          + (f", 2 channels interleaved: {functional_chan}=functional, {chan2}=alignment)" if chan2 else ")"))
    # write to a temp name and rename on success, so an interrupted conversion
    # can never be mistaken for a finished file
    # preallocate one contiguous TIFF series and fill it: writing unequal blocks
    # with TiffWriter would split the file into multiple series, and readers
    # (tifffile.imread, MATLAB) would then see only the first one
    if tr:
        from .prealign import apply_transform
        print(f"    applying pre-alignment: rotate {tr['angle_deg']:+.2f} deg, shift (dy={tr['dy']}, dx={tr['dx']})")
    tmp_path = out_path.with_name(out_path.name + ".partial")
    H, W = mmap.shape[3], mmap.shape[4]
    dst = tifffile.memmap(str(tmp_path), shape=(nframes * n_per_frame, H, W),
                          dtype=np.int16, bigtiff=True)
    for start in range(0, nframes, chunk):
        stop = min(start + chunk, nframes)
        # sbx_memmap yields NFRAMES x NPLANES x NCHAN x H x W, offset-corrected
        block = np.array(mmap[start:stop, plane, chan]) // 2
        if tr:
            block = apply_transform(block.astype(np.int16), tr["angle_deg"], tr["dy"], tr["dx"])
        if chan2:
            block2 = np.array(mmap[start:stop, plane, chan2 - 1]) // 2
            if tr:
                block2 = apply_transform(block2.astype(np.int16), tr["angle_deg"], tr["dy"], tr["dx"])
            dst[2 * start:2 * stop:2] = block.astype(np.int16)
            dst[2 * start + 1:2 * stop:2] = block2.astype(np.int16)
        else:
            dst[start:stop] = block.astype(np.int16)
    dst.flush()
    del dst
    del mmap
    tmp_path.replace(out_path)
    side.write_text(json.dumps({"run": spec.name, "frames": int(nframes), "channels": n_per_frame,
                                "transform": tr_tag}, indent=2), encoding="utf-8")
    return out_path


def convert_all(specs: list[RunSpec], out_dir: Path, functional_chan: int = 1,
                overwrite: bool = False, chan2: int | None = None,
                transforms: list[dict] | None = None) -> list[Path]:
    print("Converting sbx -> tiff:")
    transforms = transforms or [None] * len(specs)
    return [
        sbx_run_to_tiff(s, out_dir, functional_chan=functional_chan, overwrite=overwrite,
                        chan2=chan2, transform=tr)
        for s, tr in zip(specs, transforms)
    ]
