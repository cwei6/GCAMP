"""Export registered movies as per-run tiffs (AQuA2 / Fiji input).

Defaults are chosen for AQuA2/MATLAB compatibility:
  * uint16 (suite2p's registered binary is *signed* int16; some MATLAB readers
    and AQuA2 preprocessing assume unsigned data, and registration can leave
    small negative values near the edges)
  * optional temporal binning and edge cropping, because a full-rate 512x796
    x 28000-frame movie is ~23 GB, which MATLAB cannot comfortably hold.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import tifffile


def _prepare(block: np.ndarray, dtype: str, bin_t: int,
             crop: tuple[int, int, int, int]) -> np.ndarray:
    L, R, T, B = crop
    if any(crop):
        block = block[:, T:block.shape[1] - B if B else None,
                      L:block.shape[2] - R if R else None]
    if bin_t > 1:
        n = (block.shape[0] // bin_t) * bin_t
        if n == 0:
            return np.empty((0,) + block.shape[1:], dtype=dtype)
        block = block[:n].reshape(n // bin_t, bin_t, *block.shape[1:]).mean(axis=1)
    if dtype == "uint16":
        return np.clip(block, 0, 65535).astype(np.uint16)
    return block.astype(np.int16)


def export_registered_tiffs(
    plane0: Path,
    out_dir: Path,
    run_names: list[str],
    run_conditions: list[str],
    run_nframes: list[int],
    Ly: int,
    Lx: int,
    chunk_frames: int = 4000,
    overwrite: bool = False,
    dtype: str = "uint16",
    bin_t: int = 1,
    crop: tuple[int, int, int, int] = (0, 0, 0, 0),
    fs: float | None = None,
) -> list[Path]:
    """Split suite2p's registered binary back into one tiff per run.

    Output: <run>_<condition>_reg.tif — feed these to AQuA2.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    bin_path = plane0 / "data.bin"
    nbytes_frame = Ly * Lx * 2
    total = sum(run_nframes)
    if bin_path.stat().st_size < total * nbytes_frame:
        raise RuntimeError(
            f"{bin_path} smaller than expected ({bin_path.stat().st_size} bytes "
            f"< {total} frames x {nbytes_frame})")

    bin_t = max(1, int(bin_t))
    L, R, T, B = crop
    out_Ly, out_Lx = Ly - T - B, Lx - L - R
    mov = np.memmap(bin_path, dtype=np.int16, mode="r", shape=(total, Ly, Lx))
    paths, start = [], 0
    print("Exporting registered per-run tiffs (AQuA2 input):")
    if bin_t > 1 or any(crop) or dtype != "int16":
        note = f"    dtype={dtype}"
        if bin_t > 1:
            note += f", temporal bin x{bin_t}"
            if fs:
                note += f" ({fs:.2f} -> {fs / bin_t:.2f} Hz)"
        if any(crop):
            note += f", crop L{L} R{R} T{T} B{B} -> {out_Ly}x{out_Lx}"
        print(note)

    for name, cond, n in zip(run_names, run_conditions, run_nframes):
        out_path = out_dir / f"{name}_{cond}_reg.tif"
        paths.append(out_path)
        n_out = n // bin_t if bin_t > 1 else n
        gb = n_out * out_Ly * out_Lx * 2 / 1e9
        if out_path.exists() and not overwrite:
            print(f"  exists, skipping: {out_path.name}")
            start += n
            continue
        print(f"  {out_path.name}: {n_out} frames, ~{gb:.1f} GB")
        if gb > 4:
            print("    NOTE: >4 GB — written as BigTIFF. Vanilla ImageJ cannot open "
                  "it (use Fiji > Import > Bio-Formats, virtual stack). MATLAB/AQuA2 "
                  "will also need a lot of RAM; consider export.aqua_bin_t in the "
                  "config to shrink it.")
        # temporal binning must not straddle chunk edges
        step = (chunk_frames // bin_t) * bin_t if bin_t > 1 else chunk_frames
        step = max(step, bin_t)
        # preallocate one contiguous TIFF series and fill it: writing unequal
        # blocks with TiffWriter would create multiple series, and readers
        # (MATLAB/AQuA2, tifffile.imread) then see only the first one
        tmp_path = out_path.with_name(out_path.name + ".partial")
        np_dtype = np.uint16 if dtype == "uint16" else np.int16
        dst = tifffile.memmap(str(tmp_path), shape=(n_out, out_Ly, out_Lx),
                              dtype=np_dtype, bigtiff=True)
        w = 0
        for s in range(start, start + n, step):
            e = min(s + step, start + n)
            block = np.asarray(mov[s:e]).astype(np.float32) if bin_t > 1 \
                else np.asarray(mov[s:e])
            out = _prepare(block, dtype, bin_t, crop)
            if out.shape[0]:
                k = min(out.shape[0], n_out - w)
                dst[w:w + k] = out[:k]
                w += k
        dst.flush()
        del dst
        tmp_path.replace(out_path)
        # sidecar with what AQuA2 needs to know (frame rate after binning etc.)
        meta = {
            "source_run": name, "condition": cond,
            "frames_in": int(n), "frames_out": int(w),
            "bin_frames": int(bin_t),
            "fs_in_hz": float(fs) if fs else None,
            "fs_out_hz": float(fs / bin_t) if fs else None,
            "crop_LRTB_px": [int(v) for v in crop],
            "height": int(out_Ly), "width": int(out_Lx), "dtype": dtype,
            "binning": "mean of consecutive frames (same as MATLAB WriteSbxPlaneTif binT)",
        }
        out_path.with_suffix(".json").write_text(json.dumps(meta, indent=2),
                                                 encoding="utf-8")
        start += n
    del mov
    return paths
