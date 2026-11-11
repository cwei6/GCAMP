"""The behaviour csv next to a binned tiff: running (quadrature) and pupil
(DeepLabCut eye csv) averaged over the SAME bins as the movie, so row i of
<run>_<tag>_behavior.csv = recorded frame i of <run>_<tag>.tif.

unpack_bin.py writes it while unpacking; `rebuild(tif)` writes it again for a
movie that already exists (a DLC file that arrived later, a movie unpacked
before the pupil columns existed) from the numbers in the movie's .json; the
.sbx itself is not needed (only the folder and the run name).

Columns: frame, time_s, velocity_mean, still_fraction (quadrature file
present), pupil_diam_px_mean, pupil_valid_fraction (DLC eye csv present).
"""
from __future__ import annotations
import csv
import json
import re
from pathlib import Path

import numpy as np

PUPIL_COLS = ("pupil_diam_px_mean", "pupil_valid_fraction")
LOCO_COLS = ("velocity_mean", "still_fraction")


def write_behavior_csv(sbx_path: Path, N: int, n_out: int, T_total: int, fs: float,
                       out_path: Path, log=print) -> Path | None:
    """Locomotion (quad) + pupil (DLC) binned with the SAME N-frame bins as the
    tiff, so row i of the CSV = frame i of the AQuA tiff. Returns the path
    written, None when neither behaviour file exists."""
    from .loco import loco_for_run
    from .pupil import pupil_for_run
    sbx_path = Path(sbx_path)
    T = N * n_out                       # frames actually covered by the tiff
    lc = loco_for_run(sbx_path, T_total, fs, quad_glob="*quad*", still_thresh=5.0, min_still_s=5.0)
    pu = pupil_for_run(sbx_path, T_total, fs)
    if lc.get("quad_file") is None and pu is None:
        log("  behaviour: no quadrature / DLC eye file found — no behaviour csv")
        return None

    def binned(v):
        v = np.asarray(v, float)[:T].reshape(n_out, N)
        with np.errstate(all="ignore"):
            return np.nanmean(v, axis=1)

    cols = {"frame": np.arange(n_out), "time_s": (np.arange(n_out) + 0.5) * N / fs}
    if lc.get("quad_file") is not None:
        cols["velocity_mean"] = binned(lc["velocity"])
        cols["still_fraction"] = binned(lc["still"].astype(float))
    if pu is not None:
        cols["pupil_diam_px_mean"] = binned(pu["diam_px"])
        cols["pupil_valid_fraction"] = binned(pu["valid"].astype(float))
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols.keys())
        for row in zip(*cols.values()):
            w.writerow([f"{v:.4f}" if isinstance(v, float) else v for v in row])
    log(f"  wrote {Path(out_path).name} (behaviour binned like the tiff: "
        f"{', '.join(k for k in cols if k not in ('frame', 'time_s'))})")
    return Path(out_path)


def behavior_files(tif_path: Path) -> dict:
    """What a rebuild would use for this movie: its json, the .sbx, the
    quadrature and DLC files (each None when missing)."""
    from .pupil import find_pupil_file
    tif_path = Path(tif_path)
    base = tif_path.with_suffix("")
    js = Path(str(base) + ".json")
    sbx = None
    if js.exists():
        try:
            src = json.loads(js.read_text(encoding="utf-8")).get("source_sbx")
            if src and Path(src).exists():
                sbx = Path(src)
        except (OSError, ValueError):
            pass
    if sbx is None:
        cands = [p for p in tif_path.parent.glob("*.sbx") if tif_path.name.startswith(p.stem)]
        cands = cands or sorted(tif_path.parent.glob("*.sbx"))
        sbx = cands[0] if cands else None
    if sbx is None:
        # the raw file is elsewhere (archived): the quadrature / DLC readers only
        # use the folder and the run name, so a name derived from the tiff serves
        from .aqua_io import BIN_RE
        stem = BIN_RE.sub("", tif_path.stem)
        stem = re.sub(r"_ch\d+$", "", stem)
        sbx = tif_path.parent / f"{stem}.sbx"
    q = [p for p in sbx.parent.glob("*quad*") if p.is_file() and p.suffix != ".sbx"]
    q = [p for p in q if sbx.stem in p.name] or q
    return {"json": js if js.exists() else None, "sbx": sbx, "sbx_exists": sbx.exists(),
            "quad": q[0] if q else None, "dlc": find_pupil_file(sbx),
            "csv": Path(str(base) + "_behavior.csv")}


def has_pupil(csv_path: Path) -> bool:
    """Does an existing behaviour csv carry the pupil columns?"""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return False
    with open(csv_path, newline="") as fh:
        head = next(csv.reader(fh), [])
    return all(c in head for c in PUPIL_COLS)


def rebuild(tif_path: Path, log=print) -> Path | None:
    """Write <movie>_behavior.csv again from the files next to the movie,
    without touching the movie. Returns the csv path, or None when it could
    not be built (no json/sbx, or no behaviour files at all)."""
    tif_path = Path(tif_path)
    bf = behavior_files(tif_path)
    if bf["json"] is None:
        log(f"  behaviour: {tif_path.name} has no .json next to it — unpack it again to get one")
        return None
    info = json.loads(bf["json"].read_text(encoding="utf-8"))
    if bf["quad"] is None and bf["dlc"] is None:
        log(f"  behaviour: no quadrature / DLC eye file next to {tif_path.name} — no behaviour csv")
        return None
    N = int(info.get("bin_frames") or 0)
    n_out = int(info.get("frames_out") or 0)
    T_total = int(info.get("frames_in") or 0)
    fs = float(info.get("fs_in_hz") or 0)
    if not (N and n_out and fs):
        log(f"  behaviour: {bf['json'].name} lacks bin_frames/frames_out/fs_in_hz — unpack again")
        return None
    if not T_total:
        T_total = N * n_out
    return write_behavior_csv(bf["sbx"], N, n_out, T_total, fs, bf["csv"], log=log)
