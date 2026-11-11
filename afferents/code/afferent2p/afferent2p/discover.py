"""Run discovery: find .sbx runs, read their metadata, label pre/post."""
from __future__ import annotations
import os
from dataclasses import dataclass, field
from pathlib import Path

from sbxreader import sbx_get_metadata


@dataclass
class RunSpec:
    sbx_path: Path
    condition: str          # 'pre' or 'post'
    meta: dict = field(default_factory=dict)

    @property
    def name(self) -> str:
        return self.sbx_path.stem

    @property
    def nframes(self) -> int:
        return int(self.meta["num_frames"])

    @property
    def fs(self) -> float:
        return float(self.meta["frame_rate"])


def _find_sbx(data_dir: Path) -> list[Path]:
    """All .sbx files under data_dir (recursive), naturally sorted by name."""
    files = sorted(data_dir.rglob("*.sbx"), key=lambda p: p.name)
    return files


def discover_runs(data_dir: str | os.PathLike, runs_cfg) -> list[RunSpec]:
    """Build the ordered list of runs.

    runs_cfg is either the string 'auto' or a list of
    {'sbx': <filename>, 'condition': 'pre'|'post'} dicts (chronological order).
    """
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        raise FileNotFoundError(f"data_dir does not exist: {data_dir}")

    if runs_cfg in (None, "auto"):
        sbx_files = _find_sbx(data_dir)
        if not sbx_files:
            raise FileNotFoundError(f"No .sbx files found under {data_dir}")
        if len(sbx_files) == 1:
            conditions = ["pre"]
            print("Single run found — running in single-condition mode "
                  "(no pre/post comparison).")
        elif len(sbx_files) == 2:
            conditions = ["pre", "post"]
        else:
            raise ValueError(
                f"Auto-discovery found {len(sbx_files)} runs: "
                f"{[f.name for f in sbx_files]}.\n"
                "With more than 2 runs, list them explicitly in the config "
                "under experiment.runs with pre/post labels."
            )
        specs = [RunSpec(p, c) for p, c in zip(sbx_files, conditions)]
    else:
        specs = []
        for entry in runs_cfg:
            name = entry["sbx"]
            cond = str(entry["condition"]).strip()
            if not cond or any(ch in cond for ch in " /\\"):
                raise ValueError(f"bad condition label {cond!r} (no spaces or slashes)")
            matches = [p for p in _find_sbx(data_dir) if p.name == name]
            if not matches:
                raise FileNotFoundError(f"{name} not found under {data_dir}")
            specs.append(RunSpec(matches[0], cond))

    # metadata (requires the sidecar .mat next to each .sbx)
    for spec in specs:
        mat = spec.sbx_path.with_suffix(".mat")
        if not mat.exists():
            raise FileNotFoundError(
                f"Scanbox metadata file missing: {mat}\n"
                "Each .sbx needs its .mat sidecar in the same folder."
            )
        spec.meta = sbx_get_metadata(str(spec.sbx_path))

    # consistency checks
    fss = [round(s.fs, 3) for s in specs]
    sizes = [tuple(int(v) for v in s.meta["frame_size"]) for s in specs]
    chans = [int(s.meta["num_channels"]) for s in specs]
    planes = [int(s.meta["num_planes"]) for s in specs]
    if len(set(sizes)) > 1:
        raise ValueError(f"Runs differ in frame size: {sizes}")
    # runs may differ in channel count (some have a red vascular-tracer channel);
    # only the functional (green) channel is required everywhere
    if len(set(chans)) > 1:
        print(f"NOTE: runs differ in channel count {chans} — analysis uses channel 1 "
              "(green) only; red is used for registration only when ALL runs have it")
    if len(set(planes)) > 1:
        raise ValueError(f"Runs differ in plane count: {planes}")
    if len(set(fss)) > 1:
        print(f"WARNING: runs differ slightly in frame rate: {fss} (using first)")
    conds = list(dict.fromkeys(s.condition for s in specs))
    if len(conds) < 2:
        print("NOTE: only one condition present — the pipeline runs normally but "
              "no comparison will be made.")
    else:
        print(f"Conditions: {conds} (baseline = first, or analysis.baseline_condition)")

    print("Discovered runs:")
    for s in specs:
        print(
            f"  [{s.condition:>4}] {s.sbx_path.name}: {s.nframes} frames, "
            f"{s.fs:.2f} Hz, {sizes[0][0]}x{sizes[0][1]} px, "
            f"{int(s.meta['num_channels'])} chan, {planes[0]} plane(s)"
        )
    return specs
