#!/usr/bin/env python3
"""afferent2p — config-driven pre/post GCaMP afferent pipeline.

Usage:
    python run_pipeline.py config.yaml
    python run_pipeline.py config.yaml --steps convert,register
    python run_pipeline.py config.yaml --steps analyze     # e.g. after GUI curation

Steps (in order): preview, convert, register, export, loco, qc, analyze
Extra steps (not run by default):
  curate    — re-run only the ROI curation (afferent filter) with the current
              analysis.curate_* settings; writes iscell.npy + analysis/roi_curation.*
              (fast — use it to tune the filter on real data, then --steps analyze)
  manualroi — Fiji-drawn ROI extraction; see afferent2p/manual_roi.py docstring.
  compare   — compare the registered movie with a MATLAB-pipeline output:
              python run_pipeline.py <config> --steps compare --ref <matlab tif>
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

try:
    import numpy as np
    import yaml
    import sbxreader  # noqa: F401
    import suite2p    # noqa: F401
    import tifffile   # noqa: F401
    import scipy      # noqa: F401
    import matplotlib # noqa: F401
except ImportError as e:
    sys.exit(
        f"\nMissing Python package: {e.name}\n"
        "Install the pipeline's dependencies into THIS Python "
        f"({sys.executable}) with:\n\n"
        f'  "{sys.executable}" -m pip install suite2p sbxreader tifffile pyyaml scipy matplotlib\n'
    )

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.discover import discover_runs
from afferent2p.convert import convert_all
from afferent2p.register import run_suite2p, load_plane_ops
from afferent2p.export_aqua import export_registered_tiffs
from afferent2p.loco import loco_for_run
from afferent2p.pupil import pupil_for_run
from afferent2p.qc import qc_report
from afferent2p.analyze import analyze_prepost

ALL_STEPS = ["preview", "convert", "register", "export", "loco", "qc", "analyze"]


def registered_mean(plane0: Path, total: int) -> "np.ndarray | None":
    """Subsampled mean of the registered movie (for ROI/afferent maps)."""
    try:
        ops = load_plane_ops(plane0)
        Ly, Lx = int(ops["Ly"]), int(ops["Lx"])
        mov = np.memmap(plane0 / "data.bin", dtype=np.int16, mode="r", shape=(total, Ly, Lx))
        step = max(1, total // 2000)
        img = np.asarray(mov[::step]).astype(np.float32).mean(axis=0)
        del mov
        return img
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config", help="YAML config file")
    ap.add_argument("--steps", default=",".join(ALL_STEPS),
                    help=f"comma-separated subset of: {','.join(ALL_STEPS)}")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--ref", default=None,
                    help="for --steps compare: a MATLAB-pipeline projection or "
                         "registered movie tif to compare against")
    ap.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE",
                    help="override a config value, e.g. --set analysis.curate_min_autocorr=0.1 "
                         "(repeatable; numbers/true/false/auto parsed)")
    ap.add_argument("--default-curate", action="store_true",
                    help="ignore every analysis.curate_* value in the config and use the "
                         "pipeline's current defaults (configs written by older versions "
                         "carry old thresholds)")
    args = ap.parse_args()
    steps = [s.strip() for s in args.steps.split(",")]

    cfg = yaml.safe_load(Path(args.config).read_text())
    if args.default_curate:
        an = cfg.setdefault("analysis", {})
        dropped = [k for k in list(an) if k.startswith("curate_")]
        for k in dropped:
            an.pop(k)
        print(f"--default-curate: ignoring {len(dropped)} curate_* value(s) from the config")
    for item in args.set:
        key, _, val = item.partition("=")
        sec, _, k = key.partition(".")
        v = yaml.safe_load(val)
        cfg.setdefault(sec, {})[k] = v
        print(f"--set {sec}.{k} = {v!r}")
    exp = cfg["experiment"]
    data_dir = Path(exp["data_dir"])
    out_dir = Path(exp.get("output_dir") or (data_dir / "afferent2p_out"))
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "config_used.yaml").write_text(yaml.safe_dump(cfg))

    # ---- discover (always) ----
    specs = discover_runs(data_dir, exp.get("runs", "auto"))
    fs = cfg.get("imaging", {}).get("fs") or specs[0].fs
    run_names = [s.name for s in specs]
    run_conditions = [s.condition for s in specs]
    run_nframes = [s.nframes for s in specs]
    manifest = {
        "name": exp.get("name", data_dir.name),
        "fs": fs,
        "runs": [
            {"name": n, "condition": c, "nframes": f, "sbx": str(s.sbx_path)}
            for n, c, f, s in zip(run_names, run_conditions, run_nframes, specs)
        ],
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    raw_dir = out_dir / "raw_tifs"
    aqua_dir = out_dir / "aqua_input"
    qc_dir = out_dir / "qc"
    an_dir = out_dir / "analysis"
    plane0 = out_dir / "suite2p" / "plane0"

    # red (vascular tracer) channel for registration: only when EVERY run has it
    img_cfg = cfg.get("imaging", {})
    red_mode = str(img_cfg.get("register_on_red", "auto")).lower()
    all_have_red = all(int(s.meta["num_channels"]) >= 2 for s in specs)
    use_red = (red_mode == "true") or (red_mode == "auto" and all_have_red)
    if use_red and not all_have_red:
        sys.exit("register_on_red: true but not every run has a second channel")
    chan2 = 2 if use_red else None
    if use_red:
        print("Red channel present in all runs: registering on the vascular tracer (channel 2).")
    elif any(int(s.meta["num_channels"]) >= 2 for s in specs):
        print("Some runs lack the red channel: registering on green (channel 1).")
    from afferent2p.convert import raw_tiff_name
    tiffs = [raw_dir / raw_tiff_name(s, chan2) for s in specs]

    if "preview" in steps:
        from afferent2p.preview import preview_all
        preview_all(specs, out_dir / "qc_raw",
                    functional_chan=int(cfg.get("imaging", {}).get("functional_chan", 1)))

    transforms = None
    pa = cfg.get("prealign", {})
    if "convert" in steps and pa.get("enabled", True) and len(specs) > 1:
        from afferent2p.prealign import prealign_runs
        transforms = prealign_runs(specs, use_red, qc_dir,
                                   max_angle_deg=float(pa.get("max_angle_deg", 6)),
                                   min_shift_px=float(pa.get("min_shift_px", 2)),
                                   min_angle_deg=float(pa.get("min_angle_deg", 0.5)))

    if "convert" in steps:
        tiffs = convert_all(specs, raw_dir,
                            functional_chan=int(cfg.get("imaging", {}).get("functional_chan", 1)),
                            overwrite=args.overwrite, chan2=chan2, transforms=transforms)

    if "register" in steps:
        missing = [t for t in tiffs if not t.exists()]
        if missing:
            sys.exit(f"Missing raw tiffs (run the convert step first): {missing}")
        frame_size = tuple(int(v) for v in specs[0].meta["frame_size"])
        plane0 = run_suite2p(tiffs, out_dir, fs, cfg, frame_size=frame_size,
                             expected_frames=sum(run_nframes),
                             overwrite=args.overwrite, two_channel=use_red)
        # label every frame of the concatenated output with its run/condition
        from afferent2p.labels import write_frame_labels
        write_frame_labels(plane0, run_names, run_conditions, run_nframes, fs,
                           sbx_paths=[s.sbx_path for s in specs], transforms=transforms)
        print(f"  frame labels written: {plane0 / 'runs.json'}, frame_labels.csv, afferent2p_runs.mat")

    # dims for export/qc
    Ly = Lx = None
    if any(s in steps for s in ("export", "qc")):
        if not (plane0 / "data.bin").exists():
            sys.exit(f"No registered data at {plane0} — run the register step first "
                     "(python run_pipeline.py <config> --steps convert,register)")
        ops = load_plane_ops(plane0)
        Ly, Lx = int(ops["Ly"]), int(ops["Lx"])
        nbin = (plane0 / "data.bin").stat().st_size // (Ly * Lx * 2)
        if nbin != sum(run_nframes):
            sys.exit(
                f"Registered binary holds {nbin} frames but the sbx runs total "
                f"{sum(run_nframes)}.\nThis usually means suite2p ran on an "
                "incomplete raw tiff from an interrupted earlier conversion.\n"
                "Fix: re-run conversion and registration:\n"
                "  python run_pipeline.py <config> --steps convert,register,export,loco,qc,analyze\n"
                "(the convert step now detects and re-writes incomplete tiffs)")

    if "export" in steps and cfg.get("export", {}).get("aqua_tifs", True):
        ex = cfg.get("export", {})
        # binning: explicit aqua_bin_t (>1) wins; otherwise derive from the
        # target rate exactly like the MATLAB pipeline
        # (GenerateExptProjections: bin = round(scanRate / rate_target))
        bin_t = int(ex.get("aqua_bin_t", 0) or 0)
        rate = ex.get("aqua_rate_hz", 1.0)
        if bin_t < 1:
            bin_t = max(1, int(round(fs / float(rate)))) if rate else 1
        export_registered_tiffs(
            plane0, aqua_dir, run_names, run_conditions, run_nframes, Ly, Lx,
            chunk_frames=int(ex.get("chunk_frames", 4000)),
            overwrite=args.overwrite,
            dtype=str(ex.get("aqua_dtype", "uint16")),
            bin_t=bin_t,
            crop=tuple(ex.get("aqua_crop", [80, 80, 40, 40])),
            fs=fs)

    loco = None
    if "loco" in steps and cfg.get("loco", {}).get("enabled", True):
        lc = cfg.get("loco", {})
        print("Locomotion:")
        loco = [
            loco_for_run(s.sbx_path, s.nframes, fs,
                         quad_glob=lc.get("quad_glob", "*quad*"),
                         counts_per_cm=lc.get("counts_per_cm"),
                         still_thresh=float(lc.get("still_thresh", 1.0)),
                         min_still_s=float(lc.get("min_still_s", 5.0)),
                         smooth_s=float(lc.get("smooth_s", 0.5)),
                         quad_kind=lc.get("quad_kind", "cumulative"))
            for s in specs
        ]
        pc = cfg.get("pupil", {})
        if pc.get("enabled", True):
            print("Pupil:")
            for s, ld in zip(specs, loco):
                pr = pupil_for_run(s.sbx_path, s.nframes, fs,
                                   pattern=pc.get("glob", "*eye*DLC*.csv"),
                                   lik_min=float(pc.get("likelihood_min", 0.6)),
                                   max_gap_s=float(pc.get("max_gap_s", 0.5)),
                                   smooth_s=float(pc.get("smooth_s", 0.25)),
                                   min_edge_points=int(pc.get("min_edge_points", 3)))
                if pr is not None:
                    ld["pupil"] = pr["diam_px"]; ld["pupil_valid"] = pr["valid"]
                    ld["pupil_area"] = pr["area_px"]; ld["pupil_cx"] = pr["cx"]
                    ld["pupil_cy"] = pr["cy"]; ld["pupil_file"] = pr["file"]
        np.save(out_dir / "loco.npy", np.array(loco, dtype=object), allow_pickle=True)
        # flat per-frame behaviour table (frame-aligned), one CSV per run
        beh_dir = out_dir / "behavior"; beh_dir.mkdir(exist_ok=True)
        for s, ld in zip(specs, loco):
            cols = {"frame": np.arange(s.nframes), "time_s": np.arange(s.nframes) / fs,
                    "velocity": ld["velocity"], "still": ld["still"].astype(int)}
            if "pupil" in ld:
                cols["pupil_diam_px"] = ld["pupil"]; cols["pupil_valid"] = ld["pupil_valid"].astype(int)
            import csv as _csv
            with open(beh_dir / f"{s.name}_{s.condition}_behavior.csv", "w", newline="") as fh:
                w = _csv.writer(fh); w.writerow(cols.keys())
                w.writerows(zip(*[np.asarray(v) for v in cols.values()]))
    elif (out_dir / "loco.npy").exists():
        loco = list(np.load(out_dir / "loco.npy", allow_pickle=True))

    if "qc" in steps:
        qc_report(plane0, qc_dir, run_names, run_conditions, run_nframes,
                  Ly, Lx, loco=loco, fs=fs)

    if "compare" in steps:
        from afferent2p.compare import compare_to_reference
        if not args.ref:
            sys.exit("--steps compare needs --ref <MATLAB projection or registered tif>")
        if not (plane0 / "data.bin").exists():
            sys.exit("compare needs the registered movie — run the register step first.")
        ops = load_plane_ops(plane0)
        compare_to_reference(plane0, qc_dir, sum(run_nframes),
                             int(ops["Ly"]), int(ops["Lx"]), Path(args.ref))

    if "manualroi" in steps:
        from afferent2p.manual_roi import (write_roi_templates, load_rois,
                                           extract_manual_traces)
        if not (plane0 / "data.bin").exists():
            sys.exit("manualroi needs the registered movie — run the register step first.")
        ops = load_plane_ops(plane0)
        mLy, mLx = int(ops["Ly"]), int(ops["Lx"])
        roi_dir = out_dir / "manual_roi"
        mroi_cfg = cfg.get("manual_roi", {})
        masks, names = load_rois(roi_dir, mLy, mLx, explicit=mroi_cfg.get("file")) \
            if roi_dir.exists() else (None, None)
        if masks is None:
            write_roi_templates(plane0, roi_dir, sum(run_nframes), mLy, mLx)
        else:
            extract_manual_traces(plane0, roi_dir, masks, names,
                                  sum(run_nframes), mLy, mLx,
                                  ring_px=int(mroi_cfg.get("ring_px", 8)))
            if loco is None and (out_dir / "loco.npy").exists():
                loco = list(np.load(out_dir / "loco.npy", allow_pickle=True))
            analyze_prepost(roi_dir, out_dir / "analysis_manual",
                            run_conditions, run_nframes, fs, loco, cfg,
                            mean_image=registered_mean(plane0, sum(run_nframes)))

    if "curate" in steps:
        from afferent2p.analyze import curate_only
        if not (plane0 / "F.npy").exists():
            sys.exit("curate needs suite2p ROI output (F.npy) — run the register step with ROI detection first.")
        curate_only(plane0, an_dir, fs, cfg)

    if "analyze" in steps and cfg.get("analysis", {}).get("enabled", True):
        if not (plane0 / "F.npy").exists():
            print("analysis: no F.npy (ROI detection disabled or not run) — skipping")
        else:
            if not (plane0 / "runs.json").exists():
                from afferent2p.labels import write_frame_labels
                write_frame_labels(plane0, run_names, run_conditions, run_nframes, fs,
                                   sbx_paths=[s.sbx_path for s in specs])
            analyze_prepost(plane0, an_dir, run_conditions, run_nframes, fs, loco, cfg,
                            mean_image=registered_mean(plane0, sum(run_nframes)))

    print(f"\nDone. Outputs in: {out_dir}")


if __name__ == "__main__":
    main()
