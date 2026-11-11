"""Suite2p wrapper: register all runs together to one reference (+ optional ROI detection).

Registering the concatenated pre+post movie to a single reference is what makes
the pre/post comparison valid for the *same* afferents (same role as the lab's
ProcessScanboxData_concatFirst).
"""
from __future__ import annotations
from pathlib import Path

import numpy as np


def _patch_temporal_smooth():
    """Work around a suite2p 1.1.0 bug: temporal_smooth (used when
    smooth_sigma_time > 0) returns a numpy array where the torch-based
    registration expects a tensor. Wrap it to keep tensors as tensors."""
    import torch
    from suite2p.registration import utils as reg_utils
    if getattr(reg_utils.temporal_smooth, "_afferent2p_patched", False):
        return
    orig = reg_utils.temporal_smooth

    def temporal_smooth_torchsafe(data, sigma):
        if isinstance(data, torch.Tensor):
            out = orig(data.detach().cpu().numpy(), sigma)
            return torch.from_numpy(np.ascontiguousarray(out)).to(
                device=data.device, dtype=data.dtype)
        return orig(data, sigma)

    temporal_smooth_torchsafe._afferent2p_patched = True
    # patch every module that imported the name directly
    from suite2p import registration as reg_pkg
    import importlib
    for modname in ("utils", "rigid", "nonrigid", "register"):
        try:
            mod = importlib.import_module(f"suite2p.registration.{modname}")
            if hasattr(mod, "temporal_smooth"):
                mod.temporal_smooth = temporal_smooth_torchsafe
        except ImportError:
            pass


def _clamp_block(requested: int, dim: int) -> int:
    """Nonrigid blocks must be comfortably smaller than the frame."""
    b = int(requested)
    while b > max(32, dim // 2):
        b //= 2
    return max(32, b)


def run_suite2p(
    tiff_paths: list[Path],
    out_dir: Path,
    fs: float,
    cfg: dict,
    frame_size: tuple[int, int] | None = None,
    expected_frames: int | None = None,
    overwrite: bool = False,
    two_channel: bool = False,
) -> Path:
    """Run suite2p on the ordered tiff list. Returns the plane0 folder."""
    import shutil
    import suite2p
    _patch_temporal_smooth()

    # suite2p resumes from an existing output folder; if that output came from a
    # different (e.g. truncated) input, it must be cleared or it will be reused
    s2p_dir = out_dir / "suite2p"
    old_bin = s2p_dir / "plane0" / "data.bin"
    if old_bin.exists() and overwrite:
        print(f"--overwrite: deleting previous suite2p output {s2p_dir}")
        shutil.rmtree(s2p_dir)
    if old_bin.exists() and expected_frames and frame_size:
        nbin = old_bin.stat().st_size // (int(frame_size[0]) * int(frame_size[1]) * 2)
        if nbin != expected_frames:
            print(f"Existing suite2p output holds {nbin} frames but current input "
                  f"has {expected_frames} — deleting stale {s2p_dir} and re-running.")
            shutil.rmtree(s2p_dir)

    reg_cfg = cfg.get("registration", {})
    roi_cfg = cfg.get("roi_detection", {})
    img_cfg = cfg.get("imaging", {})

    settings = suite2p.default_settings()
    settings["fs"] = float(fs)
    settings["tau"] = float(img_cfg.get("tau", 0.7))
    settings["torch_device"] = reg_cfg.get("torch_device", "cpu")
    settings["run"]["do_registration"] = 1
    settings["run"]["do_detection"] = bool(roi_cfg.get("enabled", True))
    settings["run"]["do_deconvolution"] = bool(roi_cfg.get("enabled", True))
    settings["registration"]["nonrigid"] = bool(reg_cfg.get("nonrigid", True))
    # low-SNR (sparse axon) data: individual frames are too noisy to align, so
    # smooth over time for the ALIGNMENT ONLY (traces/data are not smoothed),
    # and build the reference from more frames
    settings["registration"]["smooth_sigma_time"] = float(reg_cfg.get("smooth_sigma_time", 2))
    settings["registration"]["nimg_init"] = int(reg_cfg.get("nimg_init", 1000))
    settings["registration"]["maxregshiftNR"] = float(reg_cfg.get("maxregshiftNR", 5))
    block = list(reg_cfg.get("block_size", (128, 128)))
    if frame_size is not None:
        block = [_clamp_block(block[0], int(frame_size[0])),
                 _clamp_block(block[1], int(frame_size[1]))]
    settings["registration"]["block_size"] = tuple(block)
    settings["registration"]["maxregshift"] = float(reg_cfg.get("maxregshift", 0.1))
    settings["registration"]["batch_size"] = int(reg_cfg.get("batch_size", 200))
    # --- axon-oriented detection (suite2p defaults are tuned for somata) ---
    det = settings["detection"]
    det["soma_crop"] = bool(roi_cfg.get("soma_crop", False))       # keep elongated ROIs
    det["threshold_scaling"] = float(roi_cfg.get("threshold_scaling", 0.8))  # dimmer ROIs
    det["sparsery_settings"]["spatial_scale"] = int(roi_cfg.get("spatial_scale", 1))  # smallest
    det["denoise"] = bool(roi_cfg.get("denoise", True))             # PCA denoise before detection
    det["highpass_time"] = int(roi_cfg.get("highpass_time", 100))   # frames; removes slow drift
    det["max_overlap"] = float(roi_cfg.get("max_overlap", 0.75))
    det["sparsery_settings"]["max_ROIs"] = int(roi_cfg.get("max_rois", 5000))
    # temporal binning used ONLY to find ROIs (traces are extracted at full
    # rate): binning to ~1 s pools the sparse axonal signal above the noise
    bin_s = float(roi_cfg.get("detection_bin_s", 1.0))
    if bin_s > 0:
        det["bin_size"] = max(1, int(round(bin_s * float(fs))))
    # suite2p's built-in classifier is soma-trained and rejects axons; keep
    # everything and let the analysis step auto-curate (or curate in the GUI)
    settings["classification"]["use_builtin_classifier"] = False
    settings["classification"]["preclassify"] = 0.0
    # optional: a custom classifier trained on YOUR curated afferent ROIs
    # (suite2p GUI: curate a few FOVs, then File > Classifier > add current data
    # / save classifier), which then replaces the trace-based auto-curation
    if roi_cfg.get("classifier_path"):
        settings["classification"]["classifier_path"] = str(roi_cfg["classifier_path"])
    # deconvolution baseline window (s) — shorter than default 60 s for fast sensors
    settings["dcnv_preprocess"]["win_baseline"] = float(roi_cfg.get("dcnv_win_baseline_s", 30.0))
    # Fall.mat needs ROI outputs; writing it without detection crashes suite2p
    settings["io"]["save_mat"] = bool(roi_cfg.get("save_mat", True)) and \
        bool(roi_cfg.get("enabled", True))
    # elongated axonal ROIs can have a huge median "radius", which makes
    # suite2p's lam_percentile filter allocate an enormous footprint — disable it
    settings["extraction"]["lam_percentile"] = 0.0

    # suite2p joins file_list entries onto data_path[0]; give file names relative
    # to their common parent folder
    parent = Path(tiff_paths[0]).parent
    db = {
        "data_path": [str(parent)],
        "file_list": [str(Path(p).relative_to(parent)) for p in tiff_paths],
        "input_format": "tif",
        "save_path0": str(out_dir),
        "nplanes": 1,
        "nchannels": 2 if two_channel else 1,
        "functional_chan": 1,          # interleaved tiffs: functional first
        "fs": float(fs),
    }
    # red vascular tracer = bright static structure -> align on it instead of
    # the dim GCaMP frames (the shifts are applied to both channels)
    settings["registration"]["align_by_chan2"] = bool(two_channel)
    if two_channel:
        print("Registration will be computed on channel 2 (red, structural) and applied to both.")

    print("Running suite2p (registration"
          + (" + ROI detection" if settings["run"]["do_detection"] else "")
          + f") on {len(tiff_paths)} run(s)...")
    suite2p.run_s2p(db=db, settings=settings)

    plane0 = out_dir / "suite2p" / "plane0"
    if not (plane0 / "data.bin").exists():
        raise RuntimeError(f"suite2p did not produce {plane0 / 'data.bin'}")
    return plane0


def load_plane_ops(plane0: Path) -> dict:
    for name in ("ops.npy", "settings.npy", "db.npy"):
        p = plane0 / name
        if p.exists():
            d = np.load(p, allow_pickle=True).item()
            if "Ly" in d and "Lx" in d:
                return d
    # merge db + ops-like files if needed
    merged = {}
    for name in ("db.npy", "settings.npy", "ops.npy"):
        p = plane0 / name
        if p.exists():
            merged.update(np.load(p, allow_pickle=True).item())
    if "Ly" not in merged:
        raise RuntimeError(f"Could not find Ly/Lx in {plane0}")
    return merged
