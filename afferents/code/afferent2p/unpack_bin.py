#!/usr/bin/env python3
"""unpack_bin — Scanbox .sbx -> temporally binned uint16 tif (AQuA2 input).

Reproduces the lab's UnpackSBX.m / WriteSbxPlaneTif.m 'binT' behaviour in
Python, with the bin size given in seconds (default 2 s per frame):

  * values  = 65535 - raw   (same offset correction as readSBX)
  * each output frame = MEAN of N consecutive frames, N = round(fs * seconds)
  * written as uint16, ImageJ-style TIFF (opens in plain ImageJ/Fiji/MATLAB),
    saved NEXT TO the .sbx in the recording folder — no extra output folders
  * plus <name>_meanproj.tif (mean image, for checking by eye), a small
    .json with the frame rate / bin size for AQuA2, and <name>_behavior.csv:
    running speed / still fraction (quadrature) and pupil diameter (DLC eye
    csv) averaged over the SAME bins, so row i = recorded frame i
  * then (channel 1) the global-brightness step of afferent2p/globalfix.py:
    the whole-field level is flattened and bad frames (flashes, jumps, slow
    state changes) are FLAGGED — every frame stays in the movie; the frames
    of a flagged stretch that differ clearly from their neighbours (the flash
    itself) are excluded from the statistics, the rest kept. <name>_frames.csv
    lists every frame with its class and decision, <name>_globalfix.png shows
    the traces, <name>_check_frames.tif / _check_movie.tif are ImageJ stacks
    of the flagged stretches for checking frame by frame (tab 2 of the GUI
    opens them). --keep-raw also keeps the uncorrected movie (twice the
    disk), --no-global-fix skips the step.

No registration is done here by default — AQuA2 performs its own motion
correction (one whole-pixel translation per binned frame). --register (v33,
tab 1 checkbox) registers the RAW 15.6-Hz frames of channel 1 first (rigid
translation per frame to a still-period reference, afferent2p/fullreg.py —
the legacy MATLAB CorrectExpt / RegisterSBX operation) and bins the
registered frames; the .json then carries a "registration" entry with the
shift statistics and the border crop. Decide with reg_compare.py before
switching a whole animal to it (every run must be treated the same way).

GUI:   python unpack_bin.py
CLI:   python unpack_bin.py <folder> [<folder> ...] [--sec 2] [--frames N]
       (folders may be session folders or mouse folders containing sessions)
       python unpack_bin.py <folder> --register        (register the raw frames first, channel 1)
       python unpack_bin.py <folder> --behavior-only
       (only the behaviour csv, for movies unpacked before their DLC eye csv
       existed — no re-unpacking, the movie and its flags are untouched)
"""
from __future__ import annotations
import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from sbxreader import sbx_memmap, sbx_get_metadata

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))


def write_behavior_csv(sbx_path: Path, N: int, n_out: int, T_total: int, fs: float,
                       out_path: Path, log=print):
    """Locomotion (quad) + pupil (DLC) binned with the SAME N-frame bins as the
    tiff, so row i of the CSV = frame i of the AQuA tiff (afferent2p/behavior.py)."""
    try:
        from afferent2p.behavior import write_behavior_csv as _write
    except ImportError:
        log("  behaviour: afferent2p package not found next to unpack_bin.py — skipped"); return
    return _write(sbx_path, N, n_out, T_total, fs, out_path, log=log)


def rebuild_behavior(folders, log=print) -> int:
    """--behavior-only: write <movie>_behavior.csv again for every binned movie
    under the folders (running + pupil, from the quadrature / DLC files next
    to the .sbx) without unpacking anything. For a DLC eye csv that arrived
    after the run was unpacked, or movies from before the pupil columns."""
    from afferent2p.aqua_io import find_binned_tifs
    from afferent2p.behavior import rebuild
    sessions = []
    for f in folders:                       # every folder holding a binned movie (the .sbx may be archived)
        f = Path(f)
        found = sorted({p.parent for p in f.rglob("*_bin*s.tif") if not _parked(p, f)} | set(find_sessions(f)), key=str)
        sessions.extend(s for s in found if s not in sessions)
    n = 0
    for s in sessions:
        for tif in find_binned_tifs(s):
            log(f"{tif.name}:")
            if rebuild(tif, log=log) is not None:
                n += 1
    log(f"behaviour csv written for {n} movie(s)")
    return n


# ----------------------------------------------------------------------------- core

def _parked(p: Path, root: Path) -> bool:
    """Inside a pipeline output folder or a parked _previous_<date> folder (aqua_reset.py)?"""
    try:
        from afferent2p.plan import skip_parts
        return skip_parts(p.relative_to(root).parts[:-1])
    except Exception:
        return False


def find_sessions(folder: Path) -> list[Path]:
    """A session folder contains an .sbx; a mouse folder contains session folders."""
    folder = Path(folder)
    if list(folder.glob("*.sbx")):
        return [folder]
    return sorted({p.parent for p in folder.rglob("*.sbx") if not _parked(p, folder)}, key=lambda p: str(p))


def bin_sbx(
    sbx_path: Path,
    seconds: float | None = 2.0,
    frames: int | None = None,
    chan="all",
    crop: tuple[int, int, int, int] = (0, 0, 0, 0),
    out_dir: Path | None = None,
    overwrite: bool = False,
    behavior: bool = True,
    global_fix: bool = True,
    keep_raw: bool = False,
    register: bool = False,
    log=print,
) -> Path:
    """chan: 1 = green only, 2 = red only, 'all' (default) = every channel the
    file has. Channel 1 keeps the plain name <stem>_<tag>.tif; channel 2 is
    written next to it as <stem>_<tag>_ch2.tif (+ its own _meanproj)."""
    sbx_path = Path(sbx_path)
    meta = sbx_get_metadata(str(sbx_path))
    fs = float(meta["frame_rate"])
    nchan = int(meta["num_channels"])
    if str(chan).lower() in ("all", "0"):
        chans = list(range(1, nchan + 1))
    else:
        chan = int(chan)
        if chan > nchan:
            raise ValueError(f"{sbx_path.name}: requested channel {chan} but file has {nchan}")
        chans = [chan]
    N = int(frames) if frames else max(1, int(round(fs * float(seconds))))
    sec_actual = N / fs
    out_dir = Path(out_dir) if out_dir else sbx_path.parent
    # tag WITHOUT a dot: MATLAB's fileparts() splits a name at the last dot, so a
    # tiff called ..._bin31f_1.98s.tif came out of the AQuA2 GUI as ..._bin31f_1_*
    tag = f"bin{N}f_{sec_actual:.2f}s".replace(".", "p")
    legacy_tag = f"bin{N}f_{sec_actual:.2f}s"

    def suffix(c):
        return "" if c == 1 else f"_ch{c}"

    out_paths = {c: out_dir / f"{sbx_path.stem}_{tag}{suffix(c)}.tif" for c in chans}
    proj_paths = {c: out_dir / f"{sbx_path.stem}_{tag}{suffix(c)}_meanproj.tif" for c in chans}
    info_paths = {c: out_dir / f"{sbx_path.stem}_{tag}{suffix(c)}.json" for c in chans}
    legacy_paths = {c: out_dir / f"{sbx_path.stem}_{legacy_tag}{suffix(c)}.tif" for c in chans}
    if not overwrite:
        done = [c for c in chans if out_paths[c].exists() or legacy_paths[c].exists()]
        for c in done:
            have = out_paths[c] if out_paths[c].exists() else legacy_paths[c]
            log(f"  exists, skipping: {have.name}"
                + ("  (older dotted name — fine, every tool accepts it)" if have is legacy_paths[c] else ""))
        chans = [c for c in chans if c not in done]
        if not chans:
            c0 = done[0]
            return out_paths[c0] if out_paths[c0].exists() else legacy_paths[c0]

    if register and 1 in chans:
        # channel 1: every raw frame registered to a still-period reference, then binned (afferent2p/fullreg.py);
        # the behaviour csv is written first so the reference can be taken from a still stretch
        from afferent2p.fullreg import register_bin_sbx
        mm0 = sbx_memmap(str(sbx_path)); T0 = int(mm0.shape[0]); del mm0
        if behavior:
            write_behavior_csv(sbx_path, N, T0 // N, T0, fs, out_dir / f"{sbx_path.stem}_{tag}_behavior.csv", log)
        register_bin_sbx(sbx_path, N, out_dir, chan=1, tag=tag, behavior_csv=out_dir / f"{sbx_path.stem}_{tag}_behavior.csv",
                         global_fix=global_fix, overwrite=overwrite, log=log)
        chans = [c for c in chans if c != 1]
        if not chans:
            return out_paths[1]
    mm = sbx_memmap(str(sbx_path))
    T, H, W = mm.shape[0], mm.shape[3], mm.shape[4]
    L, R, Tc, B = crop
    Hc, Wc = H - Tc - B, W - L - R
    n_out = T // N
    if n_out == 0:
        raise ValueError(f"{sbx_path.name}: only {T} frames, fewer than bin size {N}")
    gb = n_out * Hc * Wc * 2 / 1e9
    log(f"  {sbx_path.name}: {T} frames @ {fs:.2f} Hz -> {n_out} frames of "
        f"{N}-frame means ({sec_actual:.2f} s/frame, {fs / N:.3f} Hz), "
        f"{Hc}x{Wc} px, ~{gb:.2f} GB per channel; channel(s) {chans} of {nchan}")

    use_imagej = gb < 3.9          # ImageJ-style header only for classic (<4 GB) TIFF
    kw = dict(imagej=True, metadata={"axes": "TYX", "finterval": sec_actual}) \
        if use_imagej else dict(bigtiff=True)
    tmps = {c: out_paths[c].with_name(out_paths[c].name + ".partial") for c in chans}
    dsts = {c: tifffile.memmap(str(tmps[c]), shape=(n_out, Hc, Wc), dtype=np.uint16, **kw)
            for c in chans}
    acc_proj = {c: np.zeros((Hc, Wc), np.float64) for c in chans}
    bins_per_chunk = max(1, 2000 // N)         # read ~2000 raw frames at a time
    t0 = time.time()
    for b0 in range(0, n_out, bins_per_chunk):
        b1 = min(b0 + bins_per_chunk, n_out)
        for c in chans:                         # one channel at a time keeps memory as before
            raw = np.asarray(mm[b0 * N:b1 * N, 0, c - 1]).astype(np.float32)  # offset-corrected
            if any(crop):
                raw = raw[:, Tc:H - B if B else None, L:W - R if R else None]
            binned = raw.reshape(b1 - b0, N, Hc, Wc).mean(axis=1)
            dsts[c][b0:b1] = np.clip(np.round(binned), 0, 65535).astype(np.uint16)
            acc_proj[c] += binned.sum(axis=0)
        if (b0 // bins_per_chunk) % 5 == 0:
            log(f"    {b1}/{n_out} frames written ({time.time() - t0:.0f} s)")
    for c in chans:
        dsts[c].flush()
    del dsts; del mm
    for c in chans:
        tmps[c].replace(out_paths[c])
        proj = np.clip(np.round(acc_proj[c] / n_out), 0, 65535).astype(np.uint16)
        tifffile.imwrite(str(proj_paths[c]), proj, imagej=True)
        info = {
            "source_sbx": str(sbx_path), "channel": c, "channels_in_file": nchan,
            "frames_in": int(T), "fs_in_hz": fs,
            "bin_frames": N, "seconds_per_frame": sec_actual,
            "fs_out_hz": fs / N, "frames_out": int(n_out),
            "crop_LRTB_px": [int(v) for v in crop], "height": Hc, "width": Wc,
            "dtype": "uint16", "values": "mean of (65535 - raw) over each bin",
            "registration": "none (do it in AQuA2)",
        }
        info_paths[c].write_text(json.dumps(info, indent=2), encoding="utf-8")
        log(f"  wrote {out_paths[c].name}, {proj_paths[c].name}, {info_paths[c].name}")
    if behavior and not (register and (out_dir / f"{sbx_path.stem}_{tag}_behavior.csv").exists()):
        write_behavior_csv(sbx_path, N, n_out, T, fs,
                           out_dir / f"{sbx_path.stem}_{tag}_behavior.csv", log)
    if global_fix and 1 in chans:
        # flatten whole-field brightness changes and flag bad frames in the
        # functional channel before AQuA2 sees it (see afferent2p/globalfix.py)
        if keep_raw:
            raw_copy = out_dir / f"{sbx_path.stem}_{tag}_raw.tif"
            shutil.copy2(out_paths[1], raw_copy)
            log(f"  kept the uncorrected movie as {raw_copy.name} (not used by any tool)")
        try:
            from afferent2p.globalfix import apply as global_apply
            global_apply(out_paths[1], spf=sec_actual, log=log)
        except ImportError:
            log("  global correction skipped: afferent2p package not found next to unpack_bin.py")
        except Exception as e:
            log(f"  global correction FAILED ({e}) — the tiff is the plain binned movie")
    return out_paths[chans[0]]


def process_folders(folders, seconds=2.0, frames=None, chan="all", crop=(0, 0, 0, 0),
                    overwrite=False, behavior=True, global_fix=True, keep_raw=False, register=False, log=print):
    sessions = []
    for f in folders:
        sessions.extend(s for s in find_sessions(Path(f)) if s not in sessions)
    if not sessions:
        log("No .sbx files found."); return
    log(f"{len(sessions)} session(s):")
    for s in sessions:
        log(f"  {s}")
    jobs = [(s, sbx) for s in sessions for sbx in sorted(s.glob("*.sbx"))]
    sizes = sorted(sbx.stat().st_size for _, sbx in jobs)
    typical = sizes[len(sizes) // 2] if sizes else 0
    t_start = time.time()
    for k, (s, sbx) in enumerate(jobs, 1):
        log(f"\n== {s.name} ==  [{k}/{len(jobs)}]")
        if not sbx.with_suffix(".mat").exists() and typical > 0 and sbx.stat().st_size < 0.5 * typical:
            # Scanbox writes the info .mat when a recording ends properly: without it, and much smaller than the other
            # recordings, the .sbx is an aborted recording (e.g. 266 frames) — not a run; the planner skips it too
            log(f"  skipped: no {sbx.with_suffix('.mat').name} next to the .sbx and only {sbx.stat().st_size / 1e6:.0f} MB "
                f"(the others ~{typical / 1e9:.1f} GB) — an aborted recording, not a run")
            continue
        try:
            bin_sbx(sbx, seconds=seconds, frames=frames, chan=chan, crop=crop,
                    overwrite=overwrite, behavior=behavior, global_fix=global_fix,
                    keep_raw=keep_raw, register=register, log=log)
        except Exception as e:
            log(f"  FAILED: {e}")
    log(f"\nAll done ({len(jobs)} recordings, {(time.time() - t_start) / 60:.0f} min).")


# ----------------------------------------------------------------------------- GUI

def launch_gui():
    import queue, threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("unpack_bin — sbx -> binned uint16 tif for AQuA2")
            self.geometry("900x560")
            self.folders: list[Path] = []
            self.q: queue.Queue[str] = queue.Queue()
            self.running = False
            top = ttk.Frame(self); top.pack(fill="x", padx=8, pady=6)
            ttk.Button(top, text="Add folder (session or mouse)…", command=self.add).pack(side="left")
            ttk.Button(top, text="Remove selected", command=self.remove).pack(side="left", padx=6)
            ttk.Button(top, text="Clear", command=self.clear).pack(side="left")
            ttk.Button(top, text="Run", command=self.run).pack(side="right")
            lf = ttk.LabelFrame(self, text="Recordings"); lf.pack(fill="both", expand=True, padx=8)
            self.tree = ttk.Treeview(lf, columns=("sbx",), height=8)
            self.tree.heading("#0", text="Session folder"); self.tree.heading("sbx", text=".sbx")
            self.tree.column("#0", width=540); self.tree.column("sbx", width=300)
            self.tree.pack(fill="both", expand=True)
            o = ttk.LabelFrame(self, text="Options"); o.pack(fill="x", padx=8, pady=6)
            self.v = {"sec": tk.StringVar(value="2.0"), "frames": tk.StringVar(value=""),
                      "chan": tk.StringVar(value="all"), "crop": tk.StringVar(value="0,0,0,0"),
                      "overwrite": tk.BooleanVar(value=False),
                      "gfix": tk.BooleanVar(value=True)}
            def field(label, key, w):
                ttk.Label(o, text=label).pack(side="left", padx=(10, 2))
                ttk.Entry(o, textvariable=self.v[key], width=w).pack(side="left")
            field("seconds per frame", "sec", 6)
            field("or frames per bin (overrides)", "frames", 5)
            field("channel (all / 1 / 2)", "chan", 4)
            field("crop L,R,T,B", "crop", 12)
            ttk.Checkbutton(o, text="overwrite", variable=self.v["overwrite"]).pack(side="left", padx=8)
            ttk.Checkbutton(o, text="flatten global brightness + flag bad frames",
                            variable=self.v["gfix"]).pack(side="left", padx=8)
            lg = ttk.LabelFrame(self, text="Log"); lg.pack(fill="both", expand=True, padx=8, pady=(0, 8))
            self.log = tk.Text(lg, height=10, state="disabled", font=("Consolas", 9))
            self.log.pack(fill="both", expand=True)
            self.after(150, self.poll)

        def add(self):
            d = filedialog.askdirectory(title="Select a session folder or a mouse folder")
            if not d:
                return
            sessions = find_sessions(Path(d))
            if not sessions:
                messagebox.showerror("unpack_bin", f"No .sbx found under\n{d}"); return
            for s in sessions:
                if s in self.folders:
                    continue
                self.folders.append(s)
                sbx = sorted(s.glob("*.sbx"))
                self.tree.insert("", "end", iid=str(s), text=str(s),
                                 values=(", ".join(p.name for p in sbx),))

        def remove(self):
            for iid in self.tree.selection():
                self.folders = [f for f in self.folders if str(f) != iid]
                self.tree.delete(iid)

        def clear(self):
            self.folders.clear(); self.tree.delete(*self.tree.get_children())

        def run(self):
            if self.running:
                return
            if not self.folders:
                messagebox.showinfo("unpack_bin", "Add a folder first."); return
            try:
                sec = float(self.v["sec"].get())
                frames = int(self.v["frames"].get()) if self.v["frames"].get().strip() else None
                chan = self.v["chan"].get().strip() or "all"
                crop = tuple(int(x) for x in self.v["crop"].get().split(","))
                assert len(crop) == 4
            except Exception:
                messagebox.showerror("unpack_bin", "Check the option values."); return
            self.running = True
            threading.Thread(target=self.worker, args=(sec, frames, chan, crop,
                             self.v["overwrite"].get(), self.v["gfix"].get()), daemon=True).start()

        def worker(self, sec, frames, chan, crop, overwrite, gfix):
            process_folders(list(self.folders), seconds=sec, frames=frames, chan=chan,
                            crop=crop, overwrite=overwrite, global_fix=gfix,
                            log=lambda s: self.q.put(str(s) + "\n"))
            self.running = False

        def poll(self):
            try:
                while True:
                    line = self.q.get_nowait()
                    self.log.configure(state="normal"); self.log.insert("end", line)
                    self.log.see("end"); self.log.configure(state="disabled")
            except queue.Empty:
                pass
            self.after(150, self.poll)

    App().mainloop()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="*", help="session or mouse folders (none = open the GUI)")
    ap.add_argument("--sec", type=float, default=2.0, help="seconds per output frame (default 2)")
    ap.add_argument("--frames", type=int, default=None, help="frames per bin (overrides --sec)")
    ap.add_argument("--chan", default="all",
                    help="'all' (default) = every channel in the file (ch1 -> plain name, "
                         "ch2 -> _ch2), or 1 = green only, 2 = red only")
    ap.add_argument("--crop", default="0,0,0,0", help="L,R,T,B pixels to crop (default none)")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--no-global-fix", action="store_true",
                    help="do NOT flatten whole-field brightness changes / flag bad frames "
                         "(default: applied to channel 1; see afferent2p/globalfix.py)")
    ap.add_argument("--register", action="store_true",
                    help="register the raw frames of channel 1 (rigid, to a still-period reference) before binning "
                         "(afferent2p/fullreg.py); a few minutes more per run; channel 2 is binned as recorded")
    ap.add_argument("--keep-raw", action="store_true",
                    help="also keep the uncorrected binned movie as <run>_<tag>_raw.tif (twice the disk)")
    ap.add_argument("--no-behavior", action="store_true",
                    help="skip the binned locomotion/pupil csv")
    ap.add_argument("--behavior-only", action="store_true",
                    help="write only <movie>_behavior.csv (running + pupil) again for movies that already "
                         "exist — e.g. when the DLC eye csv arrived after unpacking; nothing else is touched")
    a = ap.parse_args()
    if a.behavior_only:
        if not a.folders:
            raise SystemExit("--behavior-only needs the session or mouse folder(s)")
        rebuild_behavior(a.folders)
    elif not a.folders:
        launch_gui()
    else:
        process_folders(a.folders, seconds=a.sec, frames=a.frames, chan=a.chan,
                        crop=tuple(int(x) for x in a.crop.split(",")), overwrite=a.overwrite,
                        behavior=not a.no_behavior, global_fix=not a.no_global_fix,
                        keep_raw=a.keep_raw, register=a.register)
