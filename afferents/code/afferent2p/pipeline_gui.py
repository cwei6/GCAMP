#!/usr/bin/env python3
r"""pipeline_gui — one window for the AQuA2 route.

    python pipeline_gui.py            (from the pipeline folder, in the Anaconda prompt)

Tabs in the order the steps are used (GUIDE.md):
  1 Unpack              unpack_bin.py     .sbx -> binned tiff + behaviour csv + brightness step
  2 Clean & check       aqua_clean.py     flag bad frames on existing tiffs (nothing removed); a table of
                                          every run, the flagged stretches, check stacks in ImageJ, your
                                          per-stretch decisions
  3 Detection settings  aqua_tune.py      settings by measurement: plan, sweep (MATLAB), report, apply ->
                                          aqua_opts\<experiment>_Opt.csv; optional: set by eye in the
                                          AQuA2 GUI (FOV<k>_Opt.csv); a FOV without a file gets default_Opt.csv
  4 Batch               aqua_batch.py     writes <mouse>\aqua_batch\ (which settings each FOV got, what
                                          was queued) and the MATLAB scripts
  5 MATLAB              the one-paste command run_aqua_and_cfu.m (AQuA2 batch, then CFU batch): copy it
                                          into a MATLAB window, or start MATLAB from here (matlab -batch)
  6 Curate              cfu_curate.py (its own window)
  7 Analyze             aqua_analyze.py   the comparison + the general representation, locomotion trials,
                                          heatmaps; shows summary.txt
  8 Report              aqua_report.py    the main report and the QC report (HTML)
The strip above the tabs shows where the animal is (one chip per step, read from the files on
disk: done / pending / stale / missing, the next action under it; click a chip to go to its tab,
"Details -> log" prints every step) — so the 3 -> 4 -> 5 loop can be picked up after a break.
"Run all (guided)" above the tabs runs 4 -> 3 -> 4 -> 5 -> 7 -> 8 in one go; at each MATLAB step it
shows the command to paste, copies it to the clipboard, polls the outputs every 30 s and continues
by itself (or starts MATLAB from here when that checkbox is ticked and matlab.exe is set).

Every step calls the same functions as the command-line tools, in a background
thread, with the log at the bottom. Settings are remembered in
pipeline_gui_settings.json next to this file.
"""
from __future__ import annotations
import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError:
    raise SystemExit("tkinter is not available in this Python — use the Anaconda python.exe")

SETTINGS = HERE / "pipeline_gui_settings.json"

# the CFU step's AQuA2 settings (tab 4); kept in step with aqua_batch.py (CFU_OVERLAP_DEFAULT / CFU_MIN_EVENTS_DEFAULT)
CFU_OVERLAP_DEFAULT = 0.3      # AQuA2's published overlap rule with the pipeline's threshold (AQuA2's own default: 0.5)
CFU_MIN_EVENTS_DEFAULT = 1     # every event forms a unit (AQuA2's own default: 3)

# which movies tab 3 tunes on (aqua_tune.py --per-day / --baseline-only)
DAY_RULES = ("matched sensitivity: equal recall + false-alarm rate (v36 default)",
             "one common setting for the day (most sensitive, clean on both)")
TUNE_MOVIES = ("auto: Treg = one movie per day and FOV; RearChr = first pre AND first post run per day and FOV (agreement check)",
               "one movie per condition of each experiment (same as auto)",
               "the first baseline run of each experiment only (optostim: pre only, half the sweep)")
TUNE_MOVIES_ARG = {TUNE_MOVIES[0]: "auto", TUNE_MOVIES[1]: "yes", TUNE_MOVIES[2]: "no"}


STATE_BG = {"done": "#d4edda", "partial": "#fff3cd", "stale": "#f8d7da", "pending": "#fff3cd", "missing": "#e9ecef", "na": "#f3f3f3"}
STATE_MARK = {"done": "✓", "partial": "◐", "stale": "!", "pending": "…", "missing": "○", "na": "·"}


class Stopped(BaseException):          # not an Exception: the steps' own "except Exception" must not swallow a Stop
    """Raised inside a worker when the Stop button was pressed."""


def open_path(p: Path):
    """Open a file or folder with the system viewer."""
    p = Path(p)
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(p))                      # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.Popen(["open", str(p)])
        else:
            subprocess.Popen(["xdg-open", str(p)])
    except Exception as e:
        messagebox.showerror("Open", f"Could not open {p}:\n{e}")


IMAGEJ_GUESSES = (r"C:\Fiji.app\ImageJ-win64.exe", r"D:\Fiji.app\ImageJ-win64.exe", r"C:\ImageJ\ImageJ.exe",
                  r"D:\ImageJ\ImageJ.exe", r"C:\Program Files\ImageJ\ImageJ.exe", r"C:\Program Files\Fiji.app\ImageJ-win64.exe",
                  r"~\Fiji.app\ImageJ-win64.exe", r"~\Desktop\Fiji.app\ImageJ-win64.exe", r"~\ImageJ\ImageJ.exe",
                  r"~\Downloads\Fiji.app\ImageJ-win64.exe", "/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx",
                  "/usr/local/bin/imagej", "/usr/bin/imagej")


def find_imagej() -> str:
    """A likely ImageJ / Fiji executable on this PC, or ''."""
    for g in IMAGEJ_GUESSES:
        p = Path(os.path.expanduser(g))
        if p.exists():
            return str(p)
    for name in ("ImageJ-win64.exe", "ImageJ.exe", "fiji", "imagej"):
        import shutil
        w = shutil.which(name)
        if w:
            return w
    return ""


def open_in_imagej(p: Path, imagej: str, log=print):
    """Open a tiff in ImageJ / Fiji (the executable given in the settings), or
    with the system viewer when none is set."""
    p = Path(p)
    exe = (imagej or "").strip().strip('"')
    if exe and Path(exe).exists():
        try:
            subprocess.Popen([exe, str(p)], cwd=str(Path(exe).parent))
            log(f"opened {p.name} in {Path(exe).name}")
            return
        except Exception as e:
            log(f"could not start {exe}: {e} — opening with the system viewer instead")
    else:
        log("no ImageJ executable set (top row: ImageJ) — opening with the system viewer; set the path to get "
            "the slice labels and the frame slider")
    open_path(p)


def run_status(mouse: Path) -> list[dict]:
    """One row per binned channel-1 tiff under the mouse folder: what the
    brightness step found, what is excluded, what the check said."""
    from afferent2p.aqua_io import aqua_files, find_binned_tifs
    from afferent2p.plan import skip_parts
    rows = []
    folders = sorted({p.parent for p in mouse.rglob("*.sbx")} | {p.parent for p in mouse.rglob("*_bin*s.tif")} | {mouse})
    for d in folders:
        try:
            rel = d.relative_to(mouse).parts
        except ValueError:
            rel = ()
        if skip_parts(rel):
            continue
        for tif in find_binned_tifs(d):
            base = Path(str(tif.with_suffix("")))
            row = {"tif": tif, "run": tif.stem, "frames": "", "artefact": "", "state": "", "excluded": "",
                   "movie": "", "mode": "", "corrected": "", "check": "", "stale": "", "stretches": [],
                   "png": Path(str(base) + "_globalfix.png"),
                   "sheet": Path(str(base) + "_flagged_frames.png"), "csv": Path(str(base) + "_flagged_frames.csv"),
                   "frames_csv": Path(str(base) + "_frames.csv"), "raw": Path(str(base) + "_raw.tif"),
                   "stack_frames": Path(str(base) + "_check_frames.tif"),
                   "stack_movie": Path(str(base) + "_check_movie.tif")}
            gc = None
            if Path(str(base) + ".json").exists():
                try:
                    j = json.loads(Path(str(base) + ".json").read_text())
                    gc = j.get("global_correction")
                    row["frames"] = j.get("frames_recorded") or j.get("frames_out") or ""
                except Exception:
                    pass
            if gc:
                from afferent2p.globalfix import frames_str
                row["artefact"] = gc.get("n_artefact", gc.get("n_flagged", ""))
                row["state"] = gc.get("n_state", 0 if "n_artefact" in gc else "")
                n_excl = gc.get("n_excluded", gc.get("n_flagged", ""))
                ex_list = gc.get("excluded_frames_1based")
                row["excluded"] = (f"{n_excl}: {frames_str(ex_list, max_runs=6)}" if ex_list else str(n_excl))
                row["movie"] = gc.get("frames_in_movie", "")
                row["mode"] = gc.get("bad_frames") or ("filled" if gc.get("replaced_flagged") else "kept")
                old = "n_artefact" not in gc or (bool(gc.get("stretches")) and any("core_frames" not in st for st in gc["stretches"]))
                if old:
                    row["mode"] += " (old version: re-run with force)"
                row["corrected"] = str(gc.get("written", ""))[:16].replace("T", " ")
                row["stretches"] = gc.get("stretches") or []
                if row["stretches"]:
                    n_ex = sum(1 for st in row["stretches"]
                               if st.get("excluded_frames") or (old and str(st.get("decision", st.get("default"))).startswith("ex")))
                    n_own = sum(1 for st in row["stretches"] if st.get("decision") and st.get("decision") != st.get("default"))
                    n_stack = "" if row["stack_frames"].exists() else ", no check stacks yet"
                    row["check"] = (f"{len(row['stretches'])} stretches, {n_ex} with excluded frames"
                                    + (f", {n_own} decided by you" if n_own else "")
                                    + (" — CHECK: old version, whole stretches excluded" if old else n_stack))
                elif row["artefact"] or row["state"]:
                    row["check"] = "no stretch list (old version)"
                else:
                    row["check"] = "nothing flagged"
            else:
                row["mode"] = "not corrected"
            f = aqua_files(tif)
            if f["features"].exists() or f["res"].exists():
                try:
                    newest_res = max(p.stat().st_mtime for k, p in f.items() if k in ("res", "features") and p.exists())
                    row["stale"] = "STALE (movie newer)" if tif.stat().st_mtime > newest_res + 1 else "ok"
                except Exception:
                    row["stale"] = "?"
            else:
                row["stale"] = "none"
            rows.append(row)
    return rows


class ImageViewer(tk.Toplevel):
    """A picture window with zoom and pan (v38): the file at full resolution on a scrollable canvas — Ctrl + wheel or
    the +/- buttons zoom at the cursor, drag to pan, Fit / 1:1, several files to flip through (the tuning figure and the
    'why' pictures), the path in a copyable box, Copy path / Open folder / Open file. Only the visible part is resized,
    so a 2000 x 4000 figure stays quick."""
    ZOOMS = (0.1, 0.15, 0.2, 0.25, 0.33, 0.5, 0.67, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0)

    def __init__(self, master, files, index: int = 0, title: str = "", redraw=None):
        super().__init__(master)
        self.files = [Path(f) for f in files]
        self.index = max(0, min(index, len(self.files) - 1))
        self.redraw = redraw                     # callable(path) -> path or None: make the picture again at higher resolution
        self.im = None; self.photo = None; self.zoom = 1.0; self._render_after = None
        self.title(title or self.files[self.index].name)
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry(f"{min(1400, sw - 80)}x{min(900, sh - 120)}+40+40")
        bar = ttk.Frame(self); bar.pack(fill="x", padx=6, pady=(6, 2))
        ttk.Button(bar, text="◀", width=3, command=lambda: self.step(-1)).pack(side="left")
        self.pick = ttk.Combobox(bar, values=[f.name for f in self.files], state="readonly", width=46)
        self.pick.current(self.index); self.pick.pack(side="left", padx=4)
        self.pick.bind("<<ComboboxSelected>>", lambda e: self.load(self.pick.current()))
        ttk.Button(bar, text="▶", width=3, command=lambda: self.step(1)).pack(side="left")
        ttk.Label(bar, text="   zoom").pack(side="left")
        ttk.Button(bar, text="−", width=3, command=lambda: self.zoom_step(-1)).pack(side="left", padx=(4, 0))
        self.zoom_lbl = ttk.Label(bar, text="100 %", width=6, anchor="center"); self.zoom_lbl.pack(side="left")
        ttk.Button(bar, text="+", width=3, command=lambda: self.zoom_step(1)).pack(side="left")
        ttk.Button(bar, text="Fit", command=self.fit).pack(side="left", padx=(6, 0))
        ttk.Button(bar, text="1:1", command=lambda: self.set_zoom(1.0)).pack(side="left", padx=(2, 0))
        ttk.Button(bar, text="Copy path", command=self.copy_path).pack(side="left", padx=(16, 0))
        ttk.Button(bar, text="Open folder", command=lambda: open_path(self.files[self.index].parent)).pack(side="left", padx=(4, 0))
        ttk.Button(bar, text="Open file", command=lambda: open_path(self.files[self.index])).pack(side="left", padx=(4, 0))
        if self.redraw is not None:
            ttk.Button(bar, text="Redraw at 2× resolution", command=self.do_redraw).pack(side="left", padx=(16, 0))
        prow = ttk.Frame(self); prow.pack(fill="x", padx=6, pady=(0, 2))
        ttk.Label(prow, text="file:").pack(side="left")
        self.path_var = tk.StringVar(value=str(self.files[self.index]))
        ent = ttk.Entry(prow, textvariable=self.path_var, state="readonly"); ent.pack(side="left", fill="x", expand=True, padx=4)
        ent.bind("<Double-Button-1>", lambda e: self.copy_path())
        ttk.Label(prow, text="(double-click to copy)", foreground="#666").pack(side="left")
        body = ttk.Frame(self); body.pack(fill="both", expand=True, padx=6, pady=(0, 2))
        self.canvas = tk.Canvas(body, bg="#2b2b2b", highlightthickness=0, cursor="fleur")
        vs = ttk.Scrollbar(body, orient="vertical", command=self._yview); hs = ttk.Scrollbar(body, orient="horizontal", command=self._xview)
        self.canvas.configure(yscrollcommand=vs.set, xscrollcommand=hs.set)
        vs.pack(side="right", fill="y"); hs.pack(side="bottom", fill="x"); self.canvas.pack(side="left", fill="both", expand=True)
        self.img_id = self.canvas.create_image(0, 0, anchor="nw")
        self.status = ttk.Label(self, text="", foreground="#555"); self.status.pack(anchor="w", padx=8, pady=(0, 6))
        self.canvas.bind("<Configure>", lambda e: self.schedule_render())
        self.canvas.bind("<ButtonPress-1>", lambda e: self.canvas.scan_mark(e.x, e.y))
        self.canvas.bind("<B1-Motion>", self._drag)
        self.canvas.bind("<MouseWheel>", self._wheel)                      # Windows / macOS
        self.canvas.bind("<Shift-MouseWheel>", lambda e: self._wheel(e, horizontal=True))
        self.canvas.bind("<Control-MouseWheel>", lambda e: self._wheel_zoom(e, 1 if e.delta > 0 else -1))
        self.canvas.bind("<Button-4>", lambda e: self._wheel(e, up=True)); self.canvas.bind("<Button-5>", lambda e: self._wheel(e, up=False))
        self.canvas.bind("<Control-Button-4>", lambda e: self._wheel_zoom(e, 1)); self.canvas.bind("<Control-Button-5>", lambda e: self._wheel_zoom(e, -1))
        for key, fn in (("<plus>", lambda e: self.zoom_step(1)), ("<equal>", lambda e: self.zoom_step(1)), ("<KP_Add>", lambda e: self.zoom_step(1)),
                        ("<minus>", lambda e: self.zoom_step(-1)), ("<KP_Subtract>", lambda e: self.zoom_step(-1)),
                        ("<Key-0>", lambda e: self.fit()), ("<Key-1>", lambda e: self.set_zoom(1.0)),
                        ("<Left>", lambda e: self.step(-1)), ("<Right>", lambda e: self.step(1)), ("<Escape>", lambda e: self.destroy())):
            self.bind(key, fn)
        self.canvas.focus_set()
        self.load(self.index, fit=True)

    # ---- files
    def load(self, index: int, fit: bool = True):
        from PIL import Image
        self.index = max(0, min(index, len(self.files) - 1))
        p = self.files[self.index]
        try:
            self.im = Image.open(p).convert("RGB")
        except Exception as e:
            self.status.configure(text=f"could not open {p}: {e}"); return
        self.pick.current(self.index); self.path_var.set(str(p)); self.title(p.name)
        if fit:
            self.after(50, self.fit)
        else:
            self.set_zoom(self.zoom)

    def step(self, d: int):
        if len(self.files) > 1:
            self.load((self.index + d) % len(self.files), fit=True)

    def copy_path(self):
        self.clipboard_clear(); self.clipboard_append(str(self.files[self.index]))
        self.status.configure(text=f"copied: {self.files[self.index]}")

    def do_redraw(self):
        p = self.files[self.index]
        self.status.configure(text="redrawing at higher resolution …"); self.update_idletasks()
        try:
            out = self.redraw(p)
        except Exception as e:
            self.status.configure(text=f"redraw failed: {e}"); return
        if out:
            self.load(self.index, fit=False)
            self.status.configure(text=f"redrawn: {out}")

    # ---- zoom and pan
    def fit(self):
        if self.im is None:
            return
        vw, vh = max(self.canvas.winfo_width(), 50), max(self.canvas.winfo_height(), 50)
        self.set_zoom(min(vw / self.im.width, vh / self.im.height, 1.0))

    def zoom_step(self, d: int, at=None):
        zs = self.ZOOMS
        if d > 0:
            nz = next((z for z in zs if z > self.zoom + 1e-6), zs[-1])
        else:
            nz = next((z for z in reversed(zs) if z < self.zoom - 1e-6), zs[0])
        self.set_zoom(nz, at=at)

    def set_zoom(self, z: float, at=None):
        """at = (canvas_x, canvas_y) widget coordinates to keep under the cursor."""
        if self.im is None:
            return
        z = max(0.02, min(float(z), 8.0))
        c = self.canvas
        if at is None:
            at = (c.winfo_width() / 2, c.winfo_height() / 2)
        ix = (c.canvasx(at[0])) / self.zoom; iy = (c.canvasy(at[1])) / self.zoom      # image point under the cursor
        self.zoom = z
        W, H = int(self.im.width * z), int(self.im.height * z)
        c.configure(scrollregion=(0, 0, W, H))
        # keep the same image point under the cursor
        nx, ny = ix * z - at[0], iy * z - at[1]
        c.xview_moveto(max(0.0, nx / max(W, 1))); c.yview_moveto(max(0.0, ny / max(H, 1)))
        self.zoom_lbl.configure(text=f"{100 * z:.0f} %")
        self.schedule_render()

    def _wheel(self, e, horizontal: bool = False, up=None):
        if up is None:
            up = e.delta > 0
        n = -3 if up else 3
        if horizontal or (e.state & 0x1):
            self.canvas.xview_scroll(n, "units")
        else:
            self.canvas.yview_scroll(n, "units")
        self.schedule_render()
        return "break"

    def _wheel_zoom(self, e, d: int):
        self.zoom_step(d, at=(e.x, e.y))
        return "break"

    def _drag(self, e):
        self.canvas.scan_dragto(e.x, e.y, gain=1); self.schedule_render()

    def _xview(self, *a):
        self.canvas.xview(*a); self.schedule_render()

    def _yview(self, *a):
        self.canvas.yview(*a); self.schedule_render()

    def schedule_render(self):
        if self._render_after is not None:
            try:
                self.after_cancel(self._render_after)
            except Exception:
                pass
        self._render_after = self.after(15, self.render)

    def render(self):
        """Resize only the visible part of the picture (plus a margin) at the current zoom."""
        self._render_after = None
        if self.im is None:
            return
        from PIL import Image, ImageTk
        c = self.canvas; z = self.zoom
        vw, vh = max(c.winfo_width(), 1), max(c.winfo_height(), 1)
        x0, y0 = c.canvasx(0), c.canvasy(0)
        W, H = self.im.width * z, self.im.height * z
        m = 200                                                         # margin in canvas px
        cx0, cy0 = max(0.0, x0 - m), max(0.0, y0 - m)
        cx1, cy1 = min(W, x0 + vw + m), min(H, y0 + vh + m)
        if cx1 <= cx0 or cy1 <= cy0:
            return
        ix0, iy0 = int(cx0 / z), int(cy0 / z); ix1, iy1 = int(np.ceil(cx1 / z)), int(np.ceil(cy1 / z))
        ix1 = min(ix1, self.im.width); iy1 = min(iy1, self.im.height)
        crop = self.im.crop((ix0, iy0, ix1, iy1))
        tw, th = max(1, int(round((ix1 - ix0) * z))), max(1, int(round((iy1 - iy0) * z)))
        if (tw, th) != crop.size:
            crop = crop.resize((tw, th), Image.LANCZOS if z < 1 else Image.BICUBIC)
        self.photo = ImageTk.PhotoImage(crop)
        c.itemconfigure(self.img_id, image=self.photo)
        ox = max(0.0, (vw - W) / 2); oy = max(0.0, (vh - H) / 2)        # a small picture sits in the middle
        c.coords(self.img_id, ix0 * z + ox, iy0 * z + oy)
        self.status.configure(text=f"{self.im.width} × {self.im.height} px — zoom {100 * z:.0f} %  (Ctrl + wheel zooms at the cursor, drag to pan, "
                                   f"+ / − / 0 = fit / 1 = 1:1, ◀ ▶ or the arrow keys flip through the pictures)")


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("afferent2p — AQuA2 route (v45)")
        # v45: as tall as the screen allows (up to 980 px) so that the tables of tab 2 keep their rows; the divider
        # above the log still moves
        sh = self.winfo_screenheight()
        self.geometry(f"1180x{max(700, min(980, sh - 110))}+40+20")
        self.minsize(900, 600)
        self.q: queue.Queue = queue.Queue()
        self.busy = False
        self.stop_evt = threading.Event()
        self.manual_evt = threading.Event()              # "Continue now" of a manual MATLAB step
        self.manual_force = False
        self.manual_win = None
        self._manual_done_fn = None
        self.status_rows: list[dict] = []
        saved = {}
        if SETTINGS.exists():
            try:
                saved = json.loads(SETTINGS.read_text())
            except Exception:
                saved = {}
        from afferent2p.matlab_run import find_matlab
        self.v = {
            "mouse": tk.StringVar(value=saved.get("mouse", "")),
            "aqua": tk.StringVar(value=saved.get("aqua", r"D:\AQuA2")),
            "imagej": tk.StringVar(value=saved.get("imagej") or find_imagej()),
            "matlab": tk.StringVar(value=saved.get("matlab") or find_matlab()),
            "um": tk.StringVar(value=saved.get("um", "0.792")),
            "border": tk.StringVar(value=saved.get("border", "5")),
            "sec": tk.StringVar(value=saved.get("sec", "2")),
            # unpack
            "u_overwrite": tk.BooleanVar(value=False), "u_keep_raw": tk.BooleanVar(value=False),
            "u_gfix": tk.BooleanVar(value=True), "u_behavior": tk.BooleanVar(value=True), "u_register": tk.BooleanVar(value=False),
            # clean
            "c_force": tk.BooleanVar(value=False),
            "c_report": tk.BooleanVar(value=False), "c_min_rel": tk.StringVar(value="0.03"),
            "c_shape_rel": tk.StringVar(value="0.05"), "c_shape_win": tk.StringVar(value="300"),
            "c_min_clean": tk.StringVar(value="10"), "c_win": tk.StringVar(value="60"),
            "c_max_spike": tk.StringVar(value="10"),
            "c_core_pct": tk.StringVar(value="10"), "c_core_ratio": tk.StringVar(value="3"),
            "c_stacks": tk.BooleanVar(value=True), "c_context": tk.StringVar(value="6"),
            "c_frames": tk.StringVar(value=""),
            # tune (tab 3)
            "t_days": tk.StringVar(value=""), "t_runs": tk.StringVar(value=""),
            "t_thr": tk.StringVar(value="1.5, 2, 2.5"), "t_minsize": tk.StringVar(value="20, 50, 100"),
            "t_extra": tk.StringVar(value=""), "t_ninject": tk.StringVar(value="72"),
            "t_max_false": tk.StringVar(value="0.10"), "t_min_fibre": tk.StringVar(value="0.80"),
            "t_apply": tk.BooleanVar(value=False), "t_as_root": tk.StringVar(value=saved.get("t_as_root", "")),
            "t_workers": tk.StringVar(value=saved.get("t_workers", "0")),
            "t_amps": tk.StringVar(value="0.3, 0.6, 1.0"), "t_sizes": tk.StringVar(value="20, 60, 200"),
            "t_movies": tk.StringVar(value=saved.get("t_movies") if saved.get("t_movies") in TUNE_MOVIES else TUNE_MOVIES[0]),
            "t_common": tk.BooleanVar(value=bool(saved.get("t_common", False))), "t_apply_all": tk.BooleanVar(value=False),
            "t_day_rule": tk.StringVar(value=saved.get("t_day_rule") if saved.get("t_day_rule") in DAY_RULES else DAY_RULES[0]),
            "t_tol_recall": tk.StringVar(value="0.05"),
            "t_skip_tuned": tk.BooleanVar(value=True),
            "auto_matlab": tk.BooleanVar(value=bool(saved.get("auto_matlab", False))),
            # curate
            "cur_tune": tk.StringVar(value=""),
            # batch
            "b_fov": tk.StringVar(value=""), "b_overwrite": tk.BooleanVar(value=False),
            "b_opts": tk.StringVar(value=saved.get("b_opts", str(HERE / "default_Opt.csv"))),
            "b_stim_after": tk.StringVar(value=""), "b_sec_tol": tk.StringVar(value="0.25"),
            "b_cfu_min": tk.StringVar(value=str(saved.get("b_cfu_min", CFU_MIN_EVENTS_DEFAULT))),
            "b_cfu_overlap": tk.StringVar(value=str(saved.get("b_cfu_overlap", CFU_OVERLAP_DEFAULT))),
            "b_cfu_merge": tk.StringVar(value=str(saved.get("b_cfu_merge", "0.5"))),
            "b_cfu_aqua": tk.StringVar(value=saved.get("b_cfu_aqua", "")), "b_cfu_redo": tk.BooleanVar(value=False),
            # analyze
            "a_gate": tk.StringVar(value="none"), "a_min_p": tk.StringVar(value="0"),
            "a_baseline": tk.StringVar(value=""), "a_source": tk.StringVar(value="drop-rejected-cfu"),
            "a_keep_flagged": tk.BooleanVar(value=False), "a_no_curation": tk.BooleanVar(value=False),
            "a_no_cfu": tk.BooleanVar(value=False), "a_margin": tk.StringVar(value="1"),
            "a_pupil": tk.BooleanVar(value=True), "a_maps": tk.BooleanVar(value=True), "a_tile": tk.StringVar(value="25"),
            "a_loco": tk.BooleanVar(value=True), "a_glm": tk.BooleanVar(value=True), "a_heat": tk.BooleanVar(value=True),
            "a_movies": tk.BooleanVar(value=True), "a_motion": tk.BooleanVar(value=True),
            # report
            "r_max_events": tk.StringVar(value="0"), "r_clips": tk.StringVar(value="6"),
            "r_no_clips": tk.BooleanVar(value=False), "r_open": tk.BooleanVar(value=True),
            # across animals
            "x_parent": tk.StringVar(value=saved.get("x_parent", "")),
        }
        self._build()
        self.after(150, self.poll)
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    # ------------------------------------------------------------------ layout
    def _build(self):
        top = ttk.LabelFrame(self, text="Animal and settings (shared by every step)")
        top.pack(fill="x", padx=8, pady=(6, 2))
        ttk.Label(top, text="Mouse folder").grid(row=0, column=0, sticky="w", padx=4)
        ttk.Entry(top, textvariable=self.v["mouse"], width=70).grid(row=0, column=1, columnspan=5, sticky="we", padx=4)
        ttk.Button(top, text="Browse…", command=self.pick_mouse).grid(row=0, column=6, padx=4)
        ttk.Label(top, text="AQuA2 folder").grid(row=1, column=0, sticky="w", padx=4)
        ttk.Entry(top, textvariable=self.v["aqua"], width=30).grid(row=1, column=1, sticky="w", padx=4)
        ttk.Label(top, text="µm/px").grid(row=1, column=2, sticky="e")
        ttk.Entry(top, textvariable=self.v["um"], width=7).grid(row=1, column=3, sticky="w", padx=2)
        ttk.Label(top, text="border px").grid(row=1, column=4, sticky="e")
        ttk.Entry(top, textvariable=self.v["border"], width=4).grid(row=1, column=5, sticky="w", padx=2)
        ttk.Label(top, text="s/frame").grid(row=1, column=6, sticky="e")
        ttk.Entry(top, textvariable=self.v["sec"], width=4).grid(row=1, column=7, sticky="w", padx=2)
        ttk.Label(top, text="ImageJ / Fiji").grid(row=2, column=0, sticky="w", padx=4)
        ttk.Entry(top, textvariable=self.v["imagej"], width=70).grid(row=2, column=1, columnspan=5, sticky="we", padx=4)
        ttk.Button(top, text="Browse…", command=self.pick_imagej).grid(row=2, column=6, padx=4)
        ttk.Label(top, text="(opens the check stacks of tab 2)", foreground="#555555").grid(row=2, column=7, sticky="w")
        ttk.Label(top, text="MATLAB").grid(row=3, column=0, sticky="w", padx=4)
        ttk.Entry(top, textvariable=self.v["matlab"], width=70).grid(row=3, column=1, columnspan=5, sticky="we", padx=4)
        ttk.Button(top, text="Browse…", command=self.pick_matlab).grid(row=3, column=6, padx=4)
        ttk.Label(top, text="(matlab.exe, optional: 'start MATLAB from this window'; R2019a+)", foreground="#555555").grid(row=3, column=7, sticky="w")
        top.columnconfigure(1, weight=1)
        run = ttk.Frame(self)
        run.pack(fill="x", padx=8, pady=(2, 0))
        ttk.Button(run, text="▶ Run all (guided)", command=self.run_all).pack(side="left")
        ttk.Label(run, text="settings → batch → MATLAB → analyze → reports;", foreground="#555555").pack(side="left", padx=(6, 2))
        ttk.Checkbutton(run, text="start MATLAB from here (matlab -batch) instead of pasting",
                        variable=self.v["auto_matlab"]).pack(side="left", padx=(6, 2))
        ttk.Checkbutton(run, text="skip the sweep if tuned", variable=self.v["t_skip_tuned"]).pack(side="left", padx=(8, 2))
        ttk.Button(run, text="Stop", command=self.stop_run).pack(side="left", padx=(12, 0))
        ttk.Button(run, text="Fresh start…", command=self.fresh_start_dialog).pack(side="left", padx=(12, 0))
        # where this animal is: one chip per step (click = go to the tab), the next action under it
        from afferent2p.pipeline_state import STEP_NAMES
        strip = ttk.Frame(self)
        strip.pack(fill="x", padx=8, pady=(4, 0))
        self.state_chips = []
        for i, name in enumerate(STEP_NAMES):
            lbl = tk.Label(strip, text=f"{i + 1} {name.split(' ')[0]}  ·", relief="groove", bd=1, padx=6, pady=1, bg=STATE_BG["na"],
                           cursor="hand2")
            lbl.pack(side="left", padx=(0, 3))
            lbl.bind("<Button-1>", lambda e, i=i: self.state_chip_click(i))
            self.state_chips.append(lbl)
        ttk.Button(strip, text="Refresh", width=8, command=self.refresh_state).pack(side="left", padx=(6, 0))
        ttk.Button(strip, text="Details → log", width=13, command=lambda: self.refresh_state(log_it=True)).pack(side="left", padx=(4, 0))
        self.next_lbl = ttk.Label(self, text="next: (choose a mouse folder)", wraplength=1060, justify="left", foreground="#1a3f7a")
        self.next_lbl.pack(fill="x", padx=10, pady=(1, 2))
        self.state = None
        self._state_busy = False
        # what is running right now, and how far it is (v35): fed by the [i/n] markers the steps print
        # and by the result counts of the MATLAB steps; ticks every second while a step runs
        prow = ttk.Frame(self)
        prow.pack(fill="x", padx=8, pady=(0, 2))
        self.prog_bar = ttk.Progressbar(prow, orient="horizontal", mode="determinate", length=260, maximum=100)
        self.prog_bar.pack(side="left")
        self.prog_lbl = tk.Label(prow, text="idle", anchor="w", fg="#1a5e1a", font=("Segoe UI", 9, "bold"))
        self.prog_lbl.pack(side="left", fill="x", expand=True, padx=8)
        self.prog = {"i": 0, "n": 0, "item": "", "sub": ""}
        self.run_name = ""
        self.run_t0 = 0.0
        self.run_chip = None

        # the tabs above, the log below, in a pane with a draggable divider (v36): the log never disappears on a
        # small screen and a tall tab can be given more room by dragging the divider
        self.pane = ttk.PanedWindow(self, orient="vertical")
        self.pane.pack(fill="both", expand=True, padx=8, pady=4)
        self.nb = ttk.Notebook(self.pane)
        bottom = ttk.Frame(self.pane)
        self.pane.add(self.nb, weight=4)
        self.pane.add(bottom, weight=1)
        hdr = ttk.Frame(bottom); hdr.pack(fill="x", pady=(4, 0))
        ttk.Label(hdr, text="Log").pack(side="left")
        self.state_lbl = ttk.Label(hdr, text="idle", foreground="#1a5e1a")
        self.state_lbl.pack(side="left", padx=10)
        ttk.Button(hdr, text="Clear", command=self.clear_log).pack(side="right")
        ttk.Button(hdr, text="Stop", command=self.stop_run).pack(side="right", padx=8)
        ttk.Label(hdr, text="(drag the divider above to make the log or the tab taller)", foreground="#666").pack(side="right", padx=8)
        self.log = tk.Text(bottom, height=6, state="disabled", wrap="word", font=("Consolas", 9))
        sb = ttk.Scrollbar(bottom, command=self.log.yview)
        self.log.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.log.pack(fill="both", expand=True)

        self.after(250, self._place_divider)
        self._tab_unpack()
        self._tab_clean()
        self._tab_settings()
        self._tab_batch()
        self._tab_matlab()
        self._tab_curate()
        self._tab_analyze()
        self._tab_report()

    def _tab_unpack(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="1  Unpack")
        ttk.Label(f, text="Scanbox .sbx -> binned tiff (mean of 31 frames = 1.984 s), meanproj, .json, behaviour csv "
                          "(running from the quadrature file, pupil from the DeepLabCut eye csv *eye*DLC*.csv — both must "
                          "be in the run folder), then the brightness step (level flattened, bad frames flagged — not "
                          "removed — with the QC png, frames csv and the check sheet). Runs that already have a tiff are "
                          "skipped. 'Rebuild behaviour csv' writes only the behaviour csv again for existing movies, e.g. "
                          "when the eye csv arrived after unpacking (the analysis also does this by itself).",
                  wraplength=1000, justify="left").pack(anchor="w", padx=8, pady=6)
        o = ttk.Frame(f); o.pack(anchor="w", padx=8)
        ttk.Checkbutton(o, text="overwrite existing tiffs", variable=self.v["u_overwrite"]).pack(side="left", padx=4)
        ttk.Checkbutton(o, text="brightness step (flatten + flag bad frames)", variable=self.v["u_gfix"]).pack(side="left", padx=4)
        ttk.Checkbutton(o, text="keep the uncorrected movie too (_raw.tif, twice the disk)", variable=self.v["u_keep_raw"]).pack(side="left", padx=4)
        ttk.Checkbutton(o, text="behaviour csv", variable=self.v["u_behavior"]).pack(side="left", padx=4)
        ttk.Checkbutton(o, text="register the raw frames before binning (rigid, legacy-style; +minutes per run; decide with reg_compare.py)",
                        variable=self.v["u_register"]).pack(side="left", padx=4)
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=8)
        ttk.Button(b, text="Run unpack on the mouse folder", command=self.run_unpack).pack(side="left")
        ttk.Button(b, text="Unpack one run folder…", command=self.run_unpack_one).pack(side="left", padx=8)
        ttk.Button(b, text="Rebuild behaviour csv (running + pupil) for existing movies",
                   command=self.run_behavior_only).pack(side="left", padx=8)
        lf = ttk.LabelFrame(f, text="Which run is which FOV / condition — run_labels.csv and day_labels.csv (v44)")
        lf.pack(anchor="w", fill="x", padx=8, pady=(10, 4))
        ttk.Label(lf, text="The rule: spontaneous (Treg) mouse = run 001/002/003 of a day are FOV1/2/3, condition d<date>; optostim mouse = "
                           "FOV1, FOV2, FOV1, FOV2 pre then post. When a day breaks it (a lost FOV, a shifted plane recaptured under the next "
                           "number, a second injection with pre + post recordings on one day), write the exception in <mouse>\\run_labels.csv: "
                           "one row per run — fov, condition (d<date>_pre / d<date>_post for a same-day pair), use = 0 to leave a run out, note. "
                           "'Create / complete' writes the file with every run and the rule's values (an existing file keeps its rows; new runs "
                           "are appended); edit it in Excel or Notepad, then 'Show the plan' prints what every step will use. day_labels.csv "
                           "(day,label: 260925,Baseline / 261010,D14 — a row 261010_post,… names a same-day condition) labels the days in the "
                           "analysis and the reports. The batch (tab 4) reads both; rebuild it after a change.",
                  wraplength=1000, justify="left", foreground="#444444").pack(anchor="w", padx=6, pady=(2, 4))
        rb = ttk.Frame(lf); rb.pack(anchor="w", padx=6, pady=(0, 6))
        ttk.Button(rb, text="Create / complete run_labels.csv", command=self.run_labels_make).pack(side="left")
        ttk.Button(rb, text="Open run_labels.csv", command=lambda: self.open_mouse_file("run_labels.csv")).pack(side="left", padx=6)
        ttk.Button(rb, text="Open day_labels.csv", command=lambda: self.open_mouse_file("day_labels.csv", header="day,label")).pack(side="left", padx=6)
        ttk.Button(rb, text="Show the plan (log)", command=self.run_labels_plan).pack(side="left", padx=6)

    # ------------------------------------------------------------------ tab 1: run / day labels (v44)
    def run_labels_make(self):
        m = self.mouse()
        if m is None:
            return
        import contextlib, io
        from afferent2p.plan import write_run_labels
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                p = write_run_labels(m)
        except Exception as e:
            messagebox.showerror("run_labels.csv", str(e)); return
        for line in buf.getvalue().splitlines():
            self.logln(line)
        self.logln(f"{p} — edit fov / condition / use / note, then 'Show the plan'; rebuild the batch (tab 4) afterwards")

    def open_mouse_file(self, name: str, header: str | None = None):
        m = self.mouse()
        if m is None:
            return
        p = m / name
        if not p.exists():
            if header is None:
                messagebox.showinfo(name, f"{p}\n\ndoes not exist yet — press 'Create / complete run_labels.csv' first."); return
            p.write_text(header + "\n", encoding="utf-8")
            self.logln(f"created an empty {p}")
        open_path(p)

    def run_labels_plan(self):
        m = self.mouse()
        if m is None:
            return
        import contextlib, io
        from afferent2p.plan import plan_mouse, describe
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exps = plan_mouse(m)
                txt = describe(exps)
        except Exception as e:
            messagebox.showerror("plan", str(e)); return
        for line in buf.getvalue().splitlines() + txt.splitlines():
            self.logln(line)
        if not exps:
            self.logln("no runs found")

    def _tab_clean(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="2  Clean & check")
        ttk.Label(f, text="Brightness step on existing tiffs: the level is flattened, bad frames are FLAGGED, nothing is removed. "
                          "Frames of a flagged stretch that differ clearly from the kept frames around them (>= 10 % of the frame and "
                          ">= 3 x the clean floor) are its FLAGGED FRAMES, excluded from the statistics by default; tails, mild stretches "
                          "and slow state changes are kept. CHECK BY EYE with the buttons under the table — 'Show check sheet': every "
                          "flagged stretch as pictures (each frame + its difference from the clean neighbours), in this window; 'Frame "
                          "stack' (or double-click the row): the same frames with 6 frames of context in ImageJ, slices labelled; 'QC "
                          "plot': the brightness trace. Then decide per stretch under the lower table. Force = re-flag with the current "
                          "rules; Report only = preview.",
                  wraplength=1130, justify="left").pack(anchor="w", padx=8, pady=(6, 2))
        o = ttk.Frame(f); o.pack(anchor="w", padx=8)
        ttk.Checkbutton(o, text="force (re-do corrected tiffs)", variable=self.v["c_force"]).pack(side="left")
        ttk.Checkbutton(o, text="report only (movie untouched)", variable=self.v["c_report"]).pack(side="left", padx=8)
        ttk.Checkbutton(o, text="write check stacks", variable=self.v["c_stacks"]).pack(side="left", padx=8)
        for label, key, w in (("flagged frame: min % of frame", "c_core_pct", 4), ("x clean floor", "c_core_ratio", 4),
                              ("stack context frames", "c_context", 4)):
            ttk.Label(o, text=label).pack(side="left", padx=(8, 2))
            ttk.Entry(o, textvariable=self.v[key], width=w).pack(side="left")
        o2 = ttk.Frame(f); o2.pack(anchor="w", padx=8)
        ttk.Label(o2, text="flagging rules:", foreground="#555555").pack(side="left")
        for label, key, w in (("level floor", "c_min_rel", 6), ("shape floor", "c_shape_rel", 6),
                              ("spike max s", "c_max_spike", 5), ("shape window s", "c_shape_win", 6),
                              ("level window s", "c_win", 5), ("min clean s", "c_min_clean", 5)):
            ttk.Label(o2, text=label).pack(side="left", padx=(8, 2))
            ttk.Entry(o2, textvariable=self.v[key], width=w).pack(side="left")
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=4)
        ttk.Button(b, text="Run on the mouse folder", command=self.run_clean).pack(side="left")
        ttk.Button(b, text="Run on selected run(s)", command=lambda: self.run_clean(selected=True)).pack(side="left", padx=6)
        ttk.Button(b, text="Refresh", command=self.refresh_status).pack(side="left", padx=6)
        ttk.Label(b, text="  movie actions:").pack(side="left", padx=(12, 2))
        ttk.Button(b, text="Drop excluded frames from selected movie(s)…", command=self.run_drop).pack(side="left")
        ttk.Button(b, text="Restore selected from _raw.tif…", command=self.run_restore).pack(side="left", padx=6)
        # v45: the check buttons sit directly under the recordings table and the decision buttons directly under the
        # stretch table, and both rows are packed FIRST (side=bottom) so that a short window shrinks the tables instead
        # of dropping the rows — Tk's packer leaves out the LAST-packed children when the room runs out, and from v36
        # on (Run all row, status strip, log pane) the tab was short enough for the two rows to vanish silently.
        pane = ttk.Panedwindow(f, orient="vertical"); pane.pack(fill="both", expand=True, padx=8, pady=(2, 4))
        top = ttk.Frame(pane); pane.add(top, weight=3)
        b2 = ttk.Frame(top); b2.pack(side="bottom", anchor="w", fill="x", pady=(3, 2))
        ttk.Label(b2, text="check the selected recording by eye:").pack(side="left", padx=(0, 2))
        ttk.Button(b2, text="Show check sheet", command=self.open_sheet).pack(side="left")
        ttk.Button(b2, text="Frame stack (ImageJ)", command=lambda: self.open_stack("frames")).pack(side="left", padx=(8, 0))
        ttk.Button(b2, text="Composite stack", command=lambda: self.open_stack("movie")).pack(side="left", padx=4)
        ttk.Button(b2, text="Make check stacks", command=self.make_stacks).pack(side="left", padx=4)
        ttk.Button(b2, text="QC plot", command=lambda: self.open_selected("png")).pack(side="left", padx=(12, 4))
        ttk.Button(b2, text="Run folder", command=lambda: self.open_selected("folder")).pack(side="left", padx=4)
        cols = ("run", "frames", "artefact", "state", "excluded", "movie", "mode", "corrected", "check", "stale")
        self.tree = ttk.Treeview(top, columns=cols, show="headings", height=7, selectmode="extended")
        heads = {"run": ("recording", 270), "frames": ("frames", 50), "artefact": ("artefact", 58), "state": ("state", 46),
                 "excluded": ("excluded frames", 150), "movie": ("in movie", 58), "mode": ("frames", 70),
                 "corrected": ("corrected on", 108), "check": ("stretches / check", 250), "stale": ("AQuA2 result", 95)}
        for c in cols:
            self.tree.heading(c, text=heads[c][0])
            self.tree.column(c, width=heads[c][1], anchor="w", stretch=(c in ("run", "check")))
        self.tree.tag_configure("check", foreground="#b30000")
        self.tree.tag_configure("stale", foreground="#8a5a00")
        self.tree.tag_configure("old", foreground="#555555")
        sb = ttk.Scrollbar(top, command=self.tree.yview); self.tree.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y"); self.tree.pack(fill="both", expand=True)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.fill_stretches())
        self.tree.bind("<Double-1>", lambda e: self.open_stack("frames"))
        bottom = ttk.Frame(pane); pane.add(bottom, weight=2)
        hdr = ttk.Frame(bottom); hdr.pack(fill="x")
        self.stretch_lbl = ttk.Label(hdr, text="flagged stretches of the selected recording — decide what the statistics do with each:")
        self.stretch_lbl.pack(side="left")
        b3 = ttk.Frame(bottom); b3.pack(side="bottom", anchor="w", fill="x", pady=(3, 2))
        ttk.Label(b3, text="selected stretch(es):").pack(side="left", padx=(0, 2))
        ttk.Button(b3, text="Exclude flagged frames", command=lambda: self.decide("exclude")).pack(side="left")
        ttk.Button(b3, text="Exclude whole stretch", command=lambda: self.decide("exclude-all")).pack(side="left", padx=4)
        ttk.Button(b3, text="Keep all", command=lambda: self.decide("keep")).pack(side="left")
        ttk.Button(b3, text="Back to default", command=lambda: self.decide(None)).pack(side="left", padx=4)
        ttk.Label(b3, text="   or exactly these frames:").pack(side="left", padx=(6, 2))
        ttk.Entry(b3, textvariable=self.v["c_frames"], width=13).pack(side="left")
        ttk.Button(b3, text="Exclude these frames", command=self.decide_frames).pack(side="left", padx=4)
        ttk.Label(b3, text="(e.g. 68,76-77)", foreground="#555555").pack(side="left")
        scols = ("id", "cls", "n", "sec", "level", "shape", "flagged", "verdict", "decision", "excl")
        self.stree = ttk.Treeview(bottom, columns=scols, show="headings", height=6, selectmode="extended")
        sheads = {"id": ("frames", 70), "cls": ("class", 62), "n": ("n", 34), "sec": ("s", 36), "level": ("level %", 52),
                  "shape": ("shape %", 54), "flagged": ("flagged frames", 110), "verdict": ("check: vs kept neighbours", 280),
                  "decision": ("decision", 150), "excl": ("excluded frames", 120)}
        for c in scols:
            self.stree.heading(c, text=sheads[c][0])
            self.stree.column(c, width=sheads[c][1], anchor="w", stretch=(c == "verdict"))
        self.stree.tag_configure("artefact", foreground="#b30000")
        self.stree.tag_configure("state", foreground="#a05a00")
        self.stree.tag_configure("check", background="#ffe9e9")
        ssb = ttk.Scrollbar(bottom, command=self.stree.yview); self.stree.configure(yscrollcommand=ssb.set)
        ssb.pack(side="right", fill="y"); self.stree.pack(fill="both", expand=True)

    def _tab_settings(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="3  Detection settings")
        ttk.Label(f, text="Chooses the AQuA2 detection settings (thrARScl / minSize; any other variable can join the grid) by measurement "
                          "instead of by eye (aqua_tune.py): for the chosen movies it writes a NULL copy (time series scrambled — every event "
                          "found there is false) and an INJECTED copy (synthetic transients of known amplitude / size / duration on the fibres), "
                          "a MATLAB script runs every grid setting on the three movies, and 'Report' scores each setting ON THE FIBRES (recall of "
                          "the injected transients, false-alarm share, events/min on and off the fibres) and picks the most sensitive clean one — "
                          "per day for a RearChr mouse (pre and post are both swept; v36: each gets the setting that detects with the SAME measured "
                          "recall and false-alarm rate as the other — one setting when it serves both; a day without a clean matching pair is marked "
                          "NOT COMPARABLE), one common setting per FOV for a Treg mouse (it compares its days with each other). 'Apply' writes "
                          "aqua_opts\\<experiment>[_pre|_post]_Opt.csv, which tab 4 uses from then on; a FOV without any file gets default_Opt.csv. "
                          "About 1.4 min per job ('MATLAB jobs at once' runs several). The pictures behind the choice: 'Why this setting?' after "
                          "3 Report. Details: GUIDE.md, Step 3.",
                  wraplength=1120, justify="left").pack(anchor="w", padx=8, pady=6)
        o = ttk.Frame(f); o.pack(anchor="w", padx=8, pady=2)
        ttk.Label(o, text="movies").pack(side="left", padx=(8, 2))
        ttk.Combobox(o, textvariable=self.v["t_movies"], values=TUNE_MOVIES, width=70, state="readonly").pack(side="left")
        o1 = ttk.Frame(f); o1.pack(anchor="w", padx=8, pady=2)
        for label, key, w in (("days (empty = all)", "t_days", 12), ("runs (e.g. 260921_002)", "t_runs", 14),
                              ("max false-alarm share", "t_max_false", 5), ("min on-fibre fraction", "t_min_fibre", 5)):
            ttk.Label(o1, text=label).pack(side="left", padx=(8, 2))
            ttk.Entry(o1, textvariable=self.v[key], width=w).pack(side="left")
        o2 = ttk.Frame(f); o2.pack(anchor="w", padx=8, pady=2)
        for label, key, w in (("thrARScl values", "t_thr", 16), ("minSize values", "t_minsize", 14),
                              ("other grid, e.g. sigThr=3,3.5,4", "t_extra", 22)):
            ttk.Label(o2, text=label).pack(side="left", padx=(8, 2))
            ttk.Entry(o2, textvariable=self.v[key], width=w).pack(side="left")
        o3 = ttk.Frame(f); o3.pack(anchor="w", padx=8, pady=2)
        for label, key, w in (("injected dF/F (AQuA2 reports ~half)", "t_amps", 12), ("sizes px", "t_sizes", 12),
                              ("injected transients", "t_ninject", 4), ("MATLAB jobs at once (0 = auto)", "t_workers", 3)):
            ttk.Label(o3, text=label).pack(side="left", padx=(8, 2))
            ttk.Entry(o3, textvariable=self.v[key], width=w).pack(side="left")
        o4 = ttk.Frame(f); o4.pack(anchor="w", padx=8, pady=2)
        ttk.Checkbutton(o4, text="apply: write the chosen settings", variable=self.v["t_apply"]).pack(side="left", padx=(8, 2))
        ttk.Checkbutton(o4, text="one common setting per FOV for all days (Treg: always; RearChr: per day unless ticked)", variable=self.v["t_common"]).pack(side="left", padx=(8, 2))
        ttk.Checkbutton(o4, text="trade-off if the rule is not met (Run all: always)", variable=self.v["t_apply_all"]).pack(side="left", padx=(8, 2))
        o4b = ttk.Frame(f); o4b.pack(anchor="w", padx=8, pady=2)
        ttk.Label(o4b, text="RearChr day (pre + post tuned):").pack(side="left", padx=(8, 2))
        ttk.Combobox(o4b, textvariable=self.v["t_day_rule"], values=DAY_RULES, state="readonly", width=64).pack(side="left")
        ttk.Label(o4b, text="recall gap ≤").pack(side="left", padx=(8, 2))
        ttk.Entry(o4b, textvariable=self.v["t_tol_recall"], width=5).pack(side="left")
        o5 = ttk.Frame(f); o5.pack(anchor="w", padx=8, pady=2)
        ttk.Label(o5, text="mouse folder as the MATLAB PC sees it (empty = same path; e.g. Z:\\... when MATLAB runs on another PC)").pack(side="left", padx=(8, 2))
        ttk.Entry(o5, textvariable=self.v["t_as_root"], width=36).pack(side="left")
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=(8, 2))
        ttk.Button(b, text="1  Plan", command=self.run_tune_plan).pack(side="left")
        ttk.Button(b, text="2  Copy the MATLAB command", command=self.copy_tune_matlab).pack(side="left", padx=(8, 2))
        ttk.Button(b, text="or run the sweep from here", command=self.run_tune_matlab).pack(side="left", padx=(0, 8))
        ttk.Button(b, text="Progress", command=self.check_tune_progress).pack(side="left", padx=8)
        ttk.Button(b, text="3  Report (score + recommend)", command=self.run_tune_report).pack(side="left", padx=8)
        ttk.Button(b, text="4  Apply → next: 4 Batch", command=self.rebuild_batch_tuned).pack(side="left", padx=8)
        b2 = ttk.Frame(f); b2.pack(anchor="w", padx=8, pady=(2, 6))
        ttk.Button(b2, text="Figure", command=self.show_tune_figure).pack(side="left")
        ttk.Button(b2, text="Folder", command=lambda: self.open_mouse_sub("aqua_tune")).pack(side="left", padx=8)
        ttk.Button(b2, text="Recommendation", command=lambda: self.open_mouse_sub("aqua_tune/tune_recommendation.txt")).pack(side="left")
        ttk.Button(b2, text="Why this setting? (pictures)", command=self.open_tune_why).pack(side="left", padx=8)
        ttk.Button(b2, text="Extend the grid (as suggested)", command=self.extend_grid).pack(side="left")
        ttk.Button(b2, text="Pick by eye…", command=self.open_picker).pack(side="left", padx=8)
        b3 = ttk.Frame(f); b3.pack(anchor="w", padx=8, pady=(0, 4))
        ttk.Label(b3, text="by eye in the AQuA2 GUI:", foreground="#555555").pack(side="left", padx=(0, 6))
        ttk.Button(b3, text="Copy the aqua_gui command", command=self.copy_aqua_gui_cmd).pack(side="left")
        ttk.Button(b3, text="Hand result → reference (MATLAB command)", command=self.hand_result_reference).pack(side="left", padx=8)
        ttk.Button(b3, text="List settings files", command=self.list_opts).pack(side="left")
        ttk.Label(f, text="By eye: open a recording in aqua_gui (1.984 s/frame, 0.792 µm/px, border 5; detectGlo off), Act → Temporal → Spatial until the "
                          "events look right, Detect → SaveOpts into the mouse folder as FOV1_Opt.csv (any name containing FOV1). A result SAVED from the "
                          "AQuA2 GUI (its default folder: a subfolder named after the movie) becomes the REFERENCE of that run once exported "
                          "('Hand result → reference' gives the MATLAB command): the report then says, per setting, how many of the hand-accepted "
                          "events it recovers and how many it adds. 'Extend the grid' adds what the last report suggested (grid verdict) and plans "
                          "again — finished jobs are kept, only the new combinations are swept.",
                  wraplength=1120, justify="left", foreground="#555555").pack(anchor="w", padx=8, pady=(0, 2))
        self.tune_lbl = ttk.Label(f, text="", justify="left", wraplength=1000, foreground="#333333")
        self.tune_lbl.pack(anchor="w", padx=8, pady=4)
        self.tune_box = tk.Text(f, height=5, state="disabled", font=("Consolas", 9))
        self.tune_box.pack(fill="both", expand=True, padx=8, pady=(0, 6))

    def _tab_batch(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="4  Batch")
        ttk.Label(f, text="Builds <mouse>\\aqua_batch\\ (batch.csv, parameters_for_batch.csv, run_aqua_batch.m, run_cfu_batch.m, "
                          "run_aqua_and_cfu.m, manifest.json) from the FOV settings files: aqua_opts\\<experiment>_Opt.csv (what tab 3 "
                          "writes), FOV<k>_Opt.csv, a file in the run folder, else the default below (else the pipeline's own "
                          "default_Opt.csv — so a new animal is always queued). Read the printout: which settings file each FOV got, "
                          "what was queued. Tiffs with an _AQuA2.mat are skipped unless overwrite is ticked or their settings file is "
                          "newer than the result (then they are analysed again). Next: tab 5.",
                  wraplength=1000, justify="left").pack(anchor="w", padx=8, pady=6)
        o = ttk.Frame(f); o.pack(anchor="w", padx=8)
        ttk.Label(o, text="only FOV(s), e.g. 2 or 1 3 (empty = all)").pack(side="left")
        ttk.Entry(o, textvariable=self.v["b_fov"], width=8).pack(side="left", padx=4)
        ttk.Checkbutton(o, text="overwrite (re-queue tiffs that have results)", variable=self.v["b_overwrite"]).pack(side="left", padx=10)
        ttk.Checkbutton(o, text="redo the CFU step of every result (results made with other CFU rules are redone anyway)",
                        variable=self.v["b_cfu_redo"]).pack(side="left", padx=10)
        o1 = ttk.Frame(f); o1.pack(anchor="w", padx=8, pady=2)
        ttk.Label(o1, text="settings for a FOV that has no *_Opt.csv of its own (tab 3 then chooses thrARScl / minSize by measurement)").pack(side="left")
        o1b = ttk.Frame(f); o1b.pack(anchor="w", padx=8, pady=(0, 2))
        ttk.Entry(o1b, textvariable=self.v["b_opts"], width=70).pack(side="left")
        ttk.Button(o1b, text="Browse…", command=self.pick_opts).pack(side="left", padx=6)
        o2 = ttk.Frame(f); o2.pack(anchor="w", padx=8, pady=4)
        for label, key, w in (("stim after run (RearChr, if not half)", "b_stim_after", 4), ("s/frame tolerance", "b_sec_tol", 5),
                              ("CFU min events (AQuA2 default 3)", "b_cfu_min", 4),
                              ("CFU overlap (pipeline 0.3, AQuA2 default 0.5)", "b_cfu_overlap", 5)):
            ttk.Label(o2, text=label).pack(side="left", padx=(8, 2))
            ttk.Entry(o2, textvariable=self.v[key], width=w).pack(side="left")
        ttk.Label(f, text="CFU = AQuA2's consensus functional unit: one place that produces repeated events. The pipeline applies AQuA2's "
                          "published rule: events whose weighted footprints overlap by at least 'CFU overlap' form one unit (single-linkage "
                          "clustering, cfu.CFU_tmp_function). AQuA2's own default is 0.5; the pipeline's default is 0.3 — optimised Nov 2026 on all 96 "
                          "runs of the three RearChr animals: bouton events of different extent at one spot overlap only 0.25–0.45, so at 0.5 a "
                          "spot firing seven times came out as seven one-event units; 0.3 is the threshold at which AQuA2's clustering best "
                          "reproduces same-location repeats (centres coinciding) without chaining distinct places (GUIDE.md 'CFU rules'). "
                          "'CFU min events' 1 keeps single-event places as units (AQuA2's default is 3). Changing either value makes the CFU "
                          "batch redo every result's CFU file — and the curation of tab 6 has to be redone, its decisions refer to the old units "
                          "(v43: the curation window sees that, ignores them, starts those runs clean and keeps the old file as "
                          "*_superseded_<date>.json when you save).",
                  wraplength=1000, justify="left", foreground="#444444").pack(anchor="w", padx=8, pady=(2, 0))
        o3 = ttk.Frame(f); o3.pack(anchor="w", padx=8, pady=(6, 0))
        ttk.Label(o3, text="Other AQuA2 copy for the CFU step only (normally empty = the AQuA2 folder of the top row). A copy from 28 Aug 2026 on "
                           "applies AQuA2's NEWER rules instead — watershed splitting of every cluster into basins plus a post-merge of correlated "
                           "neighbours (r ≥ 0.85 + 0.15 × the post-merge control) — which on this data triples the unit count and puts one event "
                           "into several units; the pipeline does not use them unless a copy is named here.",
                  wraplength=1000, justify="left", foreground="#666666").pack(anchor="w")
        o3b = ttk.Frame(f); o3b.pack(anchor="w", padx=8, pady=(0, 2))
        ttk.Entry(o3b, textvariable=self.v["b_cfu_aqua"], width=70).pack(side="left")
        ttk.Button(o3b, text="Browse…", command=self.pick_cfu_aqua).pack(side="left", padx=6)
        ttk.Label(o3b, text="post-merge control 0–1 (that copy only)").pack(side="left", padx=(12, 2))
        ttk.Entry(o3b, textvariable=self.v["b_cfu_merge"], width=5).pack(side="left")
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=8)
        ttk.Button(b, text="Build batch files", command=self.run_batch).pack(side="left")
        ttk.Button(b, text="Open aqua_batch folder", command=lambda: self.open_mouse_sub("aqua_batch")).pack(side="left", padx=8)
        ttk.Button(b, text="Check the overlap threshold (from the CFU results)", command=self.run_overlap_check).pack(side="left", padx=8)
        ttk.Label(f, text="The check (v42) re-asks the threshold question on THIS animal's own results after the CFU batch has run: AQuA2's pairwise "
                          "event overlaps (the batch writes them to <run>_AQuA2_cfu_pre.mat) are clustered at every threshold 0.10–0.50, the units "
                          "are compared with labels that do not use the overlap (same place = the two events' centres within half the smaller "
                          "radius; two places = each centre outside the other's extent), per day × FOV and pooled. It writes "
                          "aqua_batch\\cfu_overlap_check.txt / .csv / .png, says OK or CHECK per experiment and changes nothing — the threshold "
                          "stays the field above (one number for the dataset) unless you change it there.",
                  wraplength=1000, justify="left", foreground="#666666").pack(anchor="w", padx=8, pady=(0, 4))

    def _tab_matlab(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="5  MATLAB")
        ttk.Label(f, text="One paste: run_aqua_and_cfu.m runs the AQuA2 batch (writes _AQuA2.mat / _AQuA2_Ch1.csv next to each tiff) and then "
                          "the CFU batch (_res_cfu.mat, _events.csv, _meanimg.tif). Copy the command into a MATLAB window (any MATLAB that can "
                          "see the mouse folder; the scripts find AQuA2 on MATLAB's path when the folder in the top row is not visible there) "
                          "— or start MATLAB from here without its desktop (needs the MATLAB field in the top row; the output appears in the log and "
                          "in aqua_batch\\run_*_log.txt). Results whose CFU file is newer than the AQuA2 file are skipped, so a partial re-run "
                          "costs only that FOV. 'Check' says what has been produced so far. Next: tab 6 (curate) or tab 7 (analyze).",
                  wraplength=1000, justify="left").pack(anchor="w", padx=8, pady=6)
        self.matlab_box = tk.Text(f, height=7, font=("Consolas", 10))
        self.matlab_box.pack(fill="x", padx=8)
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=6)
        ttk.Button(b, text="Show the commands for this mouse", command=self.show_matlab).pack(side="left")
        ttk.Button(b, text="Copy the one-paste command", command=self.copy_matlab).pack(side="left", padx=8)
        ttk.Button(b, text="Check what MATLAB has produced", command=self.check_batch_done).pack(side="left", padx=8)
        b2 = ttk.Frame(f); b2.pack(anchor="w", padx=8, pady=(0, 6))
        ttk.Label(b2, text="start MATLAB from here (matlab -batch):", foreground="#555555").pack(side="left", padx=(0, 6))
        ttk.Button(b2, text="Run both now", command=lambda: self.run_batch_matlab("both")).pack(side="left")
        ttk.Button(b2, text="AQuA2 batch only", command=lambda: self.run_batch_matlab("aqua")).pack(side="left", padx=8)
        ttk.Button(b2, text="CFU batch only", command=lambda: self.run_batch_matlab("cfu")).pack(side="left", padx=8)
        ttk.Button(b2, text="Which results exist (frames table)", command=self.refresh_status_log).pack(side="left", padx=8)
        self.matlab_lbl = ttk.Label(f, text="", justify="left", wraplength=1000, foreground="#333333")
        self.matlab_lbl.pack(anchor="w", padx=8, pady=4)

    def _tab_curate(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="6  Curate")
        ttk.Label(f, text="Opens the CFU curation window for this mouse (cfu_curate.py): the mean image with the CFU outlines "
                          "(each CFU its own colour; overlapping ones never share one; the list on the left picks one of "
                          "several overlapping CFUs) next to the MOVIE with the same outlines — by default AQuA2's own REGISTERED "
                          "movie out of the _AQuA2.mat (the frames it detected on, same coordinates as the maps; 'as recorded' shows "
                          "the raw binned tiff instead) — space plays / pauses, arrows step frames, [ ] jump to the events of the "
                          "selected CFU, raw or dF/F view, scroll wheel zooms both images; move the zoomed view with right-drag, "
                          "'move view' (h) + left-drag or Ctrl+arrows. Click or pick a CFU → k real / x noise; drag a box to select several "
                          "at once (shift-drag adds, ctrl-click toggles, Esc clears) and k / x act on all of them; a + vertices + Enter "
                          "adds a region; tab next CFU; n/p runs; s saves; 'Propagate to this FOV' copies decisions to "
                          "matched CFUs of the other runs. The 'looks' column says what each CFU looks like against the fibre mask "
                          "(fibre 85 % / off-fibre 5 % / blob) — 'Select suggested' (u) picks the non-fibre ones so one x marks "
                          "them noise, f shows the mask; it is a suggestion, look before pressing x. Decisions live in "
                          "_res_cfu_curated.json next to the results. After a CFU re-run with other rules (another overlap "
                          "threshold, tab 4) the units are renumbered: decisions made for the earlier result are recognised and "
                          "ignored (the run starts clean, the title says so) and kept as *_superseded_<date>.json on the first save; "
                          "'List CFU results and decisions' says how many runs are curated under the current rules.",
                  wraplength=1000, justify="left").pack(anchor="w", padx=8, pady=6)
        o = ttk.Frame(f); o.pack(anchor="w", padx=8, pady=2)
        ttk.Label(o, text="overlay the events of tuning setting # (tab 3; empty = none)").pack(side="left", padx=(0, 4))
        ttk.Entry(o, textvariable=self.v["cur_tune"], width=4).pack(side="left")
        ttk.Label(o, text="— magenta boxes T1.. with onset and trace, read-only; the setting numbers are in tune_recommendation.txt",
                  foreground="#555555").pack(side="left", padx=6)
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=6)
        ttk.Button(b, text="Open curation window", command=self.run_curate_gui).pack(side="left")
        ttk.Button(b, text="List CFU results and decisions (log)", command=self.run_curate_list).pack(side="left", padx=8)

    def _tab_analyze(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="7  Analyze")
        ttk.Label(f, text="Reads every AQuA2 result of the mouse, converts frame numbers back to recorded frames, drops events in "
                          "excluded stretches (your decisions in tab 2) and events of CFUs marked noise (tab 6), then writes <mouse>\\aqua_analysis\\: "
                          "(a) the general representation — every session on the day's clock with running and pupil, activity by "
                          "running / resting state with event overlays, ECDFs, acute (post/pre) and chronic (day#/day 1) change of every "
                          "measure; (b) the CFU x time heatmaps (identities followed across runs, rows read from the movies) with up / down "
                          "/ unchanged CFUs; (c) the locomotion-evoked activity (Yu's trial pipeline: strict and relaxed trials, responsive "
                          "trials by three rules, per-CFU responsiveness, rest-only activity, the behaviour check, Loco2Fluor / LocoPupil2Fluor "
                          "GLMs); (d) pupil measures; (e) the comparison of each condition with its baseline — raw, at matched locomotion, "
                          "at matched pupil, behaviour-adjusted, where the activity is (common region, tile maps, fibre segments), at matched "
                          "speed, and the time course on the day's clock. Name the days' conditions in <mouse>\\day_labels.csv (day,label). "
                          "summary.txt, experiments_summary.csv, run_summary.csv, condition_summary.csv, loco_*.csv, state_rates.csv, "
                          "pupil_measures.csv, acute_chronic.csv, cfu_updown_*.csv and the figures feed tab 8.",
                  wraplength=1000, justify="left").pack(anchor="w", padx=8, pady=6)
        o = ttk.Frame(f); o.pack(anchor="w", padx=8)
        ttk.Label(o, text="gate").pack(side="left")
        ttk.Combobox(o, textvariable=self.v["a_gate"], values=("none", "still"), width=6, state="readonly").pack(side="left", padx=4)
        ttk.Label(o, text="min -log10 p").pack(side="left", padx=(8, 2))
        ttk.Entry(o, textvariable=self.v["a_min_p"], width=5).pack(side="left")
        ttk.Label(o, text="baseline condition (empty = pre / earliest day)").pack(side="left", padx=(8, 2))
        ttk.Entry(o, textvariable=self.v["a_baseline"], width=10).pack(side="left")
        ttk.Label(o, text="events").pack(side="left", padx=(8, 2))
        ttk.Combobox(o, textvariable=self.v["a_source"], values=("drop-rejected-cfu", "kept-cfu-only", "all"),
                     width=18, state="readonly").pack(side="left")
        o2 = ttk.Frame(f); o2.pack(anchor="w", padx=8, pady=4)
        ttk.Checkbutton(o2, text="keep events in excluded stretches too", variable=self.v["a_keep_flagged"]).pack(side="left", padx=4)
        ttk.Label(o2, text="margin").pack(side="left")
        ttk.Entry(o2, textvariable=self.v["a_margin"], width=3).pack(side="left", padx=(2, 8))
        ttk.Checkbutton(o2, text="ignore curation", variable=self.v["a_no_curation"]).pack(side="left", padx=4)
        ttk.Checkbutton(o2, text="skip CFU matching", variable=self.v["a_no_cfu"]).pack(side="left", padx=4)
        ttk.Checkbutton(o2, text="pupil vs activity", variable=self.v["a_pupil"]).pack(side="left", padx=4)
        ttk.Checkbutton(o2, text="maps (needs _AQuA2_events.csv)", variable=self.v["a_maps"]).pack(side="left", padx=4)
        ttk.Label(o2, text="tile um").pack(side="left")
        ttk.Entry(o2, textvariable=self.v["a_tile"], width=4).pack(side="left", padx=(2, 8))
        o3 = ttk.Frame(f); o3.pack(anchor="w", padx=8, pady=2)
        ttk.Label(o3, text="v32-v34 readouts:", foreground="#555555").pack(side="left")
        ttk.Checkbutton(o3, text="locomotion trials + the three triggers (running / pupil / rest, v34)", variable=self.v["a_loco"]).pack(side="left", padx=4)
        ttk.Checkbutton(o3, text="GLMs (the slowest part)", variable=self.v["a_glm"]).pack(side="left", padx=4)
        ttk.Checkbutton(o3, text="CFU x time heatmaps", variable=self.v["a_heat"]).pack(side="left", padx=4)
        ttk.Checkbutton(o3, text="heatmap rows from the movies", variable=self.v["a_movies"]).pack(side="left", padx=4)
        ttk.Checkbutton(o3, text="motion / registration check (v33)", variable=self.v["a_motion"]).pack(side="left", padx=4)
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=6)
        ttk.Button(b, text="Run analysis", command=self.run_analyze).pack(side="left")
        ttk.Button(b, text="Open aqua_analysis folder", command=lambda: self.open_mouse_sub("aqua_analysis")).pack(side="left", padx=8)
        ttk.Button(b, text="Show summary.txt", command=self.show_summary).pack(side="left", padx=8)
        self.summary_box = tk.Text(f, height=12, state="disabled", font=("Consolas", 9), wrap="none")
        self.summary_box.pack(fill="both", expand=True, padx=8, pady=6)

    def _tab_report(self):
        f = ttk.Frame(self.nb); self.nb.add(f, text="8  Report")
        ttk.Label(f, text="Two self-contained HTML files per animal in aqua_analysis\\ (open in a browser, click a figure to zoom, print to PDF). "
                          "THE REPORT (<mouse>_report.html, seconds): the general representation first — the sessions on the day's clock, activity "
                          "by running / resting state with event overlays, ECDFs, the acute / chronic change of every measure — then the CFU x time "
                          "heatmaps with up / down / unchanged CFUs, the locomotion-evoked activity (trials, responsive rules, responsive CFUs, "
                          "rest-only, behaviour check), pupil, coupling GLMs, the behaviour-controlled change with a verdict per experiment, the "
                          "time course, per-experiment details, settings and methods. THE QC REPORT (<mouse>_qc.html, minutes; reads the "
                          "movies): a QC table with warnings, and per recording (collapsed) the detection map with CFU outlines coloured by what "
                          "they look like, the raster with the locomotion trials, clips, a filterable gallery of EVERY detected event's montage and "
                          "the CFU table — the evidence that the detection is right. Needs tab 7.",
                  wraplength=1000, justify="left").pack(anchor="w", padx=8, pady=6)
        o = ttk.Frame(f); o.pack(anchor="w", padx=8)
        ttk.Label(o, text="QC: montages per recording (0 = every event)").pack(side="left")
        ttk.Entry(o, textvariable=self.v["r_max_events"], width=5).pack(side="left", padx=(4, 10))
        ttk.Label(o, text="clips per recording").pack(side="left")
        ttk.Entry(o, textvariable=self.v["r_clips"], width=4).pack(side="left", padx=(4, 10))
        ttk.Checkbutton(o, text="no clips", variable=self.v["r_no_clips"]).pack(side="left", padx=4)
        ttk.Checkbutton(o, text="open when done", variable=self.v["r_open"]).pack(side="left", padx=4)
        b = ttk.Frame(f); b.pack(anchor="w", padx=8, pady=8)
        ttk.Button(b, text="Build the report", command=self.run_report).pack(side="left")
        ttk.Button(b, text="Build the QC report", command=lambda: self.run_report(qc=True)).pack(side="left", padx=8)
        ttk.Button(b, text="Open the report", command=self.open_report).pack(side="left", padx=(16, 0))
        ttk.Button(b, text="Open the QC report", command=lambda: self.open_report(qc=True)).pack(side="left", padx=8)
        ttk.Button(b, text="Open aqua_analysis folder", command=lambda: self.open_mouse_sub("aqua_analysis")).pack(side="left", padx=8)
        self.report_lbl = ttk.Label(f, text="", justify="left", wraplength=1000, foreground="#333333")
        self.report_lbl.pack(anchor="w", padx=8, pady=4)
        x = ttk.LabelFrame(f, text="Across animals (v35): the stimulation protocols side by side")
        x.pack(fill="x", padx=8, pady=(10, 4))
        ttk.Label(x, text="Reads every RearChr animal's tab-7 tables under the parent folder (the folder that holds the mouse folders), "
                          "labels each day x FOV experiment with the protocol of its day (day_labels.csv: sham / 2 Hz / 5 Hz / 5 Hz + 40 ms ...) "
                          "and compares the pre -> post change per protocol: event rate (all / still / running frames), amplitude, CFU fate, "
                          "with the behavioural state next to it. No pupil step. Writes <parent>\\cross_animal\\conditions.html (+ png, csv, txt) "
                          "and lists the days that are not analysed yet.", wraplength=1000, justify="left").pack(anchor="w", padx=6, pady=(4, 2))
        xr = ttk.Frame(x); xr.pack(fill="x", padx=6, pady=(0, 6))
        ttk.Label(xr, text="Parent folder").pack(side="left")
        ttk.Entry(xr, textvariable=self.v["x_parent"], width=46).pack(side="left", padx=6)
        ttk.Button(xr, text="Browse…", command=self.pick_parent).pack(side="left")
        ttk.Button(xr, text="Compare the protocols", command=self.run_conditions).pack(side="left", padx=(16, 4))
        ttk.Button(xr, text="Open conditions.html", command=self.open_conditions).pack(side="left", padx=4)

    # ------------------------------------------------------------------ helpers
    def mouse(self) -> Path | None:
        m = self.v["mouse"].get().strip().strip('"')
        if not m:
            messagebox.showwarning("Mouse folder", "Choose the mouse folder first (the folder that holds the "
                                                   "<mouse>_<date>_<run> folders).")
            return None
        p = Path(m)
        if not p.is_dir():
            messagebox.showerror("Mouse folder", f"Not a folder:\n{p}")
            return None
        return p

    def pick_mouse(self):
        d = filedialog.askdirectory(title="Mouse folder (contains the <mouse>_<date>_<run> folders)")
        if d:
            self.v["mouse"].set(d)
            self.save_settings()
            self.refresh_status()
            self.refresh_state()

    def pick_opts(self):
        p = filedialog.askopenfilename(title="AQuA2 settings file used for FOVs without their own",
                                       filetypes=[("AQuA2 settings", "*_Opt.csv"), ("csv", "*.csv"), ("all", "*.*")])
        if p:
            self.v["b_opts"].set(p)
            self.save_settings()

    def pick_cfu_aqua(self):
        p = filedialog.askdirectory(title="AQuA2 copy for the CFU step (the folder with CFU_cmd_batch.m)")
        if p:
            self.v["b_cfu_aqua"].set(p)
            self.save_settings()

    def pick_imagej(self):
        p = filedialog.askopenfilename(title="ImageJ / Fiji executable (ImageJ.exe or ImageJ-win64.exe)",
                                       filetypes=[("executable", "*.exe"), ("all", "*.*")])
        if p:
            self.v["imagej"].set(p)
            self.save_settings()

    def pick_matlab(self):
        p = filedialog.askopenfilename(title="MATLAB executable (…\\MATLAB\\R20xx\\bin\\matlab.exe)",
                                       filetypes=[("executable", "*.exe"), ("all", "*.*")])
        if p:
            self.v["matlab"].set(p)
            self.save_settings()

    def matlab_exe(self, ask: bool = True) -> str | None:
        """The MATLAB executable from the top row, or None (with a dialog) when it is not there."""
        exe = self.v["matlab"].get().strip().strip('"')
        if exe and Path(exe).exists():
            return exe
        if ask:
            messagebox.showerror("MATLAB", "No MATLAB executable set (top row: MATLAB → Browse…, usually "
                                           "C:\\Program Files\\MATLAB\\R20xx\\bin\\matlab.exe).\n\nWithout it, copy the commands "
                                           "(tab 5 / tab 3) and run them in a MATLAB window.")
        return None

    def save_settings(self):
        try:
            SETTINGS.write_text(json.dumps({k: self.v[k].get() for k in ("mouse", "aqua", "imagej", "matlab", "um", "border", "sec",
                                                                          "t_as_root", "t_workers", "t_movies", "t_common", "b_opts", "auto_matlab",
                                                                          "x_parent", "t_day_rule", "b_cfu_min", "b_cfu_overlap", "b_cfu_merge", "b_cfu_aqua")},
                                           indent=2))
        except Exception:
            pass

    def on_close(self):
        self.save_settings()
        self.destroy()

    def logln(self, s: str):
        self.q.put(str(s) + "\n")
        self._note_progress(str(s))
        # the Stop button (v36): every step logs as it goes, so a stop request takes effect at the step's next log
        # line, inside whatever loop it is in — without the steps having to know about the button. Only from a
        # worker thread: the GUI thread logs too and must never be interrupted.
        if self.stop_evt.is_set() and threading.current_thread() is not threading.main_thread():
            raise Stopped("by the Stop button")

    # ------------------------------------------------------------------ progress (v35)
    PROG_RE = re.compile(r"\[(\d+)\s*/\s*(\d+)\]")
    STEP_CHIP = (("tune", 2), ("unpack", 0), ("behaviour", 0), ("clean", 1), ("check stacks", 1), ("restore", 1), ("drop", 1),
                 ("cfu list", 5), ("curat", 5), ("matlab", 4), ("aqua2", 4), ("cfu", 4), ("batch", 3), ("analy", 6),
                 ("report", 7), ("protocols", 7))

    def _note_progress(self, line: str):
        """A step's log line: '[i/n]' anywhere in it = item i of n is starting (the rest of the line names it)."""
        m = self.PROG_RE.search(line)
        if not m:
            return
        i, n = int(m.group(1)), int(m.group(2))
        if n <= 0:
            return
        item = " ".join((line[:m.start()] + line[m.end():]).split()).strip(" =:-")
        self.prog = {"i": i, "n": n, "item": item[:90], "sub": ""}

    def set_progress(self, done: int, total: int, text: str = ""):
        """From a worker: the count of finished items when the step cannot print [i/n] lines itself (MATLAB results)."""
        self.prog = {"i": min(done + 1, total) if total else 0, "n": total, "item": text[:90], "sub": "", "done": done}

    def _chip_for(self, name: str):
        low = name.lower()
        for key, idx in self.STEP_CHIP:
            if key in low:
                return idx
        return None

    def _highlight_chip(self, idx, on: bool):
        if idx is None or idx >= len(self.state_chips):
            return
        chip = self.state_chips[idx]
        if on:
            chip._saved_bg = chip.cget("bg")
            chip.configure(bg="#ffb347", relief="sunken")
        else:
            chip.configure(bg=getattr(chip, "_saved_bg", STATE_BG["na"]), relief="groove")

    def _tick_progress(self):
        if not self.busy:
            return
        el = time.time() - self.run_t0
        i, n, item = self.prog.get("i", 0), self.prog.get("n", 0), self.prog.get("item", "")
        done = self.prog.get("done", max(i - 1, 0))
        txt = f"running: {self.run_name} — {el / 60:.0f}:{el % 60:02.0f} elapsed"
        if n > 0:
            frac = min(done / n, 1.0)
            self.prog_bar.configure(mode="determinate"); self.prog_bar["value"] = 100 * frac
            left = (el / done) * (n - done) if done > 0 else None
            txt += f" — {'result' if 'done' in self.prog else 'item'} {min(i, n)} of {n}"
            if item:
                txt += f": {item}"
            if left is not None and done < n:
                txt += f" — about {left:.0f} s left" if left < 90 else f" — about {left / 60:.0f} min left"
        else:
            if str(self.prog_bar.cget("mode")) != "indeterminate":
                self.prog_bar.configure(mode="indeterminate"); self.prog_bar.start(80)
            txt += " (no item count for this step — see the log)"
        self.prog_lbl.configure(text=txt, fg="#b30000")
        self.after(1000, self._tick_progress)

    def _place_divider(self):
        """First layout (v36): the log keeps about seven lines, the tabs get the rest; the divider stays draggable."""
        try:
            h = self.pane.winfo_height()
            if h > 300:
                self.pane.sashpos(0, max(h - 140, 200))
        except Exception:
            pass

    def clear_log(self):
        self.log.configure(state="normal"); self.log.delete("1.0", "end"); self.log.configure(state="disabled")

    def poll(self):
        try:
            while True:
                s = self.q.get_nowait()
                if isinstance(s, tuple) and s and s[0] == "call":
                    try:
                        s[1]()
                    except Exception as e:
                        s = f"(gui) {e}\n"
                    else:
                        continue
                if s == "\x00DONE":
                    self.busy = False
                    self.state_lbl.configure(text="idle", foreground="#1a5e1a")
                    el = time.time() - self.run_t0 if self.run_t0 else 0.0
                    try:
                        self.prog_bar.stop()
                    except Exception:
                        pass
                    self.prog_bar.configure(mode="determinate"); self.prog_bar["value"] = 100 if self.prog.get("n") else 0
                    self.prog_lbl.configure(text=f"idle — last: {self.run_name} finished in {el / 60:.0f}:{el % 60:02.0f} "
                                                 f"(see the log for errors)", fg="#1a5e1a")
                    self._highlight_chip(self.run_chip, False); self.run_chip = None
                    self.refresh_status()
                    self.refresh_state(log_it=True)
                    continue
                self.log.configure(state="normal"); self.log.insert("end", s); self.log.see("end")
                self.log.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(150, self.poll)

    def start(self, name: str, fn, *args):
        if self.busy:
            messagebox.showinfo("Busy", "A step is still running — wait for it to finish (see the log) or press Stop.")
            return
        self.busy = True
        self.stop_evt.clear()
        self.state_lbl.configure(text=f"running: {name} …", foreground="#b30000")
        self.run_name, self.run_t0 = name, time.time()
        self.prog = {"i": 0, "n": 0, "item": "", "sub": ""}
        self.prog_bar.configure(mode="indeterminate"); self.prog_bar.start(80)
        self.prog_lbl.configure(text=f"running: {name} — starting …", fg="#b30000")
        self.run_chip = self._chip_for(name); self._highlight_chip(self.run_chip, True)
        self.after(1000, self._tick_progress)
        self.logln(f"\n===== {name} =====")

        def worker():
            # (plain q.put below, not logln: after a Stop every logln from this thread raises Stopped again)
            try:
                fn(*args)
            except Stopped as e:
                self.q.put(f"stopped: {e}\n")
            except SystemExit as e:                      # the command-line tools raise SystemExit with a message
                self.q.put(f"ERROR: {e}\n")
            except Exception:
                self.q.put("ERROR:\n" + traceback.format_exc() + "\n")
            finally:
                self.q.put("done.\n")
                self.q.put("\x00DONE")
        threading.Thread(target=worker, daemon=True).start()

    def set_state(self, text: str):
        """Change the running-state label from a worker thread (the chain: which of its steps runs now)."""
        def apply():
            self.state_lbl.configure(text=f"running: {text} …", foreground="#b30000")
            self.run_name = text
            self.prog = {"i": 0, "n": 0, "item": "", "sub": ""}
            self._highlight_chip(self.run_chip, False)
            self.run_chip = self._chip_for(text); self._highlight_chip(self.run_chip, True)
        self.q.put(("call", apply))

    def stop_run(self):
        """Ask the running MATLAB process / auto chain to stop (between steps; MATLAB is killed)."""
        if not self.busy:
            self.logln("nothing is running")
            return
        self.stop_evt.set(); self.manual_evt.set()
        self.q.put("stop requested — the step stops at its next log line (a recording being unpacked or analysed finishes first; "
                   "a MATLAB started from here is closed; a pasted MATLAB script keeps running in its own window)\n")
        self.prog_lbl.configure(text=f"stopping: {self.run_name} …", fg="#b36b00")

    def check_stop(self):
        if self.stop_evt.is_set():
            raise Stopped("by the Stop button")

    def fnum(self, key, default):
        s = self.v[key].get().strip()
        return type(default)(s) if s else default

    # ------------------------------------------------------------------ tabs 1-2: unpack, clean
    def run_unpack(self, folders=None):
        m = self.mouse()
        if m is None:
            return
        import unpack_bin
        kw = dict(seconds=self.fnum("sec", 2.0), overwrite=bool(self.v["u_overwrite"].get()),
                  behavior=bool(self.v["u_behavior"].get()), global_fix=bool(self.v["u_gfix"].get()),
                  keep_raw=bool(self.v["u_keep_raw"].get()), register=bool(self.v["u_register"].get()), log=self.logln)
        self.start("unpack", lambda: unpack_bin.process_folders(folders or [str(m)], **kw))

    def run_unpack_one(self):
        d = filedialog.askdirectory(title="One run folder (contains the .sbx)")
        if d:
            self.run_unpack(folders=[d])

    def run_behavior_only(self):
        m = self.mouse()
        if m is None:
            return
        import unpack_bin
        self.start("behaviour csv", lambda: unpack_bin.rebuild_behavior([str(m)], log=self.logln))

    # ------------------------------------------------------------------ step 1b
    def clean_targets(self, selected: bool):
        if selected:
            sel = self.tree.selection()
            if not sel:
                messagebox.showinfo("Select", "Select one or more recordings in the table first.")
                return None
            return [self.status_rows[int(i)]["tif"] for i in sel]
        m = self.mouse()
        return [m] if m else None

    def clean_kwargs(self):
        return dict(win_s=self.fnum("c_win", 60.0), k=4.0, min_rel=self.fnum("c_min_rel", 0.03),
                    shape_rel=self.fnum("c_shape_rel", 0.05), shape_win_s=self.fnum("c_shape_win", 300.0),
                    min_clean_s=self.fnum("c_min_clean", 10.0), max_spike_s=self.fnum("c_max_spike", 10.0),
                    core_ratio=self.fnum("c_core_ratio", 3.0), core_min_pct=self.fnum("c_core_pct", 10.0),
                    stacks=bool(self.v["c_stacks"].get()), stack_context=self.fnum("c_context", 6))

    def run_clean(self, selected: bool = False, mode: str = "keep", force: bool | None = None):
        targets = self.clean_targets(selected)
        if not targets:
            return
        import aqua_clean
        from afferent2p.globalfix import apply
        kw = self.clean_kwargs()
        force = bool(self.v["c_force"].get()) if force is None else force
        report = bool(self.v["c_report"].get()) and mode == "keep"

        def job():
            tifs = aqua_clean.collect([str(t) for t in targets])
            if not tifs:
                self.logln("no binned tiffs found (run Unpack first)")
                return
            self.logln(f"{len(tifs)} tiff(s)")
            n = 0
            for t in tifs:
                rec = apply(t, mode=mode, force=force, report_only=report, log=self.logln, **kw)
                if rec and not report:
                    n += 1
                    if mode != "keep":
                        from afferent2p.aqua_io import aqua_files
                        if any(p.exists() for k, p in aqua_files(t).items() if k in ("res", "features", "cfu")):
                            self.logln("    note: AQuA2 results exist for this tiff and were computed on the PREVIOUS "
                                       "movie — rebuild the batch with overwrite and re-run MATLAB for it")
            if not report:
                self.logln(f"done: {n} tiff(s)")
        name = {"keep": "clean & check" + (" (report only)" if report else ""),
                "drop": "drop excluded frames", "fill": "fill excluded frames"}[mode]
        self.start(name, job)

    def run_drop(self):
        rows = self.selected_rows()
        if not rows:
            return
        if not messagebox.askyesno("Drop frames",
                                   f"Leave the EXCLUDED frames out of {len(rows)} movie(s)?\n\n"
                                   "This changes the tiff on disk: AQuA2 must be re-run on it, and the frames come back "
                                   "only from _raw.tif or by re-unpacking. Keeping the frames (the default) already keeps "
                                   "their events out of the statistics."):
            return
        self.run_clean(selected=True, mode="drop", force=True)

    def run_restore(self):
        rows = self.selected_rows()
        if not rows:
            return
        missing = [r["run"] for r in rows if not r["raw"].exists()]
        if missing:
            messagebox.showinfo("No raw copy", "No _raw.tif for:\n" + "\n".join(missing)
                                + "\n\nRe-unpack those runs instead (tab 1, overwrite).")
            rows = [r for r in rows if r["raw"].exists()]
            if not rows:
                return
        if not messagebox.askyesno("Restore", f"Replace {len(rows)} movie(s) by their _raw.tif and run the brightness "
                                              "step again (frames kept)?"):
            return
        from afferent2p.globalfix import restore_raw, apply
        kw = self.clean_kwargs()
        tifs = [r["tif"] for r in rows]

        def job():
            for t in tifs:
                if restore_raw(t, log=self.logln):
                    apply(t, mode="keep", force=True, log=self.logln, **kw)
        self.start("restore from raw", job)

    # ------------------------------------------------------------------ pipeline state strip
    def refresh_state(self, log_it: bool = False):
        """Recompute where this animal is (afferent2p.pipeline_state) in a background thread and show it in the strip."""
        m = self.v["mouse"].get().strip().strip('"')
        if not m or not Path(m).is_dir() or self._state_busy:
            return
        self._state_busy = True

        def job():
            try:
                from afferent2p.pipeline_state import mouse_state
                st = mouse_state(Path(m), sec=self.fnum("sec", 2.0))
            except Exception as e:
                st = {"error": str(e)}
            self.q.put(("call", lambda: self._apply_state(st, log_it)))
        threading.Thread(target=job, daemon=True).start()

    def _apply_state(self, st: dict, log_it: bool = False):
        self._state_busy = False
        if "error" in st:
            self.next_lbl.configure(text=f"state: {st['error']}")
            return
        self.state = st
        nxt_tab = st["next"]["tab"]
        for lbl, step in zip(self.state_chips, st["steps"]):
            lbl.configure(text=f"{step['n']} {step['name'].split(' ')[0]}  {STATE_MARK.get(step['status'], '?')}",
                          bg=STATE_BG.get(step["status"], "#f3f3f3"),
                          relief="solid" if step["tab"] == nxt_tab else "groove", bd=2 if step["tab"] == nxt_tab else 1)
        self.next_lbl.configure(text="next: " + st["next"]["text"])
        if log_it:
            from afferent2p.pipeline_state import format_state
            self.logln(format_state(st))

    def state_chip_click(self, i: int):
        """A chip: go to its tab and show that step's state under the strip."""
        self.nb.select(i)
        if self.state and i < len(self.state["steps"]):
            s = self.state["steps"][i]
            self.next_lbl.configure(text=f"{s['n']} {s['name']}: {s['text']}" + (f" — {s['detail']}" if s["detail"] else "")
                                    + f"      (next: {self.state['next']['text'][:90]}{'…' if len(self.state['next']['text']) > 90 else ''})")

    def refresh_status(self):
        m = self.v["mouse"].get().strip().strip('"')
        sel_runs = {self.status_rows[int(i)]["run"] for i in self.tree.selection()} if self.status_rows else set()
        for i in self.tree.get_children():
            self.tree.delete(i)
        self.status_rows = []
        if not m or not Path(m).is_dir():
            return
        try:
            self.status_rows = run_status(Path(m))
        except Exception as e:
            self.logln(f"table: {e}")
            return
        for i, r in enumerate(self.status_rows):
            tags = ()
            if "CHECK" in str(r["check"]):
                tags = ("check",)
            elif str(r["stale"]).startswith("STALE"):
                tags = ("stale",)
            elif "old version" in str(r["mode"]):
                tags = ("old",)
            self.tree.insert("", "end", iid=str(i), values=(r["run"], r["frames"], r["artefact"], r["state"], r["excluded"],
                                                            r["movie"], r["mode"], r["corrected"], r["check"], r["stale"]),
                             tags=tags)
            if r["run"] in sel_runs:
                self.tree.selection_add(str(i))
        self.fill_stretches()

    def refresh_status_log(self):
        self.refresh_status()
        if not self.status_rows:
            self.logln("no binned tiffs under the mouse folder")
        for r in self.status_rows:
            self.logln(f"{r['run']}: AQuA2 result {r['stale']}; brightness step: {r['mode']}, "
                       f"artefact {r['artefact'] or 0}, state {r['state'] or 0}, excluded {r['excluded'] or 0}")

    def selected_rows(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("Select", "Select one or more recordings in the table first.")
            return []
        return [self.status_rows[int(i)] for i in sel]

    def selected_row(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("Select", "Select a recording in the table first.")
            return None
        return self.status_rows[int(sel[0])]

    def fill_stretches(self):
        for i in self.stree.get_children():
            self.stree.delete(i)
        sel = self.tree.selection()
        if not sel or not self.status_rows:
            return
        r = self.status_rows[int(sel[0])]
        from afferent2p.globalfix import describe_decision, frames_str, normalise_decision
        self.stretch_lbl.configure(text=f"flagged stretches of {r['run']} — select and decide (EXCLUDED frames leave the "
                                        f"statistics; frames stay in the movie). Flagged frames = differ clearly from their "
                                        f"neighbours; 'exclude' = those only:")
        for st in r["stretches"]:
            dec = st.get("decision", st.get("default", ""))
            old = "core_frames" not in st
            check = old and str(dec).startswith("ex")
            own = normalise_decision(dec) != normalise_decision(st.get("default", ""))
            tags = (st.get("class", ""),) + (("check",) if check else ())
            self.stree.insert("", "end", iid=st["id"],
                              values=(st["id"], st.get("class", ""), st.get("n_frames", ""), st.get("seconds", ""),
                                      st.get("level_pct_max", ""), st.get("shape_pct_median", ""),
                                      frames_str(st.get("core_frames") or [], max_runs=5) if not old else "(old version)",
                                      ("CHECK (old version): " if check else "") + str(st.get("verdict", "")),
                                      describe_decision(st, dec) + (" (yours)" if own else ""),
                                      frames_str(st.get("excluded_frames") or [], max_runs=5) if "excluded_frames" in st
                                      else ("whole stretch" if str(dec).startswith("ex") else "")),
                              tags=tags)

    def decide(self, decision):
        r = self.selected_row()
        if r is None:
            return
        ids = list(self.stree.selection())
        if not ids:
            messagebox.showinfo("Select", "Select one or more stretches in the lower table first.")
            return
        from afferent2p.globalfix import apply_decisions, read_decisions, write_decisions
        if decision is None:                     # back to default: drop the entries
            dec = read_decisions(r["tif"])
            for i in ids:
                dec.pop(i, None)
            write_decisions(r["tif"], dec)
            apply_decisions(r["tif"], None, log=self.logln)
        else:
            if decision == "exclude":
                empty = [st["id"] for st in r["stretches"] if st["id"] in ids and st.get("core_frames") == []]
                if empty:
                    self.logln(f"note: no flagged frames in {', '.join(empty)} — 'exclude' excludes nothing there; use "
                               f"'Exclude whole stretch' or type the frames")
            apply_decisions(r["tif"], {i: decision for i in ids}, log=self.logln)
        self.refresh_status()
        if bool(self.v["c_stacks"].get()) and r["stack_frames"].exists():
            self.logln("(the check stacks still show the previous decision in their labels — 'Make / refresh check "
                       "stacks' rewrites them)")

    def decide_frames(self):
        txt = self.v["c_frames"].get().strip()
        if not txt:
            messagebox.showinfo("Frames", "Type the recorded frame numbers to exclude, e.g. 68,76-77 (as in the stack labels).")
            return
        from afferent2p.globalfix import parse_frame_ranges
        try:
            frames = parse_frame_ranges(txt)
        except ValueError:
            messagebox.showerror("Frames", f"Could not read '{txt}': use numbers and ranges like 68,76-77")
            return
        r = self.selected_row()
        if r is None:
            return
        ids = list(self.stree.selection())
        if len(ids) != 1:
            messagebox.showinfo("Select", "Select exactly one stretch in the lower table for a frame list.")
            return
        st = next((s for s in r["stretches"] if s["id"] == ids[0]), None)
        if st is not None:
            outside = [f for f in frames if not (int(st["first"]) <= f <= int(st["last"]))]
            if outside:
                messagebox.showinfo("Frames", f"Frames {outside} are outside stretch {st['id']} and are ignored; frames "
                                              f"between stretches are never excluded.")
        self.decide("frames:" + txt)

    def open_stack(self, which: str):
        r = self.selected_row()
        if r is None:
            return
        p = r["stack_frames"] if which == "frames" else r["stack_movie"]
        if not p.exists():
            if not r["stretches"]:
                messagebox.showinfo("No stack", "Nothing is flagged in this recording — there is no check stack.")
                return
            if messagebox.askyesno("No stack yet", f"{p.name} does not exist yet (older result). Make the check stacks "
                                                   f"for this recording now (reads the movie once, ~10-30 s)?"):
                self.make_stacks(then_open=which)
            return
        open_in_imagej(p, self.v["imagej"].get(), log=self.logln)

    def make_stacks(self, then_open: str | None = None):
        rows = self.selected_rows()
        if not rows:
            return
        from afferent2p.globalfix import make_check_stacks
        ctx = self.fnum("c_context", 6)
        tifs = [r["tif"] for r in rows]

        def job():
            for t in tifs:
                make_check_stacks(t, context=ctx, log=self.logln)
            if then_open and len(rows) == 1:
                p = rows[0]["stack_frames"] if then_open == "frames" else rows[0]["stack_movie"]
                if p.exists():
                    self.q.put(("call", lambda: open_in_imagej(p, self.v["imagej"].get(), log=self.logln)))
        self.start("check stacks", job)

    def open_selected(self, what: str):
        r = self.selected_row()
        if r is None:
            return
        p = {"png": r["png"], "csv": r["csv"], "folder": r["tif"].parent}[what]
        if not Path(p).exists():
            messagebox.showinfo("Missing", f"{Path(p).name} does not exist yet — run the brightness step (or report only) first.")
            return
        open_path(p)

    def open_sheet(self):
        r = self.selected_row()
        if r is None:
            return
        if not r["sheet"].exists():
            messagebox.showinfo("No sheet", "No check sheet for this recording: nothing was flagged, or it was corrected "
                                            "by an older version (run with Force to make one).")
            return
        SheetWindow(self, r["sheet"], r["csv"], qc=r["png"])

    # ------------------------------------------------------------------ step 2
    def copy_aqua_gui_cmd(self):
        cmd = f"cd('{self.v['aqua'].get().strip()}'); aqua_gui"
        self.clipboard_clear(); self.clipboard_append(cmd)
        self.logln(f"copied: {cmd}")

    def open_mouse_sub(self, sub: str):
        m = self.mouse()
        if m is None:
            return
        p = m / sub if sub else m
        if not p.exists():
            messagebox.showinfo("Missing", f"{p} does not exist yet.")
            return
        open_path(p)

    def list_opts(self):
        m = self.mouse()
        if m is None:
            return
        from afferent2p.plan import PARK_PREFIX
        found = sorted(p for p in set(list(m.glob("*_Opt.csv")) + list((m / "aqua_opts").glob("*_Opt.csv")) + list(m.rglob("*_Opt.csv")))
                       if not any(x.startswith(PARK_PREFIX) for x in p.relative_to(m).parts))
        lines = [f"{len(found)} settings file(s):"] + [f"  {p.relative_to(m)}" for p in found]
        if not found:
            lines.append("  none yet — SaveOpts from the AQuA2 GUI into the mouse folder as FOV1_Opt.csv, FOV2_Opt.csv …")
        self.tune_box.configure(state="normal"); self.tune_box.delete("1.0", "end")
        self.tune_box.insert("end", "\n".join(lines)); self.tune_box.configure(state="disabled")
        self.logln("\n".join(lines))

    # ------------------------------------------------------------------ tab 3: detection settings
    def _tune_grid(self):
        grid = []
        thr = self.v["t_thr"].get().strip(); ms = self.v["t_minsize"].get().strip(); ex = self.v["t_extra"].get().strip()
        if thr:
            grid.append("thrARScl=" + thr.replace(" ", ""))
        if ms:
            grid.append("minSize=" + ms.replace(" ", ""))
        for part in ex.replace(" ", "").split(";"):
            if "=" in part:
                grid.append(part)
        return grid or None

    def run_tune_plan(self):
        m = self.mouse()
        if m is None:
            return
        from afferent2p import tune as T
        try:
            grid = T.parse_grid(self._tune_grid())
        except ValueError as e:
            messagebox.showerror("grid", str(e)); return
        self.save_settings()
        self.start("tune: plan", lambda: self._tune_plan_job(m, grid))

    def _set_grid_fields(self, grid: dict):
        from afferent2p import tune as T
        if "thrARScl" in grid:
            self.v["t_thr"].set(", ".join(T._fmt(v) for v in grid["thrARScl"]))
        if "minSize" in grid:
            self.v["t_minsize"].set(", ".join(T._fmt(v) for v in grid["minSize"]))
        extra = [f"{k}={','.join(T._fmt(v) for v in vals)}" for k, vals in grid.items() if k not in ("thrARScl", "minSize")]
        if extra:
            self.v["t_extra"].set("; ".join(extra))

    def extend_grid(self):
        """Add the values the last report suggested (aqua_tune\\grid_suggestion.json) to the grid fields and plan again (v37):
        finished jobs are kept, only the new combinations are swept."""
        m = self.mouse()
        if m is None:
            return
        import aqua_tune
        from afferent2p import tune as T
        sug = aqua_tune.read_grid_suggestion(m)
        if not sug or not sug.get("suggest"):
            messagebox.showinfo("Nothing to extend", "The last report suggested no extension (aqua_tune\\grid_suggestion.json is missing or "
                                                     "empty): the grid brackets the optimum of every movie, or 3 Report has not run yet.")
            return
        try:
            grid = T.parse_grid(self._tune_grid())
        except ValueError as e:
            messagebox.showerror("grid", str(e)); return
        merged = aqua_tune.merge_grid(grid, sug["suggest"])
        added = {k: [v for v in merged[k] if v not in [float(x) for x in grid.get(k, [])]] for k in merged}
        added = {k: v for k, v in added.items() if v}
        if not added:
            messagebox.showinfo("Nothing to extend", "The suggested values are already in the grid fields — press 1 Plan.")
            return
        txt = "; ".join(f"{k}: + {', '.join(T._fmt(v) for v in vals)}" for k, vals in added.items())
        why = ", ".join(f"{aqua_tune.short_run(r)} {d.get('status', '')}" for r, d in sug.get("movies", {}).items() if d.get("suggest"))
        if not messagebox.askyesno("Extend the grid", f"The last report found the grid too narrow for: {why}.\n\nAdd {txt} to the grid and plan "
                                                      f"again? Finished jobs are kept; only the new combinations are swept (then 2 Run the sweep, 3 Report)."):
            return
        self._set_grid_fields(merged)
        self.logln(f"grid extended as suggested ({txt}) — planning again")
        self.run_tune_plan()

    def open_picker(self):
        """The by-eye picker (v38): the few settings worth a look, each drawn on the real movie, with a frame slider; a click
        writes the pick as the experiment's setting."""
        m = self.mouse()
        if m is None:
            return
        if not (m / "aqua_tune" / "tune_report.csv").exists():
            messagebox.showinfo("No sweep report yet", "aqua_tune\\tune_report.csv does not exist — press '3 Report' after the MATLAB sweep first.")
            return
        try:
            from afferent2p.pick import PickWindow
            PickWindow(self, m, log=self.logln, max_false=self.fnum("t_max_false", 0.10), min_on_fibre=self.fnum("t_min_fibre", 0.80),
                       tol_recall=self.fnum("t_tol_recall", 0.05))
        except Exception as e:
            messagebox.showerror("Pick by eye", f"could not open the picker: {e}")

    def hand_result_reference(self):
        """A result saved by hand from the AQuA2 GUI -> the MATLAB command that exports its events, so that the next 3 Report
        uses it as the reference of that run (v37)."""
        m = self.mouse()
        if m is None:
            return
        path = filedialog.askopenfilename(title="The AQuA2 result saved by hand (<movie>_AQuA2.mat, usually in the subfolder named after the movie)",
                                          initialdir=str(m), filetypes=[("AQuA2 result", "*_AQuA2.mat"), ("all", "*.*")])
        if not path:
            return
        cmd = f"aqua_export_events('{path}')"
        self.clipboard_clear(); self.clipboard_append(cmd)
        self.logln(f"MATLAB (copied to the clipboard; AQuA2 and the afferent2p folder on the path): {cmd}\n"
                   "  -> writes <name>_events.csv next to the .mat; then 3 Report shows every setting against this hand result "
                   "(the run must be in the plan: the 'runs' field, then 1 Plan and the sweep, if it is not one of the tuned movies)")
        messagebox.showinfo("Hand result → reference", f"Copied to the clipboard — paste into MATLAB:\n\n{cmd}\n\nIt writes the events csv next to the "
                                                       ".mat. The next '3 Report' then scores every setting against this hand result (hand% / extra% columns, "
                                                       "white dashed circles in the pictures). If this run is not one of the tuned movies, put it in the "
                                                       "'runs' field, press 1 Plan and run the sweep for it first.")

    def _tune_plan_kwargs(self):
        days = [x.strip() for x in self.v["t_days"].get().split(",") if x.strip()] or None
        runs = [x.strip() for x in self.v["t_runs"].get().split(",") if x.strip()] or None
        as_root = self.v["t_as_root"].get().strip().strip('"') or None
        return dict(runs_filter=runs, days_filter=days, as_root=as_root,
                    per_condition=TUNE_MOVIES_ARG.get(self.v["t_movies"].get(), "auto"),
                    n_workers=int(self.fnum("t_workers", 0)), n_inject=int(self.fnum("t_ninject", 72)),
                    amplitudes=[float(x) for x in self.v["t_amps"].get().replace(" ", "").split(",") if x] or None,
                    sizes=[int(float(x)) for x in self.v["t_sizes"].get().replace(" ", "").split(",") if x] or None)

    def _tune_plan_job(self, m: Path, grid: dict):
        """aqua_tune.plan with the tab's settings; builds the batch (tab 4) first when the manifest is missing."""
        import aqua_tune
        if not (m / "aqua_batch" / "manifest.json").exists():
            self.logln("no aqua_batch\\manifest.json yet — building the batch (tab 4) first")
            self._batch_job(m, overwrite=False)
        kw = self._tune_plan_kwargs()
        aqua_tune.plan(m, self.v["aqua"].get().strip(), grid, kw["runs_filter"], False, self.fnum("um", 0.792), kw["n_inject"],
                       True, True, kw["as_root"], 0, days_filter=kw["days_filter"], n_workers=kw["n_workers"],
                       amplitudes=kw["amplitudes"], sizes=kw["sizes"], per_condition=kw["per_condition"], log=self.logln)

    def actual_file(self, rel: str) -> str:
        """The name a generated file was really written under (v40: a file locked by another program — Excel, a
        stuck reader — is written as <stem>_2<suffix> and the manifest says so): 'aqua_batch/run_aqua_and_cfu.m' ->
        'aqua_batch/run_aqua_and_cfu_2.m' when needed, else unchanged. Same for aqua_tune/run_aqua_tune.m."""
        m = self.v["mouse"].get().strip().strip('"')
        folder, _, name = rel.replace("\\", "/").rpartition("/")
        keys = {"batch.csv": "batch_csv", "parameters_for_batch.csv": "params_csv", "run_aqua_batch.m": "aqua_script",
                "run_cfu_batch.m": "cfu_script", "run_aqua_and_cfu.m": "both_script"}
        try:
            if folder == "aqua_batch" and name in keys:
                files = json.loads((Path(m) / "aqua_batch" / "manifest.json").read_text()).get("files") or {}
                return f"{folder}/{files.get(keys[name], name)}"
            if folder == "aqua_tune" and name == "run_aqua_tune.m":
                files = json.loads((Path(m) / "aqua_tune" / "tune_manifest.json").read_text()).get("files") or {}
                return f"{folder}/{files.get('script', name)}"
        except Exception:
            pass
        return rel

    def matlab_cmd_for(self, rel: str) -> str:
        """run('<mouse folder as the MATLAB PC sees it>\\<rel>') — rel like 'aqua_tune/run_aqua_tune.m'."""
        m = self.v["mouse"].get().strip().strip('"')
        root = self.v["t_as_root"].get().strip().strip('"') or m
        sep = "\\" if ("\\" in root or ":" in root[:3]) else "/"
        rel = self.actual_file(rel)
        return f"run('{root.rstrip(chr(92) + '/')}{sep}{rel.replace('/', sep)}')"

    def tune_matlab_cmd(self) -> str:
        return self.matlab_cmd_for("aqua_tune/run_aqua_tune.m")

    def copy_tune_matlab(self):
        m = self.mouse()
        if m is None:
            return
        if not (m / "aqua_tune" / "run_aqua_tune.m").exists():
            messagebox.showerror("No sweep script yet", "aqua_tune\\run_aqua_tune.m does not exist for this mouse — press "
                                                        "'1 Plan' first (it writes the movies and the script).")
            return
        cmd = self.tune_matlab_cmd()
        self.clipboard_clear(); self.clipboard_append(cmd)
        self.logln(f"copied: {cmd}")

    def _as_root_conflict(self) -> str | None:
        """'Run in MATLAB now' runs on THIS PC: the 'mouse folder as the MATLAB PC sees it' field must be empty or this path."""
        root = self.v["t_as_root"].get().strip().strip('"').rstrip("\\/")
        m = self.v["mouse"].get().strip().strip('"').rstrip("\\/")
        if root and os.path.normcase(os.path.normpath(root)) != os.path.normcase(os.path.normpath(m)):
            return (f"The tab-3 field 'mouse folder as the MATLAB PC sees it' is {root!r}, but MATLAB would run on this PC with "
                    f"{m!r}. Clear that field (same PC) or run the copied command on the other PC.")
        return None

    def run_tune_matlab(self):
        """The sweep from here: matlab -batch on run_aqua_tune.m (restartable; finished jobs are skipped)."""
        m = self.mouse()
        if m is None or self.matlab_exe() is None:
            return
        script = m / self.actual_file("aqua_tune/run_aqua_tune.m")
        if not script.exists():
            messagebox.showerror("No sweep script yet", "Press '1 Plan' first — it writes aqua_tune\\run_aqua_tune.m.")
            return
        c = self._as_root_conflict()
        if c:
            messagebox.showerror("Paths", c); return
        import aqua_tune
        try:
            pr = aqua_tune.sweep_progress(m)
        except SystemExit as e:
            messagebox.showerror("Sweep", str(e)); return
        if pr["pending"] == 0 and pr["total"]:
            messagebox.showinfo("Sweep", pr["head"] + "\n\nNothing pending — press '3 Report'.")
            return
        if not messagebox.askyesno("Run the sweep", f"{pr['head']}\n\nStart MATLAB now and run the pending jobs? This takes hours "
                                                   f"(about 1.4 min per job per worker); the GUI stays usable, Stop ends it, and a "
                                                   f"later start continues where it stopped."):
            return
        self.save_settings()
        self.start("tune: MATLAB sweep", lambda: self._matlab_job(script, "the tuning sweep"))

    def _matlab_job(self, script: Path, what: str) -> int:
        from afferent2p.matlab_run import run_script
        exe = self.v["matlab"].get().strip().strip('"')
        code = run_script(exe, script, log=self.logln, stop=self.stop_evt, name=what)
        if code == -1:
            raise Stopped("MATLAB was stopped")
        if code != 0:
            raise RuntimeError(f"MATLAB ended {what} with exit code {code} — see the lines above (and {script.stem}_log.txt)")
        return code

    def _report_job(self, m: Path, apply: bool, apply_all: bool):
        import aqua_tune
        res = aqua_tune.report(m, self.fnum("t_max_false", 0.10), self.fnum("t_min_fibre", 0.80), None, apply,
                               apply_all=apply_all, apply_common=bool(self.v["t_common"].get()),
                               day_rule="matched" if self.v["t_day_rule"].get().startswith("matched") else "common",
                               tol_recall=self.fnum("t_tol_recall", 0.05), log=self.logln)
        self.q.put(("call", lambda: self._tune_text(m / "aqua_tune" / "tune_recommendation.txt")))
        self.q.put(("call", lambda: self.show_tune_figure(quiet=True)))
        return res

    def run_tune_report(self):
        m = self.mouse()
        if m is None:
            return
        self.start("tune: report", lambda: self._report_job(m, bool(self.v["t_apply"].get()), bool(self.v["t_apply_all"].get())))

    def _tune_text(self, path: Path):
        """Show a text file (the recommendation) in the tab's box."""
        try:
            txt = Path(path).read_text(encoding="utf-8")
        except Exception as e:
            txt = f"(could not read {path}: {e})"
        self.tune_box.configure(state="normal"); self.tune_box.delete("1.0", "end")
        self.tune_box.insert("end", txt); self.tune_box.configure(state="disabled")

    def check_tune_progress(self):
        """Count the sweep's finished / failed / pending jobs from the planned folders and the MATLAB log."""
        m = self.mouse()
        if m is None:
            return
        import aqua_tune
        try:
            pr = aqua_tune.sweep_progress(m)
        except SystemExit as e:
            self.tune_lbl.configure(text=str(e)); return
        self.tune_lbl.configure(text=pr["head"])
        self.tune_box.configure(state="normal"); self.tune_box.delete("1.0", "end")
        self.tune_box.insert("end", pr["head"] + "\n" + "\n".join(pr["per_movie"]) + ("\n" + pr["log_tail"] if pr["log_tail"] else ""))
        self.tune_box.configure(state="disabled")
        self.logln(pr["head"])

    def rebuild_batch_tuned(self):
        """Tab 3 again: the runs whose FOV now has a newer settings file (what Report + apply wrote) are queued again."""
        m = self.mouse()
        if m is None:
            return
        opts = sorted((m / "aqua_opts").glob("*_Opt.csv")) if (m / "aqua_opts").is_dir() else []
        if not opts:
            messagebox.showinfo("No tuned settings", "aqua_opts\\ holds no *_Opt.csv yet — run 3 Report with 'apply' ticked first.")
            return
        if not messagebox.askyesno("Apply → Batch",
                                   f"{len(opts)} tuned settings file(s) in aqua_opts\\:\n  " + "\n  ".join(p.name for p in opts) +
                                   "\n\nBuild the batch with them now (tab 4; runs whose settings file is newer than their result are queued "
                                   "again)? Then tab 5 (MATLAB) and tab 7 (analyze)."):
            return
        self.nb.select(3)
        self.run_batch()

    # ------------------------------------------------------------------ Run all (guided)
    def run_all(self):
        """Everything from the batch to the reports in one go: tab 4 (settings file or the default) -> tab 3 plan ->
        the sweep in MATLAB -> report + apply (common setting per FOV; the trade-off when the rule is not met) ->
        tab 4 again (tuned runs queued) -> AQuA2 batch + CFU batch in MATLAB -> tab 7 -> tab 8 (both reports).
        At each MATLAB step the command is shown and copied; the chain polls the outputs and continues by itself —
        or starts MATLAB from here when 'start MATLAB from this window' is ticked. Curation (tab 6) stays manual."""
        m = self.mouse()
        if m is None:
            return
        auto_ml = bool(self.v["auto_matlab"].get())
        if auto_ml:
            if self.matlab_exe() is None:
                return
            c = self._as_root_conflict()
            if c:
                messagebox.showerror("Paths", c); return
        aqua = self.v["aqua"].get().strip().strip('"')
        if not Path(aqua).is_dir():
            self.logln(f"note: the AQuA2 folder {aqua!r} is not visible from this Python — fine when MATLAB can see it or has AQuA2 on its "
                       f"path (the generated scripts find aqua_cmd_batch by themselves)")
        from afferent2p import tune as T
        try:
            grid = T.parse_grid(self._tune_grid())
        except ValueError as e:
            messagebox.showerror("grid", str(e)); return
        n_set = 1
        for vals in grid.values():
            n_set *= len(vals)
        skip = bool(self.v["t_skip_tuned"].get())
        how = ("MATLAB is started from this window (matlab -batch)" if auto_ml else
               "at each MATLAB step the command is shown and copied to the clipboard — paste it into a MATLAB window; the chain "
               "checks the outputs every 30 s and continues by itself ('Continue now' / 'Stop' in that window)")
        if not messagebox.askyesno("Run all (guided)",
                                   f"For {m.name}:\n"
                                   f"  1. tab 4: batch from the settings files (a FOV without one gets the default)\n"
                                   f"  2. tab 3: plan — {self.v['t_movies'].get().split(':')[0]} movies, {n_set} grid settings + the current one\n"
                                   f"  3. the tuning sweep in MATLAB (hours: ~1.4 min per job per worker; restartable)\n"
                                   f"  4. report + apply: one common setting per FOV"
                                   f"{'' if self.v['t_common'].get() else ' (per day: unticked)'}; the best trade-off when no setting meets the rule\n"
                                   f"  5. tab 4 again, 6. AQuA2 batch + CFU batch in MATLAB (one paste), 7. tab 7 analysis, 8. tab 8 both reports\n"
                                   + ("  (steps 2-4 are skipped when every FOV already has a tuned aqua_opts file)\n" if skip else "")
                                   + f"\n{how}.\nStop ends it between steps. Curation (tab 6) stays manual — run tabs 7 and 8 again afterwards.\n\nStart?"):
            return
        self.save_settings()
        self.start("run all (guided)", lambda: self._auto_job(m, grid, skip))

    run_auto = run_all                                   # old name

    # ------------------------------------------------------------------ fresh start (v35)
    def fresh_start_dialog(self):
        """Park (or delete) this animal's pipeline outputs, stage by stage, and optionally run everything again
        from the raw recordings with the current version — so that no result of an older version is left in."""
        m = self.mouse()
        if m is None:
            return
        if self.busy:
            messagebox.showinfo("Busy", "A step is still running — wait for it or press Stop first.")
            return
        from afferent2p import reset as R
        win = tk.Toplevel(self); win.title(f"Fresh start — {m.name}"); win.geometry("820x560"); win.transient(self)
        ttk.Label(win, text="Everything the pipeline made for this animal is moved into <mouse>\\_previous_<date>\\ (a rename: instant, "
                            "reversible — delete that folder when the new run is through), then the steps run again from the raw "
                            "recordings with the current version. The recordings themselves, day_labels.csv, settings files saved by hand "
                            "(*_Opt.csv), hand-drawn fibre masks, reg_compare\\ and the old suite2p outputs are never touched.",
                  wraplength=790, justify="left").pack(anchor="w", padx=10, pady=(10, 6))
        v = {"unpack": tk.BooleanVar(value=True), "tune": tk.BooleanVar(value=True), "aqua": tk.BooleanVar(value=True),
             "analysis": tk.BooleanVar(value=True), "decisions": tk.BooleanVar(value=False), "delete": tk.BooleanVar(value=False)}
        f = ttk.LabelFrame(win, text="Reset which stages?"); f.pack(fill="x", padx=10, pady=4)
        for k in R.STAGES:
            ttk.Checkbutton(f, text=R.STAGE_LABEL[k], variable=v[k], command=lambda: refresh()).pack(anchor="w", padx=8, pady=1)
        g = ttk.Frame(win); g.pack(fill="x", padx=10, pady=2)
        ttk.Checkbutton(g, text="also discard my decisions (tab 2 frame decisions, tab 6 CFU curation) — default: kept, they still apply",
                        variable=v["decisions"], command=lambda: refresh()).pack(anchor="w", padx=8)
        ttk.Checkbutton(g, text="DELETE the old results instead of parking them (no way back; saves the disk space)",
                        variable=v["delete"], command=lambda: refresh()).pack(anchor="w", padx=8)
        est_lbl = ttk.Label(win, text="", justify="left", wraplength=790, foreground="#333333")
        est_lbl.pack(anchor="w", padx=12, pady=6)

        def refresh():
            stages = [k for k in R.STAGES if v[k].get()]
            try:
                found = R.collect(m, stages, bool(v["decisions"].get()))
                lines = R.describe(found, m)
                est = R.estimate_hours(m, stages)
                hours = [f"unpack {est['unpack_h']:.1f} h" if "unpack_h" in est else "",
                         f"tuning sweep ~{est['tune_h']:.0f} h of MATLAB (1 worker; /N with N workers)" if "tune_h" in est else "",
                         f"AQuA2 + CFU ~{est['aqua_h']:.1f} h of MATLAB" if "aqua_h" in est else "",
                         "analysis + reports ~0.5 h" if "analysis_h" in est else ""]
                est_lbl.configure(text=("Would " + ("delete" if v["delete"].get() else "park") + ":\n" + "\n".join(lines)
                                        + f"\n\nRe-run of {est['n_runs']} runs on {est['n_days']} day(s): " + ", ".join(h for h in hours if h)))
            except Exception as e:
                est_lbl.configure(text=f"(could not list: {e})")
        refresh()
        bar = ttk.Frame(win); bar.pack(fill="x", padx=10, pady=8)

        def go(run_after: bool):
            stages = [k for k in R.STAGES if v[k].get()]
            if not stages:
                messagebox.showinfo("Fresh start", "Tick at least one stage."); return
            what = "DELETED" if v["delete"].get() else "parked in _previous_<date>\\"
            if not messagebox.askyesno("Fresh start", f"{m.name}: the outputs of {', '.join(stages)} will be {what}"
                                                       + (", your decisions too" if v["decisions"].get() else "")
                                                       + (".\n\nThen: unpack -> tuning plan -> sweep (MATLAB) -> report + apply -> batch -> "
                                                          "AQuA2 + CFU (MATLAB) -> analysis -> reports, with the current version." if run_after else ".")
                                                       + "\n\nGo?"):
                return
            grid = None
            if run_after:
                from afferent2p import tune as T
                try:
                    grid = T.parse_grid(self._tune_grid())
                except ValueError as e:
                    messagebox.showerror("grid", str(e)); return
            opts = dict(stages=stages, discard_decisions=bool(v["decisions"].get()), delete=bool(v["delete"].get()))
            win.destroy()
            self.save_settings()

            def job():
                res = R.reset_mouse(m, log=self.logln, **opts)
                self.q.put(("call", lambda: self.refresh_state(log_it=False)))
                if not run_after:
                    return
                self.check_stop()
                if "unpack" in stages or not any(rf.glob("*_bin*s.tif") for rf in R.run_folders(m)):
                    import unpack_bin
                    self.set_state("fresh start: unpack")
                    self.logln("FRESH START — unpack (tab 1)")
                    unpack_bin.process_folders([str(m)], seconds=self.fnum("sec", 2.0), overwrite=False,
                                               behavior=bool(self.v["u_behavior"].get()), global_fix=bool(self.v["u_gfix"].get()),
                                               keep_raw=bool(self.v["u_keep_raw"].get()), register=bool(self.v["u_register"].get()),
                                               log=self.logln)
                    self.check_stop()
                self._auto_job(m, grid, skip_tuned=False)
                if res.get("park_dir"):
                    self.logln(f"The old results are in {res['park_dir']} — delete that folder when you are happy with the new run.")
            self.start("fresh start" + (" + run all" if run_after else ""), job)

        ttk.Button(bar, text="Preview (list in the log, nothing moved)",
                   command=lambda: self.start("fresh start (preview)", lambda: R.reset_mouse(
                       m, stages=[k for k in R.STAGES if v[k].get()], discard_decisions=bool(v["decisions"].get()),
                       delete=bool(v["delete"].get()), dry_run=True, log=self.logln))).pack(side="left")
        ttk.Button(bar, text="Reset only", command=lambda: go(False)).pack(side="left", padx=8)
        ttk.Button(bar, text="Reset and run everything (guided)", command=lambda: go(True)).pack(side="left", padx=8)
        ttk.Button(bar, text="Cancel", command=win.destroy).pack(side="right")

    def _matlab_step(self, script: Path, what: str, done_fn, rel: str):
        """One MATLAB step of the chain: matlab -batch from here (checkbox + matlab.exe) or the manual paste with polling."""
        exe = self.v["matlab"].get().strip().strip('"')
        if bool(self.v["auto_matlab"].get()) and exe and Path(exe).exists():
            self._matlab_job(script, what)
            st = done_fn()
            if not st.get("done"):
                self.logln(f"  NOTE: MATLAB ended but {what} is not complete: {st.get('text', '')} — continuing with what exists")
            return
        self._manual_matlab_step(what, done_fn, self.matlab_cmd_for(rel))

    def _manual_matlab_step(self, what: str, done_fn, cmd: str, poll_s: float = 30.0, stall_s: float = 1200.0):
        """Show the command, copy it, poll done_fn() every poll_s seconds until it reports done (or 'Continue now').
        When the outputs stopped changing for stall_s after some progress (MATLAB ended with results missing, e.g. a
        run whose CFU step failed), the chain goes on with what exists and says so."""
        self.manual_evt.clear(); self.manual_force = False
        self._manual_done_fn = done_fn
        self.q.put(("call", lambda: self._open_manual_window(what, cmd)))
        self.logln(f"  paste this into a MATLAB window (copied to the clipboard):\n  {cmd}\n  waiting — the outputs are checked every "
                   f"{poll_s:.0f} s; 'Continue now' in the small window goes on at once, Stop ends the chain")
        first = last = None
        t_change = time.time()
        while True:
            try:
                st = done_fn()
            except Exception as e:
                st = {"done": False, "text": f"check failed: {e}"}
            txt = str(st.get("text", ""))
            if first is None:
                first = txt
            if txt != last:
                self.logln("  " + txt); last = txt; t_change = time.time()
            self.q.put(("call", lambda t=txt: self._manual_status(t)))
            # the progress row: results present / expected (AQuA2 + CFU per tiff, or the sweep's jobs)
            pr = st.get("progress") or {}
            if pr.get("total"):
                self.set_progress(int(pr.get("done", 0)), int(pr["total"]), f"{what}: {txt}")
            elif st.get("n_queued") or st.get("n_results"):
                total = int(st.get("n_queued", 0)) + int(st.get("n_results", 0))
                self.set_progress(int(st.get("aqua_done", 0)) + int(st.get("cfu_done", 0)), total, f"{what}: {txt}")
            if st.get("done"):
                self.logln(f"  {what}: complete")
                break
            if self.manual_evt.is_set():
                self.manual_evt.clear()
                if self.manual_force:
                    self.logln(f"  continuing without waiting for {what} (Continue anyway)")
                    break
            # some progress = the status changed since the first check, or results already existed at the first check
            pr = st.get("progress") or {}
            partial = bool(st.get("aqua_done") or st.get("cfu_done") or pr.get("done"))
            if (txt != first or partial) and time.time() - t_change > stall_s:
                self.logln(f"  WARNING: the outputs of {what} have not changed for {stall_s / 60:.0f} min — MATLAB has probably finished "
                           f"with something missing ({txt}); continuing with what exists. Check the MATLAB window / aqua_batch\\run_*_log.txt.")
                break
            self.check_stop()
            self.manual_evt.wait(poll_s)
            self.check_stop()
        self.q.put(("call", self._close_manual_window))

    def _open_manual_window(self, what: str, cmd: str):
        self._close_manual_window()
        try:
            self.clipboard_clear(); self.clipboard_append(cmd)
        except Exception:
            pass
        win = tk.Toplevel(self); win.title(f"MATLAB step — {what}")
        win.geometry("760x300")
        ttk.Label(win, text=f"Paste this into a MATLAB window and press Enter ({what}). It is already on the clipboard.",
                  wraplength=720, justify="left").pack(anchor="w", padx=10, pady=(10, 4))
        box = tk.Text(win, height=3, font=("Consolas", 12), wrap="word"); box.pack(fill="x", padx=10)
        box.insert("end", cmd); box.configure(state="disabled")
        bar = ttk.Frame(win); bar.pack(fill="x", padx=10, pady=8)
        ttk.Button(bar, text="Copy again", command=lambda: (self.clipboard_clear(), self.clipboard_append(cmd))).pack(side="left")
        ttk.Button(bar, text="Continue now", command=self._manual_continue).pack(side="left", padx=8)
        ttk.Button(bar, text="Stop the chain", command=self.stop_run).pack(side="left", padx=8)
        self.manual_lbl = ttk.Label(win, text="checking the outputs every 30 s …", wraplength=720, justify="left", foreground="#333333")
        self.manual_lbl.pack(anchor="w", padx=10, pady=4)
        ttk.Label(win, text="The chain continues by itself when every expected output exists (the MATLAB window can stay open). "
                            "'Continue now' checks at once and, when something is still missing, asks before going on with what exists.",
                  wraplength=720, justify="left", foreground="#555555").pack(anchor="w", padx=10, pady=4)
        win.protocol("WM_DELETE_WINDOW", lambda: None)          # the worker closes it
        self.manual_win = win

    def _manual_status(self, text: str):
        try:
            if self.manual_win is not None and self.manual_win.winfo_exists():
                self.manual_lbl.configure(text=f"{time.strftime('%H:%M:%S')}  {text}")
        except Exception:
            pass

    def _manual_continue(self):
        try:
            st = self._manual_done_fn() if self._manual_done_fn else {"done": True, "text": ""}
        except Exception as e:
            st = {"done": False, "text": f"check failed: {e}"}
        if st.get("done") or messagebox.askyesno("Not complete", f"{st.get('text', '')}\n\nContinue anyway with what exists?"):
            self.manual_force = not st.get("done")
            self.manual_evt.set()

    def _close_manual_window(self):
        try:
            if self.manual_win is not None and self.manual_win.winfo_exists():
                self.manual_win.destroy()
        except Exception:
            pass
        self.manual_win = None

    def _auto_job(self, m: Path, grid: dict, skip_tuned: bool):
        import aqua_tune
        import aqua_batch
        t0 = time.time()
        self.logln("RUN ALL 1/8 — batch (tab 4)")
        self.set_state("run all 1/8 batch")
        st = self._batch_job(m, overwrite=False)
        self.check_stop()
        exps = []
        try:
            man = json.loads((m / "aqua_batch" / "manifest.json").read_text())
            exps = [e for e in man["experiments"] if e.get("runs")]
        except Exception:
            pass
        tuned = [e["name"] for e in exps if (m / "aqua_opts" / f"{e['name']}_Opt.csv").exists()]

        def sweep_done():
            pr = aqua_tune.sweep_progress(m)
            return {"done": pr["total"] > 0 and pr["pending"] == 0, "text": pr["head"], "progress": pr}
        if skip_tuned and exps and len(tuned) == len(exps):
            self.logln(f"RUN ALL 2-4/8 — skipped: every experiment already has a tuned settings file in aqua_opts\\ ({len(tuned)})")
        else:
            self.logln("RUN ALL 2/8 — tuning plan (tab 3, 1 Plan)")
            self.set_state("run all 2/8 tuning plan")
            self._tune_plan_job(m, grid)
            self.check_stop()
            pr = aqua_tune.sweep_progress(m)
            self.logln("RUN ALL 3/8 — the tuning sweep in MATLAB: " + pr["head"])
            if pr["pending"]:
                self.set_state(f"run all 3/8 sweep in MATLAB ({pr['pending']} jobs)")
                self._matlab_step(m / self.actual_file("aqua_tune/run_aqua_tune.m"), "the tuning sweep", sweep_done, "aqua_tune/run_aqua_tune.m")
                pr = aqua_tune.sweep_progress(m)
                self.logln("  " + pr["head"])
            self.check_stop()
            if pr["done"] == 0:
                raise RuntimeError("no sweep job finished — nothing to score (see aqua_tune\\tune_log.txt)")
            if pr["pending"] or pr["failed"]:
                self.logln(f"  WARNING: {pr['pending']} job(s) pending and {pr['failed']} failed — scoring what finished")
            self.logln("RUN ALL 4/8 — tuning report + apply (common setting per FOV; the trade-off when the rule is not met)")
            self.set_state("run all 4/8 tuning report")
            res = self._report_job(m, apply=True, apply_all=True)
            # the sweep did not bracket the optimum (v37): one more round with the grid extended as the report suggests
            # (finished jobs are kept; only the new combinations are swept), then report + apply again
            if res.get("grid_suggest"):
                from afferent2p import tune as T
                import aqua_tune as AT
                grid2 = AT.merge_grid(grid, res["grid_suggest"])
                added = {k: [v for v in grid2[k] if v not in [float(x) for x in grid.get(k, [])]] for k in grid2}
                added = {k: v for k, v in added.items() if v}
                if added:
                    self.logln("RUN ALL 4/8 — the grid only touched the optimum for " + ", ".join(f"{AT.short_run(r)} ({st_})" for r, st_ in res.get("grid_flagged", {}).items())
                               + " — extending it: " + "; ".join(f"{k} + {', '.join(T._fmt(v) for v in vals)}" for k, vals in added.items())
                               + " (one extra round: plan → sweep → report)")
                    self.q.put(("call", lambda g=grid2: self._set_grid_fields(g)))
                    self.check_stop()
                    self._tune_plan_job(m, grid2)
                    self.check_stop()
                    pr = aqua_tune.sweep_progress(m)
                    if pr["pending"]:
                        self.set_state(f"run all 3/8 sweep in MATLAB, extended grid ({pr['pending']} jobs)")
                        self._matlab_step(m / self.actual_file("aqua_tune/run_aqua_tune.m"), "the tuning sweep (extended grid)", sweep_done, "aqua_tune/run_aqua_tune.m")
                    self.check_stop()
                    self.set_state("run all 4/8 tuning report (extended grid)")
                    res = self._report_job(m, apply=True, apply_all=True)
                    if res.get("grid_suggest"):
                        self.logln("  NOTE: the report still suggests a wider grid (" + ", ".join(f"{k}={','.join(T._fmt(v) for v in vals)}" for k, vals in res["grid_suggest"].items())
                                   + ") — 'Extend the grid (as suggested)' in tab 3 goes one step further when wanted")
            if res.get("amplitudes_suggest"):
                self.logln("  NOTE: the injected transients are dimmer than this field's events (see the report) — set 'injected dF/F' to "
                           + ", ".join(f"{a:g}" for a in res["amplitudes_suggest"]) + " in tab 3, press 1 Plan and run the sweep again (only the inject jobs)")
            if not res["written"]:
                raise RuntimeError("the report applied no setting (nothing scored?) — the chain stops here")
            if res["rule_not_met"]:
                self.logln(f"  NOTE: rule not met for {', '.join(res['rule_not_met'])} — the best trade-off was written; "
                           f"look at the report / the curation window before trusting these fields")
            self.check_stop()
            self.logln("RUN ALL 5/8 — batch again with the tuned settings (tab 4)")
            self.set_state("run all 5/8 batch (tuned)")
            st = self._batch_job(m, overwrite=False)
            self.check_stop()
        bd = aqua_batch.batch_done(m)
        if st.get("queued") or not bd["done"]:
            self.logln(f"RUN ALL 6/8 — AQuA2 batch + CFU batch in MATLAB ({st.get('queued', 0)} tiff(s) queued; {bd['text']})")
            self.set_state(f"run all 6/8 AQuA2 + CFU batch in MATLAB ({st.get('queued', 0)} tiffs)")
            self._matlab_step(m / self.actual_file("aqua_batch/run_aqua_and_cfu.m"), "the AQuA2 batch + CFU batch", lambda: aqua_batch.batch_done(m),
                              "aqua_batch/run_aqua_and_cfu.m")
        else:
            self.logln("RUN ALL 6/8 — MATLAB: nothing to do (every result and CFU file is up to date)")
        self.check_stop()
        self.logln("RUN ALL 7/8 — analysis (tab 7)")
        self.set_state("run all 7/8 analysis")
        self._analyze_job(m)
        self.check_stop()
        self.logln("RUN ALL 8/8 — the report and the QC report (tab 8)")
        self.set_state("run all 8/8 reports")
        out = qc = None
        try:
            out = self._html_report_job(m, open_after=False)
        except Exception as e:
            self.logln(f"  report failed ({e}) — build it from tab 8 after checking the log")
        self.check_stop()
        try:
            qc = self._html_report_job(m, open_after=False, qc=True)
        except Exception as e:
            self.logln(f"  QC report failed ({e}) — build it from tab 8 after checking the log")
        self.logln(f"RUN ALL finished in {(time.time() - t0) / 3600:.1f} h." + (f" Report: {out}" if out else "") + (f" QC: {qc}" if qc else "")
                   + " Next: look at the QC report and tab 6 (curate the CFUs — mark noise, add missed regions), then tabs 7 and 8 again; "
                     "the tuning figure and tune_recommendation.txt say how good the detection is per FOV.")

    def _batch_job(self, m: Path, overwrite: bool | None = None) -> dict:
        """aqua_batch.build_mouse with the tab-4 settings; returns its stats; raises when nothing can run."""
        import aqua_batch
        fov_txt = self.v["b_fov"].get().replace(",", " ").split()
        st: dict = {}
        kw = dict(aqua_dir=self.v["aqua"].get().strip(), sec=self.fnum("sec", 2.0),
                  um_per_px=self.fnum("um", 0.792), border=self.fnum("border", 5),
                  stim_after_run=int(self.v["b_stim_after"].get()) if self.v["b_stim_after"].get().strip() else None,
                  cfu_min_events=int(self.fnum("b_cfu_min", CFU_MIN_EVENTS_DEFAULT)),
                  cfu_overlap=self.fnum("b_cfu_overlap", CFU_OVERLAP_DEFAULT),
                  cfu_merge=self.fnum("b_cfu_merge", 0.5), cfu_aqua_dir=self.v["b_cfu_aqua"].get().strip().strip('"') or None,
                  cfu_redo=bool(self.v["b_cfu_redo"].get()),
                  overwrite=bool(self.v["b_overwrite"].get()) if overwrite is None else overwrite, tol=self.fnum("b_sec_tol", 0.25),
                  only_fov=[int(x) for x in fov_txt] if fov_txt else None,
                  opts_default=Path(self.v["b_opts"].get().strip().strip('"')) if self.v["b_opts"].get().strip() else None,
                  log=self.logln, stats=st)
        out = aqua_batch.build_mouse(m, **kw)
        if out is None:
            raise RuntimeError(f"no runs found under {m} — unpack the recordings first (tab 1)")
        if not st.get("experiments"):
            raise RuntimeError("no experiment got a settings file — set tab 4's default to default_Opt.csv (or restore the "
                               "pipeline's default_Opt.csv) and build again")
        if not st.get("queued") and not st.get("done"):
            raise RuntimeError("nothing queued and nothing analysed: the experiments have no binned tiffs — run tab 1 (unpack) first")
        return st

    def run_batch_matlab(self, which: str):
        """Tab 5 from here: the AQuA2 batch script, the CFU script, or both, in matlab -batch."""
        m = self.mouse()
        if m is None or self.matlab_exe() is None:
            return
        c = self._as_root_conflict()
        if c:
            messagebox.showerror("Paths", c); return
        scripts = {"aqua": [(self.actual_file("aqua_batch/run_aqua_batch.m").split("/")[-1], "the AQuA2 batch")],
                   "cfu": [(self.actual_file("aqua_batch/run_cfu_batch.m").split("/")[-1], "the CFU batch")]}
        todo = scripts["aqua"] + scripts["cfu"] if which == "both" else scripts[which]
        for name, _ in todo:
            if not (m / "aqua_batch" / name).exists():
                messagebox.showerror("Batch", f"aqua_batch\\{name} does not exist — build the batch in tab 4 first.")
                return
        if which != "cfu":
            try:
                bname = "batch.csv"
                try:
                    bname = json.loads((m / "aqua_batch" / "manifest.json").read_text()).get("batch_csv", "batch.csv")
                except Exception:
                    pass
                n = sum(1 for ln in (m / "aqua_batch" / bname).read_text().splitlines()[1:] if ln.strip())
            except Exception:
                n = -1
            if n == 0:
                messagebox.showinfo("Batch", "batch.csv has no rows: every tiff already has a result (tab 4: tick overwrite to redo them)"
                                             + ("" if which == "aqua" else " — running the CFU batch only"))
                if which == "aqua":
                    return
                todo = scripts["cfu"]
            elif not messagebox.askyesno("Run in MATLAB", f"Run {n} tiff(s) through AQuA2 now? (several minutes each; Stop ends it)"):
                return

        def job():
            for name, what in todo:
                self.set_state(f"MATLAB: {what}")
                self._matlab_job(m / "aqua_batch" / name, what)
                self.check_stop()
        self.start("MATLAB: " + " + ".join(w for _, w in todo), job)

    def open_tune_why(self):
        """tune_why.html (v35): per tuned movie, what each setting finds in the real / null / injected movie, as pictures."""
        m = self.mouse()
        if m is None:
            return
        p = m / "aqua_tune" / "tune_why.html"
        if not p.exists():
            messagebox.showinfo("Missing", "aqua_tune\\tune_why.html does not exist yet — press '3 Report' after the MATLAB sweep "
                                           "(a report made before v35 has to be run again to get the pictures).")
            return
        open_path(p)

    def show_tune_figure(self, quiet: bool = False):
        """The tuning figure in a window of its own, scaled to the screen (the file itself is full size)."""
        m = self.mouse()
        if m is None:
            return
        p = m / "aqua_tune" / f"{m.name}_tune.png"
        if not p.exists():
            if not quiet:
                messagebox.showinfo("Missing", f"{p.name} does not exist yet — run Report after the MATLAB script has finished.")
            return
        self.show_pictures(p, title=f"AQuA2 tuning — {m.name}")

    def show_pictures(self, first: Path, title: str = ""):
        """The picture viewer (v38) on the tuning figure and every 'why' picture of aqua_tune\\, starting at `first`."""
        m = self.mouse()
        try:
            import importlib
            importlib.import_module("PIL.ImageTk")   # the viewer needs Pillow
        except Exception:
            open_path(first); return
        tune_dir = Path(first).parent
        files = [Path(first)]
        if m is not None:
            for pat in (f"{m.name}_tune.png", f"{m.name}_tune_detail.png", "*_why_day.png", "*_why.png"):
                for f in sorted(tune_dir.glob(pat)):
                    if f not in files:
                        files.append(f)

        def redraw(p: Path):
            if m is not None and p.name in (f"{m.name}_tune.png", f"{m.name}_tune_detail.png"):
                import aqua_tune
                aqua_tune.figure_only(m, dpi=260, log=self.logln)
                return p
            return None
        try:
            ImageViewer(self, files, index=0, title=title, redraw=redraw)
        except Exception as e:
            self.logln(f"(picture window failed: {e}) — opening the file instead")
            open_path(first)

    # ------------------------------------------------------------------ tab 4: batch
    def run_overlap_check(self):
        """v42: the overlap-threshold check on this animal's CFU results (afferent2p/cfu_overlap.py)."""
        m = self.mouse()
        if m is None:
            return
        thr = self.fnum("b_cfu_overlap", CFU_OVERLAP_DEFAULT)

        def job():
            from afferent2p import cfu_overlap
            res = cfu_overlap.check_mouse(m, chosen=thr, log=self.logln)
            png = m / "aqua_batch" / "cfu_overlap_check.png"
            if png.exists():
                self.q.put(("call", lambda: self.show_pictures(png, title=f"CFU overlap-threshold check — {m.name}")))
            n_check = sum(1 for e in res["experiments"] if e["verdict"] == "CHECK")
            self.logln(f"overlap check: pooled {res['pooled']['verdict']}, {n_check} experiment(s) say CHECK — "
                       f"details in {m / 'aqua_batch' / 'cfu_overlap_check.txt'}")
        self.start("overlap check", job)

    def run_batch(self):
        m = self.mouse()
        if m is None:
            return
        self.save_settings()

        def job():
            try:
                self._batch_job(m)
            except RuntimeError as e:
                msg = str(e)
                self.q.put(("call", lambda: messagebox.showerror("Batch", msg)))
                raise
        self.start("batch", job)

    # ------------------------------------------------------------------ tab 5: MATLAB
    def matlab_cmds(self) -> str:
        return (f"{self.matlab_cmd_for('aqua_batch/run_aqua_and_cfu.m')}\n"
                f"% = the two scripts, one after the other:\n"
                f"%   {self.matlab_cmd_for('aqua_batch/run_aqua_batch.m')}\n"
                f"%   {self.matlab_cmd_for('aqua_batch/run_cfu_batch.m')}\n")

    def show_matlab(self):
        if self.mouse() is None:
            return
        self.matlab_box.delete("1.0", "end")
        self.matlab_box.insert("end", self.matlab_cmds())
        mdir = Path(self.v["mouse"].get().strip().strip('"'))
        if not (mdir / self.actual_file("aqua_batch/run_aqua_and_cfu.m")).exists():
            self.matlab_box.insert("end", "\n(the batch files do not exist yet — build them in tab 4 first)")
        if (mdir / self.actual_file("aqua_tune/run_aqua_tune.m")).exists():
            self.matlab_box.insert("end", f"\n% tuning sweep of tab 3 (overnight, before the above):\n{self.tune_matlab_cmd()}\n")

    def copy_matlab(self):
        if self.mouse() is None:
            return
        cmd = self.matlab_cmd_for("aqua_batch/run_aqua_and_cfu.m")
        self.clipboard_clear(); self.clipboard_append(cmd)
        self.logln(f"copied: {cmd}")

    def check_batch_done(self):
        m = self.mouse()
        if m is None:
            return
        import aqua_batch
        bd = aqua_batch.batch_done(m)
        self.matlab_lbl.configure(text=("complete — " if bd["done"] else "not complete — ") + bd["text"])
        self.logln(bd["text"])

    # ------------------------------------------------------------------ tab 6: curate
    def run_curate_gui(self):
        m = self.mouse()
        if m is None:
            return
        script = HERE / "cfu_curate.py"
        try:
            cmd = [sys.executable, str(script), str(m), "--sec", str(self.fnum("sec", 2.0)), "--um-per-px", str(self.fnum("um", 0.792))]
            tune = self.v["cur_tune"].get().strip()
            if tune:
                cmd += ["--tune", str(int(float(tune)))]
            subprocess.Popen(cmd, cwd=str(HERE))
            self.logln(f"opened cfu_curate.py for {m.name} in its own window" + (f" with tuning setting {tune} overlaid" if tune else ""))
        except Exception as e:
            messagebox.showerror("cfu_curate", str(e))

    def run_curate_list(self):
        m = self.mouse()
        if m is None:
            return

        def job():
            r = subprocess.run([sys.executable, str(HERE / "cfu_curate.py"), str(m), "--list", "--sec", str(self.fnum("sec", 2.0))],
                               cwd=str(HERE), capture_output=True, text=True)
            self.logln(r.stdout.strip() or "(no output)")
            if r.stderr.strip():
                self.logln(r.stderr.strip())
        self.start("CFU list", job)

    # ------------------------------------------------------------------ tab 7: analyze
    def _analyze_job(self, m: Path):
        import aqua_analyze
        kw = dict(um_per_px=self.fnum("um", 0.792), sec=self.fnum("sec", 2.0), gate=self.v["a_gate"].get(),
                  min_p=self.fnum("a_min_p", 0.0), baseline=self.v["a_baseline"].get().strip() or None,
                  cfu=not self.v["a_no_cfu"].get(), event_source=self.v["a_source"].get(),
                  use_curation=not self.v["a_no_curation"].get(), flag_margin=self.fnum("a_margin", 1),
                  keep_flagged=bool(self.v["a_keep_flagged"].get()), pupil_analysis=bool(self.v["a_pupil"].get()),
                  maps=bool(self.v["a_maps"].get()), tile_um=self.fnum("a_tile", 25.0),
                  loco=bool(self.v["a_loco"].get()), loco_glm=bool(self.v["a_glm"].get()), heatmaps=bool(self.v["a_heat"].get()),
                  movies=bool(self.v["a_movies"].get()), motion=bool(self.v["a_motion"].get()), log=self.logln)
        aqua_analyze.analyze_mouse(m, **kw)
        self.q.put(("call", self.show_summary))

    def run_analyze(self):
        m = self.mouse()
        if m is None:
            return
        self.save_settings()
        self.start("analyze", lambda: self._analyze_job(m))

    def show_summary(self):
        m = self.v["mouse"].get().strip().strip('"')
        p = Path(m) / "aqua_analysis" / "summary.txt"
        self.summary_box.configure(state="normal"); self.summary_box.delete("1.0", "end")
        self.summary_box.insert("end", p.read_text(encoding="utf-8", errors="replace") if p.exists()
                                else "no summary.txt yet — run the analysis")
        self.summary_box.configure(state="disabled")

    # ------------------------------------------------------------------ step 8
    def report_path(self, m: Path, qc: bool = False) -> Path:
        return m / "aqua_analysis" / (f"{m.name}_qc.html" if qc else f"{m.name}_report.html")

    def _html_report_job(self, m: Path, open_after: bool = True, qc: bool = False):
        if qc:
            from afferent2p.qc_report import build_qc_report
            out = build_qc_report(m, max_events=int(self.fnum("r_max_events", 0)), clips_per_run=int(self.fnum("r_clips", 6)),
                                  no_clips=bool(self.v["r_no_clips"].get()), um_per_px=self.fnum("um", 0.792), log=self.logln)
        else:
            from afferent2p.report import build_report
            out = build_report(m, log=self.logln)
        self.q.put(("call", lambda: self.report_lbl.configure(text=f"{'QC report' if qc else 'report'}: {out}  ({out.stat().st_size / 1e6:.1f} MB)")))
        if open_after:
            self.q.put(("call", lambda: open_path(out)))
        return out

    def run_report(self, qc: bool = False):
        m = self.mouse()
        if m is None:
            return
        if not (m / "aqua_batch" / "manifest.json").exists():
            messagebox.showerror("Report", "No aqua_batch\\manifest.json — build the batch (tab 4) and run the analysis (tab 7) first.")
            return
        if not (m / "aqua_analysis" / "experiments_summary.csv").exists():
            if not messagebox.askyesno("Report", "No analysis output yet (tab 7) — most sections will be empty. Build anyway?"):
                return
        self.save_settings()
        self.start("QC report" if qc else "report", lambda: self._html_report_job(m, open_after=bool(self.v["r_open"].get()), qc=qc))

    def open_report(self, qc: bool = False):
        m = self.mouse()
        if m is None:
            return
        p = self.report_path(m, qc)
        if not p.exists():
            messagebox.showinfo("Report", f"{p.name} does not exist yet — press '{'Build the QC report' if qc else 'Build the report'}'.")
            return
        open_path(p)

    # ------------------------------------------------------------------ across animals (v35)
    def parent_dir(self) -> Path | None:
        p = self.v["x_parent"].get().strip().strip('"')
        if not p:
            m = self.v["mouse"].get().strip().strip('"')
            if m:
                p = str(Path(m).parent)
                self.v["x_parent"].set(p)
        if not p or not Path(p).is_dir():
            messagebox.showerror("Across animals", "Choose the parent folder — the folder that holds the mouse folders "
                                                   "(e.g. D:\\2photon\\Chao\\GCAMP_affrerents).")
            return None
        return Path(p)

    def pick_parent(self):
        d = filedialog.askdirectory(title="The folder that holds the mouse folders")
        if d:
            self.v["x_parent"].set(d)

    def run_conditions(self):
        p = self.parent_dir()
        if p is None:
            return
        self.save_settings()

        def job():
            import aqua_conditions
            res = aqua_conditions.run(p, log=self.logln)
            if bool(self.v["r_open"].get()) and res.get("html") and Path(res["html"]).exists():
                self.q.put(("call", lambda: open_path(Path(res["html"]))))
        self.start("protocols across animals", job)

    def open_conditions(self):
        p = self.parent_dir()
        if p is None:
            return
        h = p / "cross_animal" / "conditions.html"
        if not h.exists():
            messagebox.showinfo("Across animals", f"{h} does not exist yet — press 'Compare the protocols'.")
            return
        open_path(h)


class SheetWindow(tk.Toplevel):
    """The removed-frames sheet of one recording, with its per-frame table."""

    def __init__(self, master, png: Path, csv_path: Path, qc: Path | None = None):
        super().__init__(master)
        self.title(f"flagged frames — {png.name}")
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry(f"{min(1640, sw - 80)}x{min(980, sh - 120)}+30+30")
        top = ttk.Frame(self); top.pack(fill="x", padx=6, pady=4)
        ttk.Label(top, text=f"{png.parent.name}\\{png.name}", foreground="#555555").pack(side="left")
        ttk.Button(top, text="Open in image viewer", command=lambda: open_path(png)).pack(side="right")
        pics = [png] + ([Path(qc)] if qc and Path(qc).exists() else [])
        ttk.Button(top, text="Zoom / pan" + (" (+ QC plot)" if len(pics) > 1 else ""),
                   command=lambda: ImageViewer(self, pics, 0, title=f"check sheet — {png.parent.name}")).pack(side="right", padx=8)
        self.scale = tk.IntVar(value=1)
        ttk.Checkbutton(top, text="half size", variable=self.scale, onvalue=2, offvalue=1,
                        command=self.redraw).pack(side="right", padx=8)
        pane = ttk.Panedwindow(self, orient="vertical"); pane.pack(fill="both", expand=True, padx=6, pady=4)
        frame = ttk.Frame(pane); pane.add(frame, weight=4)
        self.canvas = tk.Canvas(frame, background="white")
        hs = ttk.Scrollbar(frame, orient="horizontal", command=self.canvas.xview)
        vs = ttk.Scrollbar(frame, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=hs.set, yscrollcommand=vs.set)
        vs.pack(side="right", fill="y"); hs.pack(side="bottom", fill="x")
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<MouseWheel>", lambda e: self.canvas.yview_scroll(-int(e.delta / 120), "units"))
        self.png = png
        self.img = None
        self.redraw()
        tbl = ttk.Frame(pane); pane.add(tbl, weight=2)
        cols = ("frame", "stretch", "cls", "reason", "level", "shape", "field", "diff", "floor", "ratio", "speed",
                "flagged", "verdict", "decision", "excluded")
        tree = ttk.Treeview(tbl, columns=cols, show="headings", height=8)
        widths = {"frame": 50, "stretch": 60, "cls": 60, "reason": 210, "level": 55, "shape": 55, "field": 60, "diff": 65,
                  "floor": 60, "ratio": 50, "speed": 50, "flagged": 55, "verdict": 200, "decision": 90, "excluded": 60}
        names = {"frame": "frame", "stretch": "stretch", "cls": "class", "reason": "why flagged", "level": "level %",
                 "shape": "shape %", "field": "field %", "diff": "diff rms %", "floor": "clean floor %", "ratio": "x floor",
                 "speed": "speed", "flagged": "flagged", "verdict": "stretch verdict", "decision": "decision",
                 "excluded": "excluded"}
        for c in cols:
            tree.heading(c, text=names[c]); tree.column(c, width=widths[c], anchor="w", stretch=(c in ("reason", "verdict")))
        sb = ttk.Scrollbar(tbl, command=tree.yview); tree.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y"); tree.pack(fill="both", expand=True)
        tree.tag_configure("artefact", foreground="#b30000")
        tree.tag_configure("state", foreground="#a05a00")
        tree.tag_configure("excluded", background="#ffe9e9")
        try:
            import csv
            for r in csv.DictReader(Path(csv_path).open(encoding="utf-8")):
                vals = (r.get("frame", ""), r.get("stretch", ""), r.get("class", ""), r.get("reason", ""),
                        r.get("level_dev_pct", ""), r.get("shape_dev_pct", ""), r.get("whole_field_change_pct", ""),
                        r.get("diff_from_neighbours_rms_pct", ""), r.get("clean_neighbour_floor_rms_pct", ""),
                        r.get("ratio", ""), r.get("speed", ""), "yes" if r.get("flagged_frame", "") == "1" else "",
                        r.get("stretch_verdict", ""), r.get("decision", ""), "yes" if r.get("excluded", "") == "1" else "")
                tree.insert("", "end", values=vals,
                            tags=(r.get("class", ""),) + (("excluded",) if r.get("excluded", "") == "1" else ()))
        except Exception as e:
            tree.insert("", "end", values=("", "", "", f"could not read {Path(csv_path).name}: {e}"))

    def redraw(self):
        try:
            img = tk.PhotoImage(file=str(self.png))
            if self.scale.get() > 1:
                img = img.subsample(self.scale.get(), self.scale.get())
        except Exception as e:
            self.canvas.delete("all")
            self.canvas.create_text(20, 20, anchor="nw", text=f"could not load {self.png.name}: {e}")
            return
        self.img = img                                   # keep a reference
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=img)
        self.canvas.configure(scrollregion=(0, 0, img.width(), img.height()))


def main():
    app = App()
    if app.v["mouse"].get().strip():
        app.after(300, app.refresh_status)
        app.after(400, app.refresh_state)
    app.mainloop()


if __name__ == "__main__":
    main()
