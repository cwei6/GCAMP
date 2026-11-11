#!/usr/bin/env python3
"""afferent2p GUI — pick one or more experiment folders, label runs pre/post, run.

Usage:  python gui.py

Needs only the Python standard library (tkinter) on top of the pipeline's own
dependencies. Each folder you add = one experiment (e.g.
...\\GCAMP_afferents\\NavRearChrTIGRE-2608-1\\NavRearChrTIGRE-2608-1_260831_001).
Runs (.sbx) are discovered inside it, including subfolders. Double-click a run
to toggle its pre/post label. "Run all" processes every queued experiment in
order and streams the log below.
"""
from __future__ import annotations
import queue
import subprocess
import sys
import threading
from pathlib import Path

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import yaml
from afferent2p.plan import plan_mouse, describe as describe_plan  # noqa: E402


# ----------------------------------------------------------------------------- helpers

def discover_sbx(folder: Path) -> list[Path]:
    return sorted(folder.rglob("*.sbx"), key=lambda p: p.name)


def default_conditions(n: int) -> list[str]:
    if n == 1:
        return ["pre"]          # single-condition mode: no comparison, just analysis
    if n == 2:
        return ["pre", "post"]
    # first half pre, second half post as a starting guess; user can toggle
    return ["pre" if i < n / 2 else "post" for i in range(n)]


def check_dependencies() -> list[str]:
    missing = []
    for mod in ("numpy", "scipy", "yaml", "tifffile", "matplotlib",
                "sbxreader", "suite2p"):
        try:
            __import__(mod)
        except ImportError:
            missing.append("pyyaml" if mod == "yaml" else mod)
    return missing


def build_config(folder: Path, runs: list[tuple[str, str]], settings: dict) -> dict:
    return {
        "experiment": {
            "name": folder.name,
            "data_dir": str(folder),
            "runs": [{"sbx": name, "condition": cond}
                     for name, cond in runs if cond != "skip"],
            "output_dir": None,
        },
        "imaging": {
            "functional_chan": int(settings["functional_chan"]),
            "register_on_red": "auto",   # use red vascular channel for alignment when all runs have it
            "nplanes": 1,
            "fs": None,
            "tau": float(settings["tau"]),
        },
        "prealign": {"enabled": True, "max_angle_deg": 6, "min_shift_px": 2, "min_angle_deg": 0.5},
        "registration": {
            "nonrigid": True,
            "block_size": [128, 128],
            "maxregshift": 0.1,
            "batch_size": 200,
            "torch_device": "cuda" if settings["use_gpu"] else "cpu",
        },
        "roi_detection": {          # axon-oriented suite2p settings
            "enabled": bool(settings["roi_detection"]),
            "soma_crop": False,
            "threshold_scaling": 0.8,
            "spatial_scale": 1,
            "denoise": True,
            "detection_bin_s": 1.0,
            "highpass_time": 100,
            "max_rois": 5000,
            "dcnv_win_baseline_s": 30,
            "save_mat": True,
        },
        "export": {
            "aqua_tifs": bool(settings["aqua_tifs"]),
            "chunk_frames": 4000,
            # uint16: ImageJ/Fiji cannot read signed 16-bit TIFFs, and MATLAB
            # preprocessing can mishandle them. Keep unless you have a reason.
            "aqua_dtype": "uint16",
            # temporal binning for the AQuA tiff: target frame rate; the bin
            # factor is round(fs / rate) exactly as the MATLAB pipeline does
            # (GenerateExptProjections). 1 Hz -> 16 frames at 15.62 Hz.
            "aqua_rate_hz": float(settings["aqua_rate_hz"]),
            "aqua_bin_t": 0,                # 0 = derive from aqua_rate_hz
            "aqua_crop": [80, 80, 40, 40],  # [L,R,T,B] px, = MATLAB projParam.edge
        },
        "loco": {
            "enabled": True,
            "quad_glob": "*quad*",
            "quad_kind": "cumulative",
            "counts_per_cm": settings["counts_per_cm"] or None,
            "still_thresh": float(settings["still_thresh"]),
            "min_still_s": 5.0,
            "smooth_s": 0.5,
        },
        "pupil": {"enabled": True, "glob": "*eye*DLC*.csv", "likelihood_min": 0.6,
                  "min_edge_points": 3, "max_gap_s": 0.5, "smooth_s": 0.25},
        "analysis": {
            "enabled": bool(settings["roi_detection"]),
            "neuropil_coeff": 0.7,
            "baseline_percentile": 10,
            "baseline_win_s": 32.0,
            "event_z": 1.0,
            "event_min_dur_s": 1.0,
            "event_min_dff": 0.05,
            "loco_gate": "still",
            "pupil_gate": "none",
            "cluster_r": "auto",        # ROIs -> afferents by activity correlation
            "cluster_r_min": 0.3,
            "cluster_smooth_s": 0.5,
            "cluster_max_dist_px": 0,
            "cluster_min_still_s": 30,
            # ROI curation (afferent filter) thresholds are NOT written here on
            # purpose: they live in afferent2p/curate.py DEFAULTS (documented in
            # config_example.yaml) so pipeline updates take effect; add
            # analysis.curate_* keys to the yaml only to override them.
            "auto_curate": True,
        },
    }


# ----------------------------------------------------------------------------- app

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("afferent2p — pre/post afferent pipeline")
        self.geometry("980x720")
        self.experiments: dict[str, dict] = {}  # folder -> {"runs": [(name, cond)]}
        self.log_q: queue.Queue[str] = queue.Queue()
        self.running = False
        self._build()
        self.after(100, self._poll_log)

    # ---- layout
    def _build(self):
        top = ttk.Frame(self); top.pack(fill="x", padx=8, pady=6)
        ttk.Button(top, text="Add mouse folder (auto-plan by paradigm)…",
                   command=self.add_mouse_planned).pack(side="left")
        ttk.Button(top, text="Add folder as ONE experiment…", command=self.add_folder).pack(side="left", padx=6)
        ttk.Button(top, text="Remove selected", command=self.remove_selected).pack(side="left", padx=6)
        ttk.Button(top, text="Run all", command=self.run_all).pack(side="right")

        mid = ttk.Frame(self); mid.pack(fill="both", expand=True, padx=8)

        # experiment list
        left = ttk.LabelFrame(mid, text="Experiments (one per session folder)")
        left.pack(side="left", fill="both", expand=True)
        self.exp_tree = ttk.Treeview(left, columns=("runs", "status"), height=8)
        self.exp_tree.heading("#0", text="Folder")
        self.exp_tree.heading("runs", text="# runs")
        self.exp_tree.heading("status", text="Status")
        self.exp_tree.column("#0", width=420)
        self.exp_tree.column("runs", width=60, anchor="center")
        self.exp_tree.column("status", width=100, anchor="center")
        self.exp_tree.pack(fill="both", expand=True)
        self.exp_tree.bind("<<TreeviewSelect>>", lambda e: self.show_runs())

        # run list for selected experiment
        right = ttk.LabelFrame(mid, text="Runs (double-click: planned label <-> skip; manual: pre / post / skip)")
        right.pack(side="left", fill="both", expand=True, padx=(8, 0))
        self.run_tree = ttk.Treeview(right, columns=("cond",), height=8)
        self.run_tree.heading("#0", text=".sbx file")
        self.run_tree.heading("cond", text="condition")
        self.run_tree.column("#0", width=320)
        self.run_tree.column("cond", width=90, anchor="center")
        self.run_tree.pack(fill="both", expand=True)
        self.run_tree.bind("<Double-1>", self.toggle_condition)

        # settings
        s = ttk.LabelFrame(self, text="Settings"); s.pack(fill="x", padx=8, pady=6)
        self.vars = {
            "functional_chan": tk.StringVar(value="1"),
            "tau": tk.StringVar(value="0.3"),
            "still_thresh": tk.StringVar(value="5.0"),
            "counts_per_cm": tk.StringVar(value=""),
            "use_gpu": tk.BooleanVar(value=False),
            "roi_detection": tk.BooleanVar(value=True),
            "aqua_tifs": tk.BooleanVar(value=True),
            "aqua_rate_hz": tk.StringVar(value="1.0"),
        }
        def field(label, key, width=6):
            ttk.Label(s, text=label).pack(side="left", padx=(10, 2))
            ttk.Entry(s, textvariable=self.vars[key], width=width).pack(side="left")
        field("Functional chan", "functional_chan", 3)
        field("tau (s)", "tau", 5)
        field("still thresh", "still_thresh", 5)
        field("counts/cm", "counts_per_cm", 7)
        ttk.Checkbutton(s, text="GPU", variable=self.vars["use_gpu"]).pack(side="left", padx=8)
        ttk.Checkbutton(s, text="ROI detection", variable=self.vars["roi_detection"]).pack(side="left")
        ttk.Checkbutton(s, text="AQuA tiffs", variable=self.vars["aqua_tifs"]).pack(side="left", padx=8)
        field("AQuA rate (Hz)", "aqua_rate_hz", 5)

        # log
        lf = ttk.LabelFrame(self, text="Log"); lf.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        self.log = tk.Text(lf, height=12, state="disabled", font=("Consolas", 9))
        self.log.pack(fill="both", expand=True)

    # ---- experiment handling
    def add_folder(self):
        folder = filedialog.askdirectory(title="Select an experiment/session folder "
                                               "(it should contain the .sbx runs)")
        if not folder:
            return
        folder = str(Path(folder))
        if folder in self.experiments:
            messagebox.showinfo("afferent2p", "That folder is already in the list.")
            return
        sbx = discover_sbx(Path(folder))
        if not sbx:
            messagebox.showerror("afferent2p", f"No .sbx files found under:\n{folder}")
            return
        missing = [p.name for p in sbx if not p.with_suffix(".mat").exists()]
        if missing:
            messagebox.showwarning("afferent2p",
                                   "These runs are missing their .mat sidecar and "
                                   "will fail:\n" + "\n".join(missing))
        conds = default_conditions(len(sbx))
        self.experiments[folder] = {"runs": [(p.name, c) for p, c in zip(sbx, conds)],
                                    "data_dir": folder, "output_dir": None, "planned": None}
        self.exp_tree.insert("", "end", iid=folder, text=folder,
                             values=(len(sbx), "queued"))
        self.exp_tree.selection_set(folder)
        self.show_runs()

    def add_mouse_planned(self):
        folder = filedialog.askdirectory(title="Select a MOUSE folder (contains the "
                                               "<mouse>_<yymmdd>_<nnn> recordings)")
        if not folder:
            return
        exps = plan_mouse(Path(folder))
        if not exps:
            messagebox.showerror("afferent2p", f"No <mouse>_<yymmdd>_<nnn>.sbx recordings found under:\n{folder}")
            return
        added = 0
        for e in exps:
            key = str(e.output_dir)
            if key in self.experiments:
                continue
            self.experiments[key] = {"runs": list(e.runs), "data_dir": str(e.data_dir),
                                     "output_dir": str(e.output_dir),
                                     "planned": {nm: lab for nm, lab in e.runs}}
            self.exp_tree.insert("", "end", iid=key, text=f"{e.name}  ({e.paradigm})",
                                 values=(len(e.runs), "queued"))
            added += 1
        self._log(f"planned {added} experiment(s) from {folder}:\n{describe_plan(exps)}\n")
        if exps:
            self.exp_tree.selection_set(str(exps[0].output_dir)); self.show_runs()

    def remove_selected(self):
        for iid in self.exp_tree.selection():
            self.experiments.pop(iid, None)
            self.exp_tree.delete(iid)
        self.show_runs()

    def _selected_exp(self) -> str | None:
        sel = self.exp_tree.selection()
        return sel[0] if sel else None

    def show_runs(self):
        self.run_tree.delete(*self.run_tree.get_children())
        exp = self._selected_exp()
        if exp and exp in self.experiments:
            for i, (name, cond) in enumerate(self.experiments[exp]["runs"]):
                self.run_tree.insert("", "end", iid=str(i), text=name, values=(cond,))

    def toggle_condition(self, _event):
        exp = self._selected_exp()
        sel = self.run_tree.selection()
        if not exp or not sel:
            return
        i = int(sel[0])
        name, cond = self.experiments[exp]["runs"][i]
        planned = self.experiments[exp].get("planned")
        if planned:   # planned label <-> skip
            newc = "skip" if cond != "skip" else planned.get(name, "pre")
        else:
            newc = {"pre": "post", "post": "skip", "skip": "pre"}.get(cond, "pre")
        self.experiments[exp]["runs"][i] = (name, newc)
        self.run_tree.item(sel[0], values=(newc,))

    # ---- running
    def run_all(self):
        if self.running:
            messagebox.showinfo("afferent2p", "Already running.")
            return
        if not self.experiments:
            messagebox.showinfo("afferent2p", "Add at least one experiment folder first.")
            return
        for folder, info in self.experiments.items():
            active = [c for _, c in info["runs"] if c != "skip"]
            if not active:
                messagebox.showerror(
                    "afferent2p",
                    f"{Path(folder).name}: every run is marked skip — nothing to do.")
                return
            if len(active) > 1 and len(set(active)) < 2:
                messagebox.showwarning(
                    "afferent2p",
                    f"{Path(folder).name}: several runs but only one condition label — "
                    "they will be pooled (single-condition mode).")
        settings = {
            "functional_chan": self.vars["functional_chan"].get() or "1",
            "tau": self.vars["tau"].get() or "0.7",
            "still_thresh": self.vars["still_thresh"].get() or "1.0",
            "counts_per_cm": (float(self.vars["counts_per_cm"].get())
                              if self.vars["counts_per_cm"].get().strip() else None),
            "use_gpu": self.vars["use_gpu"].get(),
            "roi_detection": self.vars["roi_detection"].get(),
            "aqua_tifs": self.vars["aqua_tifs"].get(),
            "aqua_rate_hz": self.vars["aqua_rate_hz"].get() or "1.0",
        }
        self.running = True
        threading.Thread(target=self._worker, args=(settings,), daemon=True).start()

    def _worker(self, settings):
        for folder, info in list(self.experiments.items()):
            self._set_status(folder, "running")
            cfg = build_config(Path(info.get("data_dir") or folder), info["runs"], settings)
            if info.get("output_dir"):
                cfg["experiment"]["output_dir"] = info["output_dir"]
                cfg["experiment"]["name"] = Path(info["output_dir"]).name
                Path(info["output_dir"]).mkdir(parents=True, exist_ok=True)
                cfg_path = Path(info["output_dir"]) / "afferent2p_config.yaml"
            else:
                cfg_path = Path(folder) / "afferent2p_config.yaml"
            cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
            self.log_q.put(f"\n===== {folder} =====\nconfig: {cfg_path}\n")
            try:
                proc = subprocess.Popen(
                    [sys.executable, str(HERE / "run_pipeline.py"), str(cfg_path)],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                    bufsize=1)
                for line in proc.stdout:
                    self.log_q.put(line)
                proc.wait()
                ok = proc.returncode == 0
            except Exception as e:
                self.log_q.put(f"ERROR: {e}\n")
                ok = False
            self._set_status(folder, "done" if ok else "FAILED")
        self.log_q.put("\nAll experiments finished.\n")
        self.running = False

    def _set_status(self, folder, status):
        def do():
            if self.exp_tree.exists(folder):
                n = self.exp_tree.item(folder, "values")[0]
                self.exp_tree.item(folder, values=(n, status))
        self.after(0, do)

    def _log(self, text: str):
        self.log_q.put(text)

    def _poll_log(self):
        try:
            while True:
                line = self.log_q.get_nowait()
                self.log.configure(state="normal")
                self.log.insert("end", line)
                self.log.see("end")
                self.log.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(150, self._poll_log)


if __name__ == "__main__":
    missing = check_dependencies()
    if missing:
        cmd = f'"{sys.executable}" -m pip install ' + " ".join(missing)
        root = tk.Tk(); root.withdraw()
        messagebox.showerror(
            "afferent2p — missing packages",
            "These Python packages are not installed in the Python running this "
            f"GUI:\n\n  {', '.join(missing)}\n\n"
            "Open a command prompt and run:\n\n" + cmd +
            "\n\n(The command has been printed to the console so you can copy it.)")
        print("\nInstall missing packages with:\n" + cmd + "\n")
        root.destroy()
        sys.exit(1)
    App().mainloop()
