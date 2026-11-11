#!/usr/bin/env python3
"""aqua_prep — unpack + register Scanbox recordings and write AQuA2-ready tiffs.

Usage:  python aqua_prep.py

Pick one or more folders. Two kinds of folder are accepted:
  * a SESSION folder that contains one recording (.sbx + .mat [+ quad]), e.g.
        D:\\2photon\\Chao\\NavTIGRETreg-2608-1\\NavTIGRETreg-2608-1_260911_001
  * a MOUSE folder that contains several session folders, e.g.
        D:\\2photon\\Chao\\NavTIGRETreg-2608-1
    -> every session folder inside it is added.

For each session: sbx -> tiff, suite2p registration (rigid + optional nonrigid,
temporally smoothed alignment for dim axon data), then the registered movie is
binned to the AQuA rate (default 1 Hz, mean of round(fs/rate) frames — the
same rule as the lab's MATLAB GenerateExptProjections), edge-cropped, and
written as uint16 tif + a .json sidecar with the frame rate. All AQuA tiffs are
also collected into one folder (<first mouse folder>\\AQuA_input by default) so
AQuA2 batch mode can be pointed at a single directory.

"Register together" mode: all selected sessions are treated as runs of the SAME
field of view (e.g. pre / post optostim) and registered to one common reference;
the first is labeled pre, the rest post (edit the generated yaml to change).
"""
from __future__ import annotations
import queue
import shutil
import subprocess
import sys
import threading
from pathlib import Path

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import yaml
from gui import check_dependencies, discover_sbx  # noqa: E402

STEPS_SINGLE = "preview,convert,register,export,loco,qc"


def classify_folder(folder: Path) -> list[Path]:
    """Return the list of session folders represented by `folder`."""
    direct = [p for p in folder.glob("*.sbx")]
    if direct:
        return [folder]
    sessions = sorted(
        {p.parent for p in folder.rglob("*.sbx")}, key=lambda p: p.name)
    return sessions


def make_config(data_dir: Path, runs: list[tuple[str, str]], opt: dict) -> dict:
    crop = [int(v) for v in opt["crop"]]
    return {
        "experiment": {
            "name": data_dir.name,
            "data_dir": str(data_dir),
            "runs": [{"sbx": n, "condition": c} for n, c in runs],
            "output_dir": None,
        },
        "imaging": {"functional_chan": int(opt["chan"]), "register_on_red": "auto", "nplanes": 1,
                    "fs": None, "tau": 0.3},
        "prealign": {"enabled": True, "max_angle_deg": 6, "min_shift_px": 2, "min_angle_deg": 0.5},
        "registration": {
            "nonrigid": not bool(opt["rigid_only"]),
            "block_size": [128, 128],
            "maxregshift": 0.1,
            "batch_size": 200,
            "smooth_sigma_time": float(opt["smooth_t"]),
            "nimg_init": 1000,
            "torch_device": "cuda" if opt["gpu"] else "cpu",
        },
        "roi_detection": {"enabled": False},      # AQuA route: no suite2p ROIs
        "export": {
            "aqua_tifs": True, "chunk_frames": 4000, "aqua_dtype": "uint16",
            "aqua_rate_hz": float(opt["rate"]), "aqua_bin_t": 0,
            "aqua_crop": crop,
        },
        "loco": {"enabled": True, "quad_glob": "*quad*", "quad_kind": "cumulative",
                 "counts_per_cm": None, "still_thresh": 5.0, "min_still_s": 5.0,
                 "smooth_s": 0.5},
        "pupil": {"enabled": True, "glob": "*eye*DLC*.csv", "likelihood_min": 0.6,
                  "min_edge_points": 3, "max_gap_s": 0.5, "smooth_s": 0.25},
        "analysis": {"enabled": False},
    }


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("aqua_prep — unpack + register + 1 Hz tiff for AQuA2")
        self.geometry("960x640")
        self.sessions: list[Path] = []
        self.log_q: queue.Queue[str] = queue.Queue()
        self.running = False
        self._build()
        self.after(100, self._poll_log)

    def _build(self):
        top = ttk.Frame(self); top.pack(fill="x", padx=8, pady=6)
        ttk.Button(top, text="Add folder (session or mouse)…",
                   command=self.add_folder).pack(side="left")
        ttk.Button(top, text="Remove selected", command=self.remove_selected).pack(side="left", padx=6)
        ttk.Button(top, text="Clear", command=self.clear).pack(side="left")
        ttk.Button(top, text="Run", command=self.run).pack(side="right")

        lf = ttk.LabelFrame(self, text="Sessions (each = one recording folder)")
        lf.pack(fill="both", expand=True, padx=8)
        self.tree = ttk.Treeview(lf, columns=("sbx", "status"), height=9)
        self.tree.heading("#0", text="Session folder")
        self.tree.heading("sbx", text=".sbx")
        self.tree.heading("status", text="Status")
        self.tree.column("#0", width=520); self.tree.column("sbx", width=260)
        self.tree.column("status", width=90, anchor="center")
        self.tree.pack(fill="both", expand=True)

        s = ttk.LabelFrame(self, text="Options"); s.pack(fill="x", padx=8, pady=6)
        self.v = {
            "rate": tk.StringVar(value="1.0"),
            "crop": tk.StringVar(value="80,80,40,40"),
            "chan": tk.StringVar(value="1"),
            "smooth_t": tk.StringVar(value="2"),
            "rigid_only": tk.BooleanVar(value=False),
            "gpu": tk.BooleanVar(value=False),
            "together": tk.BooleanVar(value=False),
            "collect": tk.StringVar(value=""),
        }
        def field(label, key, width=6):
            ttk.Label(s, text=label).pack(side="left", padx=(10, 2))
            ttk.Entry(s, textvariable=self.v[key], width=width).pack(side="left")
        field("AQuA rate (Hz)", "rate", 5)
        field("crop L,R,T,B", "crop", 12)
        field("chan", "chan", 3)
        field("align smooth (frames)", "smooth_t", 3)
        ttk.Checkbutton(s, text="rigid only (MATLAB-like)", variable=self.v["rigid_only"]).pack(side="left", padx=6)
        ttk.Checkbutton(s, text="GPU", variable=self.v["gpu"]).pack(side="left")
        ttk.Checkbutton(s, text="register all together (same FOV: pre/post)",
                        variable=self.v["together"]).pack(side="left", padx=6)

        c = ttk.Frame(self); c.pack(fill="x", padx=8)
        ttk.Label(c, text="Collect AQuA tiffs into:").pack(side="left")
        ttk.Entry(c, textvariable=self.v["collect"], width=70).pack(side="left", padx=4)
        ttk.Button(c, text="…", command=self.pick_collect).pack(side="left")
        ttk.Label(c, text="(blank = <mouse folder>\\AQuA_input)").pack(side="left", padx=6)

        lg = ttk.LabelFrame(self, text="Log"); lg.pack(fill="both", expand=True, padx=8, pady=(6, 8))
        self.log = tk.Text(lg, height=10, state="disabled", font=("Consolas", 9))
        self.log.pack(fill="both", expand=True)

    # ---- folder handling
    def add_folder(self):
        folder = filedialog.askdirectory(title="Select a session folder (one recording) "
                                               "or a mouse folder (many sessions)")
        if not folder:
            return
        sessions = classify_folder(Path(folder))
        if not sessions:
            messagebox.showerror("aqua_prep", f"No .sbx files found under:\n{folder}")
            return
        added = 0
        for sfold in sessions:
            if sfold in self.sessions:
                continue
            sbx = discover_sbx(sfold)
            if len(sbx) != 1:
                messagebox.showwarning(
                    "aqua_prep", f"{sfold.name} contains {len(sbx)} .sbx files; "
                    "expected 1 per session folder. Using the first: "
                    f"{sbx[0].name}")
            if not sbx[0].with_suffix(".mat").exists():
                messagebox.showwarning("aqua_prep", f"{sbx[0].name} is missing its .mat "
                                                    "sidecar and will fail.")
            self.sessions.append(sfold)
            self.tree.insert("", "end", iid=str(sfold), text=str(sfold),
                             values=(sbx[0].name, "queued"))
            added += 1
        self._log(f"added {added} session(s) from {folder}\n")

    def remove_selected(self):
        for iid in self.tree.selection():
            self.sessions = [s for s in self.sessions if str(s) != iid]
            self.tree.delete(iid)

    def clear(self):
        self.sessions.clear()
        self.tree.delete(*self.tree.get_children())

    def pick_collect(self):
        d = filedialog.askdirectory(title="Folder to collect all AQuA tiffs into")
        if d:
            self.v["collect"].set(d)

    # ---- run
    def run(self):
        if self.running:
            messagebox.showinfo("aqua_prep", "Already running."); return
        if not self.sessions:
            messagebox.showinfo("aqua_prep", "Add at least one folder first."); return
        try:
            crop = [int(x) for x in self.v["crop"].get().split(",")]
            assert len(crop) == 4
        except Exception:
            messagebox.showerror("aqua_prep", "crop must be four integers: L,R,T,B"); return
        opt = {
            "rate": self.v["rate"].get() or "1.0", "crop": crop,
            "chan": self.v["chan"].get() or "1",
            "smooth_t": self.v["smooth_t"].get() or "2",
            "rigid_only": self.v["rigid_only"].get(), "gpu": self.v["gpu"].get(),
        }
        together = self.v["together"].get()
        collect = self.v["collect"].get().strip()
        collect_dir = Path(collect) if collect else self.sessions[0].parent / "AQuA_input"
        self.running = True
        threading.Thread(target=self._worker, args=(opt, together, collect_dir),
                         daemon=True).start()

    def _run_config(self, cfg: dict, cfg_path: Path) -> bool:
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        self.log_q.put(f"\n===== {cfg['experiment']['data_dir']} =====\nconfig: {cfg_path}\n")
        proc = subprocess.Popen(
            [sys.executable, str(HERE / "run_pipeline.py"), str(cfg_path),
             "--steps", STEPS_SINGLE],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            self.log_q.put(line)
        proc.wait()
        return proc.returncode == 0

    def _collect(self, out_dir: Path, collect_dir: Path):
        collect_dir.mkdir(parents=True, exist_ok=True)
        aq = out_dir / "aqua_input"
        n = 0
        for f in sorted(aq.glob("*_reg.tif")) + sorted(aq.glob("*_reg.json")):
            dst = collect_dir / f.name
            if not dst.exists() or dst.stat().st_mtime < f.stat().st_mtime:
                shutil.copy2(f, dst)
                n += 1
        self.log_q.put(f"collected {n} file(s) into {collect_dir}\n")

    def _worker(self, opt, together, collect_dir):
        try:
            if together:
                parent = self.sessions[0].parent
                runs = []
                for i, sfold in enumerate(self.sessions):
                    sbx = discover_sbx(sfold)[0]
                    runs.append((sbx.name, "pre" if i == 0 else "post"))
                cfg = make_config(parent, runs, opt)
                cfg["experiment"]["name"] = parent.name + "_together"
                cfg["experiment"]["output_dir"] = str(parent / "afferent2p_out_together")
                for sfold in self.sessions:
                    self._status(sfold, "running")
                ok = self._run_config(cfg, parent / "aqua_prep_together.yaml")
                for sfold in self.sessions:
                    self._status(sfold, "done" if ok else "FAILED")
                if ok:
                    self._collect(Path(cfg["experiment"]["output_dir"]), collect_dir)
            else:
                for sfold in self.sessions:
                    self._status(sfold, "running")
                    sbx = discover_sbx(sfold)[0]
                    cfg = make_config(sfold, [(sbx.name, "pre")], opt)
                    ok = self._run_config(cfg, sfold / "aqua_prep_config.yaml")
                    self._status(sfold, "done" if ok else "FAILED")
                    if ok:
                        self._collect(sfold / "afferent2p_out", collect_dir)
        except Exception as e:
            self.log_q.put(f"ERROR: {e}\n")
        self.log_q.put(f"\nAll done. AQuA tiffs collected in: {collect_dir}\n")
        self.running = False

    def _status(self, sfold: Path, status: str):
        def do():
            iid = str(sfold)
            if self.tree.exists(iid):
                self.tree.item(iid, values=(self.tree.item(iid, "values")[0], status))
        self.after(0, do)

    def _log(self, text):
        self.log_q.put(text)

    def _poll_log(self):
        try:
            while True:
                line = self.log_q.get_nowait()
                self.log.configure(state="normal")
                self.log.insert("end", line); self.log.see("end")
                self.log.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(150, self._poll_log)


if __name__ == "__main__":
    missing = check_dependencies()
    if missing:
        cmd = f'"{sys.executable}" -m pip install ' + " ".join(missing)
        root = tk.Tk(); root.withdraw()
        messagebox.showerror("aqua_prep — missing packages",
                             f"Not installed in this Python:\n  {', '.join(missing)}\n\n"
                             "Run in a command prompt:\n\n" + cmd)
        print("\nInstall with:\n" + cmd + "\n"); root.destroy(); sys.exit(1)
    App().mainloop()
