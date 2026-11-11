"""Pick a detection setting by eye (v38).

The sweep measures every setting; the rule picks one. This window shows the few
settings worth a look — the rule's pick, the more sensitive one it rejected,
the setting in the batch now, the next cleaner one — side by side ON THE REAL
MOVIE of a tuned run: every event a setting found is a circle (green on the
fibres, orange off them), the events the same setting found in the scrambled
movie (noise) can be overlaid in red, and a frame slider plays the movie so an
event can be seen happening (an active event is filled). Zoom and pan are the
matplotlib toolbar's; the panels zoom together. "Pick this" writes the setting
as the experiment's Opt file (the whole day, or this condition only) and the
tuning json, so Step 4 queues the runs again.

Needs: aqua_tune\\tune_report.csv (3 Report), the tuned movie and its sweep
results on disk. Opened from tab 3 ("Pick by eye…") or:
    python -c "from afferent2p.pick import main; main()" <mouse folder>
"""
from __future__ import annotations

import json
import sys
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

import numpy as np

from afferent2p import tune as T

COND_COL = {"pre": "#3b7dd8", "post": "#e8702a"}


class PickWindow(tk.Toplevel):
    def __init__(self, master, mouse_dir: Path, log=print, max_false: float = 0.10, min_on_fibre: float = 0.80,
                 tol_recall: float = 0.05, min_recall: float = 0.3, run: str | None = None):
        super().__init__(master)
        import aqua_tune as AT
        self.AT = AT
        self.mouse_dir = Path(mouse_dir); self.log = log
        self.max_false, self.min_on_fibre, self.tol_recall, self.min_recall = max_false, min_on_fibre, tol_recall, min_recall
        tune_dir = self.mouse_dir / "aqua_tune"
        self.man = json.loads((tune_dir / "tune_manifest.json").read_text())
        self.grid_keys = self.man["grid_keys"]
        self.per_movie = AT.rows_from_report(tune_dir / "tune_report.csv", self.grid_keys)
        self.day_labels = AT.read_day_labels(self.mouse_dir)
        self.movies = [mv for mv in self.man["movies"] if mv["run"] in self.per_movie]
        self.by_run = {mv["run"]: mv for mv in self.movies}
        # the rule's choices, per movie and per day (as the report makes them)
        self.chosen, self.passed = {}, {}
        for r, rows in self.per_movie.items():
            best, why, passed = T.recommend(rows, max_false=max_false, min_on_fibre=min_on_fibre)
            self.chosen[r] = best; self.passed[r] = passed
        self.day_matched = {}
        exps: dict = {}
        for mv in self.movies:
            if mv.get("paradigm") != "spontaneous":
                exps.setdefault(mv["experiment"], []).append(mv["run"])
        for exp, runs in exps.items():
            if len(runs) >= 2:
                conds = {r: (self.by_run[r].get("condition") or AT.short_run(r)) for r in runs}
                res = T.recommend_matched({r: self.per_movie[r] for r in runs}, self.grid_keys, max_false=max_false,
                                          min_on_fibre=min_on_fibre, tol_recall=tol_recall, min_level=min_recall)
                self.day_matched[exp] = res + (conds,)
        self.ctx = None; self.frames = None; self.cands = []; self.axes = []; self.images = []; self.patches = []
        self.playing = False; self.lo = self.hi = None
        self.title(f"Pick a setting by eye — {self.mouse_dir.name}")
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry(f"{min(1500, sw - 60)}x{min(940, sh - 100)}+30+30")
        self._build()
        names = [self._movie_name(mv) for mv in self.movies]
        self.pick_movie["values"] = names
        start = 0
        if run:
            for i, mv in enumerate(self.movies):
                if run in mv["run"]:
                    start = i; break
        if names:
            self.pick_movie.current(start)
            self.after(100, lambda: self.load_movie(start))

    # ------------------------------------------------------------------ UI
    def _movie_name(self, mv: dict) -> str:
        dl = self.day_labels.get(self.AT.day_of(mv["run"]) or "", "")
        return f"{self.AT.short_run(mv['run'])}  FOV {mv['fov']}  [{mv.get('condition') or '-'}]" + (f"  {dl}" if dl else "")

    def _build(self):
        top = ttk.Frame(self); top.pack(fill="x", padx=8, pady=(8, 2))
        ttk.Label(top, text="movie:").pack(side="left")
        self.pick_movie = ttk.Combobox(top, state="readonly", width=52); self.pick_movie.pack(side="left", padx=4)
        self.pick_movie.bind("<<ComboboxSelected>>", lambda e: self.load_movie(self.pick_movie.current()))
        self.show_noise = tk.BooleanVar(value=False)
        ttk.Checkbutton(top, text="show the noise events (found in the scrambled movie, red)", variable=self.show_noise,
                        command=self.redraw_patches).pack(side="left", padx=(16, 4))
        self.show_ref = tk.BooleanVar(value=True)
        self.ref_check = ttk.Checkbutton(top, text="hand result (white dashed)", variable=self.show_ref, command=self.redraw_patches)
        self.ref_check.pack(side="left", padx=4)
        ttk.Label(top, text="   apply the pick to:").pack(side="left")
        self.scope = tk.StringVar(value="day")
        self.scope_day = ttk.Radiobutton(top, text="the whole day (pre and post)", variable=self.scope, value="day"); self.scope_day.pack(side="left")
        self.scope_cond = ttk.Radiobutton(top, text="this condition only", variable=self.scope, value="condition"); self.scope_cond.pack(side="left")
        row2 = ttk.Frame(self); row2.pack(fill="x", padx=8, pady=(0, 2))
        ttk.Label(row2, text="frame:").pack(side="left")
        self.frame_var = tk.IntVar(value=0)
        self.scale = tk.Scale(row2, from_=0, to=1, orient="horizontal", variable=self.frame_var, showvalue=True, length=520,
                              command=lambda v: self.set_frame(int(float(v))))
        self.scale.pack(side="left", padx=4)
        self.play_btn = ttk.Button(row2, text="▶ Play", command=self.toggle_play); self.play_btn.pack(side="left", padx=4)
        ttk.Label(row2, text="speed").pack(side="left", padx=(8, 2))
        self.speed = ttk.Combobox(row2, values=["2 fps", "5 fps", "10 fps"], state="readonly", width=7); self.speed.current(1); self.speed.pack(side="left")
        ttk.Button(row2, text="mean image", command=lambda: self.set_frame(0)).pack(side="left", padx=(12, 2))
        self.time_lbl = ttk.Label(row2, text=""); self.time_lbl.pack(side="left", padx=8)
        # the pick buttons and the status line sit at the bottom whatever the window height
        self.status = ttk.Label(self, text="", foreground="#333", wraplength=1400, justify="left")
        self.status.pack(side="bottom", anchor="w", padx=10, pady=(0, 8))
        self.btn_row = ttk.Frame(self); self.btn_row.pack(side="bottom", fill="x", padx=8, pady=(2, 2))
        body = ttk.Frame(self); body.pack(fill="both", expand=True, padx=8, pady=2)
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        self.fig = Figure(figsize=(14, 6.2), dpi=96)
        self.canvas = FigureCanvasTkAgg(self.fig, master=body)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        tb = ttk.Frame(body); tb.pack(fill="x")
        self.toolbar = NavigationToolbar2Tk(self.canvas, tb, pack_toolbar=False); self.toolbar.pack(side="left")
        ttk.Label(tb, text="   toolbar: zoom (magnifier) and pan (cross) act on every panel at once; space = play / pause", foreground="#666").pack(side="left")
        self.bind("<space>", lambda e: self.toggle_play())

    # ------------------------------------------------------------------ data
    def load_movie(self, index: int):
        if not self.movies:
            return
        self.playing = False; self.play_btn.configure(text="▶ Play")
        mv = self.movies[index]
        self.status.configure(text=f"loading {mv['run']} …"); self.update_idletasks()
        try:
            self.ctx = self.AT.movie_context(self.mouse_dir, self.man, mv, log=self.log)
        except Exception as e:
            self.status.configure(text=f"could not load this movie: {e}"); return
        try:
            self.frames = T._open(self.ctx["tif"])
        except Exception as e:
            self.frames = None; self.log(f"(movie not readable for frames: {e})")
        rm = self.ctx["raw_mean"]
        self.lo, self.hi = np.percentile(rm, [1, 99.7])
        n_fr = int(self.frames.shape[0]) if self.frames is not None else 1
        self.scale.configure(to=max(1, n_fr)); self.frame_var.set(0)
        rows = self.per_movie[mv["run"]]
        fk = T.false_key(rows)
        self.fk = fk
        self.cands = T.candidates(rows, self.chosen.get(mv["run"]), self.grid_keys, fk, self.max_false, self.min_on_fibre)
        # the day's matched choice for this movie, if it is not among the candidates yet
        dm = self.day_matched.get(mv["experiment"])
        if dm and dm[3] and dm[0] and mv["run"] in dm[0]:
            key = tuple(dm[0][mv["run"]])
            if not any(tuple(r.get(k) for k in self.grid_keys) == key for r, _, _ in self.cands):
                row = next((r for r in rows if tuple(r.get(k) for k in self.grid_keys) == key), None)
                if row is not None:
                    self.cands.insert(1, (row, "what the day uses (matched with the other condition)", "#2ecc40"))
            self.cands = self.cands[:4]
        self.scope_cond.configure(text=f"this condition only ({mv.get('condition') or 'this movie'})")
        self.ref_check.state(["!disabled"] if self.ctx.get("reference") is not None else ["disabled"])
        self.draw()
        self.set_frame(0)
        self.status.configure(text=f"{len(self.cands)} setting(s) worth a look — zoom with the magnifier, play the movie, tick the noise events; "
                                   "'Pick …' writes that setting for this experiment (then tab 4 Batch → MATLAB → 7 → 8)")

    def other_movies(self, mv: dict) -> list:
        return [m for m in self.movies if m["experiment"] == mv["experiment"] and m["run"] != mv["run"]]

    def row_at(self, run: str, key: tuple):
        return next((r for r in self.per_movie.get(run, []) if tuple(r.get(k) for k in self.grid_keys) == tuple(key)), None)

    # ------------------------------------------------------------------ drawing
    def draw(self):
        from matplotlib.patches import Circle
        self.fig.clear()
        for w in self.btn_row.winfo_children():
            w.destroy()
        n = max(1, len(self.cands))
        self.axes = []; self.images = []; self.patches = []
        mv = self.ctx["mv"]
        bg = T._bg_image(self.ctx["raw_mean"])
        off = self.ctx["offsets"].get("real") or (0, 0, None)
        off_null = self.ctx["offsets"].get("null") or off
        fm = self.ctx["fmask_raw"]
        ax0 = None
        for i, (row, tag, col) in enumerate(self.cands):
            ax = self.fig.add_subplot(1, n, i + 1, sharex=ax0, sharey=ax0)
            ax0 = ax0 or ax
            im = ax.imshow(bg, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
            try:
                ax.contour(fm.astype(float), levels=[0.5], colors=["#00d0d0"], linewidths=0.3, alpha=0.3)
            except Exception:
                pass
            det = self.ctx["read_events"](row["setting"], "real")
            null = self.ctx["read_events"](row["setting"], "null")
            pl = {"real": [], "null": [], "ref": []}
            if det is not None and len(det):
                on = T._on_mask(det, fm, off[:2])
                for (_, e), o in zip(det.iterrows(), on):
                    y = float(e["y_px"]) - 1 + off[0]; x = float(e["x_px"]) - 1 + off[1]
                    r_ = max(3.0, float(np.sqrt(max(float(e.get("n_pixels", 30)), 1) / np.pi)))
                    c = Circle((x, y), r_, fill=False, lw=0.9, color="#2ecc40" if o else "#ff9f1a", alpha=0.9)
                    ax.add_patch(c); pl["real"].append((c, int(e["t0_frame"]), int(e["t1_frame"])))
            if null is not None and len(null):
                for _, e in null.iterrows():
                    y = float(e["y_px"]) - 1 + off_null[0]; x = float(e["x_px"]) - 1 + off_null[1]
                    r_ = max(3.0, float(np.sqrt(max(float(e.get("n_pixels", 30)), 1) / np.pi)))
                    c = Circle((x, y), r_, fill=False, lw=0.9, color="#d62728", alpha=0.9, visible=False)
                    ax.add_patch(c); pl["null"].append((c, int(e["t0_frame"]), int(e["t1_frame"])))
            ref = self.ctx.get("reference")
            if ref is not None and len(ref):
                ro = self.ctx.get("ref_off") or off
                for _, e in ref.iterrows():
                    y = float(e["y_px"]) - 1 + ro[0]; x = float(e["x_px"]) - 1 + ro[1]
                    r_ = max(3.0, float(np.sqrt(max(float(e.get("n_pixels", 30)), 1) / np.pi))) + 1.5
                    c = Circle((x, y), r_, fill=False, lw=0.6, color="white", alpha=0.8, linestyle=(0, (2, 2)))
                    ax.add_patch(c); pl["ref"].append((c, int(e["t0_frame"]), int(e["t1_frame"])))
            n_on = int(row.get("n_real_on", row.get("n_real", 0)) or 0); rate = row.get("rate_real_on", row.get("rate_real", np.nan))
            title = (f"{T.setting_words(row, self.grid_keys)}\n{tag}\n"
                     f"{rate:.2f} events/min on fibres ({n_on} in {float(row.get('minutes', np.nan)):.0f} min), noise {T.noise_pct(row, self.fk):.0f} %, "
                     f"catch {100 * row.get('recall', np.nan):.0f} %, {row.get('on_fibre', np.nan):.0%} on fibres")
            others = []
            for om in self.other_movies(mv):
                orow = self.row_at(om["run"], tuple(row.get(k) for k in self.grid_keys))
                if orow is not None:
                    others.append(f"{om.get('condition') or self.AT.short_run(om['run'])} at this setting: "
                                  f"{orow.get('rate_real_on', orow.get('rate_real', np.nan)):.2f}/min, noise {T.noise_pct(orow, self.fk):.0f} %, "
                                  f"catch {100 * orow.get('recall', np.nan):.0f} %")
            if others:
                title += "\n" + "; ".join(others)
            ax.set_title(title, fontsize=8, color=col, loc="left")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color(col); sp.set_linewidth(2.0)
            self.axes.append(ax); self.images.append(im); self.patches.append(pl)
            b = ttk.Button(self.btn_row, text=f"Pick {T.setting_short(row, self.grid_keys)}: {T.setting_words(row, self.grid_keys)}",
                           command=lambda i=i: self.pick(i))
            b.pack(side="left", padx=6, expand=True)
        dl = self.day_labels.get(self.AT.day_of(mv["run"]) or "", "")
        self.fig.suptitle(f"{self.AT.short_run(mv['run'])}  FOV {mv['fov']}  [{mv.get('condition') or '-'}]" + (f"  ({dl})" if dl else "")
                          + " — the settings worth a look, each with the events it finds on the real movie (green on fibres, orange off; "
                          "an active event is filled while the movie plays)", fontsize=9, x=0.01, ha="left")
        self.fig.tight_layout(rect=(0, 0, 1, 0.95))
        self.redraw_patches()
        self.canvas.draw_idle()

    def redraw_patches(self):
        show_n = bool(self.show_noise.get()); show_r = bool(self.show_ref.get())
        for pl in self.patches:
            for c, _, _ in pl["null"]:
                c.set_visible(show_n)
            for c, _, _ in pl["ref"]:
                c.set_visible(show_r)
        self.canvas.draw_idle()

    def set_frame(self, t: int):
        """t = 0: the mean image; t >= 1: that frame (1-based, as AQuA2 counts) with the active events filled."""
        if not self.images:
            return
        if t <= 0 or self.frames is None:
            img = T._bg_image(self.ctx["raw_mean"]); self.time_lbl.configure(text="mean image")
        else:
            t = min(t, int(self.frames.shape[0]))
            fr = np.asarray(self.frames[t - 1], np.float32)
            img = np.clip((fr - self.lo) / max(self.hi - self.lo, 1e-6), 0, 1) ** 0.6
            spf = float(self.ctx["mv"].get("spf", 2.0))
            self.time_lbl.configure(text=f"frame {t} / {int(self.frames.shape[0])}  ({(t - 1) * spf / 60:.1f} min)")
        for im in self.images:
            im.set_data(img)
        for pl in self.patches:
            for kind in ("real", "null"):
                for c, t0, t1 in pl[kind]:
                    active = t > 0 and t0 <= t <= t1
                    c.set_fill(active); c.set_alpha(0.45 if active else 0.9); c.set_linewidth(2.0 if active else 0.9)
        self.canvas.draw_idle()

    def toggle_play(self):
        self.playing = not self.playing
        self.play_btn.configure(text="⏸ Pause" if self.playing else "▶ Play")
        if self.playing:
            self._tick()

    def _tick(self):
        if not self.playing or self.frames is None:
            return
        t = int(self.frame_var.get()) + 1
        if t > int(self.frames.shape[0]):
            t = 1
        self.frame_var.set(t); self.set_frame(t)
        fps = {"2 fps": 2, "5 fps": 5, "10 fps": 10}.get(self.speed.get(), 5)
        self.after(int(1000 / fps), self._tick)

    # ------------------------------------------------------------------ the pick
    def pick(self, i: int):
        row, tag, col = self.cands[i]
        mv = self.ctx["mv"]
        values = {k: row[k] for k in self.grid_keys}
        scope = self.scope.get()
        if mv.get("paradigm") == "spontaneous":
            scope = "day"
        cond = mv.get("condition")
        words = T.setting_words(row, self.grid_keys)
        what = ("every run of the day (pre and post)" if scope == "day" and mv.get("paradigm") != "spontaneous"
                else (f"the {cond} runs of this day" if scope == "condition" else "every day of this FOV"))
        if not messagebox.askyesno("Pick this setting", f"Write {words} for {what} of {mv['experiment']}?\n\n"
                                                        f"({tag}: {row.get('rate_real_on', row.get('rate_real', float('nan'))):.2f} events/min on fibres, "
                                                        f"noise {T.noise_pct(row, self.fk):.0f} %, catch {100 * row.get('recall', float('nan')):.0f} %)\n\n"
                                                        "Step 4 (Batch) then queues these runs again with the new setting.", parent=self):
            return
        try:
            paths = self.AT.apply_pick(self.mouse_dir, mv["experiment"], values, scope, cond, tag=tag, log=self.log)
        except Exception as e:
            messagebox.showerror("Pick", f"could not write the setting: {e}", parent=self); return
        self.status.configure(text=f"written: {', '.join(p.name for p in paths)} — next: tab 4 (Batch), then the MATLAB batch, 7 and 8")

    # ------------------------------------------------------------------ cleanup
    def destroy(self):
        self.playing = False
        try:
            self.frames = None
        except Exception:
            pass
        super().destroy()


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    root = tk.Tk(); root.withdraw()
    w = PickWindow(root, Path(sys.argv[1]), run=sys.argv[2] if len(sys.argv) > 2 else None)
    w.protocol("WM_DELETE_WINDOW", root.destroy)
    root.mainloop()


if __name__ == "__main__":
    main()
