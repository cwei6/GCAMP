r"""Where an animal is in the pipeline, read from the files on disk — one line
per step (1 Unpack … 8 Report) and the next thing to do. The GUI shows it as
the status strip above the tabs (refreshed after every step) and
`python pipeline_status.py <mouse>` prints it, so the 3 → 4 → 5 loop can be
picked up again after a break: nothing has to be remembered.

    st = mouse_state(mouse_dir)
    st["steps"]   [{n, name, tab, status, text, detail}]   status: done | partial | stale | pending | missing | na
    st["next"]    {tab, text}                               the first step that needs doing
    format_state(st) -> str

Rules (all from file times and counts; nothing is computed):
  1 Unpack     run folders with a .sbx vs binned tiffs
  2 Clean      tiffs with the brightness step; decisions files (optional)
  3 Settings   which settings file each experiment of the batch got; the tuning sweep: planned /
               pending jobs / done; the report (tune_recommendation.txt) newer than the sweep's
               last output; applied = aqua_opts files for every experiment
  4 Batch      manifest present; stale when a settings file, a tiff or a run folder is newer than
               the build; queued tiffs
  5 MATLAB     aqua_batch.batch_done(): results for every queued tiff, CFU file + event export
               newer than every result
  6 Curate     runs with a CFU file: curated / stale curation / not curated (optional)
  7 Analyze    aqua_analysis\summary.txt newer than every result, CFU file, curation, decision
               and day_labels.csv; v32 outputs present
  8 Report     <mouse>_report.html / _qc.html newer than summary.txt
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

from .aqua_io import aqua_files, find_binned_tifs, pick_binned_tif, read_curation, relocate, AQUA_SUFFIX

STEP_NAMES = ["Unpack", "Clean & check", "Detection settings", "Batch", "MATLAB", "Curate", "Analyze", "Report"]
RUN_RE = re.compile(r"_(\d{6})_(\d{3})$")


def _parts(path_str: str) -> list[str]:
    """Path components whatever the separator (manifests written on Windows, read anywhere)."""
    return [x for x in re.split(r"[\\/]+", str(path_str)) if x]


def _mtime(p: Path) -> float:
    try:
        return p.stat().st_mtime
    except OSError:
        return 0.0


def _ago(t: float) -> str:
    if not t:
        return "never"
    d = time.time() - t
    if d < 90:
        return "just now"
    if d < 3600:
        return f"{d / 60:.0f} min ago"
    if d < 86400 * 2:
        return f"{d / 3600:.1f} h ago"
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(t))


def _run_folders(mouse_dir: Path) -> list[Path]:
    from .plan import skip_parts
    folders = sorted({p.parent for p in mouse_dir.rglob("*.sbx")} | {p.parent for p in mouse_dir.rglob("*_bin*s.tif")})
    out = []
    for d in folders:
        try:
            rel = d.relative_to(mouse_dir).parts
        except ValueError:
            rel = ()
        if skip_parts(rel):
            continue
        out.append(d)
    return out


def _step(n: int, status: str, text: str, detail: str = "") -> dict:
    return {"n": n, "name": STEP_NAMES[n - 1], "tab": n - 1, "status": status, "text": text, "detail": detail}


def mouse_state(mouse_dir: Path, sec: float | None = 2.0) -> dict:
    mouse_dir = Path(mouse_dir)
    steps: list[dict] = []
    # ---- 1 unpack
    folders = _run_folders(mouse_dir)
    tifs: list[Path] = []
    n_sbx = 0
    missing_tif = []
    for d in folders:
        has_sbx = any(d.glob("*.sbx"))
        n_sbx += int(has_sbx)
        t = pick_binned_tif(d, sec, log=lambda *a: None) if find_binned_tifs(d) else None
        if t is not None:
            tifs.append(t)
        elif has_sbx:
            missing_tif.append(d.name)
    n_runs = max(n_sbx, len(tifs))
    if n_runs == 0:
        steps.append(_step(1, "missing", "no recordings found", "no .sbx and no binned tiff under the mouse folder"))
    elif missing_tif:
        steps.append(_step(1, "partial", f"{len(tifs)} of {n_runs} runs unpacked",
                           "no binned tiff yet: " + ", ".join(RUN_RE.search(n).group(0)[1:] if RUN_RE.search(n) else n for n in missing_tif[:8])
                           + (" …" if len(missing_tif) > 8 else "")))
    else:
        steps.append(_step(1, "done", f"{len(tifs)} runs unpacked", "newest tiff " + _ago(max((_mtime(t) for t in tifs), default=0))))
    newest_tif = max((_mtime(t) for t in tifs), default=0.0)
    # ---- 2 clean & check
    corrected = decided = 0
    for t in tifs:
        jp = t.with_suffix(".json")
        try:
            if json.loads(jp.read_text()).get("global_correction"):
                corrected += 1
        except Exception:
            pass
        if Path(str(t.with_suffix("")) + "_frames_decisions.json").exists():
            decided += 1
    if not tifs:
        steps.append(_step(2, "na", "—", "after unpacking"))
    elif corrected == len(tifs):
        steps.append(_step(2, "done", f"{corrected} of {len(tifs)} flagged" + (f", {decided} with decisions" if decided else ""),
                           "optional: look at the flagged stretches and decide (tab 2)"))
    elif corrected:
        steps.append(_step(2, "partial", f"{corrected} of {len(tifs)} flagged", "run the brightness step on the rest (tab 2) — optional"))
    else:
        steps.append(_step(2, "pending", "not run", "optional: tab 2 flags bad frames; they are excluded from the statistics"))
    # ---- manifest (needed by 3, 4, 5)
    man: dict = {}
    mp = mouse_dir / "aqua_batch" / "manifest.json"
    if mp.exists():
        try:
            man = json.loads(mp.read_text())
        except Exception:
            man = {}
    exps = [e for e in man.get("experiments", []) if e.get("runs")]
    built = float(man.get("built") or (_mtime(mp) if mp.exists() else 0.0))
    # ---- 3 detection settings
    opts_dir = mouse_dir / "aqua_opts"
    tune_dir = mouse_dir / "aqua_tune"
    src_counts = {"tuned": 0, "fov": 0, "run": 0, "default": 0, "other": 0}
    exp_opt_times = {}
    for e in exps:
        of = str(e.get("opt_file") or "")
        parts = _parts(of)
        name = parts[-1].lower() if parts else ""
        parent = parts[-2] if len(parts) > 1 else ""
        p = relocate(of, mouse_dir) if of else None
        if p is not None and p.exists():
            exp_opt_times[e["name"]] = _mtime(p)
        if "aqua_opts" in (x.lower() for x in parts):
            src_counts["tuned"] += 1
        elif name == "default_opt.csv":
            src_counts["default"] += 1
        elif "fov" in name and parent == mouse_dir.name:
            src_counts["fov"] += 1
        elif of:
            src_counts["run"] += 1
        else:
            src_counts["other"] += 1
    tuned_files = {p.stem[:-4] for p in opts_dir.glob("*_Opt.csv")} if opts_dir.is_dir() else set()
    n_tuned = sum(1 for e in exps if e["name"] in tuned_files) if exps else len(tuned_files)
    sweep = None
    sweep_last = 0.0
    if (tune_dir / "tune_manifest.json").exists():
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
            import aqua_tune
            sweep = aqua_tune.sweep_progress(mouse_dir)
        except BaseException:
            sweep = None
        try:
            import csv
            with open(tune_dir / "tune_jobs.csv", newline="") as fh:
                for r in csv.DictReader(fh):
                    sweep_last = max(sweep_last, _mtime(Path(r.get("events_csv", ""))))
        except Exception:
            pass
    rec = tune_dir / "tune_recommendation.txt"
    rec_t = _mtime(rec)
    src_txt = ", ".join(f"{v} {k}" for k, v in (("tuned", src_counts["tuned"]), ("by eye (FOV file)", src_counts["fov"]),
                                                 ("run-folder file", src_counts["run"]), ("default", src_counts["default"])) if v)
    if sweep is not None and sweep.get("total"):
        if sweep["pending"]:
            steps.append(_step(3, "pending", f"sweep {sweep['done']} of {sweep['total']} jobs done",
                               f"{sweep['pending']} pending{', ' + str(sweep['failed']) + ' failed' if sweep['failed'] else ''} — tab 3: paste the MATLAB command "
                               f"(or 'Progress'); settings in the batch now: {src_txt or 'none'}"))
        elif not rec.exists() or rec_t + 1 < sweep_last:
            steps.append(_step(3, "pending", "sweep done, report not run",
                               f"tab 3: '3 Report' (tick apply) — the sweep finished {_ago(sweep_last)}" + (f", the report is from {_ago(rec_t)}" if rec.exists() else "")))
        elif (tune_dir / "grid_suggestion.json").exists() and _mtime(tune_dir / "grid_suggestion.json") >= rec_t - 1:
            # the report found the grid too narrow (v37): the optimum lies beyond it for some movies
            try:
                import json as _json
                sug = _json.loads((tune_dir / "grid_suggestion.json").read_text(encoding="utf-8"))
                n_flag = len(sug.get("movies", {}))
                what = ", ".join(f"{k} {', '.join(str(v) for v in vals)}" for k, vals in sug.get("suggest", {}).items())
            except Exception:
                n_flag, what = 0, ""
            steps.append(_step(3, "partial", f"grid too narrow for {n_flag} movie(s)" + (f" (applied for {n_tuned})" if n_tuned else ""),
                               f"tab 3: 'Extend the grid (as suggested)' ({what}) → sweep → '3 Report'; report {_ago(rec_t)}"))
        elif exps and n_tuned < len(exps):
            steps.append(_step(3, "partial", f"report ready; applied for {n_tuned} of {len(exps)} experiments",
                               f"tab 3: '3 Report' with apply, or keep the current settings ({src_txt}); report {_ago(rec_t)}"))
        else:
            steps.append(_step(3, "done", f"tuned by measurement ({n_tuned} settings files)" if n_tuned else "report ready",
                               f"report {_ago(rec_t)}; settings in the batch: {src_txt or 'none yet (build the batch)'}"))
    elif exps:
        if n_tuned == len(exps):
            steps.append(_step(3, "done", f"tuned ({n_tuned} aqua_opts files)", src_txt))
        else:
            steps.append(_step(3, "done" if (src_counts["fov"] or src_counts["run"]) else "partial",
                               src_txt or "default settings", "not tuned by measurement (optional: tab 3 '1 Plan' → sweep → '3 Report')"))
    else:
        steps.append(_step(3, "na", "—", "build the batch first (tab 4); a FOV without a file gets default_Opt.csv"))
    # ---- 4 batch
    if not mp.exists():
        steps.append(_step(4, "missing", "not built", "tab 4: Build batch files"))
    else:
        why = []
        changed = [e["name"] for e in exps if exp_opt_times.get(e["name"], 0) > built + 1]
        if changed:
            why.append(f"{len(changed)} settings file(s) newer than the batch")
        if newest_tif > built + 1:
            why.append("tiffs newer than the batch (new recordings?)")
        if _mtime(mouse_dir / "run_labels.csv") > built + 1:
            why.append("run_labels.csv newer than the batch")
        in_man = {_parts(r["tif"])[-1] for e in exps for r in e.get("runs", []) if _parts(r["tif"])}
        extra = [t for t in tifs if t.name not in in_man]
        if extra:
            why.append(f"{len(extra)} binned tiff(s) not in the batch")
        try:
            n_queued = sum(1 for ln in (mouse_dir / "aqua_batch" / man.get("batch_csv", "batch.csv")).read_text().splitlines()[1:] if ln.strip())
        except Exception:
            n_queued = 0
        if why:
            steps.append(_step(4, "stale", "rebuild needed", "; ".join(why) + " — tab 4: Build batch files"))
        else:
            steps.append(_step(4, "done", f"{len(exps)} experiments, {n_queued} tiff(s) queued", f"built {_ago(built)}"))
    # ---- 5 MATLAB
    bd = None
    if mp.exists():
        try:
            import aqua_batch
            bd = aqua_batch.batch_done(mouse_dir)
        except BaseException:
            bd = None
    if bd is None:
        steps.append(_step(5, "na", "—", "after the batch"))
    elif bd["done"]:
        steps.append(_step(5, "done", f"{bd['n_results']} results with CFU files", bd["text"]))
    else:
        waiting = bd["pending_aqua"] + bd["pending_cfu"]
        def _short(w: str) -> str:
            m = RUN_RE.search(w.split("_bin")[0])
            return m.group(0)[1:] if m else w
        short = ", ".join(_short(w) for w in waiting[:3]) + (" …" if len(waiting) > 3 else "")
        steps.append(_step(5, "pending", f"AQuA2 {bd['aqua_done']}/{bd['n_queued']}, CFU {bd['cfu_done']}/{bd['n_results']}",
                           f"waiting for {short} — tab 5: paste run_aqua_and_cfu.m (or 'Check what MATLAB has produced')"))
    # ---- 6 curate
    n_cfu = n_cur = n_stale = 0
    newest_result = newest_cur = 0.0
    for t in tifs:
        f = aqua_files(t)
        for k in ("res", "features", "cfu"):
            newest_result = max(newest_result, _mtime(f[k]))
        if f["cfu"].exists():
            n_cfu += 1
            cu = read_curation(t)
            if cu is not None:
                cp = Path(str(t.with_suffix("")) + AQUA_SUFFIX + "_res_cfu_curated.json")
                newest_cur = max(newest_cur, _mtime(cp))
                if cu.get("stale"):
                    n_stale += 1
                else:
                    n_cur += 1
    if n_cfu == 0:
        steps.append(_step(6, "na", "—", "after the CFU batch"))
    elif n_cur == n_cfu:
        steps.append(_step(6, "done", f"{n_cur} of {n_cfu} curated", "decisions in _res_cfu_curated.json"))
    else:
        steps.append(_step(6, "partial" if n_cur else "pending", f"{n_cur} of {n_cfu} curated" + (f", {n_stale} stale" if n_stale else ""),
                           "optional but recommended: tab 6 — mark non-fibre CFUs noise ('Select suggested'), then tabs 7 and 8 again"))
    # ---- 7 analyze
    A = mouse_dir / "aqua_analysis"
    summ = A / "summary.txt"
    s_t = _mtime(summ)
    inputs = max(newest_result, newest_cur, _mtime(mouse_dir / "day_labels.csv"),
                 max((_mtime(Path(str(t.with_suffix("")) + "_frames_decisions.json")) for t in tifs), default=0.0))
    if not summ.exists():
        steps.append(_step(7, "missing", "not run", "tab 7: Run analysis"))
    elif s_t + 1 < inputs:
        steps.append(_step(7, "stale", f"older than the inputs ({_ago(s_t)})", "a result, CFU file, curation, decision or day label changed since — tab 7: Run analysis"))
    elif not (A / "loco_summary.csv").exists():
        steps.append(_step(7, "stale", f"from before v32 ({_ago(s_t)})", "tab 7: Run analysis for the general representation, heatmaps and locomotion trials"))
    else:
        steps.append(_step(7, "done", f"analysed {_ago(s_t)}", "summary.txt and the csvs/figures in aqua_analysis\\"))
    # ---- 8 report
    rep = A / f"{mouse_dir.name}_report.html"
    qc = A / f"{mouse_dir.name}_qc.html"
    r_t, q_t = _mtime(rep), _mtime(qc)
    if not summ.exists():
        steps.append(_step(8, "na", "—", "after the analysis"))
    elif not rep.exists():
        steps.append(_step(8, "missing", "not built", "tab 8: Build the report (and the QC report)"))
    elif r_t + 1 < s_t:
        steps.append(_step(8, "stale", f"older than the analysis ({_ago(r_t)})", "tab 8: Build the report" + (" and the QC report" if not qc.exists() or q_t + 1 < s_t else "")))
    elif not qc.exists():
        steps.append(_step(8, "partial", f"report {_ago(r_t)}, no QC report", "tab 8: Build the QC report"))
    elif q_t + 1 < s_t:
        steps.append(_step(8, "partial", f"report {_ago(r_t)}, QC report older than the analysis", "tab 8: Build the QC report"))
    else:
        steps.append(_step(8, "done", f"report {_ago(r_t)}, QC {_ago(q_t)}", str(rep)))
    # ---- next action: the first blocking step
    nxt = None
    order = [1, 3, 4, 5, 7, 8]                  # 2 and 6 are optional
    for n in order:
        st = steps[n - 1]
        if n == 3 and st["status"] in ("partial", "done", "na"):
            continue                              # "not applied" / "by eye" is advisory
        if st["status"] in ("missing", "stale", "pending", "partial"):
            nxt = {"tab": st["tab"], "text": f"{n} {st['name']}: {st['text']} — {st['detail']}"}
            break
    if nxt is None:
        cur = steps[5]
        if cur["status"] in ("pending", "partial"):
            nxt = {"tab": 5, "text": f"everything is up to date — optional: 6 Curate ({cur['text']}), then tabs 7 and 8 again"}
        else:
            nxt = {"tab": 7, "text": "everything is up to date — open the report (tab 8)"}
    return {"mouse": mouse_dir.name, "steps": steps, "next": nxt, "checked": time.time()}


def format_state(st: dict) -> str:
    mark = {"done": "✓", "partial": "◐", "stale": "!", "pending": "…", "missing": "○", "na": "·"}
    lines = [f"{st['mouse']} — pipeline state ({time.strftime('%H:%M', time.localtime(st['checked']))}):"]
    for s in st["steps"]:
        lines.append(f"  {mark.get(s['status'], '?')} {s['n']} {s['name']:<18} {s['text']}" + (f"   ({s['detail']})" if s["detail"] else ""))
    lines.append(f"  → next: {st['next']['text']}")
    return "\n".join(lines)
