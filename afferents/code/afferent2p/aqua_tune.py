#!/usr/bin/env python3
r"""aqua_tune — find AQuA2 detection settings by measurement, unattended.

    python aqua_tune.py <mouse folder> --plan                       # one movie per day and FOV (optostim: pre and post)
    python aqua_tune.py <mouse folder> --plan --runs 260921_002,260911_002
    python aqua_tune.py <mouse folder> --plan --grid thrARScl=1.5,2,2.5,3,4 --grid minSize=50,100,200
        -> MATLAB:  run('<mouse>\aqua_tune\run_aqua_tune.m')       # overnight; restartable
    python aqua_tune.py <mouse folder> --report                     # scores, figure, recommendation
    python aqua_tune.py <mouse folder> --report --apply             # writes the chosen Opt files

Why: "change the threshold until activity shows on the afferent" has no
criterion. For every setting in a grid this measures (afferent2p/tune.py):
  * false alarms — the setting run on a null copy of the movie (each pixel's
    time series circularly shifted by its own random offset: same pixels,
    same noise, nothing coherent) — every event found there is false; counted
    on the fibre mask, as a share of the real movie's events on the fibres;
  * recall — the setting run on the real movie plus synthetic transients of
    known amplitude (0.1 / 0.2 / 0.4 dF/F), size (80 / 200 / 400 px) and
    duration (short ~4 / long ~12 frames above half amplitude) placed on the
    fibres — the fraction recovered, corrected for chance coincidence with
    real events;
  * specificity — the share of the real movie's events on the fibres, the
    event rate on and off the fibres (off-fibre activity — autofluorescence,
    vessels — is reported, never rewarded), and how the rate changes with the
    threshold (a plateau is robust).
The recommended setting is the one with the highest recall among those whose
false-alarm share on the fibres stays at or below --max-false and whose
events are mostly on fibres. One movie per day and FOV shows
whether one setting serves every day; --report prints per movie the chosen
setting and, for a FOV recorded on several days, every day's recall and
false-alarm share at each day's chosen setting, then ONE COMMON setting for
the FOV (highest mean recall among the settings that meet the rule on every
day) — what --apply writes, as aqua_opts\<experiment>_Opt.csv:
  * spontaneous (Treg) mouse: the experiment IS the FOV across days
    (<mouse>_FOV<k>_days), so its days are compared with each other and must
    share one setting — the plan takes one movie per day per FOV
    automatically and --apply always writes the common choice;
  * optostim (RearChr) mouse: one experiment per day and FOV (pre vs post
    within the day): each day gets its own setting (a field imaged on another
    day, or a new field, is tuned on its own); with --per-day movies both a
    pre and a post run of the day are tuned and the report says whether they
    agree; what the day gets is the --day-rule (v36): 'matched' (default) =
    ONE SETTING PER CONDITION chosen so that pre and post detect the injected
    transients with the same recall (within --tol-recall) AND invent false
    events at the same rate on their null movies (the same setting when it
    serves both; aqua_opts\<experiment>_pre_Opt.csv / _post_Opt.csv +
    <experiment>_tuning.json, which aqua_analyze reads for the false-alarm-
    subtracted rates; a day without a clean matching pair is marked NOT
    COMPARABLE and gets the common setting), or 'common' = the v34/v35 rule,
    one setting common to pre and post; the day is also drawn as pictures
    (aqua_tune\<day>_FOV<k>_why_day.png, pre above post at the same
    thresholds). --common-settings forces one setting per FOV for all days.

--plan writes <mouse>\aqua_tune\: the null and injected movies and the truth
table per chosen run, parameters_for_tune.csv (one Preset column per distinct
setting: the FOV's own Opt file with the grid values replaced; setting 0 of
every movie = its current values), one batch csv per job, tune_jobs.csv and
run_aqua_tune.m — which runs the jobs one by one, exports the events and the
mean image of each, deletes the big result and skips jobs already done. At
about 1.4 min per job, the default plan for six movies (15 grid settings +
the current one, three variants = 288 jobs) is ~7 h: start it in the evening.
Re-running --plan keeps earlier movies and finished jobs; --days / --runs add
movies; --report --delete-movies frees the ~0.8 GB per movie afterwards.
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from afferent2p.aqua_io import OPT_TEMPLATE_ROWS, AQUA_SUFFIX, relocate, aqua_files            # noqa: E402
from afferent2p import tune as T                                                              # noqa: E402
from afferent2p.timecourse import day_of, read_day_labels                                     # noqa: E402


def short_run(run: str) -> str:
    """NavRearChrTIGRE-2608-1_260921_008_bin31f_1p98s -> 260921_008"""
    import re
    m = re.search(r"(\d{6}_\d{3})", run)
    return m.group(1) if m else run


def load_manifest(mouse_dir: Path) -> dict:
    p = mouse_dir / "aqua_batch" / "manifest.json"
    if not p.exists():
        raise SystemExit(f"{p} not found — run aqua_batch.py (tab 4) first; the tuning uses its presets as the base settings")
    return json.loads(p.read_text())


def choose_movies(man: dict, mouse_dir: Path, runs_filter: list[str] | None, one_per_fov: bool,
                  days_filter: list[str] | None = None, per_condition: str = "auto", log=print) -> list[dict]:
    """[{experiment, fov, run, tif, preset_values, spf, condition, paradigm}] —
    per_condition "auto" / "yes": the first run of EVERY condition of each
    experiment — for a spontaneous (Treg) mouse, whose experiments span the
    days, one movie per day and FOV so that one setting can be checked on all
    days; for an optostim mouse the first pre AND the first post run of each
    day and FOV, so that the report can say whether pre and post agree (the
    day still gets one setting). "no": only the first baseline run of every
    experiment (half the sweep for an optostim mouse). --days keeps only the
    experiments whose name contains one of the given dates (spontaneous: the
    runs whose name does), --runs names runs explicitly (any condition)."""
    presets = {p["preset"]: p["values"] for p in man.get("presets", [])}
    out, seen_fov = [], set()
    if not any(e.get("runs") for e in man["experiments"]):
        raise SystemExit("the batch has no runs: tab 4 found no settings file for any FOV — set 'settings for a FOV that "
                         "has no *_Opt.csv' to default_Opt.csv (aqua_batch.py --opts default_Opt.csv) and build the batch again")
    for e in man["experiments"]:
        runs = e.get("runs", [])
        if not runs:
            continue
        spont = e.get("paradigm") == "spontaneous"
        if days_filter and not spont and not any(d in e["name"] for d in days_filter):
            continue
        if one_per_fov and e.get("fov") in seen_fov:
            continue
        conds = [r["condition"] for r in runs]
        base = "pre" if "pre" in conds else conds[0]
        every = per_condition in ("yes", "auto")
        if runs_filter:
            picks = [r for r in runs if any(s in r["run"] for s in runs_filter)]
            if not picks:
                continue
        elif every:
            picks, seen_c = [], set()
            for r in runs:
                if r["condition"] not in seen_c:
                    picks.append(r); seen_c.add(r["condition"])
        else:
            picks = ([r for r in runs if r["condition"] == base] or runs)[:1]
        if days_filter and spont:
            picks = [r for r in picks if any(d in r["run"] for d in days_filter)]
        for r in picks:
            tif = relocate(r["tif"], mouse_dir)
            if not Path(tif).exists():
                log(f"  {r['run']}: movie not found ({tif}) — skipped"); continue
            out.append({"experiment": e["name"], "fov": e.get("fov"), "run": r["run"], "tif": str(tif),
                        "condition": r.get("condition"), "paradigm": e.get("paradigm"),
                        "preset_values": presets.get(r.get("preset"), {}), "spf": float(r.get("seconds_per_frame", 2.0))})
            seen_fov.add(e.get("fov"))
    return out


def merge_grid(grid: dict, extra: dict) -> dict:
    """The grid with the suggested values added (sorted, no duplicates); a new variable is not added."""
    out = {k: sorted(set(float(v) for v in vals)) for k, vals in grid.items()}
    for k, vals in (extra or {}).items():
        if k in out:
            out[k] = sorted(set(out[k]) | {float(v) for v in vals})
    return out


def read_grid_suggestion(mouse_dir: Path) -> dict | None:
    """aqua_tune\\grid_suggestion.json written by --report when the grid did not bracket the optimum (v37)."""
    p = Path(mouse_dir) / "aqua_tune" / "grid_suggestion.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def plan(mouse_dir: Path, aqua_dir: str, grid: dict, runs_filter, one_per_fov: bool, um_per_px: float,
         n_inject: int, do_null: bool, do_inject: bool, as_root: str | None, seed: int,
         days_filter=None, n_workers: int = 0, amplitudes=None, sizes=None, per_condition: str = "auto", extend: bool = False,
         log=print):
    from aqua_batch import _as_seen
    seen = lambda p: _as_seen(p, mouse_dir, as_root)       # noqa: E731
    man = load_manifest(mouse_dir)
    aqua_dir = aqua_dir or man.get("aqua_dir") or r"D:\AQuA2"
    movies = choose_movies(man, mouse_dir, runs_filter, one_per_fov, days_filter=days_filter, per_condition=per_condition, log=log)
    if not movies:
        raise SystemExit("no movies selected (check --runs / --days; the manifest lists the runs aqua_batch.py knows)")
    tune_dir = mouse_dir / "aqua_tune"
    tune_dir.mkdir(exist_ok=True)
    old_man: dict = {}
    old_movies = []
    if (tune_dir / "tune_manifest.json").exists():
        try:
            old_man = json.loads((tune_dir / "tune_manifest.json").read_text())
            old_movies = old_man["movies"]
        except Exception:
            old_man, old_movies = {}, []
    # the injected classes: as given, else those of the earlier plan (a changed class list means new injected movies
    # and every inject job again), else the defaults
    def _is_old(vals, olds):
        return any(tuple(round(float(v), 4) for v in vals) == tuple(round(float(o), 4) for o in old) for old in olds)
    if amplitudes and _is_old(amplitudes, T.OLD_AMPLITUDES):
        log(f"  injected amplitudes {', '.join(f'{a:g}' for a in amplitudes)} are the classes of an earlier version, which probe BELOW the real "
            f"events (reported at ~0.05-0.2 dF/F) — replaced by {', '.join(f'{a:g}' for a in T.AMPLITUDES)} (v37; a class list that is not "
            "the old default is used as given)")
        amplitudes = T.AMPLITUDES
    elif amplitudes:
        amplitudes = tuple(float(a) for a in amplitudes)
    elif old_man.get("amplitudes") and not _is_old(old_man["amplitudes"], T.OLD_AMPLITUDES):
        amplitudes = tuple(float(a) for a in old_man["amplitudes"])
    else:
        amplitudes = T.AMPLITUDES
    if sizes and _is_old(sizes, T.OLD_SIZES):
        log(f"  injected sizes {', '.join(str(int(x)) for x in sizes)} px are the classes of an earlier version — replaced by "
            f"{', '.join(str(x) for x in T.SIZES_PX)} px (v37: a small class tells whether minSize lets small events through)")
        sizes = T.SIZES_PX
    elif sizes:
        sizes = tuple(int(x) for x in sizes)
    elif old_man.get("sizes_px") and not _is_old(old_man["sizes_px"], T.OLD_SIZES):
        sizes = tuple(int(x) for x in old_man["sizes_px"])
    else:
        sizes = T.SIZES_PX
    # --extend (v37): the values the last report suggested beyond the edge of the grid are added to it
    if extend:
        sug = read_grid_suggestion(mouse_dir)
        if sug and sug.get("suggest"):
            before = dict(grid)
            grid = merge_grid(grid, sug["suggest"])
            added = {k: [v for v in grid[k] if v not in [float(x) for x in before.get(k, [])]] for k in grid}
            added = {k: v for k, v in added.items() if v}
            log("grid extended with the last report's suggestion: " + (", ".join(f"{k} += {', '.join(T._fmt(v) for v in vals)}" for k, vals in added.items()) if added else "nothing new")
                + f"  (movies: {', '.join(short_run(r) + ' ' + d.get('status', '') for r, d in sug.get('movies', {}).items() if d.get('suggest'))})")
        else:
            log("--extend: the last report suggested no extension (grid_suggestion.json missing or empty) — the grid is used as given")
    grid_keys = list(grid)
    all_presets: list[dict] = []
    jobs, movie_entries = [], []

    def add_jobs(run, variants, settings, preset_idx, mdir):
        for s, pidx in zip(settings, preset_idx):
            for vname, vpath in variants.items():
                out_dir = Path(mdir) / f"s{s['setting']:02d}" / vname
                jobs.append({"run": run, "variant": vname, "setting": s["setting"], "preset": pidx + 1,
                             "input": Path(vpath), "out_dir": out_dir,
                             "mat_path": out_dir / f"{Path(vpath).stem}{AQUA_SUFFIX}.mat",
                             "events_csv": out_dir / f"{Path(vpath).stem}{AQUA_SUFFIX}_events.csv"})

    def preset_numbers(base, settings):
        idx = []
        for s in settings:
            vals = dict(base); vals.update({k: s[k] for k in grid_keys})
            if vals not in all_presets:
                all_presets.append(vals)
            idx.append(all_presets.index(vals))
        return idx
    day_labels = read_day_labels(mouse_dir)
    log(f"tuning plan for {mouse_dir.name}: {len(movies)} movie(s), grid {', '.join(f'{k}={v}' for k, v in grid.items())}")
    for mv in movies:
        dl = day_labels.get(day_of(mv["run"]) or "", "")
        log(f"  {short_run(mv['run'])}  FOV {mv['fov']}  [{mv.get('condition', '')}]" + (f"  ({dl})" if dl else ""))
    log("  (narrow with --days / --runs; the GUI fields 'days' and 'runs')")
    for k_mv, mv in enumerate(movies, 1):
        tif = Path(mv["tif"])
        mdir = tune_dir / mv["run"]
        mdir.mkdir(exist_ok=True)
        log(f"\n{mv['run']}  ({mv['experiment']}, FOV {mv['fov']})  [{k_mv}/{len(movies)}]")
        base = dict(mv["preset_values"])
        if not base:
            log("  no preset values in the manifest for this run — using the template defaults")
            base = {var: default for _, var, _, default in OPT_TEMPLATE_ROWS if var}
            base["frameRate"] = f"{mv['spf']:.6g}"; base["spatialRes"] = f"{um_per_px:.6g}"
        settings = T.settings_from_grid(base, grid)
        # results of an earlier plan with other settings would be read under the wrong setting number: move them aside
        # (every movie folder records the settings its sNN folders belong to)
        sfile = mdir / "settings.json"
        prev = None
        if sfile.exists():
            try:
                prev = json.loads(sfile.read_text())
            except Exception:
                prev = None
        if prev != settings:                      # no settings.json = unknown provenance = stale
            # the finished results are matched to the new list by their setting VALUES (v37: whatever their old
            # number): a folder whose values are still in the grid is renamed to the new number of those values
            # and keeps its finished jobs — the MATLAB script skips them; only folders whose values left the grid
            # are moved aside. So extending the grid costs only the new combinations.
            prev_list = prev if isinstance(prev, list) else []

            def vals_of(d: dict) -> tuple:
                return tuple(sorted((x, v) for x, v in d.items() if x not in ("label", "setting")))
            new_by_vals = {vals_of(st): st["setting"] for st in settings}
            folders = sorted(d for d in mdir.glob("s[0-9][0-9]") if d.is_dir())
            moves, stale = {}, []
            for d in folders:
                k = int(d.name[1:])
                if k < len(prev_list) and vals_of(prev_list[k]) in new_by_vals:
                    moves[d] = new_by_vals[vals_of(prev_list[k])]
                else:
                    stale.append(d)
            import time
            stamp = time.strftime('%y%m%d_%H%M')
            if stale:
                keep = mdir / f"stale_{stamp}"
                keep.mkdir(exist_ok=True)
                for d in stale:
                    d.rename(keep / d.name)
            renamed = 0
            if any(d.name != f"s{n:02d}" for d, n in moves.items()):
                tmp = {d: mdir / f"_renum_{stamp}_{d.name}" for d in moves}      # two steps: numbers may swap
                for d, t in tmp.items():
                    d.rename(t)
                for d, n in moves.items():
                    tmp[d].rename(mdir / f"s{n:02d}"); renamed += d.name != f"s{n:02d}"
            log(f"  the settings list changed since the last plan — {len(moves)} finished folder(s) keep their results"
                + (f" ({renamed} renumbered)" if renamed else "")
                + (f"; {len(stale)} folder(s) with values no longer in the grid moved to stale_{stamp}\\ (delete it when done)" if stale else "")
                + ("" if prev_list else " — no settings.json: results of an older version, all redone"))
        sfile.write_text(json.dumps(settings, indent=1), encoding="utf-8")
        # presets: base values with the swept ones replaced; dedupe across movies
        preset_idx = preset_numbers(base, settings)
        # where the fibres are in the raw mean image (for the injection sites); a saved mask is reused so that a
        # deleted injected movie is regenerated with the same transients its finished jobs were scored against
        mean = T.mean_image(tif)
        saved = mdir / "fibre_mask_raw.npy"
        saved_meta = mdir / "fibre_mask_raw.json"
        saved_mask = None
        saved_version = 1
        if saved.exists():
            try:
                saved_mask = np.load(saved).astype(bool)
                if saved_mask.shape != mean.shape:
                    saved_mask = None
                if saved_meta.exists():
                    saved_version = int(json.loads(saved_meta.read_text(encoding="utf-8")).get("version", 1))
            except Exception:
                saved_mask = None
        fmask, how = T.manual_mask(tif, mouse_dir, mv["fov"], shape=mean.shape)      # a hand-drawn mask always wins
        if how.startswith("ignored"):
            log(f"  {how}")
        if fmask is None:
            if saved_mask is not None and saved_version >= T.MASK_VERSION:
                fmask, how = saved_mask, "saved mask of the earlier plan"
            else:
                fmask, how = T.injection_mask(mean, um_per_px)
                if saved_mask is not None:
                    log(f"  the saved fibre mask is of an older version (v{saved_version}) — made again (v{T.MASK_VERSION}: the fibre "
                        "network is kept whole; before, a dense network of crossing fibres could be dropped as a 'round speck')")
        mask_changed = saved_mask is not None and not np.array_equal(saved_mask, fmask)
        np.save(saved, fmask)
        saved_meta.write_text(json.dumps({"version": T.MASK_VERSION, "how": how, "px": int(fmask.sum())}), encoding="utf-8")
        try:
            import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 5.2)); ax.imshow(mean, cmap="gray", vmin=np.percentile(mean, 1), vmax=np.percentile(mean, 99.5))
            ax.contour(fmask.astype(float), levels=[0.5], colors=["#eda100"], linewidths=0.5); ax.set_axis_off()
            ax.set_title(f"{mv['run']}: injection sites are placed on this mask ({how})", fontsize=8)
            fig.tight_layout(); fig.savefig(mdir / "fibre_mask_raw.png", dpi=110); plt.close(fig)
        except Exception:
            pass
        log(f"  injection mask: {int(fmask.sum())} px ({how})")
        variants = {"real": tif}
        if do_null:
            p = mdir / f"{tif.stem}_null.tif"
            if not p.exists():
                T.make_null_movie(tif, p, seed=seed, sec_per_frame=mv["spf"], log=log)
            else:
                log(f"  null movie exists: {p.name}")
            variants["null"] = p
        truth_csv = mdir / f"{tif.stem}_inject_truth.csv"
        if do_inject:
            p = mdir / f"{tif.stem}_inject.tif"
            classes_changed = False
            if truth_csv.exists():
                try:
                    import pandas as pd
                    old = pd.read_csv(truth_csv)
                    classes_changed = (sorted(set(np.round(old["amplitude"], 4))) != sorted(set(np.round(amplitudes, 4)))
                                       or sorted(set(int(v) for v in old["size_class"])) != sorted(set(sizes)))
                except Exception:
                    classes_changed = False
            if classes_changed:
                log(f"  injected classes changed (now {', '.join(f'{a:g}' for a in amplitudes)} dF/F; {', '.join(str(v) for v in sizes)} px)")
            if (mask_changed or classes_changed) and (p.exists() or truth_csv.exists()):
                # new sites -> new injected movie; results scored against the old sites are moved aside
                import time
                old = [d for d in mdir.glob("s[0-9][0-9]/inject") if d.is_dir()]
                if old:
                    keep = mdir / f"stale_{time.strftime('%y%m%d_%H%M')}_inject"
                    keep.mkdir(exist_ok=True)
                    for d in old:
                        d.rename(keep / d.parent.name)
                log(("  the injection mask changed" if mask_changed else "  the injected classes changed") + " — the injected movie is made again"
                    + (f"; {len(old)} earlier inject result folder(s) moved aside" if old else ""))
                for q in (p, truth_csv):
                    if q.exists():
                        q.unlink()
            if not p.exists() or not truth_csv.exists():
                if fmask.sum() < 50:
                    log("  no usable fibre pixels in the mean image — no injected movie for this run (recall cannot be measured; "
                        "false alarms and the event rate still are)")
                else:
                    # the movie is made again (deleted after scoring, say): the same mask and classes give the same sites;
                    # if the sites differ from the earlier truth table, results scored against the old sites are moved aside
                    old_truth = None
                    if truth_csv.exists():
                        try:
                            import pandas as pd
                            old_truth = pd.read_csv(truth_csv)
                        except Exception:
                            old_truth = None
                    T.make_injected_movie(tif, p, truth_csv, fmask, n_sites=n_inject, amplitudes=amplitudes, sizes=sizes,
                                          seed=seed + 1, sec_per_frame=mv["spf"], log=log)
                    if old_truth is not None:
                        try:
                            import pandas as pd
                            new_truth = pd.read_csv(truth_csv)
                            cols = ["site_y", "site_x", "amplitude", "size_class", "kind", "t0_frame"]
                            same = len(old_truth) == len(new_truth) and old_truth[cols].round(4).equals(new_truth[cols].round(4))
                        except Exception:
                            same = False
                        if not same:
                            import time
                            old = [d for d in mdir.glob("s[0-9][0-9]/inject") if d.is_dir()]
                            if old:
                                keep = mdir / f"stale_{time.strftime('%y%m%d_%H%M')}_inject"
                                keep.mkdir(exist_ok=True)
                                for d in old:
                                    d.rename(keep / d.parent.name)
                                log(f"  the injected transients differ from the earlier truth table — {len(old)} earlier inject result folder(s) moved aside")
            else:
                log(f"  injected movie exists: {p.name}")
            if p.exists() and truth_csv.exists():
                variants["inject"] = p
        add_jobs(mv["run"], variants, settings, preset_idx, mdir)
        movie_entries.append({"experiment": mv["experiment"], "fov": mv["fov"], "run": mv["run"], "tif": str(tif),
                              "condition": mv.get("condition"), "paradigm": mv.get("paradigm"),
                              "spf": mv["spf"], "dir": str(mdir), "truth_csv": str(truth_csv) if "inject" in variants else None,
                              "variants": {k: str(v) for k, v in variants.items()},
                              "settings": settings, "preset": [i + 1 for i in preset_idx], "base": base})
    # movies planned earlier and not re-selected now stay in the plan (their jobs are already done or still queued) —
    # unless their run is no longer in the batch at all (left out by run_labels.csv, or removed; v44): then nothing
    # downstream would use the result, so the movie is dropped and no jobs are queued for it
    new_runs = {m["run"] for m in movie_entries}
    batch_runs = {r["run"] for e in man.get("experiments", []) for r in e.get("runs", [])}
    kept = 0
    for om in old_movies:
        if om["run"] in new_runs or not Path(om.get("dir", "")).exists():
            continue
        if batch_runs and om["run"] not in batch_runs:
            log(f"\n{om['run']}: planned earlier but no longer in the batch (run_labels.csv / removed) — dropped from the plan; "
                f"its folder aqua_tune\\{Path(om['dir']).name} can be deleted")
            continue
        if list(om.get("settings", [{}])[0].keys()) != list(movie_entries[0]["settings"][0].keys()):
            log(f"\n{om['run']}: planned earlier with other grid variables — dropped from the plan (re-plan it with --runs)")
            continue
        pidx = preset_numbers(om["base"], om["settings"])
        have = {k: v for k, v in om["variants"].items() if Path(v).exists()}
        gone = sorted(set(om["variants"]) - set(have))
        if gone:
            log(f"\n{om['run']}: its {', '.join(gone)} movie(s) were deleted — finished jobs stay in the report, nothing new is queued for them")
        add_jobs(om["run"], have, om["settings"], pidx, om["dir"])
        om = dict(om); om["preset"] = [i + 1 for i in pidx]
        movie_entries.append(om); kept += 1
    if kept:
        log(f"\n{kept} movie(s) from the earlier plan kept (finished jobs are skipped by the MATLAB script)")
    paths = T.write_tune_batch(tune_dir, jobs, all_presets, OPT_TEMPLATE_ROWS, aqua_dir, seen=seen, n_workers=n_workers)
    (tune_dir / "tune_manifest.json").write_text(json.dumps({
        "mouse": mouse_dir.name, "aqua_dir": aqua_dir, "um_per_px": um_per_px, "grid": grid, "grid_keys": grid_keys,
        "movies": movie_entries, "n_jobs": len(jobs), "amplitudes": list(amplitudes), "sizes_px": list(sizes),
        "kinds": list(T.KINDS),
        "files": {"script": Path(paths["script"]).name, "jobs_csv": Path(paths["jobs_csv"]).name,
                  "params_csv": Path(paths["params_csv"]).name}}, indent=1), encoding="utf-8")
    n_done = sum(1 for j in jobs if Path(j["events_csv"]).exists())
    log(f"\n{len(jobs)} jobs ({len(movie_entries)} movies x {len(movie_entries[0]['settings'])} settings x {len(variants)} variants; "
        f"{len(all_presets)} presets)" + (f", {n_done} done already" if n_done else "")
        + f" -> about {(len(jobs) - n_done) * 1.4 / 60:.1f} h at 1.4 min each")
    log(f"MATLAB:  run('{seen(paths['script'])}')      (restartable; then: python aqua_tune.py \"{mouse_dir}\" --report)")
    if int(n_workers) == 1:
        log("  one job at a time (--workers 0 lets MATLAB choose a parallel pool from the machine, N = that many jobs at once)")
    else:
        n_split = max(2, int(n_workers))
        log(f"  jobs at once: {'chosen from the machine' if int(n_workers) == 0 else n_workers} with the Parallel Computing Toolbox; "
            f"without it: run_aqua_tune_parallel.bat starts {n_split} MATLAB processes (or open {n_split} MATLAB windows and run "
            f"run_aqua_tune_<k>of{n_split}.m in each)")
    return tune_dir


def find_reference(tif: Path) -> tuple:
    """A hand-made AQuA2 result for this movie (v37): the AQuA2 GUI saves into a subfolder named after the movie,
    <run folder>\\<stem>\\<stem>_AQuA2.mat; exported with aqua_export_events.m it has <stem>_AQuA2_events.csv (and a
    _meanimg.tif) next to it. Also accepted: <run folder>\\<stem>_reference_events.csv. Returns (events csv or None,
    mean image tif or None, note) — the note says when a .mat exists but is not exported yet."""
    tif = Path(tif); stem = tif.stem
    sub = tif.parent / stem
    cands = [sub / f"{stem}_AQuA2_events.csv", tif.parent / f"{stem}_reference_events.csv"]
    if sub.is_dir():
        cands += sorted(sub.glob("*_AQuA2_events.csv"))
    for c in cands:
        if c.exists():
            mi = c.with_name(c.name.replace("_events.csv", "_meanimg.tif"))
            return c, (mi if mi.exists() else None), ""
    mats = sorted(sub.glob("*_AQuA2.mat")) if sub.is_dir() else []
    if mats:
        return None, None, (f"a hand-made AQuA2 result exists ({mats[0].name} in {sub.name}\\) but its events are not exported — in MATLAB: "
                            f"aqua_export_events('{mats[0]}')  (then --report again)")
    return None, None, ""


def _is_close(a, b, tol: float) -> bool:
    try:
        return len(a) == len(b) and all(abs(float(x) - float(y)) <= tol * max(float(y), 1e-9) for x, y in zip(a, b))
    except Exception:
        return False


def _tuning_row(x: dict) -> dict:
    """The numbers of one scored (movie, setting) that the analysis needs: sensitivity, contamination, rates."""
    def f(v):
        try:
            v = float(v)
            return v if np.isfinite(v) else None
        except (TypeError, ValueError):
            return None
    return {"label": x.get("label"), "recall": f(x.get("recall")), "false_share": f(x.get(T.false_key([x]))),
            "fp_rate_per_min": f(T._fp_rate(x)), "fp_rate_all_per_min": f(x.get("rate_null")),
            "rate_real_on": f(x.get("rate_real_on", x.get("rate_real"))),
            "rate_real": f(x.get("rate_real")), "on_fibre": f(x.get("on_fibre")), "n_real": x.get("n_real"), "minutes": f(x.get("minutes"))}


def _write_tuning_json(path: Path, data: dict) -> None:
    """aqua_opts\\<experiment>_tuning.json: which rule chose the day's settings and the measured sensitivity / false-alarm
    rate per condition — read by aqua_analyze.py for the false-alarm-subtracted rates and the report's note."""
    try:
        path.write_text(json.dumps(data, indent=1, default=str), encoding="utf-8")
    except Exception:
        pass


def write_why_html(tune_dir: Path, man: dict, why_pngs: dict, chosen: dict, passed_rule: dict, max_false: float, min_on_fibre: float,
                   day_labels: dict, day_pngs: dict | None = None, day_matched: dict | None = None, day_rule: str = "matched",
                   grid_keys: list | None = None, tol_recall: float = 0.05) -> Path:
    """aqua_tune\\tune_why.html: the per-movie 'why this setting' pictures with the chosen setting and the rule (v35); the
    pre/post day pictures with the matched-sensitivity verdict first (v36)."""
    import html as H
    parts = [f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>AQuA2 tuning — why this setting? ({H.escape(man['mouse'])})</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;margin:20px;color:#222;max-width:1700px}} h2{{margin-top:28px}} img{{max-width:100%;border:1px solid #ddd}}
p.rule{{background:#f4f6f8;padding:8px 10px;border-left:4px solid #2ecc40}} .k{{font-family:Consolas,monospace}}</style></head><body>
<h1>AQuA2 tuning — why this setting? <span style="font-weight:normal">{H.escape(man['mouse'])}</span></h1>
<p class="rule">For every tuned movie, one column per setting on the same movie, three rows: the <b>real</b> movie with every detected event drawn as a circle of its size
(green on the fibres, orange off them); the <b>null</b> movie — the same pixels with their time series scrambled, so every detection is a false alarm (red);
the <b>injected</b> movie — the real movie plus 72 synthetic transients of known amplitude / size / duration on the fibres (found = filled cyan, missed = hollow magenta).
The rule picks the most sensitive setting (highest recall of the injected transients) whose false alarms on the fibres stay &le; {max_false:.0%} of the real events
and whose events are &ge; {min_on_fibre:.0%} on the fibres: lower thresholds add red (false alarms) faster than green (real events); higher thresholds turn cyan into magenta (missed transients).
The chosen column is framed green; a red tag marks a setting that fails the false-alarm rule; grey tags say what the other passing settings lose. Click a picture for full size.</p>"""]
    if day_pngs:
        parts.append(f"""<h1 style="font-size:1.3em;margin-top:30px">Pre and post of each day, side by side (matched sensitivity)</h1>
<p class="rule">Same settings, same columns, the pre movie above the post movie: what one threshold finds in each. The day's rule picks <b>one setting per
condition</b> so that pre and post detect the injected transients with the same recall (within {tol_recall:.2f}) and invent false
events at the same rate on the null movies — the frame of each movie's own setting is green (pre may get one column, post another); the dashed frame is the one
common setting; a red frame fails the false-alarm rule on that movie. When no pair of clean settings matches, the day is marked NOT COMPARABLE and the common
setting is written instead.</p>""")
        for exp, png in day_pngs.items():
            dm = (day_matched or {}).get(exp)
            exp_short = exp[len(man["mouse"]) + 1:] if exp.startswith(man["mouse"] + "_") else exp
            dl = day_labels.get(exp_short.split("_")[0], "")
            head = H.escape(exp_short) + (f" — {H.escape(dl)}" if dl else "")
            if dm and dm[2] and dm[0]:
                mchoice, mrows, _, minfo, mwhy, conds = dm
                gk = grid_keys or []
                if minfo.get("same_setting"):
                    verdict = ("the same setting serves pre and post: <span class='k'>"
                               + H.escape(" ".join(f"{k}={v}" for k, v in zip(gk, next(iter(mchoice.values()))))) + "</span> — "
                               + "; ".join(f"{H.escape(str(conds[r]))}: recall {x.get('recall', float('nan')):.2f}, false alarms {T._fp_rate(x):.2f}/min" for r, x in mrows.items()))
                else:
                    verdict = ("<b>one setting per condition</b> — "
                               + "; ".join(f"{H.escape(str(conds[r]))}: <span class='k'>{H.escape(' '.join(f'{k}={v}' for k, v in zip(gk, mchoice[r])))}</span> "
                                           f"(recall {x.get('recall', float('nan')):.2f}, false alarms {T._fp_rate(x):.2f}/min)" for r, x in mrows.items())
                               + f" — recalls within {minfo.get('recall_gap', float('nan')):.2f}, false-alarm rates within {minfo.get('fp_gap', float('nan')):.2f}/min")
                if day_rule != "matched":
                    verdict += " (not applied: the day rule is 'one common setting')"
            elif dm:
                verdict = f"<span style='color:#b30000'>no matched pair — {H.escape(str(dm[4]))}</span>" + (
                    "; the day is marked NOT COMPARABLE and gets the common setting" if day_rule == "matched" else "")
            else:
                verdict = ""
            parts.append(f"<h2>{head}</h2><p>{verdict}</p><p><a href='{H.escape(png.name)}' target='_blank'><img src='{H.escape(png.name)}' alt='{H.escape(png.name)}'></a></p>")
        parts.append("<h1 style=\"font-size:1.3em;margin-top:30px\">Each tuned movie on its own</h1>")
    for mv in man["movies"]:
        run = mv["run"]
        png = why_pngs.get(run)
        if png is None:
            continue
        ch = chosen.get(run)
        dl = day_labels.get(day_of(run) or "", "")
        head = f"{H.escape(short_run(run))} — FOV {mv['fov']}" + (f" — {H.escape(dl)}" if dl else "") + (f" — {H.escape(str(mv.get('condition')))}" if mv.get("condition") else "")
        verdict = (f"chosen: <span class='k'>{H.escape(ch['label'])}</span> (recall {ch.get('recall', float('nan')):.2f}, false alarms on fibres "
                   f"{ch.get(T.false_key([ch]), float('nan')):.0%}, {ch.get('on_fibre', float('nan')):.0%} on fibres)" if ch and passed_rule.get(run)
                   else (f"no setting met the rule — best trade-off: <span class='k'>{H.escape(ch['label'])}</span>" if ch else "no setting chosen"))
        parts.append(f"<h2>{head}</h2><p>{verdict}</p><p><a href='{H.escape(png.name)}' target='_blank'><img src='{H.escape(png.name)}' alt='{H.escape(png.name)}'></a></p>")
    parts.append("<p style='color:#666'>Settings and numbers: tune_recommendation.txt and tune_report.csv in this folder; the curves: &lt;mouse&gt;_tune.png.</p></body></html>")
    out = tune_dir / "tune_why.html"
    out.write_text("\n".join(parts), encoding="utf-8")
    return out


def movie_context(mouse_dir: Path, man: dict, mv: dict, log=print) -> dict:
    """Everything the scoring, the pictures and the picker need of one tuned movie (v38 — shared by report() and the
    by-eye picker): its folder, movie, truth table, length in minutes, mean image, fibre mask in the raw frame
    (dilated), AQuA2's crop offsets per variant, a hand-made reference (if any) and `read_events(setting, variant)`."""
    import pandas as pd
    import tifffile as tf
    from scipy.ndimage import binary_dilation
    tune_dir = Path(mouse_dir) / "aqua_tune"
    mdir = Path(mv["dir"]); tif = relocate(mv["tif"], mouse_dir)
    if not mdir.exists():
        mdir = tune_dir / mv["run"]
    truth = pd.read_csv(mv["truth_csv"]) if mv.get("truth_csv") and Path(mv["truth_csv"]).exists() else None
    if truth is None and mv.get("truth_csv"):
        alt = mdir / Path(mv["truth_csv"]).name
        truth = pd.read_csv(alt) if alt.exists() else None
    fl = aqua_files(tif)
    frames = None
    if fl["json"].exists():
        try:
            j = json.loads(fl["json"].read_text()); frames = int(j.get("frames_out") or j.get("frames_recorded") or 0)
        except Exception:
            frames = None
    if not frames:
        frames = int(T._open(tif).shape[0])
    minutes = frames * float(mv["spf"]) / 60.0
    raw_mean = T.mean_image(tif)
    fmask_raw, mask_how = None, ""
    saved = mdir / "fibre_mask_raw.npy"
    if saved.exists():
        try:
            fmask_raw = np.load(saved).astype(bool); mask_how = "injection mask of the plan"
            if fmask_raw.shape != raw_mean.shape:
                fmask_raw = None
        except Exception:
            fmask_raw = None
    if fmask_raw is None:
        fmask_raw, mask_how = T.injection_mask(raw_mean, float(man["um_per_px"]))
    fmask_raw = binary_dilation(fmask_raw, iterations=3)
    offsets, crop_notes = {}, []
    for vname in mv["variants"]:
        for st in mv["settings"]:
            mi = list((mdir / f"s{st['setting']:02d}" / vname).glob("*_meanimg.tif"))
            if not mi:
                continue
            try:
                exp_mean = np.asarray(tf.imread(str(mi[0])), np.float32)
            except Exception:
                continue
            off = T.crop_offset(raw_mean, exp_mean)
            if off is None:
                log(f"  {mv['run']} ({vname}): could not locate AQuA2's crop in the raw frame — assuming it starts at (0, 0)")
                offsets[vname] = (0, 0, tuple(exp_mean.shape))
            else:
                offsets[vname] = (off[0], off[1], tuple(exp_mean.shape))
                crop_notes.append(f"{vname} {exp_mean.shape[0]}x{exp_mean.shape[1]} at ({off[0]}, {off[1]}), match {off[2]:.2f}")
            break
    ref_csv, ref_mean_tif, ref_note = find_reference(tif)
    reference, ref_off = None, None
    if ref_csv is not None:
        try:
            reference = pd.read_csv(ref_csv)
            ref_off = offsets.get("real") or (0, 0, None)
            if ref_mean_tif is not None:
                rm = np.asarray(tf.imread(str(ref_mean_tif)), np.float32)
                o = T.crop_offset(raw_mean, rm)
                ref_off = (o[0], o[1], tuple(rm.shape)) if o is not None else ref_off
            log(f"  {mv['run']}: hand-made result as the reference — {len(reference)} events ({ref_csv.name})")
        except Exception as e:
            log(f"  {mv['run']}: could not read the hand-made result {ref_csv} ({e})")
            reference = None
    elif ref_note:
        log(f"  {mv['run']}: {ref_note}")

    def read_events(setting: int, variant: str, _mdir=mdir):
        d = _mdir / f"s{int(setting):02d}" / variant
        if not d.exists():
            return None
        hits = list(d.glob("*_AQuA2_events.csv"))
        if not hits:
            return None
        try:
            return pd.read_csv(hits[0])
        except Exception:
            return None
    return {"mv": mv, "mdir": mdir, "tif": tif, "truth": truth, "minutes": minutes, "frames": frames, "raw_mean": raw_mean,
            "fmask_raw": fmask_raw, "mask_how": mask_how, "offsets": offsets, "crop_notes": crop_notes,
            "reference": reference, "ref_off": ref_off, "ref_note": ref_note, "read_events": read_events}


def report(mouse_dir: Path, max_false: float, min_on_fibre: float, target_amp: float | None, apply: bool,
           apply_all: bool = False, delete_movies: bool = False, apply_common: bool = True, day_rule: str = "matched",
           tol_recall: float = 0.05, min_recall: float = 0.3, figure_dpi: int = 150, log=print):
    import pandas as pd
    tune_dir = mouse_dir / "aqua_tune"
    verdicts: dict = {}                    # run -> grid verdict (v37)
    class_notes: dict = {}                 # run -> suggested amplitude classes when the injected ones probe below the signal
    reference_notes: dict = {}             # run -> a hand-made result that is not exported yet
    mp = tune_dir / "tune_manifest.json"
    if not mp.exists():
        raise SystemExit(f"{mp} not found — run --plan first")
    man = json.loads(mp.read_text())
    grid_keys = man["grid_keys"]
    day_labels = read_day_labels(mouse_dir)
    rows_all, per_movie, chosen, passed_rule, lines = [], {}, {}, {}, []
    why_pngs: dict = {}
    movie_ctx: dict = {}
    day_pngs: dict = {}
    lines.append(f"AQuA2 tuning report — {man['mouse']}   (grid: {', '.join(f'{k}={v}' for k, v in man['grid'].items())})")
    lines.append(f"rule: highest recall{' at ' + f'{target_amp:g}' + ' dF/F' if target_amp else ''} among the settings with "
                 f"false-alarm share <= {max_false:.0%} and on-fibre fraction >= {min_on_fibre:.0%}")
    for mv in man["movies"]:
        ctx = movie_context(mouse_dir, man, mv, log=log)
        truth, minutes = ctx["truth"], ctx["minutes"]
        raw_mean, fmask_raw, mask_how, offsets, crop_notes = ctx["raw_mean"], ctx["fmask_raw"], ctx["mask_how"], ctx["offsets"], ctx["crop_notes"]
        reference, ref_off, read_setting_events = ctx["reference"], ctx["ref_off"], ctx["read_events"]
        if ctx["ref_note"]:
            reference_notes[mv["run"]] = ctx["ref_note"]
        rows = []
        for s in mv["settings"]:
            def read(variant, _s=s):
                return read_setting_events(_s["setting"], variant)
            real, null, inj = read("real"), read("null"), read("inject")
            if real is None:
                continue
            sc = T.score_setting(real, null, inj, truth, minutes, fmask_raw=fmask_raw, offsets=offsets)
            o_real = offsets.get("real", (0, 0, None))
            if reference is not None:
                rs = T.score_reference(real, reference, o_real, ref_off)
                sc.update({k: v for k, v in rs.items() if k != "ref_hit"})
            row = {"run": mv["run"], "experiment": mv["experiment"], "fov": mv["fov"], "setting": s["setting"], "label": s["label"],
                   **{k: s[k] for k in grid_keys}, "minutes": round(minutes, 1), "crop_dy": o_real[0], "crop_dx": o_real[1], **sc}
            rows.append(row); rows_all.append(row)
        per_movie[mv["run"]] = rows
        # what the day figure (pre and post side by side, v36) needs of this movie
        movie_ctx[mv["run"]] = {"cond": mv.get("condition") or short_run(mv["run"]), "rows": rows, "raw_mean": raw_mean,
                                "fmask_raw": fmask_raw, "offsets": offsets, "truth": truth, "read_events": read_setting_events,
                                "reference": reference, "ref_off": ref_off}
        best, why, passed = T.recommend(rows, max_false=max_false, min_on_fibre=min_on_fibre, target_amp=target_amp)
        chosen[mv["run"]] = best
        passed_rule[mv["run"]] = passed
        lines.append("")
        dl = day_labels.get(day_of(mv["run"]) or "", "")
        # the sanity check as pictures: what each setting finds in the real / null / injected movie (v35)
        if rows:
            try:
                why_png = tune_dir / f"{short_run(mv['run'])}_why.png"
                T.why_figure(f"{short_run(mv['run'])}  FOV {mv['fov']}" + (f"  ({dl})" if dl else "") + (f"  [{mv.get('condition')}]" if mv.get("condition") else ""),
                             rows, best, passed, grid_keys, raw_mean, fmask_raw, offsets, truth, read_setting_events, max_false, min_on_fibre, why_png,
                             reference=reference, ref_off=ref_off)
                why_pngs[mv["run"]] = why_png
            except Exception as e:
                log(f"  {mv['run']}: the 'why this setting' figure failed ({e})")
        lines.append(f"{mv['run']}  ({mv['experiment']}{', ' + dl if dl else ''}, FOV {mv['fov']}; "
                     f"{len(rows)} of {len(mv['settings'])} settings scored"
                     + (f"; AQuA2's crop in the raw frame: {'; '.join(crop_notes)}" if crop_notes else "")
                     + f"; fibres = {mask_how}, {int(fmask_raw.sum())} px dilated)")
        if not rows:
            lines.append("  nothing to score yet — run run_aqua_tune.m")
            continue
        fk = T.false_key(rows)
        amps_lbl = [f"r{float(a):g}" for a in man.get("amplitudes", [0.1, 0.2, 0.4])]
        has_ref = any(np.isfinite(r.get("ref_recall", np.nan)) for r in rows)
        hdr = ("  setting  " + "  ".join(f"{k:>9}" for k in grid_keys)
               + "   on-fibre/min  off-fibre/min  null-on/min  false%(on)  on-fibre  recall  " + "  ".join(f"{a:>4}" for a in amps_lbl)
               + "  short  long" + ("   hand%  extra%" if has_ref else ""))
        lines.append(hdr)
        lines.append("  (events/min of the real movie on and off the fibre mask; null-on/min = null-movie events on the fibres; "
                     "false%(on) = null-on / on-fibre; recall = injected transients recovered, by amplitude and duration"
                     + ("; hand% = share of the hand-made result's events recovered, extra% = share of this setting's events the hand result does not have" if has_ref else "") + ")")
        for r in sorted(rows, key=lambda r: r["setting"]):
            def f(v, w=5, p=2):
                return f"{v:{w}.{p}f}" if isinstance(v, (int, float)) and np.isfinite(v) else f"{'-':>{w}}"
            mark = (" <- chosen" if passed_rule.get(mv["run"]) else " <- best trade-off (rule not met)") if best is r \
                else ("  (current)" if r["setting"] == 0 else "")
            lines.append(f"  {r['setting']:>7}  " + "  ".join(f"{str(r[k]):>9}" for k in grid_keys)
                         + f"   {f(r.get('rate_real_on', r['rate_real']), 12, 3)}  {f(r.get('rate_off', np.nan), 13, 3)}  "
                         + f"{f(r.get('rate_null_on', r.get('rate_null', np.nan)), 11, 3)}  "
                         + (f"{100 * r[fk]:9.0f}%" if np.isfinite(r.get(fk, np.nan)) else "         -")
                         + f"  {f(r.get('on_fibre', np.nan), 8, 2)}  {f(r.get('recall', np.nan))}  "
                         + "  ".join(f(r.get(f"recall_amp{float(a):g}", np.nan), 4) for a in man.get("amplitudes", [0.1, 0.2, 0.4])) + "  "
                         + f"{f(r.get('recall_short', np.nan), 5)}  {f(r.get('recall_long', np.nan), 4)}"
                         + ((f"   {100 * r['ref_recall']:5.0f}%  " + (f"{100 * r['ref_extra']:5.0f}%" if np.isfinite(r.get('ref_extra', np.nan)) else "     -"))
                            if has_ref and np.isfinite(r.get("ref_recall", np.nan)) else ("                " if has_ref else "")) + mark)
        if best:
            cur = next((r for r in rows if r["setting"] == 0), None)
            fw = "false-alarm share on the fibres" if fk == "false_share_on" else "false-alarm share"
            lines.append(f"  -> {'chosen' if passed_rule.get(mv['run']) else 'best trade-off (rule not met)'}: "
                         f"setting {best['setting']} ({best['label']}): recall {best.get('recall', np.nan):.2f}, "
                         f"{fw} {best.get(fk, np.nan):.0%}, on-fibre {best.get('on_fibre', np.nan):.2f}, "
                         f"{best.get('rate_real_on', best['rate_real']):.3f} events/min on the fibres   [{why}]")
            if cur and cur is not best:
                lines.append(f"     the current setting gives recall {cur.get('recall', np.nan):.2f}, {fw} "
                             f"{cur.get(fk, np.nan):.0%}, {cur.get('rate_real_on', cur['rate_real']):.3f} events/min on the fibres")
            # against the hand-made result (v37)
            if has_ref:
                rr = [r for r in rows if np.isfinite(r.get("ref_recall", np.nan))]
                n_ref = int(next((r.get("ref_n", 0) for r in rr), 0) or 0)
                def f1(r):
                    p_ = 1 - r["ref_extra"] if np.isfinite(r.get("ref_extra", np.nan)) else 0.0
                    return (2 * p_ * r["ref_recall"] / (p_ + r["ref_recall"])) if (p_ + r["ref_recall"]) > 0 else 0.0
                closest = max(rr, key=f1) if rr else None
                lines.append(f"     hand-made result ({n_ref} events): this setting recovers {best.get('ref_recall', np.nan):.0%} of them"
                             + (f" and adds {best.get('ref_extra', np.nan):.0%} the hand result does not have" if np.isfinite(best.get("ref_extra", np.nan)) else "")
                             + (f"; closest to the hand result: setting {closest['setting']} ({closest['label']}: {closest['ref_recall']:.0%} recovered, "
                                f"{closest.get('ref_extra', np.nan):.0%} extra)" if closest is not None and closest is not best else "; no other setting is closer to it"))
                if np.isfinite(best.get("ref_recall", np.nan)) and best["ref_recall"] < 0.5 and n_ref >= 10:
                    lines.append("     -> the chosen setting misses most of what was accepted by hand: the grid or the injected classes do not reach "
                                 "this field's events (see the grid verdict), or the hand result is more permissive than the false-alarm rule allows "
                                 "(its extra% at the hand setting tells)")
            # does the grid bracket the optimum, or only touch it? (v37)
            gv = T.grid_verdict(rows, best, passed, grid_keys, max_false=max_false)
            verdicts[mv["run"]] = gv
            if gv["text"]:
                lines.append(f"  -> grid: {gv['text']}")
            if int(best.get("n_real", 0) or 0) == 0 and (best.get("recall", 0) or 0) < 0.5:
                lines.append("  -> this setting finds NO event in the real movie: either the field is silent or the detection is too "
                             "strict for it (the grid verdict above); a hand-tuned AQuA2 result for this run would tell — see "
                             "'reference' in the GUIDE")
            # the injected classes against the real signal (v37): what AQuA2 reports for the brightest class it found,
            # against the real events it reports at this setting
            amps = sorted(float(a) for a in man.get("amplitudes", []))
            if amps:
                a_hi = amps[-1]
                rep_hi = best.get(f"inj_dff_reported_amp{a_hi:g}", np.nan)
                if not np.isfinite(rep_hi):
                    rep_hi = max((r.get(f"inj_dff_reported_amp{a_hi:g}", np.nan) for r in rows), default=np.nan)
                # the real signal: the hand-made result's events when there is one (what the user accepts), else the
                # chosen setting's own events when it finds enough of them to mean anything
                src_real = ""
                if reference is not None and len(reference) >= 10 and "dff_max" in reference:
                    real_dff = float(np.nanmedian(reference["dff_max"])); src_real = f"the hand result's {len(reference)} events"
                elif int(best.get("n_real", 0) or 0) >= 10:
                    real_dff = best.get("dff_median", np.nan); src_real = f"the {int(best.get('n_real', 0))} events this setting finds"
                else:
                    real_dff = np.nan
                if np.isfinite(rep_hi) and not np.isfinite(real_dff):
                    lines.append(f"  -> injected classes vs the real signal: the brightest class is reported at {rep_hi:.2f} dF/F; too few real events "
                                 f"here to say whether that is the field's range — a hand-made result of this run would tell (aqua_export_events.m)")
                if np.isfinite(rep_hi) and np.isfinite(real_dff) and rep_hi < 0.8 * real_dff:
                    factor = 1.5 * real_dff / max(rep_hi, 1e-6)
                    sug_amps = tuple(round(min(2.0, a * factor) * 20) / 20 for a in amps)
                    if _is_close(sug_amps, T.AMPLITUDES, 0.3):
                        sug_amps = T.AMPLITUDES
                    class_notes[mv["run"]] = sug_amps
                    lines.append(f"  -> the injected classes probe BELOW the real signal: the brightest class (nominal {a_hi:g} dF/F) is reported at "
                                 f"{rep_hi:.2f} dF/F, {src_real} at {real_dff:.2f} — the recall measures weaker transients than the field has; "
                                 f"plan again with --amplitudes {','.join(f'{a:g}' for a in sug_amps)} (GUI: injected dF/F)")
                elif np.isfinite(rep_hi) and np.isfinite(real_dff):
                    lines.append(f"  -> injected classes vs the real signal: brightest class reported at {rep_hi:.2f} dF/F, {src_real} at {real_dff:.2f} — ok")
            # the fibre mask: when most real events lie off it, it probably misses fibres (v37)
            n_on = best.get("n_real_on"); n_all = best.get("n_real", 0) or 0
            if n_all >= 10 and n_on is not None and np.isfinite(float(n_on)) and float(n_on) / n_all < 0.5:
                lines.append(f"  -> the fibre mask may be incomplete: {100 * (1 - float(n_on) / n_all):.0f}% of the events at this setting lie off it — "
                             f"look at aqua_tune\\{mv['run']}\\fibre_mask_raw.png; a hand-drawn <mouse>\\FOV{mv['fov']}_fibre_mask.tif overrides it")
        else:
            lines.append(f"  -> no setting chosen: {why}")
    # pre and post of one optostim day tuned separately: do they agree? (the day gets ONE setting — tuning pre
    # and post apart would confound the pre-vs-post comparison; the check tells whether the signal changed)
    day_common: dict = {}
    day_matched: dict = {}
    exp_runs: dict = {}
    for mv in man["movies"]:
        if mv.get("paradigm") != "spontaneous":
            exp_runs.setdefault(mv["experiment"], []).append(mv["run"])
    for exp, runs in exp_runs.items():
        if len(runs) < 2:
            continue
        lines.append("")
        lines.append(f"{exp}: pre and post tuned separately — do they agree?")
        conds = {}
        for mv in man["movies"]:
            if mv["run"] in runs:
                conds[mv["run"]] = mv.get("condition") or short_run(mv["run"])
        labels = {run: tuple(chosen[run][k] for k in grid_keys) for run in runs if chosen.get(run)}
        for cand_run, cand in labels.items():
            txt = []
            for run in runs:
                r = next((x for x in per_movie.get(run, []) if tuple(x[k] for k in grid_keys) == cand), None)
                if r:
                    txt.append(f"{conds[run]}: recall {r.get('recall', np.nan):.2f} / false {r.get(T.false_key([r]), np.nan):.0%} / "
                               f"{r.get('rate_real_on', r['rate_real']):.2f} ev/min on fibres")
            lines.append(f"  choice of {conds[cand_run]} ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, cand))}) -> " + "; ".join(txt))
        cand, crow, cwhy, cpass = T.recommend_common({r: per_movie.get(r, []) for r in runs}, grid_keys,
                                                     max_false=max_false, min_on_fibre=min_on_fibre, target_amp=target_amp)
        day_common[exp] = (cand, crow, cpass)
        if len(labels) == len(runs) and len(set(labels.values())) == 1:
            lines.append("  -> AGREE: pre and post recommend the same setting — one setting for the day, as it should be")
        elif len(labels) == len(runs):
            lines.append("  -> DIFFER: pre and post have different optima — the signal changed between pre and post (SNR, focus, "
                         "field) more than the settings tolerate" + (": the matched-sensitivity rule below compensates for it"
                                                                      if day_rule == "matched" else
                                                                      "; the day still gets ONE setting (below), look at the QC report's "
                                                                      "mean images / brightness of these runs before trusting the pre-vs-post change"))
        if cand:
            lines.append(f"  -> {'one common setting for the day' if cpass else 'best common trade-off for the day (rule not met)'}: "
                         + " ".join(f"{k}={v}" for k, v in zip(grid_keys, cand)) + " — "
                         + "; ".join(f"{conds[r]}: recall {x.get('recall', np.nan):.2f} / false {x.get(T.false_key([x]), np.nan):.0%} "
                                     f"({T._fp_rate(x):.2f}/min)" for r, x in crow.items()) + f"   [{cwhy}]")
        else:
            lines.append(f"  -> no common setting for the day: {cwhy}")
        # matched sensitivity (v36): one setting per movie so that pre and post detect the same injected transients with the
        # same recall AND invent false events at the same rate — the fair version of "tune pre and post apart"
        mchoice, mrows, mwhy, mpass, minfo = T.recommend_matched({r: per_movie.get(r, []) for r in runs}, grid_keys, max_false=max_false,
                                                                 min_on_fibre=min_on_fibre, tol_recall=tol_recall, target_amp=target_amp,
                                                                 min_level=min_recall)
        day_matched[exp] = (mchoice, mrows, mpass, minfo, mwhy, conds)
        if mpass:
            if minfo.get("same_setting"):
                lines.append(f"  -> matched sensitivity: the same setting serves both ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, next(iter(mchoice.values()))))}): "
                             + "; ".join(f"{conds[r]}: recall {x.get('recall', np.nan):.2f} / false {T._fp_rate(x):.2f}/min" for r, x in mrows.items())
                             + f"   [{mwhy}]")
            else:
                lines.append("  -> matched sensitivity: ONE SETTING PER CONDITION — "
                             + "; ".join(f"{conds[r]}: {' '.join(f'{k}={v}' for k, v in zip(grid_keys, mchoice[r]))} (recall {x.get('recall', np.nan):.2f}, "
                                         f"false {T._fp_rate(x):.2f}/min, {x.get('rate_real_on', x['rate_real']):.2f} ev/min on fibres)" for r, x in mrows.items())
                             + f" — recalls within {minfo['recall_gap']:.2f}, false-alarm rates within {minfo['fp_gap']:.2f}/min   [{mwhy}]")
        elif minfo.get("weak"):
            lines.append(f"  -> matched sensitivity: {mwhy}" + (" (the day is marked not comparable; the common setting above is written)"
                                                                 if day_rule == "matched" else ""))
        else:
            lines.append(f"  -> matched sensitivity: NOT POSSIBLE — {mwhy}" + (" (the day is marked not comparable; the common setting above is written)"
                                                                                 if day_rule == "matched" else ""))
        if day_rule == "matched" and mpass and not minfo.get("same_setting"):
            gets = "one setting per condition, matched sensitivity (--day-rule common chooses one common setting instead)"
        elif day_rule == "matched" and mpass:
            gets = (f"the same setting for pre and post ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, next(iter(mchoice.values()))))}): "
                    f"matched at recall {minfo.get('level', np.nan):.2f}")
        elif day_rule == "matched":
            gets = ("one common setting, because pre and post could not be matched" + (" at a useful sensitivity" if minfo.get("weak") else "")
                    + " — the day is marked NOT COMPARABLE in the analysis")
        else:
            gets = "one common setting (--day-rule matched for matched sensitivity)"
        lines.append(f"  -> the day gets: {gets}")
        # the day as pictures (v36): pre and post at the same settings side by side, the matched choice framed green
        movies = [movie_ctx[r] for r in runs if movie_ctx.get(r, {}).get("rows")]
        if len(movies) >= 2:
            try:
                dl = day_labels.get(day_of(runs[0]) or "", "")
                exp_short = exp[len(man["mouse"]) + 1:] if exp.startswith(man["mouse"] + "_") else exp
                day_png = tune_dir / f"{exp_short}_why_day.png"
                matched_by_cond = {conds[r]: tuple(st) for r, st in mchoice.items()} if (mpass and mchoice) else None
                if T.why_day_figure(f"{exp_short}" + (f"  ({dl})" if dl else ""), movies, grid_keys, matched_by_cond, cand, max_false, day_png):
                    day_pngs[exp] = day_png
            except Exception as e:
                log(f"  {exp}: the pre/post day figure failed ({e})")
    # one setting for all days of a FOV? compare each day's metrics at the other days' chosen settings
    by_fov: dict = {}
    common: dict = {}
    for mv in man["movies"]:
        by_fov.setdefault(mv["fov"], []).append(mv["run"])
    for fov, runs in by_fov.items():
        if len(runs) < 2:
            continue
        lines.append("")
        lines.append(f"FOV {fov}: one setting for every day? ({len(runs)} tuned movies)")
        labels = {}
        for run in runs:
            ch = chosen.get(run)
            if ch:
                labels[run] = tuple(ch[k] for k in grid_keys)
        if not labels:
            continue
        for cand_run, cand in labels.items():
            txt = []
            for run in runs:
                r = next((x for x in per_movie.get(run, []) if tuple(x[k] for k in grid_keys) == cand), None)
                if r:
                    txt.append(f"{short_run(run)}: recall {r.get('recall', np.nan):.2f} / false {r.get(T.false_key([r]), np.nan):.0%}")
            lines.append(f"  setting of {short_run(cand_run)} ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, cand))}) -> "
                         + "; ".join(txt))
        cand, crow, cwhy, cpass = T.recommend_common({r: per_movie.get(r, []) for r in runs}, grid_keys,
                                                     max_false=max_false, min_on_fibre=min_on_fibre, target_amp=target_amp)
        common[fov] = (cand, crow, cpass)
        if cand:
            lines.append(f"  -> {'common setting for all days' if cpass else 'best common trade-off (rule not met)'}: "
                         + " ".join(f"{k}={v}" for k, v in zip(grid_keys, cand)) + " — "
                         + "; ".join(f"{short_run(r)}: recall {x.get('recall', np.nan):.2f} / false {x.get(T.false_key([x]), np.nan):.0%}"
                                     for r, x in crow.items()) + f"   [{cwhy}]")
        else:
            lines.append(f"  -> no common setting: {cwhy}")
        lines.append("  (one setting for every day keeps the days comparable — the choice for a Treg mouse, whose comparison IS "
                     "across days, and always applied for one; an optostim mouse compares each day within itself and gets per-day "
                     "settings unless 'one common setting per FOV for all days' is ticked / --common-settings)")
    plain_block: list = []
    if rows_all:
        from afferent2p.aqua_io import write_text_free
        write_text_free(tune_dir / "tune_report.csv", pd.DataFrame(rows_all).to_csv(index=False).replace("\r\n", "\n"), log=log)
        try:
            plain_block = draw_figures(mouse_dir, man, per_movie, chosen, passed_rule, verdicts, day_matched, day_common, day_labels,
                                       max_false, min_on_fibre, figure_dpi, log=log)
            lines.append(f"\n-> tune_report.csv, {man['mouse']}_tune.png (catch vs noise), {man['mouse']}_tune_detail.png (the curves)")
        except Exception as e:
            log(f"  figure failed ({e})")
    if why_pngs or day_pngs:
        try:
            write_why_html(tune_dir, man, why_pngs, chosen, passed_rule, max_false, min_on_fibre, day_labels,
                           day_pngs=day_pngs, day_matched=day_matched, day_rule=day_rule, grid_keys=grid_keys, tol_recall=tol_recall)
            lines.append(f"-> tune_why.html (one 'why this setting' picture per tuned movie: {', '.join(p.name for p in why_pngs.values())}"
                         + (f"; pre and post side by side: {', '.join(p.name for p in day_pngs.values())}" if day_pngs else "") + ")")
        except Exception as e:
            log(f"  tune_why.html failed ({e})")
    # the grid verdict of the sweep (v37): where the grid only touches the optimum, the next sweep should go further
    sug_file = tune_dir / "grid_suggestion.json"
    suggest: dict = {}
    flagged = {r: v for r, v in verdicts.items() if v.get("status") in ("too strict", "edge", "too permissive")}
    for v in flagged.values():
        for k, vals in v.get("suggest", {}).items():
            suggest.setdefault(k, set()).update(float(x) for x in vals)
    suggest = {k: sorted(vals) for k, vals in suggest.items()}
    if flagged:
        lines.append("")
        lines.append("GRID VERDICT: the sweep did not bracket the optimum of every movie —")
        for r, v in flagged.items():
            lines.append(f"  {short_run(r)}: {v['status']} (recall {v.get('recall', np.nan):.2f})"
                         + (": sweep " + ", ".join(f"{k}={', '.join(T._fmt(x) for x in vals)}" for k, vals in v.get("suggest", {}).items()) if v.get("suggest") else ""))
        if suggest:
            lines.append("  -> extend the grid and sweep again (only the new combinations are run; finished jobs are kept): "
                         + " ".join(f"--grid {k}={','.join(T._fmt(x) for x in man['grid'].get(k, []) + [v for v in vals if v not in man['grid'].get(k, [])])}"
                                    for k, vals in suggest.items())
                         + "   — or: python aqua_tune.py <mouse> --plan --extend ... (GUI tab 3: 'Extend the grid (suggested)')")
        sug_file.write_text(json.dumps({"written": time.strftime("%Y-%m-%d %H:%M"), "grid": man["grid"], "suggest": suggest,
                                        "movies": {r: {"status": v["status"], "recall": (None if not np.isfinite(v.get("recall", np.nan)) else round(float(v["recall"]), 3)),
                                                       "suggest": v.get("suggest", {})} for r, v in flagged.items()}}, indent=1), encoding="utf-8")
    elif sug_file.exists():
        sug_file.unlink()
    # hand-made results of runs that are not in the plan (v37): worth adding to the plan
    try:
        planned = {mv["run"] for mv in man["movies"]}
        bman = load_manifest(mouse_dir)
        elsewhere = []
        for exp in bman.get("experiments", []):
            for rr in exp.get("runs", []):
                if rr["run"] in planned:
                    continue
                c, _, note = find_reference(relocate(rr["tif"], mouse_dir))
                if c is not None or note:
                    elsewhere.append(short_run(rr["run"]) + ("" if c is not None else " (not exported yet)"))
        if elsewhere or reference_notes:
            lines.append("")
            if elsewhere:
                lines.append(f"HAND-MADE RESULTS exist for runs that are not in the plan: {', '.join(elsewhere)} — add them (runs: "
                             f"{','.join(x.split(' ')[0] for x in elsewhere)}; GUI tab 3 'runs' field, then 1 Plan and the sweep) to see every "
                             "setting against what was accepted by hand")
            for r_, note in reference_notes.items():
                lines.append(f"  {short_run(r_)}: {note}")
    except Exception as e:
        log(f"  (hand-result scan failed: {e})")
    if class_notes:
        pick = max(class_notes.values(), key=lambda v: v[-1])
        lines.append("")
        lines.append(f"INJECTED CLASSES: on {len(class_notes)} movie(s) the injected transients are dimmer than the real events, so every setting "
                     f"looks insensitive and the choice among them means little — plan again with --amplitudes "
                     f"{','.join(f'{a:g}' for a in pick)} (GUI tab 3: injected dF/F), then sweep the inject jobs (the real and null results are kept)")
        result_classes = list(pick)
    else:
        result_classes = None
    result = {"tune_dir": tune_dir, "scored": len(rows_all), "written": [], "not_applied": [], "rule_not_met": [],
              "chosen": {run: (ch["label"] if ch else None) for run, ch in chosen.items()},
              "grid_suggest": suggest, "grid_flagged": {r: v["status"] for r, v in flagged.items()}, "amplitudes_suggest": result_classes}
    if apply or apply_all:
        opts_dir = mouse_dir / "aqua_opts"; opts_dir.mkdir(exist_ok=True)
        lines.append("")
        # what each experiment gets: a spontaneous mouse (its experiment spans the days) always the common choice of
        # its FOV; an optostim day the common choice only when asked (apply_common), else ONE setting for the day —
        # the choice common to its tuned movies (pre and post) when several were tuned, else the movie's own choice
        exp_movies: dict = {}
        for mv in man["movies"]:
            exp_movies.setdefault(mv["experiment"], []).append(mv)
        written = set()
        for exp, mvs in exp_movies.items():
            fov = mvs[0]["fov"]
            spont = mvs[0].get("paradigm") == "spontaneous"
            use_common = (apply_common or spont) and common.get(fov, (None,))[0] is not None
            use_day = (not use_common) and len(mvs) > 1 and day_common.get(exp, (None,))[0] is not None
            # matched sensitivity (v36): one setting per condition of the day when pre and post need different
            # settings to detect with the same measured sensitivity and false-alarm rate
            dm = day_matched.get(exp)
            use_matched = use_day and day_rule == "matched" and dm is not None and dm[2] and dm[0] is not None
            tuning_json = {"experiment": exp, "rule": "matched" if use_matched else ("common" if use_day else ("fov_common" if use_common else "single")),
                           "comparable": True, "conditions": {}, "written": time.strftime("%Y-%m-%d %H:%M")}
            if use_day:
                cand_c, crow_c, cpass_c = day_common[exp]
                tuning_json["common_setting"] = dict(zip(grid_keys, cand_c)) if cand_c else None
                tuning_json["common_rows"] = {dm[5].get(r, r) if dm else r: _tuning_row(x) for r, x in (crow_c or {}).items()}
                if dm is not None:
                    tuning_json["matched_possible"] = bool(dm[2])
                    tuning_json["matched_note"] = dm[4]
                    if not dm[2] and day_rule == "matched":
                        tuning_json["comparable"] = False
                        result.setdefault("not_comparable", []).append(exp)
            if use_matched:
                mchoice, mrows, _, minfo, mwhy, conds = dm
                tuning_json["matched"] = True
                tuning_json["level"] = minfo.get("level"); tuning_json["recall_gap"] = minfo.get("recall_gap"); tuning_json["fp_gap"] = minfo.get("fp_gap")
                same = bool(minfo.get("same_setting"))
                any_changed = False
                stale_generic = opts_dir / f"{exp}_Opt.csv"
                for run, setting in mchoice.items():
                    cond = conds.get(run, run)
                    row = mrows[run]
                    tuning_json["conditions"][cond] = dict(_tuning_row(row), run=run, setting=dict(zip(grid_keys, setting)))
                    if same:
                        continue
                    base = next((m["base"] for m in mvs if m["run"] == run), mvs[0]["base"])
                    vals = dict(base); vals.update(dict(zip(grid_keys, setting)))
                    note = (f"aqua_tune matched-sensitivity setting for the {cond} runs of {exp} ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, setting))}): "
                            f"recall {row.get('recall', np.nan):.2f}, false alarms {T._fp_rate(row):.2f}/min — the other condition is detected with the same "
                            f"sensitivity by its own file")
                    p, changed = T.write_opt_csv(opts_dir / f"{exp}_{cond}_Opt.csv", vals, OPT_TEMPLATE_ROWS, note=note)
                    any_changed = any_changed or changed
                    if changed:
                        written.add(p)
                        lines.append(f"wrote {p.name}: {note}")
                    else:
                        lines.append(f"{p.name}: unchanged (same values as before): {note}")
                if not same:
                    if stale_generic.exists():
                        stale_generic.unlink()
                        lines.append(f"removed {stale_generic.name} (superseded by the per-condition files)")
                    result["written"].append(exp)
                    _write_tuning_json(opts_dir / f"{exp}_tuning.json", tuning_json)
                    continue
                # the same setting serves both with matched sensitivity: one file for the day
                cand = next(iter(mchoice.values())); crow = mrows; cpass = True
                vals = dict(mvs[0]["base"]); vals.update(dict(zip(grid_keys, cand)))
                note = (f"aqua_tune matched-sensitivity setting of the day, the same for pre and post ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, cand))}): "
                        + "; ".join(f"{conds.get(r, r)} recall {x.get('recall', np.nan):.2f}, false {T._fp_rate(x):.2f}/min" for r, x in mrows.items()))
                ok = True
            elif use_day:
                for r, x in (day_common[exp][1] or {}).items():
                    cond = dm[5].get(r, r) if dm else r
                    tuning_json["conditions"][cond] = dict(_tuning_row(x), run=r, setting=dict(zip(grid_keys, day_common[exp][0])))
                cand, crow, cpass = day_common[exp]
                if not cpass and not apply_all:
                    lines.append(f"{exp}: not applied — no setting met the false-alarm rule on both pre and post (use --apply-all / "
                                 "'trade-off' to write the trade-off anyway)")
                    result["not_applied"].append(exp)
                    continue
                vals = dict(mvs[0]["base"]); vals.update(dict(zip(grid_keys, cand)))
                mean_rec = float(np.mean([x.get("recall", np.nan) for x in crow.values()]))
                worst_f = max(x.get(T.false_key([x]), np.nan) for x in crow.values())
                note = (f"aqua_tune setting of the day, common to its tuned runs ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, cand))}): "
                        f"mean recall {mean_rec:.2f}; worst false-alarm share {worst_f:.0%}" + ("" if cpass else "; rule not met"))
                ok = cpass
            elif use_common:
                cand, crow, cpass = common[fov]
                if not cpass and not apply_all:
                    lines.append(f"{exp}: not applied — no setting met the false-alarm rule on every day (use --apply-all to write the trade-off anyway)")
                    result["not_applied"].append(exp)
                    continue
                vals = dict(mvs[0]["base"]); vals.update(dict(zip(grid_keys, cand)))
                mean_rec = float(np.mean([x.get("recall", np.nan) for x in crow.values()]))
                worst_f = max(x.get(T.false_key([x]), np.nan) for x in crow.values())
                note = (f"aqua_tune common setting for FOV {fov}, all days ({' '.join(f'{k}={v}' for k, v in zip(grid_keys, cand))}): "
                        f"mean recall {mean_rec:.2f}; worst false-alarm share {worst_f:.0%}" + ("" if cpass else "; rule not met"))
                ok = cpass
            else:
                mv = mvs[0]
                ch = chosen.get(mv["run"])
                if not ch:
                    result["not_applied"].append(exp)
                    continue
                if not passed_rule.get(mv["run"]) and not apply_all:
                    lines.append(f"{exp}: not applied — no setting met the false-alarm rule (use --apply-all to write the trade-off setting anyway)")
                    result["not_applied"].append(exp)
                    continue
                vals = dict(mv["base"]); vals.update({k: ch[k] for k in grid_keys})
                note = (f"aqua_tune {ch['label']}: recall {ch.get('recall', float('nan')):.2f}; false-alarm share "
                        f"{ch.get(T.false_key([ch]), float('nan')):.0%}" + ("" if passed_rule.get(mv["run"]) else "; rule not met"))
                ok = bool(passed_rule.get(mv["run"]))
            p, changed = T.write_opt_csv(opts_dir / f"{exp}_Opt.csv", vals, OPT_TEMPLATE_ROWS, note=note)
            result["written"].append(exp)
            if not ok:
                result["rule_not_met"].append(exp)
            if changed:
                written.add(p)
                lines.append(f"wrote {p.name}: {note}")
            else:
                lines.append(f"{p.name}: unchanged (same values as before — its runs are not queued again): {note}")
            # one setting for the day: per-condition files of an earlier matched apply would override it — remove them
            for old_pc in opts_dir.glob(f"{exp}_*_Opt.csv"):
                if old_pc.name != f"{exp}_Opt.csv" and old_pc.stem[len(exp) + 1:-4] in {"pre", "post"} | set(tuning_json["conditions"]):
                    old_pc.unlink()
                    lines.append(f"removed {old_pc.name} (one setting for the day now)")
            if use_day or use_matched:
                _write_tuning_json(opts_dir / f"{exp}_tuning.json", tuning_json)
        if written:
            lines.append(f"-> {len(written)} settings file(s) changed: rebuild the batch (tab 4; runs whose settings file is newer than "
                         "their result are queued again) and run MATLAB (tab 5) again to use them")
        elif result["written"]:
            lines.append("-> no settings file changed: nothing needs to be analysed again")
    if delete_movies and rows_all:
        n = 0
        for mv in man["movies"]:
            for k in ("null", "inject"):
                p = Path(mv["variants"].get(k, ""))
                if p.name and p.exists():
                    p.unlink(); n += 1
        lines.append(f"deleted {n} null/injected movie(s) (the truth tables and the event exports stay; --plan writes them again)")
    if plain_block:
        lines = lines[:2] + [""] + plain_block + ["", "DETAILS"] + lines[2:]
    txt = "\n".join(lines) + "\n"
    (tune_dir / "tune_recommendation.txt").write_text(txt, encoding="utf-8")
    log(txt)
    return result


def apply_pick(mouse_dir: Path, exp: str, values: dict, scope: str = "day", cond: str | None = None, tag: str = "", log=print) -> list:
    """Write a setting picked by eye in the picker (v38) as the experiment's Opt file(s) and tuning json. scope 'day' =
    one file for every run of the experiment (aqua_opts\\<exp>_Opt.csv; per-condition files removed); 'condition' =
    aqua_opts\\<exp>_<cond>_Opt.csv for that condition only (the other condition keeps its file, or gets a copy of the
    day's file so that it is not left to the default). The tuning json records rule 'picked' with the measured numbers
    of every movie of the experiment at this setting, so Step 7 still reports recall and false alarms. Returns the
    paths written."""
    mouse_dir = Path(mouse_dir)
    tune_dir = mouse_dir / "aqua_tune"
    man = json.loads((tune_dir / "tune_manifest.json").read_text())
    grid_keys = man["grid_keys"]
    mvs = [mv for mv in man["movies"] if mv["experiment"] == exp]
    if not mvs:
        raise ValueError(f"{exp} is not in the tuning plan")
    per_movie = rows_from_report(tune_dir / "tune_report.csv", grid_keys) if (tune_dir / "tune_report.csv").exists() else {}
    key = tuple(T._fmt(values[k]) if np.isfinite(T._num(values[k])) else str(values[k]) for k in grid_keys)
    rows_at = {}
    for mv in mvs:
        row = next((r for r in per_movie.get(mv["run"], []) if tuple(r.get(k) for k in grid_keys) == key), None)
        if row is not None:
            rows_at[mv["run"]] = row
    label = " ".join(f"{k}={v}" for k, v in zip(grid_keys, key))
    numbers = "; ".join(f"{mv.get('condition') or short_run(mv['run'])}: recall {rows_at[mv['run']].get('recall', np.nan):.2f}, "
                        f"false alarms {T._fp_rate(rows_at[mv['run']]):.2f}/min, {rows_at[mv['run']].get('rate_real_on', rows_at[mv['run']].get('rate_real', np.nan)):.2f} ev/min on fibres"
                        for mv in mvs if mv["run"] in rows_at)
    note = f"picked by eye in the tuning picker ({tag}) for {exp} [{scope}{(' ' + cond) if scope == 'condition' and cond else ''}] ({label}): {numbers}"
    opts_dir = mouse_dir / "aqua_opts"; opts_dir.mkdir(exist_ok=True)
    base = next((mv["base"] for mv in mvs if mv.get("condition") == cond), mvs[0]["base"])
    vals = dict(base); vals.update(dict(zip(grid_keys, key)))
    written = []
    tj_path = opts_dir / f"{exp}_tuning.json"
    tj = {}
    if tj_path.exists():
        try:
            tj = json.loads(tj_path.read_text(encoding="utf-8"))
        except Exception:
            tj = {}
    conds_json = dict(tj.get("conditions") or {})
    if scope == "condition" and cond and mvs[0].get("paradigm") != "spontaneous":
        generic = opts_dir / f"{exp}_Opt.csv"
        for mv in mvs:
            oc = mv.get("condition")
            if oc and oc != cond and not (opts_dir / f"{exp}_{oc}_Opt.csv").exists() and generic.exists():
                import shutil
                shutil.copy2(generic, opts_dir / f"{exp}_{oc}_Opt.csv")
                written.append(opts_dir / f"{exp}_{oc}_Opt.csv")
                log(f"  {exp}_{oc}_Opt.csv: a copy of the day's file, so that the {oc} runs keep their setting")
        p, changed = T.write_opt_csv(opts_dir / f"{exp}_{cond}_Opt.csv", vals, OPT_TEMPLATE_ROWS, note=note, skip_unchanged=False)
        written.append(p)
        if generic.exists():
            generic.unlink()
            log(f"  removed {generic.name} (per-condition files now)")
        for mv in mvs:
            if mv.get("condition") == cond and mv["run"] in rows_at:
                conds_json[cond] = dict(_tuning_row(rows_at[mv["run"]]), run=mv["run"], setting=dict(zip(grid_keys, key)), picked=True)
        rule = "picked"
    else:
        for mv in mvs:
            oc = mv.get("condition")
            pc = opts_dir / f"{exp}_{oc}_Opt.csv" if oc else None
            if pc is not None and pc.exists():
                pc.unlink(); log(f"  removed {pc.name} (one setting for the day now)")
        p, changed = T.write_opt_csv(opts_dir / f"{exp}_Opt.csv", vals, OPT_TEMPLATE_ROWS, note=note, skip_unchanged=False)
        written.append(p)
        conds_json = {}
        for mv in mvs:
            if mv["run"] in rows_at:
                conds_json[mv.get("condition") or short_run(mv["run"])] = dict(_tuning_row(rows_at[mv["run"]]), run=mv["run"],
                                                                                setting=dict(zip(grid_keys, key)), picked=True)
        rule = "picked"
    recs = [c.get("recall") for c in conds_json.values() if c.get("recall") is not None]
    tj.update({"experiment": exp, "rule": rule, "picked_by": "eye", "picked_tag": tag, "comparable": True, "conditions": conds_json,
               "written": time.strftime("%Y-%m-%d %H:%M"), "matched_note": note,
               "level": (min(recs) if recs else None), "recall_gap": ((max(recs) - min(recs)) if len(recs) >= 2 else None)})
    _write_tuning_json(tj_path, tj)
    log(f"wrote {', '.join(p.name for p in written)}: {note}")
    return written


def summary_groups(man: dict, day_labels: dict) -> tuple:
    """[(panel title, experiment, [runs])] — the pre and post movies of a day together, other movies on their own —
    and {run: condition}."""
    groups, conds, by_exp = [], {}, {}
    for mv in man["movies"]:
        conds[mv["run"]] = mv.get("condition") or short_run(mv["run"])
        by_exp.setdefault(mv["experiment"], []).append(mv)
    for exp, mvs in by_exp.items():
        dl = day_labels.get(day_of(mvs[0]["run"]) or "", "")
        exp_short = exp[len(man["mouse"]) + 1:] if exp.startswith(man["mouse"] + "_") else exp
        if mvs[0].get("paradigm") != "spontaneous" and len(mvs) >= 2:
            groups.append((f"{exp_short}" + (f"  ({dl})" if dl else "") + "  — pre and post", exp, [m["run"] for m in mvs]))
        else:
            for m in mvs:
                dl_m = day_labels.get(day_of(m["run"]) or "", "")
                groups.append((f"{short_run(m['run'])}  FOV {m['fov']}" + (f"  ({dl_m})" if dl_m else ""), exp, [m["run"]]))
    return groups, conds


def draw_figures(mouse_dir: Path, man: dict, per_movie: dict, chosen: dict, passed_rule: dict, verdicts: dict, day_matched: dict,
                 day_common: dict, day_labels: dict, max_false: float, min_on_fibre: float, dpi: int, log=print) -> list:
    """<mouse>_tune.png (the summary: catch vs noise, plain words) and <mouse>_tune_detail.png (the curves) — v38.
    Returns the lines for the report text (the plain-words block)."""
    tune_dir = Path(mouse_dir) / "aqua_tune"
    grid_keys = man["grid_keys"]
    groups, conds = summary_groups(man, day_labels)
    classes = (", ".join(f"{float(a):g}" for a in man.get("amplitudes", [])) + " dF/F x "
               + ", ".join(str(int(x)) for x in man.get("sizes_px", [])) + " px x short/long")
    out_lines = []
    try:
        T.tune_summary_figure(groups, per_movie, chosen, passed_rule, grid_keys, max_false, min_on_fibre, day_matched, conds,
                              tune_dir / f"{man['mouse']}_tune.png", dpi=dpi, classes=classes, day_common=day_common)
    except Exception as e:
        log(f"  summary figure failed ({e})")
    try:
        titles = {}
        for mv in man["movies"]:
            dl = day_labels.get(day_of(mv["run"]) or "", "")
            titles[mv["run"]] = f"{short_run(mv['run'])}  FOV {mv['fov']}" + (f"  ({dl})" if dl else "") \
                + ("" if passed_rule.get(mv["run"]) or not chosen.get(mv["run"]) else "  [rule not met]")
        notes = {}
        for run, v in verdicts.items():
            if v.get("status") in ("too strict", "edge", "too permissive", "weak"):
                sug = ", ".join(f"{k}={T._fmt(x)}" for k, vals in v.get("suggest", {}).items() for x in vals)
                notes[run] = f"grid: {v['status']}" + (f" — sweep {sug}" if sug else "")
        T.tune_figure(per_movie, tune_dir / f"{man['mouse']}_tune_detail.png", grid_keys,
                      [s_ for mv in man["movies"] for s_ in mv["settings"]], chosen, titles=titles, dpi=dpi, notes=notes)
    except Exception as e:
        log(f"  detail figure failed ({e})")
    # the plain-words block
    out_lines.append("IN PLAIN WORDS")
    for title, exp, runs in groups:
        rows_all = [r for run in runs for r in per_movie.get(run, [])]
        if not rows_all:
            out_lines.append(f"  {title}: not scored yet"); continue
        fk = T.false_key(rows_all)
        out_lines.append(f"  {title}")
        dm = day_matched.get(exp) if len(runs) >= 2 else None
        if dm is not None:
            dc = (day_common or {}).get(exp)
            for w in T.plain_day_words(dm[5], chosen, dm[0], dm[1], dm[2], dm[3], dm[4], grid_keys, fk,
                                       common=dc[0] if dc else None, crow=dc[1] if dc else None):
                out_lines.append(f"    {w}")
        for run in runs:
            if per_movie.get(run) and chosen.get(run) is not None:
                w = T.plain_movie_words(per_movie[run], chosen[run], passed_rule.get(run, False), grid_keys, fk, max_false, min_on_fibre)
                out_lines.append(f"    {conds.get(run, '') or run}: {w[0]}")
                out_lines += [f"      {x}" for x in w[1:]]
    out_lines.append(f"  (pictures: {man['mouse']}_tune.png — catch vs noise per setting; {man['mouse']}_tune_detail.png — the curves; "
                     "tune_why.html — the events on the image)")
    return out_lines


def figure_only(mouse_dir: Path, dpi: int = 260, max_false: float = 0.10, min_on_fibre: float = 0.80, tol_recall: float = 0.05,
                min_recall: float = 0.3, log=print):
    """<mouse>_tune.png and <mouse>_tune_detail.png again from tune_report.csv, without rescoring (v38: the GUI's
    'Redraw at 2x resolution'; the day rules are recomputed from the scored rows)."""
    tune_dir = Path(mouse_dir) / "aqua_tune"
    rep = tune_dir / "tune_report.csv"
    if not rep.exists():
        raise SystemExit(f"{rep} not found — run --report first")
    man = json.loads((tune_dir / "tune_manifest.json").read_text())
    grid_keys = man["grid_keys"]
    per_movie = rows_from_report(rep, grid_keys)
    order = [mv["run"] for mv in man["movies"] if mv["run"] in per_movie]
    per_movie = {r: per_movie[r] for r in order} if order else per_movie
    day_labels = read_day_labels(Path(mouse_dir))
    chosen, passed_rule, verdicts = {}, {}, {}
    for run, rows in per_movie.items():
        best, why, passed = T.recommend(rows, max_false=max_false, min_on_fibre=min_on_fibre)
        chosen[run] = best; passed_rule[run] = passed
        if best:
            verdicts[run] = T.grid_verdict(rows, best, passed, grid_keys, max_false=max_false)
    day_matched, day_common = {}, {}
    exp_runs: dict = {}
    for mv in man["movies"]:
        if mv.get("paradigm") != "spontaneous" and mv["run"] in per_movie:
            exp_runs.setdefault(mv["experiment"], []).append(mv["run"])
    for exp, runs in exp_runs.items():
        if len(runs) < 2:
            continue
        conds = {mv["run"]: (mv.get("condition") or short_run(mv["run"])) for mv in man["movies"] if mv["run"] in runs}
        cand, crow, cwhy, cpass = T.recommend_common({r: per_movie[r] for r in runs}, grid_keys, max_false=max_false, min_on_fibre=min_on_fibre)
        day_common[exp] = (cand, crow, cpass)
        mchoice, mrows, mwhy, mpass, minfo = T.recommend_matched({r: per_movie[r] for r in runs}, grid_keys, max_false=max_false,
                                                                 min_on_fibre=min_on_fibre, tol_recall=tol_recall, min_level=min_recall)
        day_matched[exp] = (mchoice, mrows, mpass, minfo, mwhy, conds)
    draw_figures(mouse_dir, man, per_movie, chosen, passed_rule, verdicts, day_matched, day_common, day_labels, max_false, min_on_fibre, dpi, log=log)
    out = tune_dir / f"{man['mouse']}_tune.png"
    log(f"-> {out} and {man['mouse']}_tune_detail.png at {dpi} dpi")
    return out


def rows_from_report(rep: Path, grid_keys: list) -> dict:
    """tune_report.csv -> {run: [row dicts]} with the grid values as the strings the scorer uses."""
    import pandas as pd
    df = pd.read_csv(rep)
    per_movie: dict = {}
    for run, g in df.groupby("run", sort=False):
        rows = []
        for _, r in g.iterrows():
            d = {k: (v if not isinstance(v, float) or np.isfinite(v) else np.nan) for k, v in r.items()}
            for k in grid_keys:
                d[k] = T._fmt(d[k]) if np.isfinite(T._num(d[k])) else str(d[k])
            rows.append(d)
        per_movie[run] = rows
    return per_movie


def sweep_progress(mouse_dir: Path) -> dict:
    """How far the MATLAB sweep is, from the planned folders: {done, failed,
    pending, total, per_movie: [str], min_per_job: float|None, log_tail: str,
    head: str}; raises SystemExit when there is no plan."""
    tune_dir = Path(mouse_dir) / "aqua_tune"
    mp = tune_dir / "tune_manifest.json"
    if not mp.exists():
        raise SystemExit("no plan yet — run --plan (tab 3, 1 Plan) first")
    man = json.loads(mp.read_text())
    lines, done, failed, pending, times = [], 0, 0, 0, []
    for mv in man["movies"]:
        mdir = Path(mv["dir"])
        if not mdir.exists():
            mdir = tune_dir / mv["run"]
        d = f = p = 0
        for st in mv["settings"]:
            for vname, vpath in mv["variants"].items():
                out = mdir / f"s{st['setting']:02d}" / vname
                ev = out / f"{Path(vpath).stem}_AQuA2_events.csv"
                if ev.exists():
                    d += 1; times.append(ev.stat().st_mtime)
                elif (out / "FAILED.txt").exists():
                    f += 1
                else:
                    p += 1
        done += d; failed += f; pending += p
        lines.append(f"  {mv['run']}: {d} done, {f} failed, {p} pending of {d + f + p}")
    total = done + failed + pending
    per_job = None
    if len(times) >= 3:
        gaps = np.diff(np.sort(times)) / 60.0
        gaps = gaps[(gaps > 0.2) & (gaps < 30)]
        if gaps.size >= 2:
            per_job = float(np.median(gaps))
    head = f"{done} of {total} jobs done, {failed} failed, {pending} pending"
    if per_job and pending:
        head += f" — about {per_job:.1f} min per job so far, ~{pending * per_job / 60:.1f} h left"
    elif pending == 0 and total:
        head += " — all done: press 3 Report" if not failed else " — finished (failed jobs are reported without a score)"
    tail = ""
    log_p = tune_dir / "tune_log.txt"
    if log_p.exists():
        try:
            last = [ln for ln in log_p.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()][-3:]
            tail = "MATLAB log, last lines:\n  " + "\n  ".join(last)
        except Exception:
            pass
    return {"done": done, "failed": failed, "pending": pending, "total": total, "per_movie": lines,
            "min_per_job": per_job, "log_tail": tail, "head": head}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mouse", type=Path, help="mouse folder (with aqua_batch\\manifest.json from aqua_batch.py)")
    ap.add_argument("--plan", action="store_true", help="write the movies, the grid and the MATLAB script")
    ap.add_argument("--report", action="store_true", help="score the finished jobs and recommend a setting per movie")
    ap.add_argument("--progress", action="store_true", help="how many sweep jobs are done / failed / pending")
    ap.add_argument("--figure-only", action="store_true", help="draw <mouse>_tune.png again from tune_report.csv (no rescoring), e.g. at another --dpi")
    ap.add_argument("--dpi", type=int, default=None, help="resolution of <mouse>_tune.png (report: 150; --figure-only: 260)")
    ap.add_argument("--apply", action="store_true", help="with --report: write the chosen settings as aqua_opts\\<experiment>_Opt.csv")
    ap.add_argument("--apply-all", action="store_true", help="like --apply, but also write the trade-off settings of movies where no setting met the false-alarm rule")
    ap.add_argument("--delete-movies", action="store_true", help="with --report: delete the null and injected movies (~0.8 GB per tuned run) once scored")
    ap.add_argument("--per-day", dest="per_condition", action="store_const", const="yes", default="auto",
                    help="plan: one movie per condition of every experiment (the default: Treg = every day; optostim = pre and post of every day)")
    ap.add_argument("--baseline-only", dest="per_condition", action="store_const", const="no", help="plan: only the first baseline run per experiment")
    ap.add_argument("--day-rule", choices=("matched", "common"), default="matched",
                    help="optostim day with a pre and a post movie tuned: 'matched' (default, v36) = one setting per condition "
                         "chosen so that pre and post detect the injected transients with the same recall and invent false events "
                         "at the same rate (the same setting when it serves both); 'common' = one setting for the day")
    ap.add_argument("--tol-recall", type=float, default=0.05, help="matched rule: largest allowed recall gap between pre and post")
    ap.add_argument("--min-recall", type=float, default=0.3, help="matched rule (v37): a pair matched below this common sensitivity is 'too weak' "
                                                                  "(the day is marked not comparable; the grid verdict says how to extend the grid)")
    ap.add_argument("--extend", action="store_true", help="plan (v37): add the values the last report suggested (aqua_tune\\grid_suggestion.json) "
                                                          "to the grid; finished jobs are kept, only the new combinations are swept")
    ap.add_argument("--common-settings", dest="apply_common", action="store_true", default=False,
                    help="apply: one setting per FOV for ALL days of an optostim mouse too (a spontaneous mouse always gets one)")
    ap.add_argument("--per-day-settings", dest="apply_common", action="store_false",
                    help="apply: each optostim day its own setting (the default since v34; kept for old command lines)")
    ap.add_argument("--runs", default=None, help="comma-separated run name parts to tune (default: the first baseline run of every day and FOV)")
    ap.add_argument("--days", default=None, help="comma-separated dates (e.g. 260911,260921): only these days' experiments")
    ap.add_argument("--one-per-fov", action="store_true", help="only the first day of each FOV")
    ap.add_argument("--grid", action="append", default=None, metavar="VAR=v1,v2,...",
                    help="a swept variable and its values (repeatable); default thrARScl=1.5,2,2.5,3 and minSize=20,50,100 (v37)")
    ap.add_argument("--aqua-dir", default=None, help="AQuA2 folder (default: from the batch manifest)")
    ap.add_argument("--um-per-px", type=float, default=0.792)
    ap.add_argument("--n-inject", type=int, default=72, help="injected transients per movie (default 72: 4 per amplitude x size x duration class)")
    ap.add_argument("--amplitudes", default=None, help="injected amplitudes, dF/F of the pixel's baseline (default 0.1,0.2,0.4; AQuA2 reports roughly half, "
                                                       "because it grows the footprint — 0.3,0.6,1.0 brackets events it reports at 0.15-0.5)")
    ap.add_argument("--sizes", default=None, help="injected footprint sizes in pixels (default 30,80,200 since v37; an existing plan keeps its own "
                                                   "classes unless given — a changed list means new injected movies and every inject job again)")
    ap.add_argument("--no-null", action="store_true", help="skip the null movie (no false-alarm measurement)")
    ap.add_argument("--no-inject", action="store_true", help="skip the injected movie (no recall measurement)")
    ap.add_argument("--as-root", default=None, help=r"write MATLAB paths with the mouse folder at this root (e.g. Z:\2photon\...\<mouse>)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0, help="MATLAB jobs at once: 0 = let MATLAB choose from cores and memory (default), "
                                                          "1 = one at a time, N = that many (parallel pool, or N MATLAB processes via the .bat)")
    ap.add_argument("--max-false", type=float, default=0.10, help="accepted false-alarm share (null events / real events; default 0.10)")
    ap.add_argument("--min-on-fibre", type=float, default=0.80, help="required on-fibre fraction of the real movie's events (default 0.80)")
    ap.add_argument("--target-amp", type=float, default=None, help="rank by recall at this injected amplitude (default: all amplitudes)")
    a = ap.parse_args()
    if a.apply_all:
        a.report = True
    if not a.plan and not a.report and not a.progress and not a.figure_only:
        ap.error("say what to do: --plan, --progress, --report or --figure-only")
    if a.figure_only:
        figure_only(a.mouse, dpi=a.dpi or 260, max_false=a.max_false, min_on_fibre=a.min_on_fibre)
    if a.progress:
        pr = sweep_progress(a.mouse)
        print(pr["head"]); print("\n".join(pr["per_movie"]))
        if pr["log_tail"]:
            print(pr["log_tail"])
    if a.plan:
        plan(a.mouse, a.aqua_dir, T.parse_grid(a.grid), [s.strip() for s in a.runs.split(",")] if a.runs else None,
             a.one_per_fov, a.um_per_px, a.n_inject, not a.no_null, not a.no_inject, a.as_root, a.seed,
             amplitudes=[float(x) for x in a.amplitudes.split(",") if x.strip()] if a.amplitudes else None,
             sizes=[int(float(x)) for x in a.sizes.split(",") if x.strip()] if a.sizes else None, per_condition=a.per_condition,
             days_filter=[s.strip() for s in a.days.split(",")] if a.days else None, n_workers=a.workers, extend=a.extend)
    if a.report:
        report(a.mouse, a.max_false, a.min_on_fibre, a.target_amp, a.apply, apply_all=a.apply_all, delete_movies=a.delete_movies,
               apply_common=a.apply_common, day_rule=a.day_rule, tol_recall=a.tol_recall, min_recall=a.min_recall, figure_dpi=a.dpi or 150)


if __name__ == "__main__":
    main()
