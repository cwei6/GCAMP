"""Experiment planner: turn a mouse folder's recordings into FOV-level experiments.

File naming (Scanbox): <mouse>_<yymmdd>_<nnn>.sbx, run numbers sequential
within a day. Two lab paradigms:

  spontaneous  (animal name WITHOUT 'RearChr'):
      each day: n_fov recordings (001 = FOV1, 002 = FOV2, ...); recorded on
      several days (e.g. day -1, 5, 7, 9). One experiment per FOV, its runs =
      that FOV's recording on each day; conditions = the dates (earliest =
      baseline). Same FOV across days is registered to one reference.

  optostim     (animal name WITH 'RearChr'):
      one day: FOV1, FOV2, FOV1, FOV2 (pre) -> 1 h stimulation (no recording)
      -> FOV1, FOV2, FOV1, FOV2 (post). Runs 001..008 in order, so runs in the
      first half are pre, second half post, FOV = position within the block
      modulo n_fov. One experiment per (day, FOV) with conditions pre / post.

The planner only reads file names (and optionally frame counts); nothing is
moved. Each experiment gets its own output folder:
      <mouse folder>\\afferent2p_out\\<experiment name>\\

  run_labels.csv (v44, optional, in the mouse folder): when the run number
  does not say which FOV a recording is (a lost FOV, a shifted plane that was
  recaptured with the next number, a second injection with pre + post runs on
  one day), one row per run that deviates:
      run,fov,condition,use,note
      NavTIGRE-2609-1_261010_006,1,d261010_post,1,recapture of FOV1
      NavTIGRE-2609-1_261010_004,,,0,shifted plane
  run = the run name (or just <yymmdd>_<nnn>); fov = the FOV it belongs to
  (empty = unchanged); condition = its label (empty = the default: d<date> for
  a spontaneous mouse, pre / post for optostim; a same-day pre / post pair is
  written d<date>_pre / d<date>_post, the date first so the day's label still
  applies); use = 0 drops the run from every step. Rows for runs that follow
  the rule are not needed (write_run_labels() writes a complete template).
"""
from __future__ import annotations
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

RUN_RE = re.compile(r"^(?P<mouse>.+?)_(?P<date>\d{6})_(?P<run>\d{3})$")

# folders created by this pipeline / earlier routes: never scan them for raw runs
SKIP_DIRS = {"afferent2p_out", "aqua_batch", "aqua_analysis", "aqua_opts",
             "aqua_input", "aqua_analysis_manual", "suite2p", "raw_tifs",
             "qc", "qc_raw", "analysis", "analysis_manual", "behavior", "reg_compare", "aqua_tune"}
PARK_PREFIX = "_previous_"          # aqua_reset.py parks old outputs in <mouse>\_previous_<date>\ — never a place to scan


def skip_parts(parts) -> bool:
    """True when a path (its parts relative to the mouse folder) lies in a pipeline output folder or a parked one."""
    return any(str(x).lower() in SKIP_DIRS or str(x).lower().startswith(PARK_PREFIX) for x in parts)


@dataclass
class Experiment:
    name: str
    data_dir: Path
    runs: list[tuple[str, str]]          # (sbx file name, condition label)
    paradigm: str
    fov: int
    output_dir: Path
    notes: list[str] = field(default_factory=list)


def parse_runs(mouse_dir: Path):
    """[(sbx_path, mouse, date, run)] for every .sbx under the mouse folder.

    Pipeline output folders are skipped, and duplicated run names (e.g. a copy
    of a recording in a 'test_run' folder) are reported and used once — the copy
    whose parent folder is named like the run wins, since that is the layout the
    microscope writes.
    """
    found = {}
    cands = [p for p in sorted(mouse_dir.rglob("*.sbx"), key=lambda p: p.name) if not skip_parts(p.relative_to(mouse_dir).parts[:-1])]

    def _size(p):
        try:
            return p.stat().st_size
        except OSError:
            return 0
    sizes = sorted(_size(p) for p in cands)
    typical = sizes[len(sizes) // 2] if sizes else 0
    for p in cands:
        m = RUN_RE.match(p.stem)
        if not m:
            print(f"  plan: skipping {p.name} (name is not <mouse>_<yymmdd>_<nnn>)")
            continue
        if not p.with_suffix(".mat").exists() and typical > 0 and _size(p) < 0.5 * typical:
            # Scanbox writes the info .mat when a recording ends properly: a .sbx without it that is much smaller
            # than the others is an aborted recording (261016_005 of 2608-1: 434 MB / 266 frames, re-recorded as 051)
            print(f"  plan: skipping {p.name} — no {p.with_suffix('.mat').name} next to it and only {_size(p) / 1e6:.0f} MB "
                  f"(the others ~{typical / 1e9:.1f} GB): an aborted recording, not a run — rename or remove the folder if it is one")
            continue
        prev = found.get(p.stem)
        if prev is not None:
            better = p if p.parent.name == p.stem else prev
            skipped = prev if better is p else p
            print(f"  plan: {p.stem} found twice — using {better.parent.name}\\, "
                  f"ignoring the copy in {skipped.parent.name}\\")
            found[p.stem] = better
            continue
        found[p.stem] = p
    out = []
    for stem, p in sorted(found.items()):
        m = RUN_RE.match(stem)
        out.append((p, m["mouse"], m["date"], int(m["run"])))
    return out


def guess_paradigm(mouse_name: str) -> str:
    return "optostim" if "rearchr" in mouse_name.lower() else "spontaneous"


RUN_LABELS = "run_labels.csv"
RUN_KEY_RE = re.compile(r"(\d{6})_(\d{3})")


def run_key(name: str) -> str | None:
    """'<yymmdd>_<nnn>' out of a run name / stem / file name (None when it has none)."""
    m = RUN_KEY_RE.search(str(name))
    return f"{m.group(1)}_{m.group(2)}" if m else None


def read_run_labels(mouse_dir: Path) -> dict:
    """{'<yymmdd>_<nnn>': {fov: int|None, condition: str|None, use: bool, note: str}} from
    <mouse>\\run_labels.csv — empty when there is no file. Bad rows are reported and skipped."""
    p = Path(mouse_dir) / RUN_LABELS
    out = {}
    if not p.exists():
        return out
    import csv
    try:
        with open(p, newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
    except Exception as e:
        print(f"  plan: {p.name} could not be read ({e}) — ignored")
        return out
    for i, r in enumerate(rows, start=2):
        r = {str(k).strip().lower(): (v.strip() if isinstance(v, str) else "") for k, v in r.items() if k is not None}
        key = run_key(r.get("run", ""))
        if not key:
            if any(r.values()):
                print(f"  plan: {p.name} line {i}: no run name in '{r.get('run', '')}' — skipped")
            continue
        fov = None
        if r.get("fov"):
            try:
                fov = int(float(r["fov"]))
            except ValueError:
                print(f"  plan: {p.name} line {i}: fov '{r['fov']}' is not a number — skipped"); continue
        use_txt = r.get("use", "")
        use = use_txt.lower() not in ("0", "no", "false", "n", "exclude", "x") if use_txt else True
        cond = r.get("condition") or None
        if cond and not re.match(r"^[A-Za-z0-9_+-]+$", cond):
            print(f"  plan: {p.name} line {i}: condition '{cond}' has characters other than letters, digits, _ + - — skipped"); continue
        out[key] = {"fov": fov, "condition": cond, "use": use, "note": r.get("note", "")}
    return out


def write_run_labels(mouse_dir: Path, paradigm: str = "auto", n_fov: int | None = None,
                     stim_after_run: int | None = None) -> Path:
    """Create <mouse>\\run_labels.csv with every run and its current FOV / condition (the rule's
    or the file's own), or complete an existing file with the runs it does not list. Nothing
    else is changed, so the file can be edited in Excel / Notepad and re-completed after a new day."""
    mouse_dir = Path(mouse_dir)
    p = mouse_dir / RUN_LABELS
    existing = read_run_labels(mouse_dir)
    exps = plan_mouse(mouse_dir, paradigm=paradigm, n_fov=n_fov, stim_after_run=stim_after_run)
    seen = {}
    for e in exps:
        for nm, cond in e.runs:
            k = run_key(nm)
            if k:
                seen[k] = (Path(nm).stem, e.fov, cond)
    # excluded runs are in no experiment: keep their rows
    import csv
    rows = []
    if p.exists():
        with open(p, newline="", encoding="utf-8-sig") as fh:
            rows = [dict(r) for r in csv.DictReader(fh)]
    listed = {run_key(r.get("run", "")) for r in rows}
    added = 0
    for (pth, _, date, run) in parse_runs(mouse_dir):
        k = f"{date}_{run:03d}"
        if k in listed:
            continue
        if k in seen:
            stem, fov, cond = seen[k]
            rows.append({"run": stem, "fov": fov, "condition": cond, "use": 1, "note": ""})
        else:
            rows.append({"run": pth.stem, "fov": "", "condition": "", "use": 0, "note": "not in any experiment"})
        added += 1
    rows.sort(key=lambda r: run_key(r.get("run", "")) or "")
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["run", "fov", "condition", "use", "note"], extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})
    print(f"  plan: {p.name}: {len(rows)} runs listed ({added} added; {len(existing)} were there)")
    return p


def plan_mouse(mouse_dir: Path, paradigm: str = "auto", n_fov: int | None = None,
               stim_after_run: int | None = None) -> list[Experiment]:
    mouse_dir = Path(mouse_dir)
    recs = parse_runs(mouse_dir)
    if not recs:
        return []
    mouse = recs[0][1]
    if paradigm == "auto":
        paradigm = guess_paradigm(mouse)
    by_date = defaultdict(list)
    for p, _, date, run in recs:
        by_date[date].append((run, p))
    for d in by_date:
        by_date[d].sort()
        nums = [r for r, _ in by_date[d]]
        if paradigm == "optostim" and nums != list(range(1, len(nums) + 1)):
            # run numbers out of sequence (a re-take after an aborted run gets the next free number,
            # e.g. 051 standing in for 005): the recording order is the file order, not the number order
            try:
                by_date[d].sort(key=lambda rp: rp[1].stat().st_mtime)
                print(f"  plan: {mouse} {d}: run numbers {nums} are not 1..{len(nums)} — ordered by recording time instead: "
                      + ", ".join(f"{r:03d}" for r, _ in by_date[d]))
            except OSError:
                pass
    dates = sorted(by_date)
    exps: list[Experiment] = []
    labels = read_run_labels(mouse_dir)
    applied: list[str] = []            # what run_labels.csv changed, for the notes

    def override(date: str, run: int, fov: int, cond: str):
        """(fov, cond, keep) after run_labels.csv."""
        lab = labels.get(f"{date}_{run:03d}")
        if lab is None:
            return fov, cond, True
        why = f" ({lab['note']})" if lab.get("note") else ""
        if not lab["use"]:
            applied.append(f"{date}_{run:03d} excluded{why}")
            return fov, cond, False
        f2 = lab["fov"] if lab["fov"] is not None else fov
        c2 = lab["condition"] or cond
        if (f2, c2) != (fov, cond):
            applied.append(f"{date}_{run:03d} -> FOV{f2} {c2}{why}")
        return f2, c2, True

    if paradigm == "spontaneous":
        # FOV k = k-th recording of each day, unless run_labels.csv says otherwise
        counts = [len(by_date[d]) for d in dates]
        if len(set(counts)) > 1 and not labels:
            print(f"  plan: {mouse}: days have different recording counts {dict(zip(dates, counts))}"
                  " — FOVs matched by run number; missing days are simply absent (run_labels.csv can say which run is which FOV)")
        per_fov: dict[int, list] = defaultdict(list)
        for d in dates:
            for run, p in by_date[d]:
                fov, cond, keep = override(d, run, run, f"d{d}")
                if keep:
                    per_fov[fov].append((d, run, p, cond))
        nf = n_fov or (max(per_fov) if per_fov else 0)
        for k in range(1, nf + 1):
            recs_k = sorted(per_fov.get(k, []), key=lambda t: (t[0], t[1]))
            runs = [(p.name, cond) for _, _, p, cond in recs_k]
            if not runs:
                continue
            name = f"{mouse}_FOV{k}_days"
            e = Experiment(name=name, data_dir=mouse_dir, runs=runs, paradigm=paradigm,
                           fov=k, output_dir=mouse_dir / "afferent2p_out" / name)
            conds = list(dict.fromkeys(c for _, c in runs))
            e.notes.append(f"baseline = {conds[0]} (this FOV's earliest condition); conditions = {', '.join(conds)}")
            if len(runs) == 1:
                e.notes.append("only one day found — single-condition analysis")
            own = [a for a in applied if any(a.startswith(f"{d}_{r:03d}") for d, r, _, _ in recs_k)]
            if own:
                e.notes.append("run_labels.csv: " + "; ".join(own))
            exps.append(e)
        dropped = [a for a in applied if " excluded" in a]
        if dropped and exps:
            exps[0].notes.append("run_labels.csv excluded (whole mouse): " + "; ".join(a.replace(" excluded", "") for a in dropped))

    elif paradigm == "optostim":
        nf = n_fov or 2
        for d in dates:
            runs_d = by_date[d]
            n = len(runs_d)
            # protocol: FOV1, FOV2, FOV1, FOV2 (2 x nf runs) -> stim -> same again,
            # so runs <= split are pre. A short day (fewer than 4 x nf runs) is
            # taken as "the pre block, possibly followed by a partial post block"
            # rather than split in half — e.g. 2 runs = FOV1 pre + FOV2 pre, no post.
            if stim_after_run:
                split = stim_after_run
            elif n <= 4 * nf:
                split = min(n, 2 * nf)
            else:
                split = n // 2
            if stim_after_run is None and n != 4 * nf and not labels:
                print(f"  plan: {mouse} {d}: {n} runs (protocol expects {4 * nf} = 2 x {nf} FOVs pre + post) — "
                      f"assuming runs 1-{split} pre, {'runs ' + str(split + 1) + '-' + str(n) + ' post' if split < n else 'NO post runs'}; "
                      f"set stim_after_run (or run_labels.csv) if that is wrong")
            per_fov = defaultdict(list)
            pre_block = [(r, p) for r, p in runs_d if r <= split]
            post_block = [(r, p) for r, p in runs_d if r > split]
            for block, label in ((pre_block, "pre"), (post_block, "post")):
                for i, (r, p) in enumerate(block):
                    fov, cond, keep = override(d, r, i % nf + 1, label)
                    if keep:
                        per_fov[fov].append((p.name, cond, r))
            for k in sorted(per_fov):
                runs = [(nm, lab) for nm, lab, _ in per_fov[k]]
                name = f"{mouse}_{d}_FOV{k}"
                e = Experiment(name=name, data_dir=mouse_dir, runs=runs, paradigm=paradigm,
                               fov=k, output_dir=mouse_dir / "afferent2p_out" / name)
                e.notes.append("runs " + ", ".join(f"{r:03d}={lab}" for _, lab, r in per_fov[k]))
                if not any(l == "post" for l in (lab for _, lab, _ in per_fov[k])):
                    e.notes.append("no post runs yet (stimulation not done?) — single-condition")
                own = [a for a in applied if any(a.startswith(f"{d}_{r:03d}") for _, _, r in per_fov[k])]
                if own:
                    e.notes.append("run_labels.csv: " + "; ".join(own))
                exps.append(e)
    else:
        raise ValueError(f"unknown paradigm {paradigm!r}")
    if applied:
        print(f"  plan: {mouse}: run_labels.csv applied — " + "; ".join(applied))
    return exps


def describe(exps: list[Experiment]) -> str:
    lines = []
    for e in exps:
        lines.append(f"{e.name}  [{e.paradigm}, FOV {e.fov}]  -> {e.output_dir}")
        for nm, lab in e.runs:
            lines.append(f"    {lab:>8}  {nm}")
        for n in e.notes:
            lines.append(f"    note: {n}")
    return "\n".join(lines)


def experiment_config(e: Experiment, template: dict) -> dict:
    """Fill a config template (everything but `experiment:`) for this experiment."""
    cfg = dict(template)
    cfg["experiment"] = {
        "name": e.name,
        "data_dir": str(e.data_dir),
        "runs": [{"sbx": nm, "condition": lab} for nm, lab in e.runs],
        "output_dir": str(e.output_dir),
        "paradigm": e.paradigm,
        "fov": e.fov,
    }
    return cfg
