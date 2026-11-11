# Treg-paradigm animals through the pipeline (Nov 3–5, 2026) — plan, state, how to continue

The user (Nov 3): "help me process the Treg data" → all five Treg-paradigm animals on samson, full route to the curation point (tuning sweep with one common setting per FOV across days, AQuA2, CFU at overlap 0.3), run remotely where possible. Progress of the wider pipeline work: `claude/v34-progress.md`; the paradigm notes of Oct 21: v34-progress "Treg paradigm".

## The animals (D:\2photon\Chao\GCAMP_affrerents)

| animal | days (label) | runs / day | state before Nov 3 |
|---|---|---|---|
| NavTIGRE-2609-1 (control) | 260925 Baseline, 261001 D5, 261003 D7, 261005 D9, 261010 D14 (second injection: pre + post) | 3, 3, 3, 3, 6 | unpacked + frame-checked only |
| NavTIGRETreg-2609-1 | same days | 3, 3, 3, 3, 5 | unpacked + frame-checked only |
| NavTIGRETreg-2609-2 | 260925 Baseline, 261001 D5, 261003 D7, 261005 D9 | 4 (FOV1–4, regular) | unpacked + frame-checked only |
| NavTIGRETreg-2608-1 | 260911 Baseline, DTx 260912, 260917 D5, 260919 D7, 260921 D9 | 3, 4, 4, 3 | AQuA2 + CFU of Oct 11 (old pipeline, old-class sweep); FOV Opt files in the root |
| NavTIGRE-2608-1 | 260911 Baseline only (no injection yet) | 3 | aqua_batch of Oct 11, nothing else |

Run → FOV exceptions (from the user, Nov 3), now in each animal's `run_labels.csv`:
- NavTIGRE-2609-1 261010: 001 FOV1 pre, 002 FOV2 pre, 003 lost FOV (bone growth) → out, 004 FOV1 post but shifted plane → out, 005 FOV2 post, 006 = recapture of FOV1 = FOV1 post. FOV3 has no D14.
- NavTIGRETreg-2609-1 261010: 001 FOV1 pre, 002 FOV2 pre, 003 first FOV1 post attempt, shifted → out, 004 FOV2 post, 005 FOV1 post; FOV3 lost that day.
- NavTIGRETreg-2608-1: run 004 of 260917 and 260919 is a field that matches none of the baseline FOVs (mean-image cross-correlation ≤ 0.31 vs 0.55–0.75 for the matched runs) → left out (user's choice: a FOV without a baseline day).
- The other animals/days follow the rule (run k = FOV k).
`day_labels.csv` per animal as in the table; the D14 pre/post conditions are `d261010_pre` / `d261010_post` with their own label rows.

## What was built — v44 (delivered Nov 3: working copy, `afferent2p_v44\`, zip; md5-verified)

- `afferent2p/plan.py`: `read_run_labels`, `write_run_labels` (creates / completes `<mouse>\run_labels.csv` with every run), `run_key`; `plan_mouse` applies the file in both paradigms (fov, condition, use = 0 → the run is in no experiment; notes list what changed). Conditions of a same-day pair: `d<date>_pre` / `d<date>_post` (date first → the day's label applies, recording order kept; the analysis takes the baseline = first condition in run order).
- `afferent2p/timecourse.py`: `cond_day`, `cond_label` (exact key `261010_post` in day_labels, else the day's label + suffix); `report._cond_name` and `general.py` labels use it.
- GUI tab 1: block "Which run is which FOV / condition": Create / complete run_labels.csv, Open run_labels.csv, Open day_labels.csv, Show the plan (log); title v44. `plan.py --run-labels`. pipeline_state: batch stale when run_labels.csv is newer; aqua_reset keeps run_labels.csv. GUIDE "Which run is which FOV: run_labels.csv and day_labels.csv (v44)", README "Since v44".
- Nov 4 patch: `aqua_tune.plan` drops movies of an earlier plan whose run is no longer in the batch (left out by run_labels.csv). Found because the NavTIGRETreg-2608-1 plan of Nov 3 had kept the two excluded run-004 movies of its Oct plan (harmless: the report scored them under a FOV4 experiment the batch does not have and wrote an unused `FOV4_days_Opt.csv`; wasted sweep time only).
- Nov 5 patch: `cfu_overlap.check_mouse` re-roots the manifest's tif paths with `relocate()` (it silently found no recordings when the folder was read from another machine). Working copy + v44 package refreshed each time.
- Tested on synthetic folders with the five animals' real run lists; the tuning's movie choice sweeps both D14 conditions of FOV1/2 (16 movies for the 2609-1 animals).

## The driver: `D:\2photon\Chao\GCAMP_affrerents\treg_prepare.py` (three pastes for all five animals)

Run from the Administrator Anaconda Prompt: 1. `python treg_prepare.py` (plan printout, batch with `default_Opt.csv`, tuning plan with the v37 grid / classes, one movie per condition and FOV; writes `treg_sweeps_A.m` / `treg_sweeps_B.m` for two MATLAB instances). 2. `python treg_prepare.py --after` (tuning report + apply: common setting per FOV, apply_all, max false 0.10, on-fibre 0.80; batch with the tuned settings; writes `treg_aqua_all.m`). 3. `python treg_prepare.py --check`. The GUI can take over at any stage (same files); `--mice <names>` limits it.

## What happened (Nov 3–5)

- Nov 3 ~21:40 samson time: the user ran stage 1 — 1,791 sweep jobs planned (432 / 432 / 432 / 414 / 81).
- Nov 4 03:36 / 03:37 (remote, screen on): the two sweep drivers typed into the idle MATLAB R2021b windows; ~60 jobs/h per driver, 6 workers each. All done without errors: NavTIGRE-2609-1 7.0 h, NavTIGRETreg-2608-1 3.7 h (234 of its jobs were done from the October sweep), NavTIGRE-2608-1 1.4 h, NavTIGRETreg-2609-1 7.5 h, NavTIGRETreg-2609-2 7.4 h (18:33).
- The samson screen stayed locked afterwards and the PC was shut down overnight; nothing was lost (the sweeps were complete). Nov 5 16:08 samson time: the user ran stage 2 from the admin prompt — reports + apply + batches for all five, `treg_aqua_all.m` written 16:15. 16:32: `run('…\treg_aqua_all.m')` typed into the freshly opened MATLAB R2021b (remote); AQuA2 for the 59 queued movies + the CFU step for all 63 runs finished 18:14 samson time.
- **Settings applied (one common setting per FOV across days; thrARScl / minSize):** NavTIGRE-2609-1 FOV1 1.5/50, FOV2 2/20, FOV3 2/20. NavTIGRETreg-2609-1 FOV1 2.5/20 (worst day false-alarm share 15 %), FOV2 2.5/20, FOV3 2/20. NavTIGRETreg-2609-2 FOV1 2/50, FOV2 2/20 (rule not met: 64 % false alarms on D9), FOV3 2/20 (rule not met: 68 % on D9), FOV4 2.5/50. NavTIGRETreg-2608-1 FOV1 2.5/20, FOV2 2/20 (rule not met: 39 % on D5), FOV3 unchanged from October (its AQuA2 results stand; CFU redone), plus an unused FOV4 file. NavTIGRE-2608-1 FOV1 2.5/20, FOV2 2/50, FOV3 2/20. Many fields at the grid edge (minSize 20 / thr 1.5) → the grid-extension round (`grid_suggestion.json` per animal) is available later. Reports: `<mouse>\aqua_tune\tune_recommendation.txt`, `<mouse>_tune.png`, `tune_why.html`.
- **Results (overlap 0.3, min events 1, all verified from the files):** NavTIGRE-2609-1 16 runs, 172 units / 229 events (FOV1 goes 9 → 6 → 1 → 0 units across Baseline…D9, 45 units after the D14 injection; FOV2 2 → 13 → 9 → 2, 11 pre / 27 post on D14; FOV3 18 → 4 → 21 → 3); NavTIGRETreg-2609-1 16 runs, 549 / 659 (FOV3 80 → 121 → 70 → 5; FOV2 D14 post 137 units); NavTIGRETreg-2609-2 16 runs, 481 / 664 (FOV2 69 → 49 → 87 → 39, FOV3 10 → 16 → 89 → 50, FOV4 0 → 0 → 7 → 2); NavTIGRETreg-2608-1 12 runs, 610 / 916 (FOV2 38 → 102 → 102 → 98); NavTIGRE-2608-1 3 runs, 70 / 79. Raw unit counts, uncurated, under the rule-not-met caveats above.
- **Overlap check (committed to each animal's `aqua_batch\`):** 0.3 OK in every judged experiment — pooled best 0.25–0.3 with 0.3 at 96–100 % of the best (NavTIGRE-2609-1 1 OK / 2 too few pairs; NavTIGRETreg-2609-1 2 / 1; NavTIGRETreg-2609-2 3 / 1; NavTIGRETreg-2608-1 3 / 0; NavTIGRE-2608-1 too few pairs everywhere).
- **Curation status: 0 of 63 runs curated** → the five animals are at the curation point (tab 6 per animal), then tab 7 analysis and tab 8 reports.

## Cloud / scratch

Working tree /home/claude/afferent2p = v44 + the Nov 4/5 patches; `scratchpad/bedtreg/<mouse>/` = symlink beds of the staged results (the check ran there); check outputs `/mnt/user-data/outputs/treg_check_nov5/`; the driver's source `/mnt/user-data/outputs/treg_driver/treg_prepare.py`; label files `/mnt/user-data/outputs/treg_labels/`; staged under `/mnt/user-data/uploads/GCAMP_affrerents/`: all 63 `_res_cfu.mat` + `_cfu_pre.mat`, the five manifests, tune logs, tune manifests, tune recommendations and `_tune.png` figures.

## Open for later

- The cross-animal comparison for Treg (control vs Treg by day) is not built — `aqua_conditions.py` compares protocols of RearChr animals; a Treg version (genotype × day) is the natural next step once the per-animal analyses exist.
- The rule-not-met FOVs (NavTIGRETreg-2609-2 FOV2/FOV3, NavTIGRETreg-2608-1 FOV2): the analysis's false-alarm-corrected rates and the curation will matter there; a stricter setting for those FOVs or a grid-extension round are the alternatives.
- Leftover files: `NavTIGRETreg-2608-1\aqua_opts\NavTIGRETreg-2608-1_FOV4_days_Opt.csv` and the two run-004 movie folders under its `aqua_tune\` (unused; can be deleted).
