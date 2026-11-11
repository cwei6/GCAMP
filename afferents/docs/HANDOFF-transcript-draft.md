> Draft built from the transcript of Cowork session 'Afferent activity analysis pipeline' on 2026-11-09,
> checked against the session's own cumulative summary of 2026-11-09 22:06 UTC and its project notes as quoted in the transcript;
> the session's own export (branch cowork-export/afferents) will be merged in and supersedes it where they differ.

# afferents: HANDOFF

Status as of **2026-11-09T22:32:17Z**. Times are UTC unless marked "samson". samson's local time is UTC-5.

---

## 1. What this project is

The project uses **two-photon GCaMP imaging of afferent fibres in awake mice**. Recordings are Scanbox `.sbx` files at 15.62 Hz. In the user's words, the afferents express jGCaMP8s and co-express ReaChR-mCitrine in the RearChr animals. Only the green channel carries the signal; ReaChR adds no red signal. Most runs also have a red vascular tracer in channel 2. In the Treg animals the imaged afferents are mostly near the superior sagittal sinus.

There are two paradigms, told apart by the animal name:

| Paradigm | Animals | Recording layout | Question |
|---|---|---|---|
| **RearChr optostim** (name contains `RearChr`) | NavRearChrTIGRE-2608-1, -2608-2, -2609-1 | Runs 001–008 per day. Pre: FOV1, FOV2, FOV1, FOV2 (15 min each). Then 1 h of optostimulation. Post: FOV1, FOV2, FOV1, FOV2. | Does optostimulation change afferent activity? Pre vs post within each day. Days are different protocols (sham, 2 Hz or 5 Hz, 3.3 mW or max power, 40 ms pulse width). |
| **Treg / spontaneous** (no `RearChr`; `NavTIGRETreg-*` = Treg, `NavTIGRE-*` = controls) | NavTIGRE-2609-1, NavTIGRETreg-2609-1, NavTIGRETreg-2609-2, NavTIGRETreg-2608-1, NavTIGRE-2608-1 | 15 min per FOV, 3 (or 4) FOVs per day. Days: Baseline, then D5, D7, D9 after a DTx injection. The 2609 cohort also has D14 with a second injection, recorded before and after it. | Chronic changes across days after DTx. |

The user's main concerns:

- The FOV cannot be kept exactly the same across sessions. For optostim the objective is lifted and put back between pre and post.
- Afferents do not respond at the same location every time.
- **Locomotion** is a major driver of spontaneous activity. The analysis therefore uses:
  - locomotion from `_quadrature.mat`;
  - pupil from the DeepLabCut eye csv;
  - three triggers: running onsets, pupil dilations and rest.

**Method.** Each recording is binned to 2 s/frame (31 frames at 15.62 Hz = 1.984 s). AQuA2 then detects events, and AQuA2's CFU step groups events into units. The user curates the CFUs by hand, and then the statistics and HTML reports are made. All of this runs through the `afferent2p` pipeline: a Python Tk GUI plus MATLAB. Claude wrote the pipeline and keeps it up to date.

---

## 2. Where things live

### 2.1 Data on samson

Data root: **`D:\2photon\Chao\GCAMP_affrerents\`**. Note the spelling: "affrerents".

| Path | What |
|---|---|
| `<root>\<mouse>\` | Animal folders: `NavRearChrTIGRE-2608-1`, `NavRearChrTIGRE-2608-2`, `NavRearChrTIGRE-2609-1`, `NavTIGRE-2609-1`, `NavTIGRETreg-2609-1`, `NavTIGRETreg-2609-2`, `NavTIGRETreg-2608-1`, `NavTIGRE-2608-1` |
| `<mouse>\<mouse>_<yymmdd>_<nnn>\` | One run folder. Raw files: `.sbx`, `.mat` (Scanbox info), `_quadrature.mat` (locomotion), DLC eye csv (e.g. `<run>_eye_nativeDLC_HrnetW48_Pupil_9P_WB_260402shuffle1_snapshot_best-70_filtered.csv`, 9-point pupil, point 9 = centre) |
| `<run>_bin31f_1p98s.tif` | Binned movie: 31 frames at 15.62 Hz = 1.98 s per frame. Older files use `_bin31f_1.98s` (NavTIGRE-2608-1, older NavTIGRETreg-2608-1). Match both spellings. |
| `<run>_bin31f_1p98s_AQuA2.mat` | AQuA2 result (contains the registered movie) |
| `<run>_bin31f_1p98s_AQuA2_res_cfu.mat` | CFU result. Every file now has rules "overlap ≥ 0.3, min events 1". |
| `<run>_bin31f_1p98s_AQuA2_cfu_pre.mat` | v42+. Field `pairs` = AQuA2's pairwise overlap table `s_t0`; plus per event `centre_yx`, `n_pixels`, `weight_sum`, `t0_t1`, and `max_counts`, `sz`, `overlap_thr`, `min_num_evt`, `aqua_dir`, `made` |
| `_AQuA2_Ch1.csv`, `_AQuA2_Ch1_curves.xlsx`, `_AQuA2_eventmap.tif`, `_AQuA2_events.csv`, `_AQuA2_meanimg.tif` | Other AQuA2 exports per run |
| `_AQuA2_res_cfu_curated.json` | Curation sidecar: rejected/added units, notes, and since v43 `cfu_rules`, `cfu_result_mtime` (plus `aqua_result_mtime` and the `written` stamp) |
| `_res_cfu_curated_superseded_<YYYY-MM-DD_HHMM>.json` | Stale decisions archived by v43 (a counter suffix is added if the name exists). Never deleted automatically; no `aqua_reset` pattern touches them. |
| `_flagged_frames.png` / `.csv`, `_frames.csv`, `_globalfix.png`, `<run>.json`, `_check_frames.tif`, `_check_movie.tif` | Clean-step outputs: check sheet, flagged-frame table, per-frame table (tab 7 reads excluded frames from it), QC plot (brightness trace), run metadata, ImageJ check stacks (frame stack / composite stack with traces) |
| `<run>_bin31f_1p98s_frames_decisions.json` | The user's per-stretch clean-step decisions (tab 2). Kept by Fresh start; take effect in tab 7. |
| `_raw.tif`, `_meanproj.tif`, `_ch2.tif` / `_ch2_meanproj.tif` / `_ch2.json`, `_behavior.csv`, `_motion.csv` / `.json` | Other per-run files: uncorrected copy (tab 2 "Restore from _raw.tif…"), mean projection, red channel 2, behaviour table, motion/registration check |
| `<mouse>\aqua_batch\` | `batch.csv`, `parameters_for_batch.csv`, `run_aqua_batch.m`, `run_cfu_batch.m`, `run_aqua_and_cfu.m`, `manifest.json`, `cfu_overlap_check.txt/.csv/.png`; RearChr also has `run_cfu_batch_v41.m` |
| `<mouse>\aqua_tune\` | Tuning sweep: `tune_manifest.json`, `tune_jobs.csv`, `run_aqua_tune.m`, `tune_log.txt`, `tune_recommendation.txt`, `<mouse>_tune.png` (catch vs noise), `<mouse>_tune_detail.png`, `<day_run>_why.png`, `tune_why.html`, `grid_suggestion.json`, one folder per movie (settings s00..s08, null/inject tifs) |
| `<mouse>\FOV<k>_Opt.csv` | The user's own AQuA2 settings per FOV, saved from the AQuA2 GUI (Sep convention: one file per FOV in the animal root, applied to that FOV on all days). Outranked by `aqua_opts\`. NavTIGRETreg-2608-1 still has such October-era files. |
| `<mouse>\aqua_opts\` | Applied AQuA2 settings. RearChr: `<exp>_pre_Opt.csv` / `<exp>_post_Opt.csv` + `<exp>_tuning.json`. Treg: `<mouse>_FOVk_days_Opt.csv`. |
| `<mouse>\run_labels.csv` | v44 run overrides. Columns: `run,fov,condition,use,note`. `run` = full run name (e.g. `NavTIGRE-2609-1_261010_003`; read as key `<yymmdd>_<nnn>`); empty `fov`/`condition` keep the default rule; `use` = 0/no/false/n/exclude/x leaves the run out (the row is kept); condition characters `[A-Za-z0-9_+-]`. Present for NavTIGRE-2609-1, NavTIGRETreg-2609-1 and NavTIGRETreg-2608-1 (contents in section 6.2). Tab 1 "Create / complete" adds every run with the rule's values. |
| `<mouse>\day_labels.csv` | `day,label`. Rows like `261010_pre` / `261010_post` label same-day pairs (else the day's label + suffix is used). Written Nov 3 for all five Treg animals: the 2609 animals `260925 Baseline`, `261001 D5`, `261003 D7`, `261005 D9`, `261010 D14 (second injection)`, `261010_pre D14 before the second injection`, `261010_post D14 after the second injection` (NavTIGRETreg-2609-2 only the first four); NavTIGRETreg-2608-1 `260911 Baseline`, `260917 D5`, `260919 D7`, `260921 D9`; NavTIGRE-2608-1 `260911 Baseline (no injection yet)`. |
| `<mouse>\aqua_analysis\` | Tab 7/8 outputs: `summary.txt`, `experiments_summary.csv`, `run_summary.csv`, `events_<exp>.csv`, figures (e.g. `<mouse>_rates.png`), the self-contained HTML report `<mouse>_report.html` (10–60 MB, v31 naming) and, since v32, a separate QC report. 2608-1 also has an Oct 11 `NavRearChrTIGRE-2608-1_report_preview.html`. Not re-run since the re-tuning. |
| `<mouse>\*FOV_across_days*.pdf` / FOV png | "FOV across days" figure for the PI (rows = FOV, columns = day). Memory says it was delivered Oct 27 "for all four animals" without naming them; the listing confirms `NavTIGRE-2608-1_FOV_across_days.pdf` (Oct 27 samson). Some Treg folders also hold `FOV1.png`/`FOV2.png`/`FOV3.png` dated Oct 10. |
| `<mouse>\_previous_<date>\` | Outputs parked by "Fresh start" (`aqua_reset.py`) |
| `NavRearChrTIGRE-2608-1\reg_compare\` | Registration comparison (`reg_compare_verdict.md/.png`) |
| `<root>\Claude outputs\` | Figures made by Claude: `cfu_rules_nov2\` (`260911_001_cfu_rules_compared.png`, `overlap_threshold_all_runs.png`, `overlap_threshold_optimisation.png`) and `curation_preview_nov3\` (`NavRearChrTIGRE-2608-1_curation_preview_0p3.png`, `NavRearChrTIGRE-2608-1_260911_001_before_after.png`) |
| `<parent>\cross_animal\conditions.html` | Oct 20 cross-animal protocol comparison (`aqua_conditions.py`, tab 8 "Across animals"). `<parent>` is the folder that holds the RearChr animal folders, i.e. the data root. Predates the re-tuning and the 0.3 CFUs. |
| `<root>\run_cfu_v41_all.m` + `run_cfu_v41_all_log.txt` | Nov 2 CFU re-run at 0.3 for the 96 RearChr runs |
| `<root>\treg_prepare.py`, `treg_sweeps_A.m`, `treg_sweeps_B.m`, `treg_aqua_all.m` | Treg driver and the MATLAB drivers it generated |
| `<root>\AQuA2_2026-10\` | Newer AQuA2 release (commit 8251f74). Installed but **not used** (see section 5). |
| `<root>\lagacy_matlab\` | Copy of the legacy lab MATLAB repo (spelled "lagacy"). Registration is in `ProcessScanboxData_concatFirst.m`. |
| `<root>\_DLC_batch_logs\`, `<root>\DLC_9point_selected_conditions_log_20260915_*.csv` | DLC batch logs (belong to the deeplabcut project) |
| `<root>\test_run_NavRearChrTIGRE-2608-1\` | Nov 2 CFU-rule test run |

### 2.2 Code

| Where | What |
|---|---|
| `D:\2photon\Chao\GCAMP_affrerents\afferent2p\` | **The working copy.** The user wants fixes written here directly. It is v45 (plus the v44 patches). |
| `D:\2photon\Chao\GCAMP_affrerents\afferent2p_vNN\` + `afferent2p_vNN.zip` | Versioned packages. Zips `afferent2p_v22.zip` … `afferent2p_v45.zip`; folders only from `afferent2p_v29\` on (root listing of Nov 3 plus the v44/v45 deliveries). The working copy has been the place for fixes since Oct 6 (before that the user copied from `afferent2p_v30`). |
| `D:\2photon\Chao\GCAMP_affrerents\afferent2p\pipeline_gui_settings.json` | **The user's GUI settings. Never deliver or overwrite it.** As of Nov 9 it holds: mouse `D:/2photon/Chao/GCAMP_affrerents/NavTIGRETreg-2608-1`, aqua `H:\MATLAB\AQuA2`, imagej `D:\ImageJ\ImageJ.exe`, matlab `C:\Program Files\MATLAB\R2021b\bin\matlab.exe`, um 0.792, border 5, sec 2, `b_opts` `…\afferent2p\default_Opt.csv`, `b_cfu_overlap` 0.3, `b_cfu_min` 1, `b_cfu_merge` 0.5, `b_cfu_aqua` '', auto_matlab false, t_common false, t_day_rule "matched sensitivity …". |
| Cowork container (not reachable from Claude Code) | Working tree `/home/claude/afferent2p/` = v45 as delivered (its `pipeline_gui_settings.json` is a local test copy, never delivered); packages `/mnt/user-data/outputs/afferent2p_vNN/` + zip; `treg_prepare.py` source `/mnt/user-data/outputs/treg_driver/`; label-file sources `/mnt/user-data/outputs/treg_labels/`; files staged from samson under `/mnt/user-data/uploads/GCAMP_affrerents/`; test beds in the session scratchpad (`bed43/`, `bed44/`, `bed44b/`, `bedtreg/`, `bed45/`, `tab2/` render scripts). The export should come from there. |
| This repo, `afferents/code/` | Does not exist yet; no pipeline code is in this repo. The session's export should bring the code here. |

The package has 73 files.

Top level (32 files):

- docs: `ANALYSIS_PLAN.md`, `GUIDE.md`, `README.md`, `REMOTE_PC.md`
- step scripts: `aqua_analyze.py`, `aqua_batch.py`, `aqua_clean.py`, `aqua_conditions.py`, `aqua_export_events.m`, `aqua_prep.py`, `aqua_report.py`, `aqua_reset.py`, `aqua_tune.py`
- CFU tools: `cfu_curate.py`, `cfu_overlap_check.py`
- config: `config_example.yaml`, `default_Opt.csv`, `requirements.txt`
- GUIs and runners: `gui.py`, `pipeline_gui.py`, `pipeline_status.py`, `plan.py`, `run_gui.bat`, `run_pipeline.py`, `setup_remote.bat`, `unpack_bin.py`
- registration and other tools: `debug_dump.py`, `legacy_register_run.m`, `patch_aqua2_toolbar.m`, `reg_compare.py`, `reg_demo.py`, `remap_batch.py`

Package `afferent2p\` (41 modules):

`__init__`, `adjust`, `analyze`, `aqua_io`, `behavior`, `cfu_overlap`, `cfu_shape`, `compare`, `convert`, `curate`, `discover`, `export_aqua`, `fullreg`, `general`, `globalfix`, `heatmaps`, `labels`, `loco`, `locotrials`, `manual_roi`, `matlab_run`, `motion`, `pick`, `pipeline_state`, `plan`, `prealign`, `preview`, `pupil`, `pupil_coupling`, `qc`, `qc_report`, `register`, `report`, `reset`, `runbundle`, `spatial`, `statematch`, `stats`, `timecourse`, `triggers`, `tune`.

### 2.3 Project docs and memory (in the claude.ai project, not in this repo)

claude.ai project `01a05a2b-e5f1-712c-9f39-83e753250429`, read with the Projects tool in Cowork. A Claude Code session cannot read these unless they are exported.

| Doc | Content |
|---|---|
| `claude/v34-progress.md` | **Master continuation notes**, "v34 … v45 work in progress (Oct 19 – Nov 9, 2026)". Written because a session may run out of usage mid-work; the session updates it at every milestone (standing instruction). One section per version (what changed, why, tests, delivery) plus "Notes for continuation": cloud working tree and scratch beds, the delivery recipe, open ideas, what is still open. |
| `claude/treg-processing-nov3.md` | Treg animals (rewritten Nov 5): the animals and their state before Nov 3, run → FOV exceptions, v44, the `treg_prepare.py` driver, what happened Nov 3–5, settings applied, results, overlap-check verdicts, open points |
| `claude/cfu-rules-overlap-threshold-nov2.md` | CFU rule and the reasoning for 0.3: the user's question, AQuA2's published vs newer rules, the 96-run optimisation, "Per recording?", what v41/v42 built, a v43 paragraph |
| `claude/pipeline-v31-treg-auto-report.md` | v31 Treg route, settings priority, Auto chain, MATLAB from the GUI, the HTML report's layout, scale-aware statistics (older tab numbering) |
| `claude/registration-review-and-motion-check.md` | Registration comparison; decision: AQuA2 alone |
| `claude/protocol-comparison-oct20.md` | First cross-animal protocol comparison (before re-tuning) |
| `claude/aqua-tune-v37-why-the-sweep-was-blind.md`, `claude/aqua-tune-2609-1-first-v37-sweep-read.md`, `claude/aqua-tune-settings-by-measurement.md`, `claude/aqua-tune-first-real-sweep-results.md` | Tuning design and readings |
| `claude/2609-1-three-days-quick-read-oct27.md` | Quick uncurated read of 2609-1 (261006 / 261014 / 261025) |
| `claude/locomotion-normalisation-and-timecourse.md`, `claude/todo-report-redesign-oct12.md`, `claude/pipeline-v32-report-redesign.md`, `claude/navrearchr-2608-1-run-status.md` | Earlier analysis and report design |
| Memory `/projects/01a05a2b-e5f1-712c-9f39-83e753250429/areas/gcamp-afferent-analysis.md` (last write seen: Nov 3, version `ef03320a8980`, about 34 of 48 KB used) | Dated "[stated]" lines since Sep 2026: the user's requests, decisions and what was delivered for each (the main source for the version history in section 4) |
| Memory `…/areas/treg-project.md` (last write seen: Nov 5, version `734874af0707`) | Treg cohort: animals, recording days and day labels, DTx dates, run → FOV exceptions, the "FOV across days" PDF request, the Nov 5 processing state |

Already converted into this repo:

- `afferents/docs/artifacts/afferent2p-how-to-run-aqua2-pipeline.md` (Sep 22 guide)
- `afferents/docs/artifacts/navrearchrtigre-2608-1-pre-vs-post-optostim.md` (Oct 3 report; predates re-tuning and the 0.3 CFUs)
- `afferents/docs/artifacts/mouse-trigeminal-fiber-planner.md`

The pipeline's own `GUIDE.md` and `README.md` (in the package) are the most detailed user docs. GUIDE sections added recently: Step 2 paragraph (v45); Step 3 "Does the sweep reach this field at all?", "A result tuned by hand as the reference", "The tuning at a glance", "Pick by eye" (v37–v39); troubleshooting rows for PermissionError / Administrators, matlab -batch cannot see H:, diary permission denied (v40); "CFU rules (v41) — what a unit is, and why the overlap threshold is 0.3" after Step 4, with "Does 0.3 fit this animal?" (v42); Step 6 "After a CFU re-run (v43)"; "Which run is which FOV: run_labels.csv and day_labels.csv (v44)" and a note in "Adding a day". README has a "Since vNN" note per version.

---

## 3. Environment and how the user runs things

**samson**:

- Shared Windows PC. The user's Windows account is `cwei6`, and other lab users log in too. One of them runs large MATLAB jobs, which slows the sweeps.
- Local time UTC-5.
- The connected folders in Cowork were `D:\2photon\Chao\GCAMP_affrerents` and one run folder of 2608-1 260911_001.

Software:

| Tool | Path / notes |
|---|---|
| Python | Anaconda at `C:\ProgramData\anaconda3` (the user says they use base Anaconda). On Nov 9 the running pipeline GUI process was `C:\ProgramData\anaconda3\envs\deeplabcut\python.exe`. A `deeplabcut` env exists; the transcript does not say which env the user activates. Also installed: Python 3.13 (user install) and Python 3.8. |
| MATLAB | **R2021b** `C:\Program Files\MATLAB\R2021b\bin\matlab.exe` (R2022b and R2025a also installed) |
| AQuA2 | **`H:\MATLAB\AQuA2`** (mid-August 2026 copy, published CFU rules) on the mapped H: drive. UNC path: `\\its\ifs\users\cwei6\MATLAB\AQuA2`. Not used: `GCAMP_affrerents\AQuA2_2026-10`, `H:\MATLAB\AQuA-master`, and an old manifest's `H:\MATLAB\AQuA2-main`. |
| ImageJ | `D:\ImageJ\ImageJ.exe` |

How the user runs things:

- **Pipeline GUI:**
  - Open an **Administrator Anaconda Prompt** and start `pipeline_gui.py` from `D:\2photon\Chao\GCAMP_affrerents\afferent2p\`. The exact command line is not recorded in this period.
  - Or double-click `afferent2p\run_gui.bat`, which activates base Anaconda (since v40).
  - The GUI has tabs 1–8; see section 4.
  - A window started before a code update keeps running the old code in memory. **Restart the GUI after every delivery.** v45 needs a restart.
- **MATLAB steps:** the user pastes `run('<path>.m')` into an interactive MATLAB R2021b window. Several MATLAB windows may be open at once. Commands typed into a busy window queue up and run after it.
- **Elevation:**
  - Most outputs were written by elevated processes, so they belong to BUILTIN\Administrators. A non-elevated window, MATLAB or bridge cannot overwrite them (PermissionError). v40's `write_text_free` then falls back to a new name such as `…_2.csv`.
  - `manifest.json` is written with `write_text_free(..., allow_alt=False)`, so a Python step that rebuilds a batch (e.g. `treg_prepare.py --after`, tab 4) fails instead of writing a second name. Such steps must run in the user's **Administrator** Anaconda Prompt; Claude cannot type into that prompt remotely (it is click-only for computer use), which is why the user ran `treg_prepare.py` and `--after` themselves.
  - Elevated processes cannot see the mapped H: drive; use the UNC `addpath`.
  - Oct 27 workaround to get a second elevated GUI without the user: set the ImageJ field to `afferent2p\run_gui.bat`, double-click a recording in tab 2 (this runs the "ImageJ" command, i.e. the bat, with the window's Administrator token), then restore the ImageJ field.
  - Suggested one-time fix (not known to be applied): run `icacls "D:\2photon\Chao\GCAMP_affrerents" /grant cwei6:(OI)(CI)F /T` from an Administrator prompt.

**Claude Code access:**

- A Claude Code cloud session **cannot reach samson** unless the user shares the folder through Remote Control in the Claude desktop app on samson (Settings → Claude Code → Remote Control). Then use `list_computer_folders` / `request_computer_folder`.
- Otherwise, write or fix code in this repo and give the user exact commands or pastes to run on samson.
- The Cowork session used remote-device tools: `device_list_dir`, `device_stage_files` and `device_commit_files` (at most about 50 files per call), plus computer-use to type into MATLAB. Those tools are not available in Claude Code. During long MATLAB jobs it scheduled its own check-ins (`send_later`), reported progress rather than acting beyond the request, and stopped polling after repeated "device not connected" rounds (Nov 5) to wait for the user.

**Delivery convention.** The user relies on this:

1. Write the changed files into the working copy `afferent2p\`.
2. Write the full package to `afferent2p_vNN\` plus `afferent2p_vNN.zip`, excluding `__pycache__` and `pipeline_gui_settings.json`.
3. Check the delivered files against the originals (md5 or equivalent).
4. Bump the GUI title (`afferent2p — AQuA2 route (vNN)`) and add a README "Since vNN" note (and the GUIDE paragraph for the change).
5. Update the progress doc `claude/v34-progress.md` (in Claude Code: this HANDOFF.md, per the repo `CLAUDE.md`) and tell the user to restart the GUI.

Small fixes between versions were patched into the working copy **and** the current `afferent2p_vNN\` + zip (Nov 4 `aqua_tune.py`, Nov 5 `afferent2p/cfu_overlap.py`), without a new version number.

---

## 4. Pipeline and tools

The GUI tab numbering below is the one in use since v32. The pipeline-state strip shows "1 … 5 MATLAB ✓ … next: 6 Curate / 7 Analyze".

| # | Step (GUI tab) | What it does | How to run |
|---|---|---|---|
| 1 | Unpack + plan (tab 1) | Unpack `.sbx` and bin to 2 s/frame (31 frames → `_bin31f_1p98s.tif`). The planner (`afferent2p/plan.py` `plan_mouse`) builds experiments. RearChr: day × FOV, pre/post. Treg: `<mouse>_FOV<k>_days` with conditions `d<date>` in run order; the **baseline is the first condition** (the earliest day; for a same-day pair `d<date>_pre` comes before `d<date>_post`). Run → FOV overrides come from `run_labels.csv`, labels from `day_labels.csv` (`timecourse.cond_day` / `cond_label`). Everything downstream (batch, tuning, analysis, curation, reports) reads the experiments from `aqua_batch\manifest.json`, which the batch builds from the planner. Aborted `.sbx` files (no info `.mat`) are skipped, and runs are ordered by time (e.g. 2608-1 261016: 051 replaces the aborted 005). | Tab 1 has a "Which run is which FOV / condition" block: Create/complete run_labels.csv, Open run_labels.csv, Open day_labels.csv, Show the plan. CLI: `python plan.py <mouse> --run-labels`. |
| 2 | Clean & check (tab 2) | Brightness step on the binned tiffs (`afferent2p/globalfix.py`): the level is flattened and artefact frames (e.g. whole-FOV flashes) are flagged and grouped into stretches. Frames that differ clearly from the clean neighbours (≥ 10 % of the frame and ≥ 3 × the clean floor) are the stretch's flagged frames, excluded by default; tails, mild stretches and slow state changes are kept. **Frames stay in the movie.** Per stretch the user decides "exclude flagged frames", "exclude whole stretch", "keep all", "back to default" or types exact frames (e.g. `68,76-77`); decisions go to `_frames_decisions.json` and take effect in tab 7, where only AQuA2 events in excluded frames are left out. | Buttons: *Run on the mouse folder*, *Run on selected run(s)*, *Refresh*, *Drop excluded frames…*, *Restore from _raw.tif…*; options force (re-flag with current rules), report only (preview), write check stacks, flagged-frame and flagging-rule parameters. v45 has a row under the recordings table, "check the selected recording by eye": *Show check sheet* (with *Zoom / pan (+ QC plot)*), *Frame stack (ImageJ)* (also double-click the row), *Composite stack*, *Make check stacks*, *QC plot*, *Run folder*. A second row under the stretch table holds the decision buttons. After changing a decision, *Make check stacks* refreshes the stack labels. The handlers have been unchanged since v29. |
| 3 | Tune (tab 3) | `aqua_tune.py` runs an AQuA2 settings grid on one movie per condition × FOV. Each setting runs on a real, a null (scrambled) and an injected copy of the movie. The report gives false-alarm share, recall of injected transients and on-fibre fraction (counted on the fibre mask; off-fibre reported separately). Rules: `tune.recommend_matched` (RearChr, per condition) and `tune.recommend_common` (one setting per FOV across days: highest mean recall whose false-alarm share is ≤ max-false on every day, else the best trade-off flagged "rule not met"); since v38/v39 `is_clean()` allows one null event of grace (`GRACE_NULL_EVENTS`, only when there are ≥ 4 real events). Apply writes `aqua_opts\…_Opt.csv`. Settings priority: `aqua_opts\<exp>_Opt.csv` > `FOV<k>_Opt.csv` > Opt file in the run folder > tab default > `default_Opt.csv`. | `aqua_tune.py --plan [--extend]`, then paste `run('<mouse>\aqua_tune\run_aqua_tune.m')` into MATLAB. Progress: `aqua_tune.py --progress` / `aqua_tune.sweep_progress(<mouse>)`. Then tab 3 *Report* (apply ticked) and *Apply*. Options: `--day-rule matched\|common`, `--tol-recall`, `--min-recall 0.3`, `--common-settings`, `--figure-only --dpi`. Tools: *Why this setting? (pictures)*, *Pick by eye* (`afferent2p/pick.py`), a hand-tuned result as reference (`aqua_export_events.m`). |
| 4 | Batch files (tab 4) | `aqua_batch.py` builds `aqua_batch\` (scripts + `manifest.json` with a `cfu` block). The CFU defaults are overlap 0.3 and min events 1. Tab 4 also has "Check the overlap threshold (from the CFU results)". | Tab 4 *Build batch files* (leave overwrite unticked unless re-detection is wanted). CLI: `python aqua_batch.py … --cfu-overlap 0.3 --cfu-min-events 1 [--cfu-redo] [--cfu-aqua <dir>] [--cfu-merge]`. |
| 5 | MATLAB: AQuA2 + CFU | Runs AQuA2 per queued run, then the CFU step. The step writes `_res_cfu.mat` and `_cfu_pre.mat`, and removes other AQuA2 copies from the path for the run. A run is re-queued automatically when its settings file is newer than its `_AQuA2.mat` (the shipped `default_Opt.csv` is exempt). | Paste `run('D:\2photon\Chao\GCAMP_affrerents\<mouse>\aqua_batch\run_aqua_and_cfu.m')` into MATLAB. |
| 4b | Overlap-threshold check | `afferent2p/cfu_overlap.py` checks the CFU overlap threshold per experiment and pooled. It labels event pairs by space only (same place = weighted centres within 0.5 × the smaller equivalent radius, r = √(n_pixels/π); two places = each centre outside the other's extent), clusters them by union-find single linkage at 0.10–0.50 in steps of 0.025, and computes F1 / TPR / FPR per experiment (day × FOV) and pooled. Verdicts: OK (≥ 90 % of the experiment's best F1), CHECK, or too few pairs (< 30 same-place pairs). Missing or stale `_cfu_pre.mat` tables are listed. | `python cfu_overlap_check.py <mouse> [--thr] [--pos-radius] [--min-pairs] [--ok-ratio]`, or the tab-4 button. Writes `aqua_batch\cfu_overlap_check.txt/.csv/.png`. |
| 6 | Curate (tab 6) | `cfu_curate.py` window: mean image and AQuA2's registered movie (v33), zoom/pan, reject or keep units. A shape rule (`cfu_shape.py`) suggests non-fibre units. Decisions go to `_res_cfu_curated.json`. Stale decisions are ignored (v43): the run list shows "(old decisions ignored)", the image title says "earlier decisions (other CFU result) ignored — curate afresh", and a message at start lists such runs only when stale files hold decisions. The image title also names the rule that made the result ("AQuA2 rule: overlap ≥ 0.3, min events 1"). | Tab 6 *Open curation window* / *List CFU results and decisions*. CLI: `python cfu_curate.py <mouse> [--list [--suggest] \| --reject 3,7 \| --keep-all \| --propagate --iou 0.3 \| --fps 15]`. `--list` ends with "N of M runs … curated under the current CFU rules". |
| 7 | Analyze (tab 7) | `aqua_analyze.py`: event rates per condition with locomotion matching, the three-trigger comparison (running / pupil / rest), detection columns and false-positive-corrected ratios, and a day time course. | Tab 7 |
| 8 | Reports (tab 8) | `aqua_report.py` makes the main HTML report and a QC report per animal. Comparisons are rate difference, exact Poisson intervals on the ratio, pseudo-count log ratio and modulation index, on shared axes. Report section 9 states the CFU rule and the overlap-check verdict. "Compare the protocols" (`aqua_conditions.py`) compares RearChr animals. | Tab 8 |
| – | Fresh start | `aqua_reset.py` parks outputs per stage in `<mouse>\_previous_<date>\` (or deletes them) and keeps decisions, `run_labels.csv` and `day_labels.csv`. | GUI "Fresh start…" |
| – | Status | `pipeline_status.py` / the GUI pipeline-state strip | |
| – | Treg driver | `D:\2photon\Chao\GCAMP_affrerents\treg_prepare.py` runs all five Treg animals. `python treg_prepare.py`: plan + batch + tuning plan, writes `treg_sweeps_A.m` / `treg_sweeps_B.m`. `--after`: report + apply (common per FOV) + batch, writes `treg_aqua_all.m`. `--check`: overlap check + status per animal. `--mice <names>` limits the animals. Constants: `AQUA_DIR = H:\MATLAB\AQuA2`, µm/px 0.792, border 5, 2 s, max false 0.10, on-fibre 0.80, tol_recall 0.05, 72 injection sites; v37 grid and classes (`T.DEFAULT_GRID`, `T.AMPLITUDES`, `T.SIZES_PX`), `per_condition="auto"`; `--after` calls `aqua_tune.report(..., apply_all=True, apply_common=True, day_rule="matched")`. Sweep groups (two MATLAB instances, about equal work): A = NavTIGRE-2609-1, NavTIGRETreg-2608-1, NavTIGRE-2608-1; B = NavTIGRETreg-2609-1, NavTIGRETreg-2609-2. Same files as the GUI, so the GUI can take over at any stage. | Run from an Administrator Anaconda Prompt in `D:\2photon\Chao\GCAMP_affrerents`. MATLAB: `run('D:\2photon\Chao\GCAMP_affrerents\treg_sweeps_A.m')` (and `_B`), then `run('D:\2photon\Chao\GCAMP_affrerents\treg_aqua_all.m')`. |
| – | Registration tools | `reg_compare.py`, `reg_demo.py`, `afferent2p/fullreg.py`, `legacy_register_run.m` (legacy lab registration wrapper) | Kept for reference; AQuA2's own registration is used. |

Fixed imaging parameters: 2 s/frame (1.98 s), **0.792 µm/px**, **border exclusion 5**.

### Version history (from `claude/v34-progress.md` and the memory file, as read in the transcript)

| Version (date 2026) | What changed | Why |
|---|---|---|
| ≤ v17 (Sep) | suite2p route with a per-ROI afferent filter | Dropped: about 4,000 ROIs per FOV, mostly non-afferent; the v17 filter rejected almost all. Sep 19: AQuA2 became the primary route. |
| v22–v28 (zips Sep 22–29) | Brightness / flash step (`globalfix`) that corrected or removed artefact frames | A whole-FOV flash (2608-1 260911_005 frame ~73) was detected as events. v24 removed 66 frames where one was an artefact, and v25–v28 stayed far too broad, so the user set the rule: keep frames, flag per stretch, decide by eye, keep `_raw.tif`. |
| v29 (zip Oct 5) | Tab-2 check stacks, check sheet and per-stretch decisions in their current form (handlers unchanged since v29) | The user wanted a frame stack in ImageJ, not a video, and an explicit decision before any drop. |
| Oct 2–5 (version not recorded) | CFU min events 1 (Oct 2); curation window zoom, movie panel and distinct colours for overlapping CFUs (asked Oct 5) | Every event should form a unit; overlapping CFUs were hard to tell apart and select. |
| v30 (Oct 5–10) | Matched-speed (locomotion) comparison and day time course; `aqua_tune.py`: AQuA2 grid sweep on real, null (scrambled) and injected copies, metrics on the fibre mask | Locomotion drives activity and post runs sit at different delays after stimulation; tuning by eye was arbitrary and slow, and off-fibre autofluorescence must not drive the choice. |
| v31 (Oct 11) | Treg route (one common setting per FOV, `recommend_common`), `default_Opt.csv` fallback, MATLAB from the GUI, Auto chain, HTML report, scale-aware comparison (rate difference, exact Poisson intervals, pseudo-count log ratio, modulation index) | "run_aqua_tune.m not found" on the first Treg animal; the user objected that raw numbers made the sham day look like an effect. |
| v32 (Oct 13) | Report redesign (general representation first, CFU × time heatmaps, locomotion-evoked activity with three responsive-trial rules), separate QC report, guided "Run all" with MATLAB pastes, tabs renumbered 1–8, pipeline-state strip + `pipeline_status.py` | The user's Oct 12 list; Auto always stopped on "AQuA2 folder cannot be located"; the user loses track of each animal's state. |
| v33 (Oct 14–15) | AQuA2's registered movie in curation, reports and heatmaps; panning; per-run motion check; `reg_compare.py`, `fullreg.py`, `legacy_register_run.m` | The field shifts with locomotion; registration had to be compared with the lab's legacy code. |
| v34 (Oct 19) | Day labels; RearChr tuning per day × FOV (pre + post movie); three-trigger comparison (`triggers.py`, `stats.py`); registration decision: AQuA2 alone | The objective is lifted for stimulation, so pre and post fields can differ; days differ anyway. |
| v35 (Oct 20) | `reg_demo.py`, `aqua_conditions.py` + tab 8 "Across animals", progress row, planner skips aborted `.sbx` and orders runs by time, "why this setting?" pictures, Fresh start | User requests (movie proof of drift, quick protocol comparison, a universal fresh start). |
| v36 (Oct 21) | Matched-sensitivity rule, per-condition `aqua_opts`, detection columns + fp-corrected ratio, Stop fixed, tabs/log in a PanedWindow, window 1180×900 | Different pre/post settings must not create false effects. (The taller top strip from here on hid the tab-2 rows until v45.) |
| v37 (Oct 23) | Injected classes 0.3/0.6/1.0 dF/F × 20/60/200 px, grid 1.5/2/2.5 × 20/50/100, fibre-mask fix, grid verdict + `--extend`, hand-tuned result as reference (`aqua_export_events.m`), TOO WEAK | The v30–v36 probes were dimmer than real events, so the sweeps were blind (hand-tuned 2609-1 261014_008 found events where the auto setting found none). |
| v38 / v39 (Oct 26) | Picture viewer with zoom and *Redraw at 2×*; one null event of grace; catch-vs-noise summary figure; *Pick by eye* | The tuning figures were hard to read; the user wanted to pick by eye among the best candidates. |
| v40 (Oct 27) | `write_text_free`, sweep diary fallback, `run_gui.bat` activates base Anaconda | Administrator-owned files and an H: drive invisible to elevated processes. |
| v41 / v42 (Nov 2–3) | Published AQuA2 CFU rule at overlap 0.3 / min events 1, other AQuA2 copies removed from the path; `_cfu_pre.mat` and the overlap-threshold check | See section 5. |
| v43 (Nov 3) | Stale curation decisions recognised and archived | Old 0.5 decisions would have been applied to the 0.3 units. |
| v44 (Nov 3) + Nov 4/5 patches | `run_labels.csv`; `aqua_tune.plan` drops runs no longer in the batch; `cfu_overlap` relocates paths | Treg days that break run k = FOV k. |
| v45 (Nov 9) | Tab-2 check and decision rows always visible; screen-sized window | The rows had been silently dropped by Tk since v36. |

---

## 5. Key decisions and parameters

| Decision | Value | Why |
|---|---|---|
| Analysis route | AQuA2 (events + CFUs). The suite2p route was left out (Sep 2026). | suite2p gave thousands of non-afferent ROIs; AQuA2 works well on these afferents. |
| Binning / scale | 31 frames → 1.98 s/frame, 0.792 µm/px, border 5 | Colleague's AQuA route; the user's AQuA2 convention. |
| Artefact frames | Kept in the movie. The user decides per stretch; only events in excluded frames are dropped from the statistics. | The user rejected automatic broad exclusion (v24–v28 removed far too much). A flash may be a real locomotion response. |
| Registration (Oct 19) | AQuA2 alone (variant A) | Three-way comparison with the legacy lab code (`registration-review-and-motion-check.md`). |
| Curation and evidence views (v33) | Show AQuA2's registered movie (from `_AQuA2.mat`), not the raw one | The field shifts with locomotion and was hard to track in the raw movie. |
| Settings by measurement (Oct 6, v30) | AQuA2 settings chosen by a sweep (real / null / injected), scored on the fibre mask, not by eye | The user found by-eye tuning arbitrary, slow and possibly different per day. |
| Comparisons on one scale (Oct 11, v31) | Rate difference, exact Poisson interval on the ratio, pseudo-count log ratio, modulation index, a "few events" flag, shared y axes | The user objected that raw numbers made the sham day look like an effect and high-baseline runs incomparable. Normalising by fibre area was mentioned, not built. |
| Tuning grid (v37) | thrARScl 1.5/2/2.5 × minSize 20/50/100. Injected 0.3/0.6/1.0 dF/F × 20/60/200 px at 72 sites. Variants real/null/inject (27 jobs per movie). | The v30–v36 probe transients were dimmer than real events, so the sweeps were "blind". |
| RearChr setting rule (v36) | **Matched sensitivity** per day × FOV: equal recall and equal on-fibre null false-alarm rate pre vs post. A day is NOT COMPARABLE when no clean pair exists, and TOO WEAK when recall is low. tol_recall 0.05 is a judgement call. | Different pre/post settings must not create false differences; false alarms are subtracted in the analysis. |
| Treg setting rule | **One common setting per FOV across all days** (apply_all). Max false-alarm share 0.10 on every day; on-fibre ≥ 0.80, skipped when no setting meets it on every day; otherwise best trade-off. | The Treg comparison is across days. |
| CFU rule (Nov 2, v41) | AQuA2's **published** overlap clustering: single linkage on duration-weighted footprint overlap. **Overlap threshold 0.3** (AQuA2 default 0.5), **min events 1** (AQuA2 default 3). | The user's principle: no custom merge strategy, follow the cited tool's rules and tune only its documented setting. 0.3 maximises F1 on all 96 RearChr runs (pooled 0.3; per animal 0.275/0.275/0.3). At 0.3, 74 % of same-place pairs are one unit and 8 % of two-place pairs. At 0.2, chains join 44 % of two-place pairs. |
| Newer AQuA2 release (Oct 2026, basin splitting + post-merge) | **Not adopted.** Installed at `AQuA2_2026-10` but unused; the CFU script removes it from the path. | On 260911_001 it turned 53 events into 154 units, with one event in up to 13 units. |
| Per-recording threshold | **No.** One threshold per dataset plus an automatic check (the user chose "Automatic check only"). | Unestimable per recording in 75 of 96 runs. Pre and post must share the unit definition. |
| Robustness | Run the final pre/post statistics once more at 0.5, uncurated | The user wants this. Event rates do not depend on the threshold; unit counts and unit-matched measures do. |
| Stale curation (v43) | `aqua_io.curation_stale_reason`: decisions are stale if `_AQuA2.mat` is newer than the recorded `aqua_result_mtime`, or `_res_cfu.mat` newer than `cfu_result_mtime` (both falling back to `written`, with 60 s tolerance), or the recorded `cfu_rules` differ from the CFU file's. Stale decisions are ignored (the run opens fresh) and, on the next save, archived as `_superseded_<date>.json` if they held decisions or notes. | Old 0.5 decisions would otherwise be applied silently to the renumbered 0.3 units; the CFU batch never rewrites `_AQuA2.mat`. |
| Run labels (v44) | `run_labels.csv` per animal; same-day pairs are `d<date>_pre` / `d<date>_post` | Lost or shifted FOVs and the D14 pre/post pair break the run k = FOV k rule. |
| Treg exclusions | See section 6.2 | The user's calls on FOV identity. NavTIGRETreg-2608-1 run 004 of 260917/260919: mean-image cross-correlation with the baseline FOVs ≤ 0.31 (0.55–0.75 for matched runs), so a field without a baseline day. |
| Treg sweeps in two MATLAB instances (groups A/B) | About 1,800 jobs | Roughly 20 h on one 6-worker instance. The running NavTIGRETreg-2608-1 sweep was not re-planned after the run-004 bug: only time was lost. |
| Treg processing scope (Nov 3) | All five Treg-paradigm animals, "Full route to curation (Recommended)": tuning sweep with one common setting per FOV across days → AQuA2 → CFU at 0.3 / 1 → ready for curation. The user also answered "Yes, go ahead" to Claude using the samson desktop. | The user's AskUserQuestion answers of 2026-11-03 20:39. |
| Delivery | Working copy + `afferent2p_vNN\` + zip; settings file never delivered | The user's standing instruction. |

---

## 6. Current status (as of 2026-11-09T22:32:17Z)

**Code:** afferent2p **v45** is on samson: the working copy, `afferent2p_v45\` and the zip. The working copy also carries the v44 patches (Nov 4 `aqua_tune.plan`, Nov 5 `cfu_overlap` relocate). `treg_prepare.py` carries the Nov 7 tStart patch. The user has **not yet restarted the GUI** or confirmed the v45 tab-2 fix. Unlike v43 and v44 (md5 checked by staging back), the v45 delivery was checked only by a size listing of the working copy.

### 6.1 RearChr optostim animals (96 runs)

All three animals were re-tuned with v40 rules (Oct 28–29), and AQuA2 + CFU were run. The CFU step was re-run at 0.3 / 1 on Nov 2 (19:43–20:17 samson time) and all 96 files were verified. The overlap check is OK, with a pooled best of 0.3 in each animal. **0 of 96 runs are curated under the 0.3 rules.** The only two sidecars (2608-1 261016_001 and 2609-1 261014_008) are empty and stale. Analysis (tab 7) and reports (tab 8) have **not** been re-run since the re-tuning; the Oct 20/27 comparisons predate it.

| Animal | Days (protocol) | Runs | Units at 0.3 (at 0.5) | Overlap check |
|---|---|---|---|---|
| NavRearChrTIGRE-2608-1 | 260831 test, no stim (not analysed); 260911 2 Hz 3.3 mW; 260917 sham; 260921 2 Hz max power; 261002 5 Hz max power; 261016 5 Hz max power + 40 ms (new FOV) | 40 | 1,076 (1,451); 1,489 events; 93 suggested non-fibre | 4 OK, 6 too few pairs |
| NavRearChrTIGRE-2608-2 | 260915 5 Hz 3.3 mW; 260922 sham; 261001 5 Hz max power; 261018 5 Hz max power + 40 ms (new FOV) | 32 | 1,200 (1,796) | 4 / 4 OK |
| NavRearChrTIGRE-2609-1 | 261006 sham; 261014 5 Hz max power + 40 ms; 261025 5 Hz max power (previous width, same FOVs as 261014) | 24 | 289 (450) | 3 / 3 OK |
| Total | | 96 | 2,565 (3,697); 3,860 events | |

For every RearChr day, FOV1 = runs 001/003 (pre) and 005/007 (post), and FOV2 = 002/004 (pre) and 006/008 (post).

2608-1 at 0.3, per day × FOV (units (at 0.5) / events / suggested non-fibre):

| Day | FOV1 | FOV2 |
|---|---|---|
| 260911 | 72 (81) / 82 / 14 | 32 (33) / 33 / 6 |
| 260917 | 3 (3) / 3 / 1 | 3 (3) / 3 / 0 |
| 260921 | 128 (157) / 161 / 10 | 373 (493) / 502 / 41 |
| 261002 | 278 (418) / 431 / 11 | 88 (148) / 159 / 7 |
| 261016 | 33 (39) / 39 / 0 | 66 (76) / 76 / 3 |

In 260911_001, 53 events form 47 units at 0.3. The only multi-event units are {1,2,5,50,51,52} and {6,41}.

Tuning settings applied with v40 (thrARScl/minSize; "a vs b" = pre vs post):

| Animal | Matched pairs | NOT COMPARABLE / TOO WEAK |
|---|---|---|
| 2608-1 | 260911 FOV2 (2/20 vs 2/50); 260917 FOV1 (1.5/50 vs 2.5/50); 260917 FOV2 (2.5/50 both); 260921 FOV2 (2/20 both); 261002 FOV1 (2/20 both) | 260911 FOV1 (2/20, too weak); 260921 FOV1 (2/20); 261002 FOV2 (2.5/20); 261016 FOV1 (2/50, too weak); 261016 FOV2 (2.5/20) |
| 2608-2 | 260915 FOV2 (2.5/20 vs 2/50); 261018 FOV1 (2.5/20 both); 261018 FOV2 (2/50 vs 1.5/100) | 260915 FOV1 (2/100, too weak); 260922 FOV1 (1.5/100, too weak); 260922 FOV2 (1.5/50); 261001 FOV1 (2/100); 261001 FOV2 (2/20, too weak; 40 ev/min pre, noisy FOV) |
| 2609-1 | 261006 FOV1 2/50 both; 261006 FOV2 2/20 both; 261014 FOV2 2/20 both | 261014 FOV1 2.5/50 NOT COMPARABLE; 261025 FOV1 2/50 and FOV2 2/50, both TOO WEAK trade-offs. The three empty FOV1 runs of 261014 have no events. |

The batch manifests on samson are still v40 (no `cfu` block). Rebuilding the batch in a v42+ window with overwrite unticked refreshes them without re-detecting.

### 6.2 Treg-paradigm animals (63 runs)

Tuning sweeps (Nov 4), stage 2 (`--after`, run by the user from the Administrator prompt on Nov 5 at 16:08 samson time after a PC shutdown) and AQuA2 + CFU at 0.3 / 1 (`treg_aqua_all.m` typed into MATLAB by Claude at 16:32, finished 18:14 samson time: AQuA2 on 59 queued movies, CFU on all 63 runs) are all done, and the overlap check was run. **All five animals are at the curation point: 0 of 63 runs curated.**

| Animal | Days (labels) | Runs in batch | Exclusions (run_labels.csv) | Settings per FOV (thr/minSize) | Units / events |
|---|---|---|---|---|---|
| NavTIGRE-2609-1 (control) | 260925 Baseline; DTx 260926; 261001 D5; 261003 D7; 261005 D9; 261010 D14 second injection (`d261010_pre` / `d261010_post`) | 16 | 261010: 003 out (lost FOV, bone growth), 004 out (FOV1 post, shifted plane); 005 = FOV2 post, 006 = FOV1 post. FOV3 has no D14. | FOV1 1.5/50 (mean recall 0.47, worst false 7 %); FOV2 2/20 (0.56, 7 %); FOV3 2/20 (0.57, 0 %) | 172 / 229 |
| NavTIGRETreg-2609-1 | same days | 16 | 261010: 003 out (shifted); 004 = FOV2 post, 005 = FOV1 post; FOV3 lost that day | FOV1 2.5/20 (mean recall 0.43, worst-day false 15 %); FOV2 2.5/20 (0.46, 10 %); FOV3 2/20 (0.58, 8 %) | 549 / 659 |
| NavTIGRETreg-2609-2 | 260925, 261001, 261003, 261005 (no 261010) | 16 | none (runs 001–004 = FOV1–4 every day) | FOV1 2/50 (0.43, 7 %); FOV2 2/20 (**rule not met**: false-alarm share 42 / 11 / 12 / 64 % on Baseline / D5 / D7 / D9); FOV3 2/20 (**rule not met**: 0 / 20 / 9 / 68 %); FOV4 2.5/50 (mean recall 0.26, worst 12 %) | 481 / 664 |
| NavTIGRETreg-2608-1 | 260911 Baseline; DTx 260912; 260917 D5; 260919 D7; 260921 D9 | 12 | run 004 of 260917 and 260919 out (field with no baseline) | FOV1 2.5/20 (0.38, 7 %); FOV2 2/20 (**rule not met**: false-alarm share 21 / 39 / 29 / 2 % on Baseline / D5 / D7 / D9); FOV3 kept its October values (AQuA2 not redone for its 4 runs, CFU only; 8 of 12 runs queued for AQuA2) | 610 / 916 |
| NavTIGRE-2608-1 | 260911 Baseline only (no injection yet) | 3 | none | FOV1 2.5/20 (recall 0.29), FOV2 2/50 (the existing setting, in the grid; recall 0.40), FOV3 2/20 (0.41); false alarms 0 % | 70 / 79 |

Per-run units/events at 0.3, uncurated (from the Nov 5 summary of all 63 `_res_cfu.mat`; D14 runs: NavTIGRE-2609-1 pre 001/002, post 006 = FOV1 and 005 = FOV2; NavTIGRETreg-2609-1 pre 001/002, post 005 = FOV1 and 004 = FOV2):

| NavTIGRE-2609-1 | Baseline | D5 | D7 | D9 | D14 pre | D14 post |
|---|---|---|---|---|---|---|
| FOV1 | 9/16 | 6/11 | 1/1 | 0/0 | 1/1 | 45/52 |
| FOV2 | 2/2 | 13/28 | 9/11 | 2/4 | 11/16 | 27/36 |
| FOV3 | 18/21 | 4/4 | 21/23 | 3/3 | – | – |

| NavTIGRETreg-2609-1 | Baseline | D5 | D7 | D9 | D14 pre | D14 post |
|---|---|---|---|---|---|---|
| FOV1 | 3/3 | 0/0 | 8/14 | 23/25 | 11/14 | 29/33 |
| FOV2 | 5/6 | 1/1 | 1/1 | 52/62 | 3/3 | 137/162 |
| FOV3 | 80/100 | 121/132 | 70/97 | 5/6 | – | – |

| NavTIGRETreg-2609-2 | Baseline | D5 | D7 | D9 |
|---|---|---|---|---|
| FOV1 | 2/2 | 4/4 | 50/74 | 7/9 |
| FOV2 | 69/98 | 49/91 | 87/123 | 39/52 |
| FOV3 | 10/11 | 16/20 | 89/106 | 50/64 |
| FOV4 | 0/0 | 0/0 | 7/8 | 2/2 |

| NavTIGRETreg-2608-1 | Baseline 260911 | D5 260917 | D7 260919 | D9 260921 |
|---|---|---|---|---|
| FOV1 | 9/9 | 91/142 | 54/85 | 20/33 |
| FOV2 | 38/49 | 102/146 | 102/131 | 98/136 |
| FOV3 | 3/7 | 52/89 | 23/58 | 18/31 |

NavTIGRE-2608-1 (Baseline 260911 only): FOV1 7/7, FOV2 9/11, FOV3 54/61.

Overlap check (run by the session on Nov 5 after the `relocate()` fix; outputs in each `<mouse>\aqua_batch\cfu_overlap_check.*`):

| Animal | Pooled best (F1) | At 0.3 | Per experiment |
|---|---|---|---|
| NavTIGRE-2609-1 | 0.25 (0.88) | F1 0.86 = 99 %, OK | FOV2 OK (best 0.225); FOV1, FOV3 too few pairs |
| NavTIGRETreg-2609-1 | 0.25 (0.87) | F1 0.84 = 96 %, OK | FOV2 OK (best 0.3), FOV3 OK (best 0.25, 93 %); FOV1 too few pairs |
| NavTIGRETreg-2609-2 | 0.25 (0.94) | F1 0.91 = 97 %, OK | FOV1 OK (best 0.225, 90 %, borderline), FOV2 OK, FOV3 OK (91 %); FOV4 too few pairs |
| NavTIGRETreg-2608-1 | 0.3 (0.89) | 100 %, OK | FOV1, FOV2, FOV3 all OK |
| NavTIGRE-2608-1 | **cannot judge**: too few pairs pooled (13 same-place / 54 two-place) and in all three FOVs | – | – |

So 0.3 is OK in every experiment that can judge it. The session's Nov 5 message said "the pooled optimum at 0.25–0.3 everywhere"; that is true for four animals, but NavTIGRE-2608-1 has too few pairs to judge at all.

`treg_aqua_all.m` does not need to be run again unless a day is added.

Sweep facts:

| Animal | Group | Movies | Jobs | Time |
|---|---|---|---|---|
| NavTIGRE-2609-1 | A | 16 | 432 | 7.0 h |
| NavTIGRETreg-2609-1 | B | 16 | 432 | 7.5 h |
| NavTIGRETreg-2609-2 | B | 16 | 432 | 7.4 h |
| NavTIGRETreg-2608-1 | A | 14 | 414 | 3.7 h (234 jobs reused from October; 90 jobs were for the excluded run-004 movies) |
| NavTIGRE-2608-1 | A | 3 | 81 | 1.4 h |

Stage 1 was run by the user on Nov 3 (about 21:41 samson time, 1,791 jobs). Claude typed both drivers into two idle MATLAB R2021b windows on Nov 4 (03:36 / 03:37 samson time; 6 workers each). Throughput was about 60–63 jobs/h per driver, about 6.5 min per job per worker with two instances running.

`run_labels.csv` contents (written Nov 3):

- NavTIGRE-2609-1: `261010_001` FOV1 `d261010_pre`; `_002` FOV2 `d261010_pre`; `_003` use 0 (lost FOV, bone growth); `_004` use 0 (FOV1 post, plane shifted, replaced by 006); `_005` FOV2 `d261010_post`; `_006` FOV1 `d261010_post` (recapture).
- NavTIGRETreg-2609-1: `261010_001` FOV1 pre; `_002` FOV2 pre; `_003` use 0 (first FOV1 post attempt, shifted, replaced by 005); `_004` FOV2 post; `_005` FOV1 post (FOV3 lost that day).
- NavTIGRETreg-2608-1: `260917_004` and `260919_004` use 0 (a field without a baseline day).

Leftovers: NavTIGRETreg-2608-1 has an unused `FOV4_days_Opt.csv` and `aqua_tune` folders for the excluded run-004 movies. They can be deleted.

---

## 7. In flight and next steps

- **Nothing is running.** No MATLAB jobs and no scheduled check-ins as of 2026-11-09 22:32.
- **Pending in the Cowork session:** the user's answer to the "GitHub export" question:
  - export and push to branch `cowork-export/afferents` of `cwei6/GCAMP` (artifact fallback);
  - export without pushing;
  - skip.

  When the export arrives, merge it here. It supersedes this draft, and its code should go to `afferents/code/`.
- **User's next steps:**
  1. Restart the pipeline GUI to get v45. Check that tab 2 shows the check-the-recording row and the decision row.
  2. **Curate** in tab 6, one animal at a time (set the mouse field to each GCAMP_affrerents animal): the three RearChr animals (0/96) and the five Treg animals (0/63).
  3. For each animal: tab 7 analysis, then tab 8 main + QC reports.
  4. RearChr: tab 8 "Compare the protocols" across the three animals.
  5. Run the final statistics once more at **0.5, uncurated**, as a robustness check. A report block showing 0.5 vs 0.3 side by side is an idea, not built.
  6. Optional for Treg: `python treg_prepare.py --check` (status + overlap check).
- **Before tab 7 on NavRearChrTIGRE-2609-1:** on Oct 27 the pupil (DLC) data of day 261025 was not analysed yet and the quick read left the pupil out. Check that the run folders now hold the eye csv.
- **To build later:** a cross-animal **Treg comparison** (control vs Treg by day, genotype × day). `aqua_conditions.py` only compares RearChr protocols.
- **Open ideas noted by the session, not requested and not built** (`claude/v34-progress.md`): the matched rule on per-class recall curves; the hand-tuned setting as an extra grid point when a hand `Opt` exists in the run folder; injected size classes adapted to the field; make `run_aqua_batch.m` robust to a MATLAB that cannot see H: (explicit `addpath` / UNC); route apply and analysis outputs through `write_text_free`; a CFU-threshold robustness block in the report (0.5 vs 0.3 side by side); a per-animal curation preview figure from the GUI (the Nov 3 figure came from a scratch script, `curation_preview_2608-1.py`).
- **Optional:** a grid-extension round. The grid verdict says many fields sit at the edge (thr 1.5 / minSize 20).
  - RearChr: suggested thrARScl 1.25 and minSize 10; 2608-2 also 3.0 / 200. About 10 h of MATLAB.
  - Treg: minSize 10 suggested.
  - The user has not decided.

---

## 8. Open questions for the user

1. GitHub export: push, zip only, or skip?
2. Is the v45 tab-2 fix confirmed after restarting the GUI?
3. Run a grid-extension round (RearChr about 10 h; Treg fields at the grid edge)?
4. Rule-not-met Treg FOVs (no common setting keeps the on-fibre false-alarm share ≤ 10 % on every day): NavTIGRETreg-2609-2 FOV2 (42 / 11 / 12 / 64 % on Baseline / D5 / D7 / D9) and FOV3 (0 / 20 / 9 / 68 %), and NavTIGRETreg-2608-1 FOV2 (21 / 39 / 29 / 2 %). Accept the best trade-off and use false-alarm-corrected rates, re-tune, or exclude? Note also that NavTIGRETreg-2609-1 FOV1 (worst day 15 %) and NavTIGRETreg-2609-2 FOV4 (worst 12 %) were not labelled "rule not met" in the `--after` output although their worst-day share is above 10 %; the transcript does not explain why. A likely cause, not checked: the one-null-event grace of `is_clean()` (v38/v39, `GRACE_NULL_EVENTS`, when there are ≥ 4 real events), which can pass a setting whose single null event is more than 10 % of a small count.
5. Which responsive-trial rule to use for locomotion-evoked activity? Three rules are implemented; the lab has not chosen.
6. Curate or exclude the vessel-band CFUs of 2608-1 260921_004/006?
7. Should the motion-flagged post events of 2609-1 261025 (3 of 4 post runs had motion warnings) be reviewed?
8. Delete the stray files? They must be deleted from an Administrator prompt:
   - `NavRearChrTIGRE-2608-1\aqua_batch\parameters_for_batch_2.csv` (+ `run_aqua_batch_2.m`)
   - `NavRearChrTIGRE-2608-2\aqua_batch\parameters_for_batch_2.csv`
   - other `*_2.*` / `*~2.*` files
   - NavTIGRETreg-2608-1's unused `FOV4_days_Opt.csv` and run-004 `aqua_tune` folders
   - `NavRearChrTIGRE-2608-2\aqua_tune\run_aqua_tune_b.m` and `tune_log_20261027_1807.txt` from the Oct 27 night sweep
9. Normalise rates by fibre area per FOV? Mentioned in v31; not done.
10. Two unresolved observations from the v45 work:
    - In the smoke test, "Exclude whole stretch" did not change the displayed decision of stretch 25–31. The test bed was read-only symlinks, so this may be a test artefact; it was not checked.
    - Sizes differ between the samson working copy and the cloud tree: `legacy_register_run.m` (31,286 vs 33,902 bytes) and `reg_demo.py` (21,481 vs 21,893 bytes). `afferent2p_v45\` on samson holds the cloud versions. The session's final progress doc (v45 section) says the samson files are "CRLF copies … content identical, left alone", but no check of this appears in the transcript, and the samson copies are the *smaller* ones (a CRLF copy of an LF file would be larger). Treat it as unverified; compare the files before relying on either.
11. The repo's `sessions/docs/artifacts/work-todo-2026-11-02-project-items.md` lists an open to-do: "update the recording note for 2608-1 and 2608-2 … the 1101 recordings are 5 Hz, 10 ms width". This session's transcript does not mention a 261101 day, so it is unknown whether those recordings are on samson or in the pipeline.

---

## 9. Working with this user

- A neuroscientist, not a coder: Claude writes and maintains all the pipeline code; the user runs it and makes the scientific calls.
- Writes short, informal messages with typos ("2508-1" for 2608-1, "lastest", "anacoda"). Interpret them sensibly and **state the interpretation**.
- Wants a clear yes/no on "is that okay?" questions, plus a plain explanation of what happened and what to do now (a workaround first, then the fix).
- Follows **published, citable tool rules** (AQuA2's) rather than home-made strategies. Tunes only documented settings.
- Wants parameters justified on **all recordings**, not one example, and wants robustness checks. Prefers one dataset-wide parameter plus an automatic check over adaptive per-recording settings. Checks key numbers with a colleague.
- Makes the scientific calls on FOV identity and exclusions. Ask when a file-name rule does not hold.
- Wants to **see** things: image summaries, before/after figures, check sheets, zoomable HTML reports rather than embedded doc figures, and a "FOV across days" PDF per animal for the PI.
- Multitasks and loses track of pipeline state. Values status summaries with counts, and the pipeline-state strip. Keep the progress doc current at every milestone, because a session may run out of usage mid-work.
- Wants to resume without revisiting resolved problems. No stale decisions or leftovers should be carried over silently.
- Runs elevated Python steps themselves and pastes MATLAB commands. Accepted Claude typing into MATLAB remotely when the screen was free, but ask first: samson is shared, and the user once said "Wait until I say the desktop is free". Give **exact commands and paths**.
- When re-running, wants a clean overwrite of earlier results ("a universal start", hence Fresh start), not a mix of old and new outputs.
- Figures: the session used the palette `#2a78d6 #eb6834 #1baf7a #eda100 #e87ba4 #008300 #6250d6 #e34948`, no dual axes, short two-line panel titles, and a caption that states any interpretation (e.g. "2508-1" read as 2608-1).
- Wants code fixes written into the samson working copy, with versioned packages, and the settings file never touched.
- Output convention since Sep: per-run outputs are uint16 tiffs that open in ImageJ, saved next to the `.sbx` in the run folder, with no extra aggregation folders.
- Curates by hand before analysis and reports are re-run.

---

## 10. Pitfalls and known issues

- **Spelling:** the data root is `GCAMP_affrerents` and the legacy repo copy is `lagacy_matlab`. Binned-file suffixes come in two forms: `_bin31f_1p98s` and `_bin31f_1.98s`.
- **Restart the GUI after a delivery.** A running window keeps old code in memory. The curation window and `cfu_curate.py --list` are subprocesses, so they pick up new files at once.
- **Administrator-owned files:** non-elevated processes cannot overwrite them; `write_text_free` falls back to `…_2.csv`. Elevated processes cannot see H:, so `matlab -batch` from an elevated GUI gives "AQuA2 not found". Use `addpath('\\its\ifs\users\cwei6\MATLAB\AQuA2')` or paste into the user's MATLAB.
- **Two AQuA2 copies on the MATLAB path** merge their `+cfu` packages. The v41+ scripts remove other copies (`which('aqua_cmd_batch','-all')` + `rmpath`).
- **CFU re-runs do not touch `_AQuA2.mat`** (`whetherUpdateRes=false`). Anything keyed on `_AQuA2.mat` dates misses CFU changes; v43 fixed this for curation.
- **Manifest paths are absolute Windows paths.** Code reading a moved or staged copy must use `aqua_io.relocate()` (the Nov 5 bug in `cfu_overlap`). `relocate()` also needs the binned `.tif` to exist.
- **Tuning plans reuse earlier movies.** Fixed Nov 4 to drop runs that are no longer in the batch. Check `aqua_tune` when `run_labels.csv` changes.
- **Batch scripts clear the MATLAB workspace.** Keep driver state in root appdata (`setappdata(0,…)`), not in variables (the Nov 7 tStart error).
- **Tk packer:** widgets packed last are dropped silently when space runs out (the tab-2 rows were hidden from v36 to v44). Pack must-show rows first with `side="bottom"`.
- **tkinter in a Linux container:** run `apt-get update` first (without it `python3-tk` had "no installation candidate"), then `apt-get install -y python3-tk tk`, and use `/usr/bin/python3.12` (Tk 8.6); `/usr/local/bin/python3` has no tkinter. Render with `xvfb-run -a -s '-screen 0 1920x1080x24'`. Installing only `numpy pillow` was enough to import the GUI module; a big pip install timed out.
- **Remote work on samson (Cowork, or Claude Code if a computer is attached):**
  - The screen lock gives an empty capture ("desktopCapturer returned no screen sources"); access grants expire after 30 min idle.
  - The bridge can be down for hours, and samson may be shut down overnight.
  - Apps resolve by executable path, not name ("Python"/"Matlab" by name were "not installed"): `c:\programdata\anaconda3\python.exe` appears as "Jupyter Notebook" (full tier), `c:\program files\matlab\r2021b\bin\win64\matlab.exe` as "Matlab". File Explorer and the Anaconda Prompt (`cmd.exe`) are click-only, so no typing into prompts. Type MATLAB commands in chunks of at most 13 characters, then Return; commands typed into a busy MATLAB window queue. Taskbar thumbnails do not show in screenshots.
  - `device_bash` runs in a Linux VM, not on Windows: it cannot start Windows Python or MATLAB.
  - Recursive listing is capped at 2000 entries; staging and committing take at most 50 files per call, and `device_commit_files` only from `/mnt/user-data/outputs`.
  - Commit from a fresh staged folder (a reused staged path once uploaded stale content) and verify by staging back + md5.
- **Windows paths on Linux:** `pathlib.Path(r"D:\…").parent.name` is empty on Linux (the first `treg_prepare.write_driver` bug, fixed by passing `(mouse, script)` pairs). Use `PureWindowsPath` or pass names explicitly. Test beds need placeholder tifs + symlinks because `relocate()` checks the tif.
- **`pipeline_state` API:** `mouse_state(mouse_dir, sec)` and `format_state(...)` (the first `treg_prepare.py` imported non-existent `pipeline_state` / `state_text`).
- **The tune-log time estimate is wrong:** it says 1.4 min per job; the real figure is about 6.5 min per job per worker with two instances.
- **Settings precedence:** `aqua_opts\<exp>_Opt.csv` outranks everything else, and NavTIGRETreg-2608-1 still has October-era files. A run whose settings file is newer than its `_AQuA2.mat` is re-queued.
- **Statistics caveats:** RearChr days marked NOT COMPARABLE or TOO WEAK, and Treg FOVs where the rule was not met, need false-alarm-corrected rates. The Oct 20/27 cross-animal results predate the re-tuning and the 0.3 CFUs.
- Check the settings file's `mouse` field before running a step. On Nov 3 it pointed at `D:/2photon/Chao/GCAMP_visual_exp/NavRearChr-2609-1`; on Nov 9 it pointed at `D:/2photon/Chao/GCAMP_affrerents/NavTIGRETreg-2608-1`.

---

## 11. Sources

- Cowork session **"Afferent activity analysis pipeline"**, `session_01QvBfQetv2TTEfBpfqGr53o`. It was still running at fetch time; transcript up to event `de17a2f4-5101-46df-a716-72d3de4acc6f` (2026-11-09T22:32:17Z).
- `sessions/afferent-pipeline/README.md`: session facts and sync cursor.
- `sessions/afferent-pipeline/compaction-summary-2026-11-09.md`: **the newest summary**, the session's own cumulative summary written at the 2026-11-09 22:06 UTC compaction (covers Nov 3 to the user's 22:03 tab-2 question; v45 came after it and is covered by the digest).
- `sessions/afferent-pipeline/compaction-summary.md`: the session's own summary at the 2026-11-03 18:58 UTC compaction (work up to v42 and the CFU re-run at 0.3).
- `sessions/afferent-pipeline/user-messages.md`: the user's messages, Nov 3–9, verbatim.
- `sessions/afferent-pipeline/digest.md`: chronological digest, Nov 3–9.
- The project docs `claude/v34-progress.md` (Nov 3 version as read, plus the v44/v45 sections and continuation notes the session wrote on Nov 9), `claude/treg-processing-nov3.md` (Nov 5), `claude/pipeline-v31-treg-auto-report.md`, and the memory files, as read or written inside the transcript.
- The session's own export (branch `cowork-export/afferents`), when it arrives, is authoritative.
