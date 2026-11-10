> **Newer source (10 Nov):** the user's own handoff from the old "Bahvior" claude.ai project is in [`docs/claude-project-handoff/`](docs/claude-project-handoff/).
> It covers all six behavior chats and the current samson `CLAUDE.md`, and wins where it differs from this draft (e.g. the Treg openness reference uses both the ambient and the dark habituation blocks).

> Draft built from the transcript of Cowork session 'DLC batch video analysis GUI' on 2026-11-09; the session's own export (branch cowork-export/dlc-behavior) will be merged in and supersedes it where they differ.
> At capture the session had not yet written its own `behavior/HANDOFF.md`. This draft draws on its summary and on the samson files it read on 9 Nov (`CLAUDE.md`, `README.md`, `docs\project-history.md`, settings and paradigm files).

# HANDOFF: behavior (squint / eyelid DLC pipeline, Treg paradigm)

## 1. What this project is

- **Squinting behaviour:**
  1. mouse eye video (Bonsai AVI);
  2. DeepLabCut **8-point eyelid** tracking;
  3. eye-openness timeline;
  4. detection of **severe squint events**, plus blinks and face **wipes**.
  - It runs across several paradigms in parallel: standard optogenetic stimulation, NTG drug, Treg ambient/dark, and LED-intensity tests.
- **MATLAB → Python port:** the analysis used to be MATLAB (`dlc_eyelid_area_auto_by_day.m`, a multi-folder batch runner and a cross-animal aggregator exporting a Prism workbook), which required every recording renamed `D1_habituation`, `D1_baseline`, …. It is being replaced by a Python package plus GUI (`python -m squint`) that works on raw Bonsai file names **without renaming**.
- **Goal (from samson's `CLAUDE.md`):** make squint and wipe characterisation more convincing: rework the normalisation, add a visual atlas to add/remove events, and improve detection.
- **What this Cowork session did for the project** (~1–5 Oct 2026):
  - built `DLC_batch_gui_squint.py`, the batch DLC GUI for the Treg videos;
  - fixed the Treg paradigm block meanings;
  - added `day_zero` experiment-day numbering;
  - added automatic per-animal paradigm recognition and multi-condition runs (per `claude/squint_pipeline_notes.md`, in both the DLC GUI and the `squint` package);
  - repaired `squint\cli.py` (a duplicate subparser and a garbled merge).
- **Rest of the package:** v2 analysis, review tool, clips, aggregate and legacy check. The samson docs credit these to the former "behavior" Claude project; this draft takes their state from `README.md` / `CLAUDE.md`.
- **Also in scope (from `CLAUDE.md`, not worked on in this session):** a **CPP** (conditioned place preference) setup at an early stage.
  - The chamber is recorded with a Logitech HD 720p USB webcam via a webcam app on a Samsung Galaxy tablet; no audio.
  - Analysis has not started. The chamber and its 3D-printed floors are covered in a separate 3D-printer project.
  - CPP folder on samson (per `CLAUDE.md`): `D:\Behavior\CPP_recorder`.

## 2. Where things live

### Data on samson

| What | Path / pattern |
|---|---|
| Raw squint recordings | dated folders in `D:\Behavior\Squinting` (e.g. `2025_10_1`, `2026_2_14`; also `Filr3`, `logitechCam`, `Sample_videos`) |
| Organised data | `D:\Behavior\Squinting\compiled\<condition>\<animal>\...` |
| Treg (raw names, not renamed) | `D:\Behavior\Squinting\compiled\Treg\<animal>\<yymmdd>\<animal>_<yymmdd>_<camera>_<block>.avi`, e.g. `NavTIGRE-2608-1_260915_F2_0`, `NavTIGRETreg-2609-1_260924_F3_Right_3`. Animals: NavTIGRE-2608-1, NavTIGRE-2609-1, NavTIGRETreg-2608-1, NavTIGRETreg-2609-1, NavTIGRETreg-2609-2. Holds `paradigm.yaml` (the reference for Treg) and `_DLC_batch_logs\DLC_batch_squint_log_20261001_142653.csv`, `..._20261005_155358.csv`. |
| NTG | `compiled\NTG\` (10mg, 1mg, NavWT-2606-3, NavWT-2606-4, plot-poses). Never renamed; sidecars may sit beside the condition folder. |
| Other conditions | e.g. `compiled\Bone-400-2hz\NavRearChr-2511-2-15ms` (the validation animal), other `Bone-*`, opto, and `compiled\High_led_test` (LED-intensity tests) |
| Renamed legacy files | `D1_baseline_1_<animal>_<date>_<cam>_<block>.avi`; DLC CSVs like `D1_stim_2_NavRearChr-2511-2DLC_..._filtered.csv` |
| Bonsai sidecars | `StartTime*` / `StopTime` (used for the true frame rate), `Quad_*`, `Optostim_*` |

### DLC models (in `D:\Behavior`)

| Model | Use |
|---|---|
| `D:\Behavior\Squint_8P_ispi\config.yaml` | **In use**: named in both settings files, `README.md` and `CLAUDE.md` |
| `D:\Behavior\Squint_8P_Large-Rodrigo-2025-09-11\config.yaml` | The default constant in the Claude-project copy of `DLC_batch_gui_squint.py` (1 Oct). Whether the samson copy still defaults to it is unknown. |
| `Squint_2cam`, `Squint_8P_contra`, `Squint_8P_Flir3`, `Squint_8P_webcam`, `Pupil` | Other model folders; check what each is before use. The pupil models for the 2-photon videos belong to `deeplabcut/`. |

### Code (`D:\Behavior\Chao\squint_pipeline\`, shared with the deeplabcut project)

| File | Size, mtime (UTC) | Role |
|---|---|---|
| `DLC_batch_gui_squint.py` | 96,584 B, 2026-10-06 01:52 | Squint DLC batch GUI/CLI, written by this session |
| `DLC_batch_gui_squint_settings.json` | 537 B, 2026-10-06 18:47 | Its saved settings (see §5) |
| `squint\` package | — | `__init__` 108, `__main__` 50, `cli` 12,997, `clips` 21,515, `config` 14,525, `detect` 13,298, `discover` 19,501, `dlc` 8,682, `features` 9,084, `gui` 66,852, `legacy_compare` 4,896, `naming` 6,338, `plots` 6,791, `review` 5,019, `run` 15,432, `session` 15,442 B. `config`, `discover`, `gui`, `plots`, `run` and `session` last changed 2026-11-03 02:05 UTC. |
| `paradigms\` | — | `standard_opto.yaml` (640 B), `treg_ambient_dark.yaml` (1,882 B), `ntg_drug.yaml` (321 B) |
| `params.yaml` | 1,202 B, 2026-09-28 | Default analysis parameters |
| `README.md` | 15,743 B, 2026-11-03 02:05 | **The full manual**: read it before changing anything |
| `requirements.txt` | 235 B | Python ≥ 3.9 deps |
| `squint_gui.bat` | 67 B, 2026-11-05 21:03 | `cd /d "%~dp0"`, `python -m squint gui`, pause on error |
| `squint_dlc_settings.json` | 374 B, 2026-11-05 21:07 | The app's DLC-step settings (see §5) |
| `CLAUDE.md`, `docs\project-history.md`, `docs\switch-prompts.txt` | 9,250 / 4,178 / 5,356 B, all 2026-11-07 01:44 | Handover docs for a local Claude Code project in this folder. They refer to `docs\claude_notes\behavior\`, but that folder **did not exist** on 9 Nov 22:26 UTC (`docs\` held only these files and `dlc-pupil-history.md`) |

**Not yet in this repo:** the code will arrive with the session's export (`cowork-export/dlc-behavior`: `behavior/code/squint_pipeline/`, `behavior/code/config_examples/`).

### Outputs

- **DLC outputs next to each video:**
  - `<video>DLC_<model>_snapshot_<snap>.h5/.csv`;
  - `..._filtered.h5` / `..._filtered.csv` (what `squint` reads);
  - `..._p60_labeled.mp4`, or clips in `_sanity_clips\`;
  - run log `<root>\_DLC_batch_logs\DLC_batch_squint_log_<stamp>.csv`;
  - archives in `_DLC_archive_<timestamp>\`.
- **Analysis in `<animal>\analysis_v2\`:**
  - `<animal>_D<n>_timeline.csv`, `_norm_binned.csv`, `_events.csv`, `_phase_summary.csv`, `_qc.png`, `_shape.png`;
  - `*_block_summary.csv`;
  - `<animal>_phase_summary_all_days.csv`, `<animal>_events_all_days.csv`;
  - `<animal>_review_status.csv`;
  - `clips\` (per-event MP4s, `<animal>_gallery.html`) and `<animal>_D<n>_atlas.png`.
- **Condition level:** `<condition>\analysis_v2\squint_summary_long.csv`, `..._events_long.csv`, `..._prism.xlsx`. With `--all`, also `<root>\analysis_v2\squint_summary_all_conditions_*`.
- **`discover`** writes `manifest.csv` into the folder.

### Project docs and notes

- **Claude project `01a0a5ec-134b-74fc-98e5-f88816f02cec`:**
  - `claude/squint_pipeline_notes.md` (2026-10-06, paradigm recognition and multi-condition runs);
  - `claude/DLC_batch_gui_squint.py` (copy from 2026-10-01; may be older than samson's);
  - memory area `areas/squint-pipeline.md`.
- **samson:** `README.md`, `CLAUDE.md`, `docs\project-history.md`. (`docs\claude_notes\behavior\`, named in `CLAUDE.md`, was absent on 9 Nov: Step 2a of `switch-prompts.txt`, which would create it, had evidently not been run.)

## 3. Environment and how the user runs things

- **samson:** Windows, local time UTC-5, no GPU. samson is shared: touch nothing outside the user's folders.
- **Squint app:** install once with `cd D:\Behavior\Chao\squint_pipeline` then `pip install -r requirements.txt` (Python ≥ 3.9). Then run `python -m squint gui` or double-click `squint_gui.bat`. The bat uses whatever `python` is on PATH; which env that is on samson is not recorded. `squint\__pycache__` on samson holds `cpython-312` files from 6–13 Oct and `cpython-310` files from 3 Nov 02:10 and 5 Nov 21:06 UTC, so the package has been run under Python 3.12 and, since Nov, Python 3.10.
- **DLC step:** needs the DeepLabCut env. The app launches `DLC_batch_gui_squint.py` as a subprocess with the env's python, as remembered in `squint_dlc_settings.json`: `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe` (**Detect** button, or browse once).
  - samson's `CLAUDE.md` names `C:\Users\cwei6\anaconda3\envs\deeplabcut` instead. Check which exists.
- **Standalone DLC GUI:** run from the DLC env:
  ```
  python DLC_batch_gui_squint.py
  python DLC_batch_gui_squint.py --no-gui --root "D:\Behavior\Squinting\compiled\Treg" [--list] [--only ...] [--exclude ...] [--dry-run] [--force-reanalyze] [--labeled full|clip|none]
  ```
- **Labeling:** frame labelling and GPU work run on another PC that maps samson's D: as `Z:` (DLC 3.0.0rc13, `python -m deeplabcut`). Don't work on the same folder from both PCs at once.
- **Reaching samson from Claude Code:** a Claude Code cloud session **cannot reach samson** unless the user shares the folder through Remote Control in the Claude desktop app on samson (Settings → Claude Code → Remote Control); then use `list_computer_folders` / `request_computer_folder`. Otherwise edit code here and give the user commands to run.
- **Local Claude Code on samson:** `docs\switch-prompts.txt` plans a local Claude Code project in `D:\Behavior\Chao\squint_pipeline` with git for code and notes only. Whether it exists is unknown.

## 4. Pipeline and tools

All of this is also in the desktop app (`python -m squint gui`). Its tabs are: **0 DLC tracking**, **1 Discover**, **2 Parameters**, **3 Run**, **4 Results**, **5 Review**, **6 Clips & gallery**, **7 Aggregate**, **8 Legacy check**.

1. **DLC tracking** (tab 0, `python -m squint dlc <folder> [--only ...] [--dry-run] [--then-analyze]`, or `DLC_batch_gui_squint.py`):
   - The constants below are from the 1 Oct Claude-project copy of `DLC_batch_gui_squint.py`; check them against the samson copy (6 Oct).
   - The scanner walks up to `MAX_SCAN_DEPTH = 5` levels. It lists videos (`.avi/.mp4/.mkv`) whose stem matches `^(?P<animal>.+?)_(?P<date>\d{6})_?(?P<camera>F\d+(?:_[A-Za-z]+)?)_(?P<block>\d+)$`.
   - It skips folders starting with `_`, `.` or `old`, and folders containing `analysis`, `plot-poses`, `labeled`, `dlc-models`, `training-datasets` or `evaluation-results`.
   - For each video: `analyze_videos` → `filterpredictions` (median, window 5) → `create_labeled_video`.
   - **Resume:** each stage is skipped if its output exists.
   - **Force re-analysis** archives old outputs to `_DLC_archive_<timestamp>`.
   - The AVI goes to DLC unchanged (no conversion or rescaling), and outputs land next to the video.
   - With *then run the squint analysis* ticked, discovery and analysis follow automatically.
2. **Discover:** `python -m squint discover <condition folder>` writes `manifest.csv`. Edit `include`, `phase` or `day` there, then run with `--manifest manifest.csv`. `python -m squint conditions <folder> --all` lists which paradigm each condition will use.
3. **Parameters:** `python -m squint params --out params.yaml`, then pass `--params params.yaml`.
4. **Run:**
   - `python -m squint run <folder> [--days 1 4] [--mode v2|legacy]`;
   - `python -m squint run D:\Behavior\Squinting\compiled --all` for every condition.
5. **Review:**
   - Tab 5: event list + trace + synced video frame. Click to move the cursor, drag to add an event, keep/discard, trim, set the type, add notes. Keys `1`–`4` set the type (squint_strong, squint_mild, blink, wipe).
   - Each animal-day has a review status (○ / ◐ / ✓).
6. **Clips:**
   - `python -m squint clips <animal or condition> [--days 1] [--types squint_strong wipe] [--pad 5] [--max-per-type 20]` makes per-event MP4s, an atlas and an HTML gallery.
   - Untick false positives in the gallery, then **Export review csv** and run `python -m squint apply-review <folder>`.
7. **Aggregate:** `python -m squint aggregate <condition folder>` writes the long tables and the Prism workbook. Events with `keep = False` are excluded.
8. **Legacy check:** `python -m squint compare <animal folder> --day 3` compares against an existing MATLAB result.

## 5. Key decisions and parameters

| Decision / parameter | Reason |
|---|---|
| **Paradigm resolved per animal folder:** (1) the nearest `paradigm.yaml` at or above it (complete, or `use: treg_ambient_dark` + overrides + `day_zero`); (2) otherwise a guess from the library by folder keyword (treg / ntg / opto / hz), then blocks per day; (3) else `standard_opto`. Guesses are flagged `(guessed)` / `(default!)` and recorded in the manifest (`paradigm`, `paradigm_source`). | The user runs several paradigms in parallel and wants them processed together and recognised automatically. |
| **Treg paradigm** (`treg_ambient_dark`): one day = 4 × 30 min blocks, video `_0`…`_3` = ambient habituation, ambient recording, dark habituation, dark recording. These map onto the slots `habituation / baseline / stim / post` with `phase_labels` `ambient_hab / ambient_rec / dark_hab / dark_rec`. | The user's correction (Oct). Reusing the standard slot names keeps all downstream code working. |
| **3-block first days** (`NavTIGRE-2608-1/260915`, `NavTIGRETreg-2608-1/260915`) → `[habituation, baseline, post]`. The lone dark block counts as the dark **recording**. | The user: the first runs had only one 30 min dark block, really a habituation, "but we can treat it as the recording for now". |
| **Experiment day** = `day_numbers["animal/yymmdd"]` > days since `day_zero[animal]` > `day_schedule` > rank of the date. Negative days = baseline. | The user's day scheme. Rank-based days are flagged `rank!` because they are not experiment days. |
| **Treg `day_zero`:** NavTIGRE-2608-1 and NavTIGRETreg-2608-1 = **260912** (no baseline; 260915 = D3, 260918 = D6, 260922 = D10, 260929 = D17). NavTIGRE-2609-1, NavTIGRETreg-2609-1 and NavTIGRETreg-2609-2 = **260926** (260924 = D-2 baseline, 260929 = D3). | The user: the first cohort starts at D3; for the others the first day is a baseline 1–2 days before D0 and the second is D3. |
| **Treg eye-openness reference** = ambient habituation + ambient recording of the same day. The dark blocks are expressed relative to ambient. | Decided in the behavior project (`CLAUDE.md` says "both the ambient and the dark habituation blocks"; `README.md` says ambient hab + rec). **This conflict is unresolved.** |
| **Standard (opto) paradigm:** 30 min habituation, 30 min baseline, 1 h stim, 30 min post (`blocks: [habituation, baseline, stim, stim, post]`). It applies to all conditions except NTG and Treg. | From the user. |
| **Squint DLC settings** (`DLC_batch_gui_squint_settings.json`): config `D:\Behavior\Squint_8P_ispi\config.yaml`, data_root `...\compiled\Treg`, shuffle null (newest trained), snapshot `best`, trainingsetindex 0, raw_csv true, median filter window 5, labeled_mode `full`, clip_seconds 120, label_filtered false, pcutoff/dotsize from config.yaml (pcutoff 0.6 → `_p60_labeled.mp4`), skeleton on, skip_completed on. | Matches how the NTG / Bone-* conditions were processed. |
| **App DLC step** (`squint_dlc_settings.json`): python `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe`, same model, snapshot best, filter on, labeled_mode **clip** (120 s), skip completed, then_analyze true. | Small, fast sanity clips when run from the app. |
| The model `config.yaml` is never edited. A temporary `_runtime_<project>_<stamp>_config.yaml` (local `project_path`, chosen `snapshotindex`) is written next to it and deleted after a fully successful run. | The original config's `project_path` points to a network path. |
| **Duplicate handling:** a raw copy of a recording that also exists renamed (`D1_phase_...`) is unticked. The copy that already has DLC outputs is kept, else the renamed one. **One DLC output per recording:** the newest snapshot wins, the rest go to `dropped_alternatives`. | MATLAB's glob double-counted 28 recordings in `Bone-400-2hz` and concatenated a duplicate `best-60` file. |
| **Discovery hygiene:** stubs < 20 s are excluded; blocks < 5 min are flagged `short`. Frame rate comes from the Bonsai sidecars (n_frames / wall time); otherwise the nominal **15.49 fps**, with `possible_frame_drop` when a block is 50–97 % of the expected length. | A recording-PC bug sometimes drops the frame rate. |
| **v2 analysis** (default): (1) per-point rejection at likelihood < 0.25 (`lik_floor`) or a jump > 6 MAD (`jump_mad_k`) from the 1 s rolling median; gaps ≤ 0.5 s interpolated; up to 2 points filled from the block's median shape; a frame is valid with ≥ 6 tracked points. (2) Openness = polygon area / (canthal width)², width = 5 s rolling median. (3) Reference = 95th percentile over `[habituation, baseline]` of the **same day**. (4) Squint candidates kept from 1 s (`squint_min_dur_s`); strong / mild thresholds 0.5 / 0.7 (or `baseline_pct`); the old 30 s cut only in `*_ge30s` columns (`squint_report_min_dur_s: 30`); blinks (< 0.5 s) are a separate channel; wipes from the legacy detector (`wipe_n`) and an occlusion detector (`wipe_v2_n`, ≥ 4 points lost). | Eye points shift and jump, so single point pairs are unreliable. The 30 s minimum was an arbitrary MATLAB choice. Per-day references because camera angle differs between days (0.47–0.56). |
| **`--mode legacy`** reproduces MATLAB: convex-hull area, p99 whole-day reference, all-8-point likelihood gates. | Validated on `NavRearChr-2511-2-15ms` D2–D5: areas equal to 1e-11 px, phase means to 1e-6, identical squint/wipe counts. |

## 6. Current status (as of 2026-11-09T22:30:48.301883Z)

This session's summary records no squint work after the Treg / paradigm-recognition work of 1–5 Oct; on 9 Nov it only read the files. The status below is from samson's `CLAUDE.md` / `README.md` as read that day.

| Item | Status |
|---|---|
| Python squint pipeline | Built: discovery without renaming, v2 analysis, squint/blink/wipe detectors, review tool, clips + gallery, aggregate to Prism, legacy mode (verified on NavRearChr-2511-2-15ms D2–D5; D1 differs only because MATLAB had a duplicate `best-60` post file). |
| Treg (`compiled\Treg`, 5 animals, 66 recordings) | DLC done via `DLC_batch_gui_squint.py` (logs from 2026-10-01 and 2026-10-05). Discovery reproduces the batch log's day and phase for all 66 recordings. Days come out as baseline (−2), D3, D6, D10, D17. Day 260915 recorded at only 10–13 fps (frame drops), so event times in those blocks may be off by seconds. **Whether `analysis_v2` has been run and reviewed for Treg is unknown.** |
| NTG, Bone-*, opto conditions | DLC-processed earlier (the squint DLC GUI copies "the way the NTG / Bone-* conditions were processed"). MATLAB results exist at least for `Bone-400-2hz` (used for the legacy check). Per-condition Python (`analysis_v2`) status was not checked in this session. |
| `High_led_test` | Tests whether the set LED intensity is enough to trigger squinting: light-only pilots (1500 / 2000 lux, 2–3 h) and single-day NTG injection tests (30 min hab + 30 min baseline, injection, 30 min wait, 4 × 30 min post, then a 5000 lux + dark test or 2 × max power ~6000–8000 lux). Analysis status unknown. |
| CPP | Setup stage; no analysis. |
| Settings changes after this session | `squint_gui.bat` and `squint_dlc_settings.json` were re-saved 2026-11-05 ~21:03–21:07 UTC, not by this session. Their content on 9 Nov was identical to what the session read on 3 Nov, and the squint app was started at 21:06 (new `cpython-310` .pyc files). |

## 7. In flight and next steps

- **Running jobs or scheduled check-ins:** none for this project.
- **The Cowork session's export** (branch `cowork-export/dlc-behavior` of `cwei6/GCAMP`, or the private artifact "dlc-behavior export") was being assembled at capture. `behavior/code/squint_pipeline/`, `behavior/code/config_examples/` and an empty `behavior/docs/` existed; **`behavior/HANDOFF.md` was not written yet**. Merge it when it arrives; it supersedes this draft.
- **Next analysis steps** (`README.md` "Next" and `CLAUDE.md`):
  1. A hand-scored validation clip → precision/recall of the squint and wipe detectors.
  2. The duration distribution of squint candidates across animals, to justify a report cut (replacing the arbitrary 30 s).
  3. `CLAUDE.md` says to check `docs\claude_notes\behavior\` for anything newer, but that folder did not exist on samson on 9 Nov. Ask the user whether Step 2a of `switch-prompts.txt` (copy the old "behavior" project's notes there) is still planned.
- **Housekeeping:** compare the samson `DLC_batch_gui_squint.py` (96,584 B, 6 Oct) with the Claude-project copy (1 Oct) and keep the newer. Confirm the default model constant.

## 8. Open questions for the user

1. Treg openness reference: ambient hab + ambient rec (`README.md`), or ambient + dark habituation (`CLAUDE.md`)?
2. The single dark block on 260915 is treated as the dark recording "for now". Keep that, or relabel it as dark habituation?
3. Is `Squint_8P_ispi` the only current squint model? What are `Squint_8P_Large-Rodrigo-2025-09-11`, `Squint_2cam` and the others still for?
4. New animals or conditions since 5 Oct that need a `day_zero` / `paradigm.yaml` entry? Without one, days fall back to `rank!`.
5. Has the Treg `analysis_v2` been run and reviewed? Is a combined all-conditions aggregate wanted?
6. Which python/env runs `squint_gui.bat` on samson, and which DLC env path is canonical (`C:\ProgramData\anaconda3\envs\DEEPLABCUT` vs `C:\Users\cwei6\anaconda3\envs\deeplabcut`)?
7. Is the local Claude Code project on samson (`squint_pipeline`, with git) in use? Which copy of the code is canonical going forward?

## 9. Working with this user

- **Messages:** short, informal, with typos ("not quiet right", "qeustion"). Infer intent; corrections come as plain statements of how the experiment really runs.
- **Not a coder:** a neuroscientist. Explain in plain language. Before anything non-trivial, give 2–3 options with trade-offs and a recommendation, then wait. This is from the user's `CLAUDE.md`.
- **Parallel paradigms:** wants the tools to recognise paradigm and day by themselves ("can the pipeline/gui automatically recognize what paradigm it is doing? Because I have those running in parallel and might need to process them together").
- **Tools:** wants GUIs that expose every function and also have a console mode; no renaming of raw recordings.
- **Results:** wants things to look at (traces with synced video frames, labeled videos, figures, galleries) and paths they can copy.
- **Code changes:** edit code in place in `squint_pipeline` (no zips or copies). Record decisions, dated, in the history file.
- **Never delete, rename or overwrite** videos, models or earlier results without asking.

## 10. Pitfalls and known issues

- **SNAPSHOT_MISMATCH:** outputs tagged `best-NN` are skipped when the Snapshot setting is an integer (and vice versa), unless the setting is changed or Force re-analysis is used.
- **Never edit a DLC `config.yaml` directly** (and never with regex). Use the runtime-copy approach.
- **Animals without `day_zero`/`day_numbers`** get rank days (`rank!`), not experiment days.
- **Frame drops:** day 260915 (Treg first run) recorded at 10–13 fps. Sidecar timestamps correct the time base, but drops are bursty.
- **Point tracking:** eye points shift, mistrack and jump across recordings. `eyelid_2` / `eyelid_8` have shown 60 px excursions, so prefer the 8-point area over single point pairs (EAR is exported but not primary).
- **Per-day references** vary with camera angle (0.47–0.56 on the validation animal), so pooling across days is off by default.
- **Bugs fixed earlier in this session (Sept–Oct; the summary doesn't always say which file):**
  - Treeview iids in a tkinter GUI;
  - `robust_smooth` crashed on short traces;
  - `squint\cli.py` had a duplicate subparser and a garbled merge (rewritten cleanly).
- **Remote-device bridge** (Cowork, probably also Remote Control):
  - writes need the current `expectedMtimeMs` (mtime drift rejected commits);
  - files > 400 MB can't be staged;
  - verify the size/md5 on samson after any write.
- **samson's docs:** `CLAUDE.md` is also the deeplabcut project's file; edit only the squint/CPP parts when working on behavior.

## 11. Sources

- **Cowork session** `session_01TXsqNvay4sMGNwGK4XyuVe`, "DLC batch video analysis GUI", read to event `b3880823-a8ba-421c-9f2f-0cd6c89bec65` (2026-11-09 22:30:48 UTC).
- **`sessions/dlc-batch-gui/`:**
  - `README.md`: session facts and sync cursor;
  - `compaction-summary.md`: the Treg/paradigm work up to 3 Nov, with the user's original messages;
  - `user-messages.md`;
  - `digest.md`.
- **Samson files the session read on 2026-11-09** (content summarised here, not copied): `D:\Behavior\Chao\squint_pipeline\CLAUDE.md`, `README.md`, `docs\project-history.md`, `docs\switch-prompts.txt`, `paradigms\treg_ambient_dark.yaml`, `DLC_batch_gui_squint_settings.json`, `squint_dlc_settings.json`, `squint_gui.bat`.
- **Claude-project files the session read on 2026-11-09:** `claude/squint_pipeline_notes.md`, `claude/DLC_batch_gui_squint.py`.
