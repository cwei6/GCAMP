# Two-camera pupil DLC batch, 3–4 Nov 2026 (samson, CPU only)

Started 2026-11-03 05:08 (samson local time, EST) with `pupil_run_5hz40ms.bat`.
Claude launched it from File Explorer through computer use, so the run was **not elevated**.
The .bat runs three steps one after another:

1. `DLC_batch_gui.py --no-gui --root "D:\2photon\Chao\GCAMP_visual_exp\5hz 40ms"`
2. `... --root "D:\2photon\Chao\GCAMP_visual_exp\visualstim_test\NavRearChr-2609-1" --only 261020`
3. `... --root (same) --only 261014_001 --force-reanalyze`

Logs:
- Whole run: `5hz 40ms\_DLC_batch_logs\console_run.txt` (~8 MB, mostly tqdm progress).
- Step 1: `5hz 40ms\_DLC_batch_logs\DLC_batch_log_20261103_050813.csv`.
- Steps 2 and 3: `visualstim_test\NavRearChr-2609-1\_DLC_batch_logs\DLC_batch_log_20261103_235538.csv` and `..._20261104_021737.csv`.

## Camera check, done before the run

`DLC_batch_gui.py --no-gui --root <folder> --frames <dir>` saved 2 frames per camera video into `5hz 40ms\_camera_frames\`.

| Recording (NavRearChr-2609-1) | `_eye.mj2` | `_face.mj2` | Models used |
|---|---|---|---|
| visualstim_test 261014_001, 261020_001, 261020_002 (single-camera test runs) | angled side view | – | eye → angled (via `pupil_models.yaml`, `eye: angled`) |
| 5hz 40ms 261027_001 … 261031_002 (7 recordings) | frontal close-up, pupil ~25–50 px | angled side view, pupil ~9 × 4–5 px | eye → frontal, face → angled (default rule) |

The user expected the first few two-camera recordings to have the cameras the other way round. They don't: all seven 5hz 40ms recordings follow the default rule.
The recordings with the angled camera saved as `_eye` are the single-camera test runs in visualstim_test. They are covered by the `pupil_models.yaml` in `visualstim_test\NavRearChr-2609-1`.

The short NavWT-2609-1 10/20 two-camera tests (visualstim_test) have `_eye` = angled and `_face` = frontal close-up. They were not analysed.

## Results

Values are the median likelihood, with the % of frames that have any point < 0.6 in brackets. Every video below has status SUCCESS unless noted.

| Recording | eye → frontal (Pupil_9P_WB_260402) | face → angled (Pupil_9P_angled_261016) |
|---|---|---|
| 261027_001 | 0.91 (8 %) | 0.85 (4 %) |
| 261027_002 | 0.89 (81 %) | 0.93 (0 %) |
| 261029_001 | 0.88 (8 %) | 0.93 (0 %) |
| 261029_002 | 0.81 (56 %) | 0.90 (2 %) |
| 261029_003 | 0.86 (26 %) | 0.89 (2 %) |
| 261031_001 | 0.89 (19 %) | 0.89 (0 %) |
| 261031_002 | 0.85 (51 %) | 0.89 (2 %) |

Test runs (eye → angled):

| Recording | Result |
|---|---|
| 261020_001 | 0.98 (0 %) |
| 261020_002 | 0.89 (36 %). Clean for the first ~11 min; in the last ~6 min all 9 points drop together, probably the eye closing or the pupil reflection disappearing. Check the labeled video before using that part. |
| 261014_001 | **ARCHIVE_FAILED, not re-analysed** (see below) |

Checks that were done:
- **Frames:** raw predicted points were plotted on native frames for 261027_001 (eye and face), 261027_002 eye and 261029_002 eye.
- **Point placement:** the edge points sit on the pupil boundary in both views.
- **Why some eye videos have a high "% low":** it comes from the lower/right edge points (pupil_5 … pupil_8). In those sessions that part of the pupil is under the lower lid (261027_002) or low-contrast against the iris (261029_002). The confident points are still placed correctly.
- **Recommendation for pupil size:** fit an ellipse to the points with likelihood ≥ 0.6 instead of using all 8.
- **Angled view:** the pupil is small (about 9 × 4–5 px), so diameters from that camera are coarse.

## 261014_001: archive step failed

`--force-reanalyze` first moves the old DLC outputs into `_DLC_archive_<stamp>`. Windows refused to delete the old raw H5 after copying it (`PermissionError: [WinError 5] Access is denied`).

- **Most likely cause:** those outputs were written on 10/16 from the user's *Administrator* Anaconda Prompt, while this run was not elevated.
- **What was left behind:** the originals (frontal model, median likelihood 0.13, useless for this angled view) are all still in place. `261014_001\_DLC_archive_20261104_021737\` holds only a copy of the raw H5.
- **Runtime config:** because that model had a failure, its runtime config was kept: `Pupil_9P_angled_261016\_runtime_Pupil_9P_angled_261016_20261104_021737_config.yaml`. It can be deleted.
- **Fix in the code:** `archive_dlc_outputs` in `DLC_batch_gui.py` now uses `os.replace` (plain rename). It reports a clear error and leaves no duplicate or empty archive folder.
- **To redo:** run from the Administrator prompt:
  `python DLC_batch_gui.py --no-gui --root "D:\2photon\Chao\GCAMP_visual_exp\visualstim_test\NavRearChr-2609-1" --only 261014_001 --force-reanalyze`
- **Alternative:** angled-model outputs for this video already exist in `Pupil_9P_angled_261016\_test_output\` (step3 of model building).

**Correction:** `CLAUDE.md` and `docs\dlc-pupil-history.md` on samson (written 5–7 Nov from the run's plan) say 261014_001 was re-analysed and its old outputs archived. That is wrong; see above.

## Speed

HRNet-W48 on samson's CPU:
- About 26–34 frames/s when the machine is otherwise idle (Oct 16/18 runs).
- About 7–9 frames/s during this run, with other jobs running (MATLAB, Jupyter and the squint GUI were open), dropping to ~1 frame/s for a while in the evening.
- At about 8 frames/s, a 30 000-frame (16.7 min) video takes about 65 min, plus about 2 min for conversion, filtering and the labeled video.
