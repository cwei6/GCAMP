> Draft built from the transcript of Cowork session 'DLC batch video analysis GUI' on 2026-11-09; the session's own export (branch cowork-export/dlc-behavior) will be merged in and supersedes it where they differ.
> Sources and session records: `sessions/dlc-batch-gui/`; the session's own pupil handoff, rebuilt from its transcript, is in `sessions/dlc-batch-gui/recovered/HANDOFF.md`.

# HANDOFF: deeplabcut (pupil DeepLabCut for the 2-photon eye / face videos)

## 1. What this project is

- **Goal:** track the pupil in the camera videos recorded during the user's 2-photon GCaMP experiments (trigeminal afferents in `GCAMP_affrerents`, visual stimulation in `GCAMP_visual_exp`).
- **Model:** a 9-point DeepLabCut model, with 8 points on the pupil edge plus `pupil_9` at the centre, so the radius survives partial occlusion by the lid.
- **Use downstream:** the filtered DLC CSVs feed the afferent and visual-exp imaging pipelines, where the pupil serves as an arousal control and as a trigger. Those pipelines belong to `afferents/` and `gcamp/`.
- **Main tool:** `DLC_batch_gui.py`, a tkinter GUI with a console mode. For every raw camera video (`*.mj2`) it:
  - converts MJ2 → MP4 at native size;
  - runs DLC `analyze_videos` → `filterpredictions` → filtered CSV → labeled sanity video;
  - reports a likelihood QC.
- **Two models since Nov 2026,** chosen by camera:
  - the colleague's **frontal** model `Pupil_9P_WB_260402` ("Rodrigo's model");
  - an **angled** model `Pupil_9P_angled_261016`, built in Oct 2026 for the new angled camera of the visual-stim setup.
- **History:**
  - Older script `DLC_9point_selected_conditions.py` (frontal model only). Its logs in the `GCAMP_visual_exp` root date from 10–13 Aug 2026; the Claude-project copy is dated 2026-09-15.
  - First `DLC_batch_gui.py` on 25 Sep 2026.
  - Two-camera rebuild on 2–3 Nov 2026.
  - First real two-model batch on 3–4 Nov 2026.

## 2. Where things live

### Data on samson (Windows lab server, drive D:)

| What | Path / pattern |
|---|---|
| Afferent experiments (one frontal eye camera) | `D:\2photon\Chao\GCAMP_affrerents\<animal>\<animal>_<YYMMDD>_<NNN>\<animal>_<YYMMDD>_<NNN>_eye.mj2`. Example: `D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\NavRearChrTIGRE-2608-1_260911_001\NavRearChrTIGRE-2608-1_260911_001_eye.mj2`. The folder name really is spelled `affrerents`. The recording folder also holds `.sbx`, `.mat` and `_quadrature.mat` (imaging data, not this project's). Older data used condition → animal → recording. |
| Visual-stim experiments | `D:\2photon\Chao\GCAMP_visual_exp\<condition>\<animal>\<recording>\`. Condition folders seen on 3 Nov: `15ms test`, `2hz test`, `2hz_wheellock`, `50ms test`, `5hz`, `5hz 40ms`, `old-200umfiber`, `old-200umfiber_2`, plus `visualstim_test\`. The root also holds analysis folders (`_gcamp_stats`, `gcamp_pipeline`, `All_test_D1-5`, …; gcamp project), the angled model folder, and the old batch logs `DLC_9point_selected_conditions_log_20260810_*` … `_20260813_*`, `DLC_9point_upscaled_mj2_log_20260731_*` / `_20260807_*` and `DLC_new_model_rerun_log_20260729_*`. |
| Two-camera recordings (from 10/27) | `GCAMP_visual_exp\5hz 40ms\NavRearChr-2609-1\NavRearChr-2609-1_<YYMMDD>_<NNN>\` with `<rec>_eye.mj2` (**frontal close-up**, pupil ~25–52 px) and `<rec>_face.mj2` (**angled side view**, pupil ~9 × 4–5 px). Recordings: 261027_001, 261027_002, 261029_001, 261029_002, 261029_003, 261031_001, 261031_002. Each MJ2 is ~435–570 MB. |
| Single-camera test runs | `GCAMP_visual_exp\visualstim_test\NavRearChr-2609-1\`: 261014_001 (32,500 frames, the angled model's training video), 261020_001, 261020_002. The **angled** camera was saved as `_eye`. |
| Short two-camera tests, cameras reversed | `GCAMP_visual_exp\visualstim_test\NavWT-2609-1\` 261020_001, 261020_002: `_eye` = angled, `_face` = frontal close-up. |
| Camera-check frames | `GCAMP_visual_exp\5hz 40ms\_camera_frames\<rec>_<eye\|face>_<sec>s.png` (34 PNGs: 2 per camera video, at ~250 s and ~600 s) plus `frames_log.txt` |

**Video format:**
- Current: 240×180 grey JPEG2000 MJ2 at 30 fps, ~30,000 frames (~16.7 min).
- Old (afferents): 160×112, tagged 30 fps but truly 15.62 fps. These were upscaled to 640×448 by the old script.

**Recording schedule (two-camera visual stim):** a recording every other day (10/27, 10/29, 10/31), with a test run one week before.

### Models

| Alias | config.yaml | Facts |
|---|---|---|
| frontal | `D:\2photon\Rodrigo\DLC_Models-RN\Pupil_9P_WB_260402\config.yaml` | The colleague's folder, **read-only by the user's rule**. HRNet-W48, shuffle 1, snapshot `...\dlc-models-pytorch\iteration-0\Pupil_9P_WB_260402-trainset95shuffle1\train\snapshot-best-070.pt`. Scorer `DLC_HrnetW48_Pupil_9P_WB_260402shuffle1_snapshot_best-70`. |
| angled | `D:\2photon\Chao\GCAMP_visual_exp\Pupil_9P_angled_261016\config.yaml` | Created 2026-10-16. HRNet-W48, shuffle 1, best snapshot `snapshot-best-090.pt`, initialised from the frontal model's best-070. 79 labeled frames from one animal/session (261014_001). Train RMSE 1.04 px, test RMSE 1.06 px. Whole-video median likelihood 0.96–1.00 on the training video. Scorer `..._Pupil_9P_angled261016shuffle1_snapshot_best-90`. The folder also has `scripts\`, `_qc\`, `_test_output\`, `labeled-data\`, `training-datasets\`, `videos\`, `README.md`. |

**Point convention (both models):** `pupil_1` = left-most edge point, then clockwise on screen (2 upper-left, 3 top, 4 upper-right, 5 right, 6 lower-right, 7 bottom, 8 lower-left); `pupil_9` = centre.

### Outputs (next to the MJ2, in the recording folder)

- **Model input:** `<rec>_<eye|face>_native.mp4`.
- **Frontal outputs:** `<rec>_eye_nativeDLC_HrnetW48_Pupil_9P_WB_260402shuffle1_snapshot_best-70{.h5, _filtered.h5, _filtered.csv, _full.pickle, _meta.pickle, _filtered_p10_labeled.mp4}`.
- **Angled outputs:** `<rec>_<face|eye>_nativeDLC_HrnetW48_Pupil_9P_angled261016shuffle1_snapshot_best-90{...same...}`.
- **Run logs:** `<scanned root>\_DLC_batch_logs\DLC_batch_log_<YYYYmmdd_HHMMSS>.csv` (stamp in samson local time). When started with the launcher there is also `console_run.txt`.
- **Archives:** forced re-analysis moves old outputs to `<recording>\_DLC_archive_<stamp>\`. Folders starting with `_` or `.` are never scanned.
- **Old script's log:** `<run_root>\DLC_9point_selected_conditions_log_<stamp>.csv`.

### Code (on samson, `D:\Behavior\Chao\squint_pipeline\`, shared with the behavior project)

| File | Version / facts |
|---|---|
| `DLC_batch_gui.py` | **Current: 2026-11-04 08:48:05 UTC**, 115,749 B, md5 `60cfe6672c6427edd68f38546c0607a9` (adds the archive fix). Identical to the Cowork working copy on 2026-11-09. Earlier versions: 2026-11-03 10:05 UTC (114,929 B, md5 `af6fc789…`, adds `--frames`); 2026-11-03 02:56 UTC (113,334 B, two-camera); 2026-09-25 21:58 UTC (82,483 B, frontal only). |
| `DLC_batch_gui_settings.json` | The user's saved GUI settings (598 B): `config_path` = frontal, `data_root` = `D:\2photon\Chao\GCAMP_affrerents`, dotsize 3, pcutoff 0.1, median window 5, `converted_suffix` `_native`, `h264_crf` 18. No angled key, so the default angled path is used. |
| `pupil_check_frames.bat` | 986 B, CRLF. Runs `--frames` for "5hz 40ms" and `visualstim_test\NavRearChr-2609-1`. |
| `pupil_run_5hz40ms.bat` | 1,585 B. The 3–4 Nov batch launcher (see §4). |
| `D:\2photon\Chao\GCAMP_visual_exp\visualstim_test\NavRearChr-2609-1\pupil_models.yaml` | 174 B, `eye: angled` |
| `D:\2photon\Chao\GCAMP_visual_exp\Pupil_9P_angled_261016\scripts\` | `angled_pupil_autolabel.py` (23,519 B, md5 `b8e214fdbc03d47b750c1c6ec25dae72`), `step1_review_labels.py`, `step2_train.py`, `step3_test_model.py`, `qc_dlc_output.py` |
| `DLC_9point_selected_conditions.py` | Older frontal-only batch with a condition picker and pause/resume. Known only as the Claude project file `claude/DLC_9point_selected_conditions.py` (created 2026-09-15). **Where it is on samson is unknown**; samson's `CLAUDE.md` says it is not in `squint_pipeline`. |

**Not yet in this repo:**
- The code will arrive with the session's export (`cowork-export/dlc-behavior`: `deeplabcut/code/`).
- Rebuilt from the transcript: `sessions/dlc-batch-gui/recovered/pupil_models.example.yaml`.

### Project docs and notes

- **Claude project `01a0a5ec-134b-74fc-98e5-f88816f02cec` (claude.ai):**
  - `claude/DLC_batch_gui.py` (doc `4aa6c869-1cdd-4e5c-bfc8-444c98b6363b`, updated 2026-11-04);
  - `claude/angled_pupil_model_notes.md` (doc `5219e46a-1556-43d9-a7c3-dda973ef433b`, updated 2026-11-04);
  - `claude/DLC_9point_selected_conditions.py`;
  - `claude/angled_pupil_autolabel.py`;
  - memory area `areas/dlc-pupil-pipeline.md`.
- **On samson (written ~5–7 Nov for a local Claude Code move, not by this session; all three saved 2026-11-07 01:44 UTC):**
  - `D:\Behavior\Chao\squint_pipeline\CLAUDE.md`;
  - `docs\dlc-pupil-history.md`;
  - `docs\switch-prompts.txt`.
  - `CLAUDE.md` and the history point to `docs\claude_notes\deeplabcut\` (including an `angled_pupil_model_notes.md` there), but **that folder did not exist on 9 Nov 22:26 UTC**: `docs\` held only the three files above. Step 2b of `switch-prompts.txt`, which would create it, had evidently not been run. The angled-model notes exist only in the Claude project and in `recovered/` here.
- **In this repo:** `sessions/dlc-batch-gui/recovered/angled_pupil_model_notes.md` (full notes, 4 Nov) and `recovered/batch_2026-11-03_results.md`.

## 3. Environment and how the user runs things

- **samson:**
  - Windows (win32 x64), local time UTC-5, **no GPU**, so DLC inference runs on the CPU. napari does not work there.
  - Also installed: MATLAB R2021b, ImageJ, Jupyter, Python 3.8 and 3.13.
  - samson is shared: touch nothing outside the user's folders.
- **DLC environment:** DeepLabCut 3.0.0rc13 (PyTorch); the 4 Nov console log shows "Loading DLC 3.0.0rc13" and "GPU: none detected".
  - The launchers use `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe`, with FFmpeg/ffprobe from `...\envs\DEEPLABCUT\Library\bin`. This worked for the whole 3–4 Nov run.
  - samson's `CLAUDE.md` names `C:\Users\cwei6\anaconda3\envs\deeplabcut`. **Which is canonical is unknown**; check before relying on either.
- **Labeling PC:** frame labelling and GPU work run on another PC. It sees samson's D: as `Z:` and starts the DLC GUI with `python -m deeplabcut` (DLC 3.0.0rc13). Don't work on the same folder from both PCs at once.
- **How the user runs the pupil batch:** from an *Administrator: Anaconda Prompt* with the deeplabcut env active:
  ```
  cd /d D:\Behavior\Chao\squint_pipeline
  python DLC_batch_gui.py
  ```
  Settings are remembered in `DLC_batch_gui_settings.json`.
- **Speed (HRNet-W48, CPU):**
  - 26–34 fps on an idle samson (the 16/18 Oct runs);
  - 7–9 fps during the 3–4 Nov run, with MATLAB, Jupyter and the squint GUI open (~65 min per 30,000-frame video);
  - once ~1 fps.
- **Reaching samson from Claude Code:** a Claude Code cloud session **cannot reach samson** unless the user shares the folder through Remote Control in the Claude desktop app on samson (Settings → Claude Code → Remote Control). If they have, use `list_computer_folders` / `request_computer_folder` to attach it. Otherwise, write or fix code here and give the user copy-paste commands to run on samson.
- **Local Claude Code on samson:** `docs\switch-prompts.txt` on samson plans a local Claude Code project in `D:\Behavior\Chao\squint_pipeline` (desktop app → Code tab → Local). Whether it was set up, and whether that folder is under git, is unknown.

## 4. Pipeline and tools

### Steps for a new batch

1. **Check the camera views first** for any new setup or recording type.
   - Console: `python DLC_batch_gui.py --no-gui --root "<folder>" --frames "<outdir>"` saves 2 PNGs per camera video, named `<rec>_<camera>_<sec>s.png` and labelled with the camera → model mapping, then exits.
   - GUI: highlight recordings and click **Show frames**.
2. **If a folder deviates from the default rule,** put a `pupil_models.yaml` in the condition, animal or recording folder (the nearest one wins):
   ```yaml
   eye: frontal            # *_eye.mj2 without a *_face.mj2
   eye_with_face: frontal  # *_eye.mj2 when a *_face.mj2 exists
   face: angled            # *_face.mj2
   # values: frontal | angled | path to a config.yaml (or its folder; relative to this file)
   ```
   Rows whose model came from a yaml show "*" in the Model column. Typos give a MODEL_MAP_ERROR row.
3. **Scan the specific condition folder,** not all of `GCAMP_visual_exp`, so that older eye-only recordings are not queued. Use `--list` to print the table without running anything.
4. **Run** the GUI or the console. Per video:
   1. find the raw MJ2;
   2. FFmpeg → native H.264 MP4 (CRF 18), verifying geometry and frame count;
   3. write a runtime config copy `_runtime_<model>_<stamp>_config.yaml` **in the model folder**, with `project_path` and `pose_config_path` corrected (deleted after a clean run);
   4. `analyze_videos` → `filterpredictions` (median, window 5) → filtered CSV (made from the filtered H5 if needed) → `create_labeled_video` (filtered, pcutoff 0.1);
   5. likelihood QC from the raw H5;
   6. a row in the run CSV.

   Recordings are processed one at a time and are resumable: completed stages are skipped.
5. **QC** (see §5): read the Lik. column and the % low frames. For flagged videos, look at which points are low and at the frames or labeled MP4.

### Console flags

```
python DLC_batch_gui.py --no-gui --root <dir>
    [--config <frontal config.yaml>] [--angled-config <angled config.yaml>]
    [--shuffle N] [--angled-shuffle N] [--camera eye|face|both]
    [--only a,b] [--exclude a,b]                 # substring filters on the path
    [--include-completed] [--force-reconvert] [--force-reanalyze] [--dry-run]
    [--no-filter] [--no-labeled] [--pcutoff X] [--dotsize N] [--device D]
    [--list] [--frames DIR]
```

- `--root` overrides the saved `data_root`.
- Without `--config`, the frontal model comes from the saved settings and the angled model from `DEFAULT_ANGLED_CONFIG_PATH`.

**GUI:**
- **Top rows:** frontal and angled config, each with a shuffle setting (blank/auto = newest trained); data folder with Scan.
- **Columns:** sel, folder, recording, camera, model, source, mp4, h5, filt, csv, lab, Lik., status.
- **Buttons:** Select all / none / incomplete, Eye only, Face only, Show frames, Open folder, Play labeled video.
- **Options:** filter, window, labeled video, pcutoff, dot size, skeleton, skip completed, force re-convert, **Force re-analysis** (archives old outputs), dry run, device.

**Statuses and stages seen:**
- SUCCESS, ALREADY_COMPLETE, NEED_ANALYSIS (and other NEED_* stages), COMPLETE;
- MODEL_MISMATCH (outputs from the other model; skipped unless forced);
- MODEL_PREFLIGHT_FAILED (that model's config/shuffle/runtime config failed; other models still run);
- MODEL_MAP_ERROR (bad yaml); NO_MJ2; ARCHIVE_FAILED.

**Run-CSV columns include:** recording, camera, model, model_source, model_config, existing_model, status, likelihood_median, frames_low_likelihood, error, archived_files.

### Launchers (double-click in Explorer)

They were written because Claude, working through computer use, could double-click but not type into terminals. Both hard-code the env path, set `PYTHONUTF8=1` / `PYTHONIOENCODING=utf-8`, and `cd` to their own folder.

- **`pupil_check_frames.bat`:**
  - runs `--frames` on `%DATA%\5hz 40ms` and `%DATA%\visualstim_test\NavRearChr-2609-1`;
  - writes to `5hz 40ms\_camera_frames\` with `frames_log.txt`.
- **`pupil_run_5hz40ms.bat`:**
  - relaunches itself minimized (`start "pupil DLC batch" /min cmd /k`);
  - logs everything to `GCAMP_visual_exp\5hz 40ms\_DLC_batch_logs\console_run.txt`;
  - runs, in order:
    1. `--no-gui --root "%DATA%\5hz 40ms"`;
    2. `--root "%DATA%\visualstim_test\NavRearChr-2609-1" --only 261020`;
    3. same root `--only 261014_001 --force-reanalyze`.
  - Milestone lines: `[start]`, `[5hz 40ms finished, exit code N]`, `[test runs 10/20 finished, exit code N]`, `[test run 10/14 finished, exit code N]`, `[all done]`.
  - **Re-running it resumes step 1 but always redoes step 3.**
  - The 3–4 Nov run ended `[test run 10/14 finished, exit code 1]` / `[all done]` at 2026-11-04 02:17:52 local (the exit code 1 is the 261014_001 ARCHIVE_FAILED).

**Monitoring `console_run.txt` (~0.5 MB per video, ~8 MB total):**
- Split on `\r` and `\n`.
- Progress lines contain `/30000 [` and `it/s]`.
- Headers look like `[N/14] <animal>/<rec>   eye camera -> frontal model`.
- Each video ends with `Likelihood check: median X; frames with a point < 0.6: Y %`.

### Building or extending the angled model (`Pupil_9P_angled_261016\scripts\`)

1. `angled_pupil_autolabel.py --video "<...>_eye_native.mp4" --project "D:\2photon\Chao\GCAMP_visual_exp\Pupil_9P_angled_261016"`
   - Options: `--scorer`, `--n-kmeans 60`, `--n-uniform 12`, `--n-hard 8`, `--seed 0`, `--copy-video`, `--no-config-update`, `--dry-run`.
   - Picks frames and auto-labels them from a bright-pupil ellipse fit.
2. `step1_review_labels.py`: napari review, on the labeling PC.
3. `step2_train.py [--shuffle 2]`: trains from the frontal snapshot.
   - Passes `WeightInitialization(..., with_decoder=False)`, then patches `pytorch_config.yaml` to `with_decoder: true`.
   - Calls the PyTorch `train_network` API directly with 150 epochs, batch 8, milestones [100, 130].
   - Trained on samson's CPU in 57 min.
4. `step3_test_model.py`: whole-video check.
5. `qc_dlc_output.py`: GOOD / CHECK / RETRAIN health check for any DLC output.

### Legacy: `DLC_9point_selected_conditions.py` (frontal only)

- **Picker:** a tkinter condition picker that opens at `D:\2photon\Chao\GCAMP_affrerents` and can navigate up.
- **`PREPROCESS_MODE`:** `native` (default; output `<rec>_eye_native.mp4`) or `upscale`, for the old 160×112 videos. Upscale scales to 640×448 bicubic and retimes to 15.62 fps (output `<rec>_eye_640x448.mp4`).
- **Pause/resume:** P/R keys or a `PAUSE_DLC.flag` file beside the script.
- **Exit codes:** 0 OK, 1 failures, 2 no MJ2, 130 Ctrl-C.

## 5. Key decisions and parameters

| Decision / parameter | Reason |
|---|---|
| Feed current videos at **native resolution**: no rescale, reframe or retime (MJ2 → H.264 MP4, CRF 18 from the `h264_crf` setting; the old script also used `-vsync passthrough`). | The user's request on 25 Sep: new videos don't need rescaling. MJ2 must still be converted because DLC can't read it. |
| The raw MJ2 is the only source; it is never modified, moved or deleted. | Protects raw data. Old outputs made from unconverted `*_eye.mp4` were wrong input. |
| **Camera → model:** `_eye` → frontal, with or without `_face`; `_face` → angled (`CAMERA_MODEL_DEFAULTS = {eye: frontal, eye_with_face: frontal, face: angled}`). | The user's explicit rule (~2 Nov). Confirmed from frames for all 7 "5hz 40ms" recordings on 3 Nov. |
| **Never run the frontal model on angled views.** | It gives median likelihood 0.02–0.5 there (0.13 on 261014_001). |
| Per-folder `pupil_models.yaml` overrides (nearest wins). Unknown or empty values raise an error. | Exceptions (test runs, reversed cameras) without hard-coding, and never a silently wrong model. |
| `visualstim_test\NavRearChr-2609-1\pupil_models.yaml` = `eye: angled` (written by Claude on 3 Nov). | The frames proved these single-camera test runs are the angled camera saved as `_eye`. |
| MODEL_MISMATCH skip; Force re-analysis **archives** old outputs (now with `os.replace`), never deletes. | Never mix models in one folder; protect earlier results; no duplicate on a refused move. |
| The model's `config.yaml` is never edited; a runtime copy is used. Shuffle auto = newest trained snapshot. | `project_path` is machine-specific (DLC rewrites it on the labeling PC, which sees the data as `Z:`), so the copy sets `project_path` / `pose_config_path` for samson. The old script edits the copy by a ruamel round-trip and validates it with DLC's strict loader, so a broken YAML fails early. |
| Filter: median, window 5. Labeled video: pcutoff 0.1, dotsize 3 (the user's saved settings). The filtered CSV is made from the filtered H5 with pandas (`ensure_filtered_csv`; the old script always did this). | Downstream pipelines read the classic scorer/bodyparts/coords CSV header, and DLC sometimes fails to write the CSV. |
| **QC:** Lik. = median of per-point median likelihoods from the raw H5, flagged "!" when < 0.6. Also report % frames with any point < 0.6. | A quick check that the model fits the camera view. A high % with correct visible points is a data property (lid, iris contrast), not a failure. |
| Pupil size (recommended, **not yet implemented**): fit an ellipse to points with likelihood ≥ 0.6. In the angled view use the horizontal diameter (pupil_1–pupil_5) or the ellipse's major axis. | Lower/right points 5–8 are often low on the close-up eye camera. The angled view's vertical axis is foreshortened. |
| Blinks: mask frames where any edge point < 0.6, plus ±4 frames. | The angled model keeps tracking a half-covered pupil with high likelihood. |
| If `qc_dlc_output.py` says CHECK/RETRAIN on a new angled recording: add ~30 frames (autolabel → review → `step2_train.py --shuffle 2`). | The angled model was trained on one animal/session only. |
| Scan one condition folder at a time. | Avoids queueing old eye-only recordings. |

## 6. Current status (as of 2026-11-09T22:30:48.301883Z)

| Dataset | Status |
|---|---|
| `GCAMP_visual_exp\5hz 40ms\NavRearChr-2609-1` (7 recordings × eye + face) | **14/14 SUCCESS**, finished 2026-11-03 23:55 local. Run CSV `5hz 40ms\_DLC_batch_logs\DLC_batch_log_20261103_050813.csv`. Table below. |
| `visualstim_test\NavRearChr-2609-1` 261020_001 | SUCCESS, angled (via yaml), median 0.98, ~0 % low. Run CSV for both 10/20 runs: `visualstim_test\NavRearChr-2609-1\_DLC_batch_logs\DLC_batch_log_20261103_235538.csv`. |
| `visualstim_test\NavRearChr-2609-1` 261020_002 | SUCCESS, angled, median 0.89 (36 % low). The last ~6 min are bad: all 9 points drop together (eye closing or reflection lost?). **Check the labeled video before using that stretch.** |
| `visualstim_test\NavRearChr-2609-1` 261014_001 | **Not re-analysed.** The forced redo on 2026-11-04 02:17 local failed with WinError 5 (Access is denied) while archiving. It still has only the old 10/16 frontal outputs (median 0.13, unusable; that run's log is `visualstim_test\_DLC_batch_logs\DLC_batch_log_20261016_175228.csv`), shown as MODEL_MISMATCH. The failed redo's run CSV is `visualstim_test\NavRearChr-2609-1\_DLC_batch_logs\DLC_batch_log_20261104_021737.csv` (status ARCHIVE_FAILED). Angled-model outputs for it exist in `Pupil_9P_angled_261016\_test_output\`. Leftovers: `...\NavRearChr-2609-1_261014_001\_DLC_archive_20261104_021737\` (a duplicate of one .h5) and `Pupil_9P_angled_261016\_runtime_Pupil_9P_angled_261016_20261104_021737_config.yaml`. |
| `visualstim_test\NavWT-2609-1` 261020_001/002 | Not analysed. Cameras reversed; no yaml there. |
| `GCAMP_affrerents\NavRearChrTIGRE-2608-1` | Earlier recordings were analysed with the frontal model before Nov (per the session notes; not re-checked). `NavRearChrTIGRE-2608-1_261101_001` (eye-only, 258 MB MJ2) had **no DLC outputs** when listed on 3 Nov. |
| `GCAMP_affrerents\NavRearChrTIGRE-2609-1` day 261025 | Listed as pending in samson's `CLAUDE.md`; not run in this session. |
| Older visual-exp conditions (`2hz test`, `2hz_wheellock`, `5hz`, `50ms test`, `old-200umfiber`) | Analysed with the frontal model before Nov (per memory; "works well"); not re-checked here. The batch logs in the `GCAMP_visual_exp` root date from 29 Jul to 13 Aug 2026 (`DLC_new_model_rerun_log_*`, `DLC_9point_upscaled_mj2_log_*`, `DLC_9point_selected_conditions_log_*`). |

"5hz 40ms" results: median likelihood (% frames with any point < 0.6):

| Recording | eye → frontal | face → angled |
|---|---|---|
| 261027_001 | 0.91 (8 %) | 0.85 (4 %) |
| 261027_002 | 0.89 (81 %) | 0.93 (0 %) |
| 261029_001 | 0.88 (8 %) | 0.93 (0 %) |
| 261029_002 | 0.81 (56 %) | 0.90 (2 %) |
| 261029_003 | 0.86 (26 %) | 0.89 (2 %) |
| 261031_001 | 0.89 (19 %) | 0.89 (0 %) |
| 261031_002 | 0.85 (51 %) | 0.89 (2 %) |

**Low eye frames come from points 5–8:**
- 261027_002: per-point medians p6 0.42, p7 0.36, p8 0.59; the lower lid covers the edge.
- 261029_002: p5–p7 ~0.62; low contrast against the iris.
- In both, the confident points sit on the pupil edge (checked on frames).

**Code and docs:**
- `DLC_batch_gui.py` on samson is the 4 Nov version.
- samson's `CLAUDE.md` and `docs\dlc-pupil-history.md` **wrongly say 261014_001 was re-analysed**. Correct them. They also send readers to `docs\claude_notes\deeplabcut\angled_pupil_model_notes.md`, which does not exist on samson (§2).

## 7. In flight and next steps

- **Running jobs:** none. **Scheduled check-ins:** none pending; all 7 `send_later` check-ins of the 3–4 Nov run have fired, the last at 2026-11-04 08:45 UTC.
- **The Cowork session is still running.** At capture it was writing its own export: branch `cowork-export/dlc-behavior` of `cwei6/GCAMP`, or the private artifact "dlc-behavior export" if the push is refused. Its `gh` token was invalid. When the export arrives, merge it here: it supersedes this draft. Then re-sync from the cursor in `sessions/dlc-batch-gui/README.md`.
- **Waiting on the user:**
  - Redo 261014_001 with the angled model, from the **Administrator** Anaconda Prompt:
    ```
    cd /d D:\Behavior\Chao\squint_pipeline
    python DLC_batch_gui.py --no-gui --root "D:\2photon\Chao\GCAMP_visual_exp\visualstim_test\NavRearChr-2609-1" --only 261014_001 --force-reanalyze
    ```
    Or use the GUI with Force re-analysis ticked.
  - Approve deleting the two leftovers (§6).
  - Decide on 261020_002's last ~6 min.
- **Unfinished work:**
  - Pupil-size extraction (ellipse fit on points ≥ 0.6; horizontal diameter for the angled view) is recommended but not written.
  - `qc_dlc_output.py` has not been run on the new "5hz 40ms" face videos (only the likelihood QC was used). Run it, because the angled model was trained on one animal.
  - Before Rodrigo's folder is made read-only (as samson's `CLAUDE.md` / switch plan want), change `make_runtime_config` to write the runtime copy somewhere local (e.g. `squint_pipeline\_runtime\`). Otherwise every frontal video fails preflight.
  - Generalise the two .bat launchers (hard-coded env path and folders) before reuse.
  - Fix the 261014_001 statement in samson's docs, and their pointer to the missing `docs\claude_notes\deeplabcut\` (§2).
- **New data to run:** `GCAMP_affrerents\NavRearChrTIGRE-2608-1_261101_001`, `NavRearChrTIGRE-2609-1` day 261025, and any two-camera days after 10/31. Check the camera views first.

## 8. Open questions for the user

1. 261014_001: redo it elevated, or use the angled outputs already in `Pupil_9P_angled_261016\_test_output\`?
2. 261020_002: drop the last ~6 min, or inspect the labeled video first?
3. Analyse the NavWT-2609-1 10/20 two-camera tests? They need a `pupil_models.yaml` there with `eye_with_face: angled` and `face: frontal`.
4. OK to delete `_DLC_archive_20261104_021737` (in 261014_001) and `_runtime_Pupil_9P_angled_261016_20261104_021737_config.yaml`?
5. The user never replied to Claude's reading of "for the first few recordings, it is the opposite way around". The frames show the 10/27–10/31 recordings follow the rule, and only the test runs are reversed. Does the user agree?
6. Is the frontal model precise enough on the new close-up eye camera, where the lower edge is often low-confidence? Or should close-up frames be added to a model?
7. Which DLC env is canonical: `C:\ProgramData\anaconda3\envs\DEEPLABCUT` (used by the launchers) or `C:\Users\cwei6\anaconda3\envs\deeplabcut` (samson's `CLAUDE.md`)?
8. What made samson drop to ~1 fps on the evening of 3 Nov? Should MATLAB, Jupyter and the squint GUI be closed during DLC runs?
9. Where does `DLC_9point_selected_conditions.py` live on samson, and is it still needed? Are any old 160×112 recordings still waiting (they would need upscale mode)?
10. Going forward, which copy of the code is canonical: samson's `squint_pipeline` folder (and its local Claude Code / git, if set up) or this repo's `deeplabcut/code/`?

## 9. Working with this user

- **Messages:** short, informal, with typos (e.g. "recoridng", "qeustion"). Infer the intent; ask only when it really matters.
- **Not a coder:** a neuroscientist who describes goals; Claude writes and runs the code. Explain in plain language, with tables of per-video numbers.
- **Options first:** for anything non-trivial, give 2–3 options with trade-offs and a recommendation, then wait. This is from the user's `CLAUDE.md` on samson.
- **Results to look at:** labeled frames, camera contact sheets, labeled videos, figures, plus full Windows paths and copy-paste commands for the Administrator Anaconda Prompt.
- **Keep existing behaviour:** "if a animal folder only has the previous eye video, then nothing need to be changed".
- **Unattended runs:** happy to let Claude run long jobs on samson while away ("I am away from the pc, can you help me check and run the analysis?"; "set up the computer use/access request first? I am Going to rest").
  - Wants short status updates (progress, speed, ETA, QC).
  - Messages only on failure or a big timing change, plus milestone summaries.
  - Asks briefly: "How is the status now?"
  - Give times in samson local time (UTC-5).
- **Never delete, rename or overwrite** videos, models or earlier results without asking. Rodrigo's folder is read-only. Edit code in place (no zips or new copies).
- **Wants tools that recognise things themselves,** because several experiments run in parallel.

## 10. Pitfalls and known issues

- **Windows permissions:** files written from an Administrator prompt may not be movable by a non-elevated run (WinError 5). Run forced re-analysis from the same kind of prompt. The archive now uses `os.replace` with a clear error. The old `shutil.move` left a duplicate.
- **Runtime config written into the model folder:** this breaks if Rodrigo's folder becomes read-only (§7).
- **CPU-only and shared:** speed swings from 34 to ~1 fps. Check what else is running before promising an ETA.
- **DLC output naming:** DLC drops the leading zero of the snapshot (`snapshot-best-070.pt` → `best-70`). Output ownership is matched by the model token, lower-case `Task + date` (`pupil_9p_wb_260402`, `pupil_9p_angled261016`).
- **What counts as COMPLETE:** raw H5, filtered H5, filtered CSV, and a labeled MP4 of matching geometry. Deleting only the labeled MP4 redoes only that step.
- **Large logs:** `console_run.txt` is large; parse it, don't print it.
- **Restarting the launcher:** re-running `pupil_run_5hz40ms.bat` always redoes 261014_001 (`--force-reanalyze`).
- **Old fps tags:** old 160×112 MJ2s report 30 fps, but the true rate is 15.62.
- **Remote-device bridge** (Cowork, and probably Remote Control):
  - 400 MB per-file staging limit, so MJ2s can't be pulled; extract frames on samson with `--frames`.
  - A commit from a reused staging path once wrote stale content. Verify the size, mtime and md5 on samson after every write.
  - Writes need the current `expectedMtimeMs`.
  - A recursive listing of `GCAMP_visual_exp` truncates at 2,000 entries; list subfolders.
  - A file held open by a running process can be read but not written.
- **Computer use on samson:** Command Prompt and File Explorer were click-tier only (no typing), hence the .bat launchers. "Windows Terminal" and "Python" are not valid app names.
- **Testing:** the GUI needs tkinter. The Cowork container had it only for `/usr/bin/python3.12`. The 3 Nov samson run was the first real inference with the two-model script; earlier tests used a fake DLC module.
  - That test kit (`fake_dlc3/deeplabcut/__init__.py`, driven by env vars `FAKE_LOW_LIK`, `FAKE_DLC_CALLS`, `FAKE_DLC_FAIL`, `FAKE_DLC_SLEEP`; mock trees built by `build_mock_pupil.py`; a Python 3.10 venv for compile checks) lived only in the Cowork container's scratchpad. It is **not** in this repo or the export, so rebuild a stub `deeplabcut` package on `PYTHONPATH` to test without DLC.
  - Regression pattern the session used: `PYTHONPATH=<stub> python DLC_batch_gui.py --no-gui --root <mock 2photon tree> --config <mock frontal config> --angled-config <mock angled config> [--only 261027_001 --force-reanalyze]`. Keep the code compiling on Python 3.10.
- **Recording folders hold imaging outputs too** (`*_bin31f_1p98s_behavior.csv`, `_check_frames.tif`, AQuA2 outputs, ...). They belong to the gcamp/afferents pipelines; don't touch them.

## 11. Sources

- **Cowork session** `session_01TXsqNvay4sMGNwGK4XyuVe`, "DLC batch video analysis GUI". The transcript was read to event `b3880823-a8ba-421c-9f2f-0cd6c89bec65` (2026-11-09 22:30:48 UTC).
- **`sessions/dlc-batch-gui/`:**
  - `README.md`: session facts and sync cursor;
  - `compaction-summary.md`: the session's own summary, 25 Sep – 3 Nov;
  - `user-messages.md`: the user's messages since then;
  - `digest.md`: dated digest, 3–9 Nov;
  - `recovered/`: the session's own `HANDOFF.md`, `batch_2026-11-03_results.md`, `angled_pupil_model_notes.md` and `pupil_models.example.yaml`, rebuilt from the transcript (see `recovered/RECOVERED.md`).
