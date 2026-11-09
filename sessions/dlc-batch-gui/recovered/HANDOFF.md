# HANDOFF: pupil DeepLabCut (2-photon eye / face videos)

Snapshot exported 9 Nov 2026 from the Cowork session that built and ran this tooling (Sept–Nov 2026).
It is for a fresh Claude Code session with none of that context.
Code here is the version on samson on 9 Nov. `DLC_batch_gui.py` matches `D:\Behavior\Chao\squint_pipeline\DLC_batch_gui.py` byte for byte.

## 1. Goal

Track the pupil with a 9-point DeepLabCut model in the camera videos of the user's 2-photon GCaMP recordings:
- **Points:** 8 edge points plus the centre (`pupil_9`), so the radius survives partial occlusion.
- **Where the output goes:** the filtered CSVs feed the imaging pipelines (pupil as an arousal control and as a trigger).
- **The main tool:** `code/DLC_batch_gui.py`, a tkinter GUI plus console mode. For every raw camera video it:
  - converts MJ2 → MP4 at native size (no rescale, reframe or retime; this was decided);
  - runs DLC `analyze_videos` → `filterpredictions` (median, window 5) → filtered CSV → `create_labeled_video` as a visual sanity check;
  - reports a likelihood QC.
- **Two models since Nov 2026:** a frontal one (Rodrigo's) and an angled one (built in Oct 2026). Each video gets the model that matches its camera.

## 2. Data layout on samson (data server, drive D:)

- **Afferent experiments:** `D:\2photon\Chao\GCAMP_affrerents\<animal>\<animal>_<YYMMDD>_<NNN>\<rec>_eye.mj2`.
  - Only one camera (frontal), so they use the frontal model.
  - The same folder holds `.sbx`, `.mat`, `_quadrature.mat` and `suite2p\`.
  - Older data used condition → animal → recording.
- **Visual-stim experiments:** `D:\2photon\Chao\GCAMP_visual_exp\<condition>\<animal>\<recording>\`.
  - Conditions include `2hz test`, `2hz_wheellock`, `5hz`, `50ms test`, `old-200umfiber`, `5hz 40ms`, plus `visualstim_test\`.
  - **Two-camera setup (from 10/27, `5hz 40ms\NavRearChr-2609-1`):** each recording has `<rec>_eye.mj2`, a **frontal close-up** (the "previous eye camera") → frontal model. It also has `<rec>_face.mj2`, the new **angled side view** → angled model.
  - **Schedule:** recordings every other day (10/27, 10/29, 10/31), with a test run one week before.
  - **Test runs (`visualstim_test\NavRearChr-2609-1`, 261014_001, 261020_001, 261020_002):** only the angled camera, saved as `_eye`. A `pupil_models.yaml` with `eye: angled` sits in that animal folder.
  - **`visualstim_test\NavWT-2609-1` 261020_001/002:** short two-camera tests with the cameras the other way round (`_eye` angled, `_face` frontal). Not analysed.
- **Videos:** current videos are 240×180 grey JPEG2000 MJ2 at 30 fps (~30 000 frames ≈ 16.7 min, 400–570 MB each). The old videos were 160×112, tagged 30 fps but really 15.62 fps, and were upscaled to 640×448. The older script `DLC_9point_selected_conditions.py` still has that "upscale" mode; it is not exported here and remains in the Claude project docs.
- **Outputs (next to the MJ2):** `<rec>_<cam>_native.mp4`, then `<rec>_<cam>_nativeDLC_HrnetW48_<Task><date>shuffle1_snapshot_best-NN{.h5,_filtered.h5,_filtered.csv,_full.pickle,_meta.pickle,_filtered_p10_labeled.mp4}`.
- **Logs:** `<scanned root>\_DLC_batch_logs\DLC_batch_log_<stamp>.csv`, plus `console_run.txt` when started with the launcher.
- **Archives:** forced re-analysis moves old outputs to `<recording>\_DLC_archive_<stamp>\`. Folders starting with `_` are never scanned.

## 3. Models and the camera → model rule

| Alias | config.yaml | Notes |
|---|---|---|
| frontal | `D:\2photon\Rodrigo\DLC_Models-RN\Pupil_9P_WB_260402\config.yaml` | Colleague's (Rodrigo's) folder. Treat it as read-only. HRNet-W48, shuffle 1, best-070. Scorer `..._Pupil_9P_WB_260402shuffle1_snapshot_best-70` |
| angled | `D:\2photon\Chao\GCAMP_visual_exp\Pupil_9P_angled_261016\config.yaml` | Built Oct 2026 from the frontal model's weights. Trained on one animal and one session (261014_001). Test RMSE 1.06 px. Scorer `..._Pupil_9P_angled261016shuffle1_snapshot_best-90` |

- **Points:** `pupil_1` is the left-most edge point, then clockwise on screen (2 upper-left … 5 right … 8 lower-left); `pupil_9` is the centre.
- **Default rule:**
  - `*_eye.mj2` alone → frontal.
  - `*_eye.mj2` with a `*_face.mj2` next to it → frontal.
  - `*_face.mj2` → angled.
- **Overrides:** a `pupil_models.yaml` (see `code/pupil_models.example.yaml`) in a condition, animal or recording folder overrides the rule. The nearest one wins; the keys are `eye`, `eye_with_face` and `face`. Overridden rows show "*" in the GUI.
- **Never use the frontal model on angled videos.** It gives median likelihoods of 0.02–0.5 there; on 261014_001 it gave 0.13.
- **MODEL_MISMATCH:** output names contain the model token, so the GUI can tell which model made existing outputs. A video whose outputs come from the other model is shown as MODEL_MISMATCH and skipped. "Force re-analysis" archives those outputs and redoes the video.

## 4. Environment

- **samson:**
  - Windows 10, no GPU, so DLC inference runs on the CPU.
  - napari doesn't work there, so frame labelling is done on the other PC, which sees samson's D: as `Z:`, with `python -m deeplabcut` (DLC 3.0.0rc13).
  - samson is shared: touch nothing outside the user's folders.
- **DLC environment:** DeepLabCut 3.0.0rc13 (PyTorch) in a conda env.
  - The launchers call `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe` directly and put its folders on PATH. This worked for the whole 3–4 Nov run, and the squint app remembers the same path in `squint_dlc_settings.json`.
  - The new `CLAUDE.md` on samson names `C:\Users\cwei6\anaconda3\envs\deeplabcut`. Check which one exists before relying on it.
  - FFmpeg and ffprobe come from the env's `Library\bin`.
- **How the user runs things:** from an *Administrator: Anaconda Prompt*, with `(deeplabcut)` active, in `D:\Behavior\Chao\squint_pipeline`: `python DLC_batch_gui.py`. The GUI remembers its settings in `DLC_batch_gui_settings.json` (see `code/DLC_batch_gui_settings.example.json`).
- **Speed (HRNet-W48 on CPU):**
  - About 26–34 frames/s on an otherwise idle samson.
  - About 7–9 frames/s when other jobs run (suite2p, MATLAB/AQuA, Jupyter); it once fell to ~1 frame/s.
  - Budget about 65 min per 30 000-frame video when samson is busy.

## 5. How to run

**GUI:** `python DLC_batch_gui.py`.
- **Top rows:** the frontal and angled `config.yaml`, each with a shuffle setting (auto = newest trained), and the data folder with Scan.
- **Table columns:** Camera, Model ("*" = from a `pupil_models.yaml`), MJ2 source, MP4/H5/Filt/CSV/Labeled ticks, Lik. ("!" = median below 0.6) and Status.
- **Buttons:** Select all / none / incomplete, Eye only, Face only, Show frames (a frame of every camera of the highlighted recordings, labelled with its model), Open folder, Play labeled video.
- **Options:** filter, window, labeled video, pcutoff, dot size, skeleton, skip completed, force re-convert, force re-analysis (archives old outputs), dry run, device.

**Console:** `python DLC_batch_gui.py --no-gui --root <folder> [options]`
- `--list`: scan and print the table, run nothing.
- `--frames <dir>`: save 2 PNG frames per camera video (to check which camera is which), then stop.
- `--camera eye|face|both`
- `--only a,b` / `--exclude a,b`: substring filters on the path.
- `--config` / `--angled-config` / `--shuffle` / `--angled-shuffle`
- `--include-completed`, `--force-reconvert`, `--force-reanalyze`, `--dry-run`, `--no-filter`, `--no-labeled`, `--pcutoff`, `--dotsize`, `--device`

The scanner:
- walks up to `MAX_SCAN_DEPTH` levels;
- skips folders starting with `_` or `.`, skip-listed folders and DLC project folders;
- finds `<recording>\<recording>_eye.mj2` / `_face.mj2`.

**Launchers (double-click):**
- `code/pupil_check_frames.bat`: runs `--frames` for "5hz 40ms" and the NavRearChr-2609-1 test runs.
- `code/pupil_run_5hz40ms.bat`: the 3–4 Nov run. It relaunches itself minimized and writes everything to `5hz 40ms\_DLC_batch_logs\console_run.txt`.

Both hard-code the env path and should be generalised before reuse. They were used because Claude, working through computer use, could double-click files but not type into terminals.

**Building the angled model:** `code/angled_model_scripts/`
- `angled_pupil_autolabel.py`: picks frames and auto-labels them from a bright-pupil ellipse fit.
- `step1_review_labels.py`: napari review.
- `step2_train.py`: trains initialised from the frontal model (with_decoder trick).
- `step3_test_model.py`
- `qc_dlc_output.py`: GOOD / CHECK / RETRAIN health check.

The recipe and history are in `docs/angled_pupil_model_notes.md`.

## 6. QC rules

- Look at the Lik. column, the median of the per-point median likelihoods from the raw H5. Below 0.6 the model probably doesn't fit that camera view.
- Also look at the % of frames with any point < 0.6, then check *which* points are low (per-point likelihood from the H5) and look at frames. A high % with correct placement on the visible edge is a data property (lid, contrast, blink), not a model failure.
- With a new setup or a new recording type, check the camera views first (`--frames` or Show frames) before spending CPU hours.
- **Blinks:** mask frames where any edge point is < 0.6, plus ±4 frames.
- **Pupil size:**
  - Fit an ellipse to the points with likelihood ≥ 0.6.
  - In the angled view, use the horizontal diameter (pupil_1–pupil_5) or the ellipse's major axis, because the vertical axis is foreshortened.
- If `qc_dlc_output.py` says CHECK/RETRAIN for a new angled recording: add ~30 frames from it (autolabel → review → `step2_train.py --shuffle 2`).

## 7. Status on 9 Nov 2026

The full table is in `docs/batch_2026-11-03_results.md`.
- **5hz 40ms (7 recordings × 2 cameras):** 14/14 SUCCESS.
  - Eye (frontal): medians 0.81–0.91. Many low frames in 261027_002 (81 %), 261029_002 (56 %) and 261031_002 (51 %), from the lower-edge points; placement is OK.
  - Face (angled): 0.85–0.93, ≤ 4 % low.
- **Test runs:**
  - 261020_001 → 0.98.
  - 261020_002 → 0.89, but its last ~6 min are lost (all points drop together).
- **261014_001: not re-analysed (ARCHIVE_FAILED, WinError 5).** It still has only the bad frontal outputs. `CLAUDE.md` / `docs\dlc-pupil-history.md` on samson wrongly say it was redone. Please correct them.
- **Leftovers (delete only after the user says so):**
  - `Pupil_9P_angled_261016\_runtime_Pupil_9P_angled_261016_20261104_021737_config.yaml`
  - `visualstim_test\NavRearChr-2609-1\NavRearChr-2609-1_261014_001\_DLC_archive_20261104_021737\`, which holds a copy of one H5.
- **Earlier batches:** the GCAMP_affrerents and the older visual_exp conditions were analysed with the frontal model before Nov; see the `DLC_9point_selected_conditions_log_*.csv` and `_DLC_batch_logs` in those roots.

## 8. Open questions for the user

1. 261014_001: redo it from the Administrator prompt (command in the results doc), or use the angled outputs in the model's `_test_output\`?
2. 261020_002: drop or inspect the last ~6 min?
3. Analyse the NavWT-2609-1 10/20 two-camera tests? They need `eye_with_face: angled` and `face: frontal`.
4. OK to delete the two leftovers above?
5. New recordings to run, e.g. afferent NavRearChrTIGRE-2609-1 day 261025 (listed as pending in `CLAUDE.md`), and any new two-camera days after 10/31.
6. The frontal model works on the new close-up eye camera, but the lower edge is often low-confidence. Is that precise enough, or should close-up frames be added to a model?
7. Which DLC env path is canonical (section 4)?

## 9. How the user likes to work

- Not a coder: they describe goals and ask questions; Claude writes and runs the code. Explain in plain language.
- For anything non-trivial, give 2–3 options with trade-offs and a recommendation, then wait for their choice.
- Show results they can look at (labeled frames, figures, tables) and give paths they can copy.
- Edit code in place in `D:\Behavior\Chao\squint_pipeline`: no zips, no new copies.
- Never delete, rename or overwrite videos, models or earlier results without asking. Rodrigo's folder is read-only.
- They like GUIs that also have a console mode, long runs that go unattended, and short status updates with result tables.

## 10. Known pitfalls

- **Rodrigo's folder being read-only will break the frontal model.** `DLC_batch_gui.py` (like `DLC_batch_gui_squint.py`) writes a temporary `_runtime_<model>_<stamp>_config.yaml` *into the model folder*, with `project_path` corrected for samson, and deletes it after a clean run. It does this in Rodrigo's folder too. The new `CLAUDE.md` asks for that folder to be read-only, and the switch plan wants model folders protected; with that in place, every frontal video fails at preflight (MODEL_PREFLIGHT_FAILED). Before making it read-only, change `make_runtime_config` to write the runtime copy somewhere local (e.g. `squint_pipeline\_runtime\`, with `project_path` pointing at the model folder). Then test that DLC 3.0.0rc13 still finds the snapshots.
- **Windows permissions:** files written from the Administrator prompt may not be movable or deletable by a non-elevated run (WinError 5). Run forced re-analysis from the same kind of prompt.
- **CPU contention:** check what else is running before promising an end time.
- **DLC naming:** DLC drops the leading zero of the snapshot (`snapshot-best-070.pt` → `best-70` in output names). Outputs are matched by the model token, lower-case `Task + date`.
- **What counts as complete:** a video is COMPLETE only with raw H5, filtered H5, filtered CSV and a labeled MP4 of the right geometry. Deleting the labeled MP4 makes it redo only that step.
- **Large logs:** tqdm progress makes `console_run.txt` large (~0.5 MB per video). Parse it; don't print it.
- **Staging limit (Cowork only):** with the device bridge, files over 400 MB (every MJ2) can't be pulled, so extract frames on samson with `--frames` instead.
- **Stale file writes (Cowork only):** writing a file back from a reused staging path once wrote stale content. Verify the size on samson after any write.
