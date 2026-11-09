# Angled-camera 9-point pupil model (Pupil_9P_angled_261016)

New DLC 3 (PyTorch) project created 2026-10-16 at
`D:\2photon\Chao\GCAMP_visual_exp\Pupil_9P_angled_261016\` for the visual-stim
setup where the eye camera looks at the eye from an angle. First recording of
that setup: `GCAMP_visual_exp\visualstim_test\NavRearChr-2609-1\NavRearChr-2609-1_261014_001`
(240x180 @ 30 fps, 32 500 frames).

## Status (2026-10-20): trained and validated — in use
* Labels reviewed in napari-deeplabcut on the remote PC (DLC 3.0.0rc13 GUI,
  project reached via network share; DLC rewrites `project_path` per machine).
  Final CollectedData_Chao: 79 frames, all 9 points (closed-eye frame left empty).
* Trained on the data server (no GPU, CPU, 57 min): shuffle 1, hrnet_w48,
  initialised from `Pupil_9P_WB_260402` snapshot-best-070 (backbone + heads),
  150 epochs, batch 8. Best snapshot `snapshot-best-090.pt`.
  Evaluation: train RMSE 1.04 px, test RMSE 1.06 px, test mAP 89 / mAR 92.5.
* Whole-video check (step3): all 9 points median likelihood 0.96–1.00; 35 frames
  (0.11 %) below 0.6 = the 24 blinks; 99.6 % of non-blink frames >= 0.85;
  DLC pupil_1–pupil_5 diameter vs independent blob tracker r = 0.98, constant
  offset +1.4 px (labels sit on the outer edge, blob tracker on the half-max).
  Frame-to-frame jitter of the raw DLC centre ~0.9 px (median filter in the
  batch pipeline reduces it). Verdict: good; use via DLC_batch_gui.py.
* Trained on one animal / one session only — run `scripts\qc_dlc_output.py`
  on each new recording; if it says CHECK/RETRAIN, add ~30 frames of that
  recording (autolabel -> review -> `step2_train.py --shuffle 2`).
* Blink handling: model keeps tracking a half-covered pupil with high
  likelihood, so mask frames with any edge point < 0.6 plus +/-4 frames.
* rc13 quirks seen: `deeplabcut.train_network` wrapper uses `displayiters`
  (script now calls the PyTorch API directly); evaluation image plotting warns
  "DataFrame reshape failed ... 2 individuals" for single-animal projects —
  cosmetic, metrics unaffected.

## Why the frontal model fails there
* The pupil is a BRIGHT (retro-reflecting) disk in these videos; the 9-point
  model = 8 edge points + centre (`pupil_9`).
* Angled view: pupil is a small flat ellipse (~20 x 9 px, long axis tilted ~+10°)
  right under the upper eyelid. Frontal model `Pupil_9P_WB_260402`
  (HrnetW48, shuffle1, snapshot_best-70) gives median likelihoods 0.02–0.5 per
  point on this video -> unusable.

## Point convention (verified against the frontal project's labels)
pupil_1 = left-most edge point, then clockwise on screen (2 upper-left, 3 top,
4 upper-right, 5 right, 6 lower-right, 7 bottom, 8 lower-left); pupil_9 = centre.

## What the project contains
* `config.yaml` — Task Pupil_9P_angled, scorer Chao, date 261016, hrnet_w48,
  TrainingFraction 0.9, skeleton = octagon 1-2-...-8-1, snapshotindex -1
  (DLC lists the "best" snapshot last, so -1 = best, like the old project).
* `labeled-data\NavRearChr-2609-1_261014_001_eye_native\` — frames + labels
  (72 auto-labeled by bright-blob ellipse fit, 8 blink frames hand-labeled).
* `_qc\*_autolabel_01..04.png` contact sheets; `_qc\*_blobtrack.csv` = independent
  per-frame pupil ellipse (cx, cy, major, minor, tilt) from the blob tracker.
* `_test_output\` — step3 results on the training video (labeled mp4, h5/csv,
  check plot).
* `scripts\angled_pupil_autolabel.py` (frame picking + auto labels for any new
  angled video), `step1_review_labels.py` (napari GUI), `step2_train.py`,
  `step3_test_model.py`, `qc_dlc_output.py` (health check of a DLC output:
  GOOD / CHECK / RETRAIN), `README.md`.

## Training recipe (step2_train.py; validated on DLC 3.0.2 and run on 3.0.0rc13)
1. `create_training_dataset(config, Shuffles=[1], net_type="hrnet_w48",
   weight_init=WeightInitialization(snapshot_path=<old snapshot-best-070.pt>, with_decoder=False))`
2. patch `dlc-models-pytorch\...\train\pytorch_config.yaml`
   `train_settings.weight_init` -> `with_decoder: true`,
   `conversion_array: [0..8]`, `bodyparts: [pupil_1..pupil_9]`
   (with_decoder=True directly in create_training_dataset is reserved for
   SuperAnimal fine-tuning in DLC 3.0.x; PoseModel.build then loads backbone +
   heads from the old snapshot — verified: all 1834 tensors identical).
3. `pose_estimation_pytorch.apis.train_network(config, shuffle=1, epochs=150,
   batch_size=8, save_epochs=25, display_iters=50,
   pytorch_cfg_updates={"runner.scheduler.params.milestones": [100, 130]})`
4. `evaluate_network(config, Shuffles=[1], plotting=True)`.
Old project (`D:\2photon\Rodrigo\DLC_Models-RN\Pupil_9P_WB_260402`) is read
directly by the script (also found via a network share path);
`--include-old-labels` merges its labeled frames.

## Usage afterwards: DLC_batch_gui.py (two-camera version, 2026-11-02)
`D:\Behavior\Chao\squint_pipeline\DLC_batch_gui.py` holds both models
("Frontal model config.yaml" + "Angled model config.yaml" rows) and picks one per video:
* `*_eye.mj2`, with or without a `*_face.mj2` next to it -> frontal `Pupil_9P_WB_260402`
  (eye-only recordings behave exactly as before).
* `*_face.mj2` (new two-camera setup, "5hz 40ms" from 10/27 on) -> angled `Pupil_9P_angled_261016`.
* Per-folder override: `pupil_models.yaml` (nearest at or above a recording wins), keys
  `eye`, `eye_with_face`, `face`; values `frontal`, `angled` or a config path.
* Outputs carry the model in their name, e.g.
  `..._face_nativeDLC_HrnetW48_Pupil_9P_angled261016shuffle1_snapshot_best-90...`
  (DLC drops the leading zero of snapshot-best-090). Outputs of the other model ->
  MODEL_MISMATCH, left untouched unless "Force re-analysis" (archives them first).
* "Show frames" = one frame per camera of the highlighted recordings, labelled with
  the model; "Lik." column = median likelihood (flagged when < 0.6).

### Camera check and first batch (2026-11-03/04, run on samson, CPU ~5-8 fps)
Frames of every camera (`DLC_batch_gui.py --no-gui --root <folder> --frames <dir>`,
saved in `5hz 40ms\_camera_frames`) showed:
* `5hz 40ms\NavRearChr-2609-1` 261027_001 ... 261031_002 (7 recordings): `_eye` = frontal
  close-up (pupil ~25-50 px), `_face` = angled side view (pupil ~9 x 4-5 px) -> the
  default rule is right, no pupil_models.yaml there.
* `visualstim_test\NavRearChr-2609-1` 261014_001, 261020_001, 261020_002 = single-camera
  test runs, angled camera saved as `_eye` -> `pupil_models.yaml` with `eye: angled` in
  that animal folder. (NavWT-2609-1 261020_00x, the short two-camera tests, have
  `_eye` angled and `_face` frontal close-up; not analysed.)
Results (median likelihood, % frames with any point < 0.6), all SUCCESS:

| recording | eye -> frontal | face -> angled |
|---|---|---|
| 261027_001 | 0.91 (8 %) | 0.85 (4 %) |
| 261027_002 | 0.89 (81 %) | 0.93 (0 %) |
| 261029_001 | 0.88 (8 %) | 0.93 (0 %) |
| 261029_002 | 0.81 (56 %) | 0.90 (2 %) |
| 261029_003 | 0.86 (26 %) | 0.89 (2 %) |
| 261031_001 | 0.89 (19 %) | 0.89 (0 %) |
| 261031_002 | 0.85 (51 %) | 0.89 (2 %) |

Test runs with the angled model: 261020_001 0.98 (0 %); 261020_002 0.89 (36 %: clean for
the first ~11 min, then all points drop together in the last ~6 min).
* The frontal model works on the new close-up eye camera: points sit on the pupil edge.
  High "% low" on some eye sessions = lower/right edge points (5-8) under the lower lid
  or low-contrast against the iris; use an ellipse fit to points >= 0.6.
* 261014_001 was NOT redone: archiving its old frontal outputs failed with WinError 5
  (Access is denied). They were most likely written from an Administrator prompt on
  10/16 and the batch ran non-elevated. Redo it from the Administrator prompt with Force
  re-analysis. The angled outputs in the model's `_test_output\` exist already.
  `archive_dlc_outputs` now renames (os.replace) instead of copy+delete, so a refused
  move leaves no duplicate in the archive.

## Analysis note for the angled view
Pupil size: use the horizontal diameter (pupil_1–pupil_5) or an ellipse fitted
to the 8 edge points (major axis), not the mean centre-to-edge radius — the
vertical axis is foreshortened. Low likelihood on pupil_2/3/4 = (partial) blink.
