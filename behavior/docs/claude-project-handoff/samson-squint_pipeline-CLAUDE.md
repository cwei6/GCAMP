# Behavior + DeepLabCut (squint pipeline and pupil DLC)

This folder holds two lines of work that used to be two Claude projects ("behavior" and "deeplabcut"):
1. **Squinting behavior:** eye video → DeepLabCut 8-point eyelid tracking → eye-shape analysis → severe squinting events (and wiping events). Being ported from MATLAB to a Python GUI. Plus a CPP setup at an early stage.
2. **Pupil DLC for the 2-photon recordings:** batch DeepLabCut 9-point pupil tracking of the eye/face videos from the GCaMP afferent and GCaMP visual experiments (9th point = centre, so the radius survives partial occlusion). The CSVs feed both imaging pipelines.

## How to work with me
- I'm a neuroscientist, not a coder: I describe goals and ask questions, you write and run the code.
- Explain in plain language. Before building anything non-trivial, give me 2–3 options with trade-offs and your recommendation, then wait for my choice.
- Show results I can look at: traces with synced video frames, labeled videos, figures. Give file paths I can copy.
- Write code changes directly into this folder (no zips or new copies); git keeps the history.
- Never delete, rename or overwrite videos, models or earlier results without asking.
- When we decide something, add it (dated) to the matching history file and update "Where we left off" below.

## Folders
- This folder = `D:\Behavior\Chao\squint_pipeline`:
  - Squint pipeline: package `squint\`, GUI via `squint_gui.bat` or `python -m squint gui` (tabs 0 DLC tracking, 1 Discover, 2 Parameters, 3 Run, 4 Results, 5 Review, 6 Clips & gallery, 7 Aggregate, 8 Legacy check). `README.md` is its full manual (commands, modes, outputs, parameters) — read it before changing anything. Paradigm files in `paradigms\` (`standard_opto`, `treg_ambient_dark`, `ntg_drug`); defaults in `params.yaml`.
  - Squint DLC step: `DLC_batch_gui_squint.py`, run by the app as a subprocess with the deeplabcut env's python (remembered in `squint_dlc_settings.json`; on 5 Nov set to `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe`), using the model's snapshot "best" (= best-120).
  - Pupil DLC: `DLC_batch_gui.py` (GUI batch: converts `.mj2` → `.mp4`, runs DLC, makes a labeled video as a sanity check; settings in `DLC_batch_gui_settings.json`), launchers `pupil_run_5hz40ms.bat` and `pupil_check_frames.bat`. The older `DLC_9point_selected_conditions.py` isn't in this folder (its logs are in the GCAMP data roots).
- Raw squint recordings: dated folders in `D:\Behavior\Squinting` (e.g. `2025_10_1`, `2026_2_14`; also `Filr3`, `logitechCam`, `Sample_videos`).
- Squint data (organised): `D:\Behavior\Squinting\compiled` → condition → animal. `compiled\Treg` (animal / date / camera_index, not renamed) holds the DLC batch log, `paradigm.yaml` (reference for Treg paradigms) and the config file; `compiled\NTG` recordings were noted as never renamed, yet on 10 Nov `NTG\1mg` and `NTG\10mg` hold `D1_<phase>_k_…` videos tracked with the older snapshot best-60, with MATLAB outputs (sidecars in `NTG\<animal>\<date>`) — check with me; `compiled\High_led_test` = LED-intensity tests.
- Squint DLC model (8-point eye): `D:\Behavior\Squint_8P_ispi\` (Squint_8P_LargeSep11, iteration-1 shuffle2, ResNet-50; retrained 28–30 Sep 2026; evaluated 1 Oct: best-120, test RMSE 7.7 px, mAP 97.3; Treg tracked with best-120). Other model folders in `D:\Behavior` (`Squint_2cam`, `Squint_8P_contra`, `Squint_8P_Flir3`, `Squint_8P_webcam`, `Squint_8P_Large-Rodrigo-2025-09-11`, `Pupil`) — check what each is before using it. CPP: `D:\Behavior\CPP_recorder`.
- Pupil models:
  - Frontal eye model: `D:\2photon\Rodrigo\DLC_Models-RN\Pupil_9P_WB_260402` — a colleague's (Rodrigo's) folder: read only, never write there.
  - Angled-camera model: `D:\2photon\Chao\GCAMP_visual_exp\Pupil_9P_angled_261016` (auto-label / review / train / test / QC scripts); use it via `DLC_batch_gui.py` with that folder's `config.yaml`.
- Pupil data roots (the folder picker defaults to the current root but can navigate up):
  - `D:\2photon\Chao\GCAMP_affrerents` → animal → `<animal>_<YYMMDD>_<NNN>` → `<rec>_eye.mj2`
  - `D:\2photon\Chao\GCAMP_visual_exp` → condition (`2hz test`, `2hz_wheellock`, `5hz`, `50ms test`, `old-200umfiber`, `5hz 40ms`) → animal → recording; plus `visualstim_test\`
  - Output CSV next to the video, e.g. `<rec>_eye_nativeDLC_HrnetW48_Pupil_9P_WB_260402shuffle1_snapshot_best-70_filtered.csv`; batch logs in `<condition>\_DLC_batch_logs\`.
- MATLAB squint original: `dlc_eyelid_area_auto_by_day` + a multi-folder batch runner + a cross-animal aggregator exporting a Prism master workbook. It needed recordings renamed `D1_habituation`, `D1_baseline`, `D1_stim`, `D1_post`, … — the Python version must not need renaming.
- Notes: `docs\project-history.md` (squint + CPP), `docs\dlc-pupil-history.md` (pupil DLC), `docs\claude_notes\behavior\` and `docs\claude_notes\deeplabcut\` (detailed notes from the two old projects). `docs\claude_notes\behavior\SESSIONS.md` lists the behavior project's chats: links, what each did, what is and isn't on disk.

## Machines and environment
- DeepLabCut 3.0 (PyTorch) in conda env `deeplabcut` (`C:\Users\cwei6\anaconda3\envs\deeplabcut`; activate by full path if the name doesn't resolve).
- samson (data server, D:) has no GPU and can't open napari, so DLC frame labeling runs on the other PC, which sees samson's D: as `Z:`; there the DLC GUI starts with `python -m deeplabcut` (DLC 3.0.0rc13). Batch runs slow down when other heavy jobs run on samson.
- samson is shared: touch nothing outside my folders.

## Squinting: paradigms (several run in parallel — the pipeline should recognise each recording's paradigm itself)
- Standard (opto), all conditions except NTG and Treg: 30 min habituation, 30 min baseline, 1 h stimulation, 30 min post. Recordings ~30 min (some longer). A rare recording-PC bug drops the frame rate, so a 30 min recording can play shorter.
- Treg: 30 min ambient habituation, 30 min ambient test, 30 min dark habituation, 30 min dark recording. First runs had only one 30 min dark block (really a habituation) — treat it as the dark recording for now. Days: NavTIGRE-2608-1 and NavTIGRETreg-2608-1 started at D3 (no baseline day); all others have a baseline day 1–2 days before D0, then D3.
- High_led_test: light-only pilots (1500 / 2000 lux, 2–3 h) and single-day NTG tests (30 min hab + 30 min baseline, injection, 30 min wait, 4 × 30 min post, then a 5000 lux + dark test or 2 × max power ~6000–8000 lux).

## Squinting: rules already decided
- The pipeline is a GUI exposing every function, including the DLC step (recognises the paradigm and processes videos from raw).
- Treg eye-openness reference uses both the ambient and the dark habituation blocks.
- The 30 s minimum squint duration in the MATLAB pipeline was an arbitrary choice from test runs — open to change, with justification.
- Tracking caution: eye points shift, mistrack and jump across recordings, so metrics on single point pairs may be unreliable.
- Event review tool = desktop Python app: trace + synced video frame, add/remove events, writes a curated events CSV.
- v2 analysis choices are explained in `README.md`: squint candidates kept from 1 s (the old 30 s cut only in the `*_ge30s` columns), thresholds 0.5 / 0.7 or per-animal percentiles, a separate blink channel, two wipe detectors, openness reference = 95th percentile of the same day's reference phases.
- Goal: make squint and wipe characterisation more convincing (normalisation rework, visual atlas to add/remove events, better detection).

## Pupil DLC: rules already decided
- Feed current videos at native resolution as they are: no rescaling, reframing or retiming. (Old videos were 160×112, tagged 30 fps but truly 15.62 fps, upscaled to 640×448; current ones are ~240×180 at 30 fps.)
- Which model: if a recording also has `<rec>_face.mj2`, the face video is the new angled camera → angled model, and `<rec>_eye.mj2` is the previous eye camera → Rodrigo's frontal model. Folders with only the eye video stay as before. Recordings with only the angled camera are test runs. Never use the frontal model on angled videos.
- Two-camera visual-stim paradigm: recordings every other day, with a test run one week before.

## CPP (setup stage)
Records the chamber with a Logitech HD 720p USB webcam connected to a Samsung Galaxy tablet running a webcam app (the tablet's Chrome only sees its built-in cameras). No audio needed. Analysis not started.
- Decision 12 Oct: the Android app "USB Camera" (ShenYao) at 1280×720 MJPG 30 fps, with audio, auto-delete and auto-segment off; camera overhead ~55–60 cm above the floor, centred over the corridor (the first test frame was 320×240 and filmed from the front).
- `D:\Behavior\CPP_recorder` = browser recorder built 6–12 Oct (auto-named files, timed stop, crash-safe, CSV log). It works on a PC but can't see the webcam on this tablet, and was never tested with the real webcam. Details: `docs\claude_notes\behavior\cpp-recorder.md`.

## Where we left off
- Squint: the Python pipeline is built — discovery without renaming, v2 analysis (per-point cleaning, scale-corrected openness, per-day reference), squint / blink / wipe detectors, the review tool, clips + gallery, aggregate to a Prism workbook, and a legacy mode that reproduces the MATLAB numbers (checked on NavRearChr-2511-2-15ms D2–D5). Treg videos went through the DLC batch (log + `paradigm.yaml` in `compiled\Treg`); Treg days come out as baseline (−2), D3, D6, D10, D17; day 260915 recorded at only 10–13 fps.
- Squint next (from the README): a hand-scored validation clip → precision/recall of the squint and wipe detectors; the squint-duration distribution across animals to justify a report cut. Check `docs\claude_notes\behavior\` for anything newer.
- Squint on disk (checked 10 Nov): Treg NavTIGRE-2609-1 days 261006 (D10) and 261014 (D18) are raw video only, with no DLC yet (other animals not checked); its D-2/D3/D6 are analysed and event review started 13 Oct. `compiled\Treg\NavTIGRE-2609-1\261002\analysis_v2\` (`261002_D1_*`) looks like a stray run on one date folder — ask before removing.
- Other behavior work (see `SESSIONS.md`): von Frey analyses, hindpaw NTG 10 mg/kg (4 Nov) and facial 0.1 mg/kg NTG scored with SUDO (11 Nov; rules and findings in `docs\claude_notes\behavior\facial-von-frey-0.1NTG.md`). Their scripts and outputs exist only in the old chat. Also a lux/photophobia literature summary (29 Sep) and a SUDO vs Chaplan method review (10 Nov).
- Pupil (4 Nov): the two-camera batch on `GCAMP_visual_exp\5hz 40ms` finished 14/14 (eye medians ~0.81–0.91, face ~0.85–0.93; low lower-edge points on some eye videos are a lid/contrast effect, not a failure). Then test runs `visualstim_test\NavRearChr-2609-1` 261020_001/_002 and 261014_001 (re-analysed with `--force-reanalyze`, old outputs in `_DLC_archive_*`) — results in `docs\claude_notes\deeplabcut\angled_pupil_model_notes.md`.
- Pending: pupil for new recordings (e.g. afferent NavRearChrTIGRE-2609-1 day 261025).
