# Project history — squinting pipeline + CPP (was the "behavior" Claude project)

Copied from the Claude project notes on 5 Nov 2026, when the project moved to Claude Code. Detailed notes from earlier sessions are in `docs\claude_notes\behavior\`.

## Squinting pipeline
- Workflow: record mouse eye video → DeepLabCut tracks 8 points around the eye → analyse the shape; the main output is severe squinting events.
- Analysis was in MATLAB (`dlc_eyelid_area_auto_by_day` + a multi-folder batch runner + a cross-animal aggregator that exports a Prism master workbook); being reworked as a Python pipeline.
- Data root `D:\Behavior\Squinting\compiled`, one subfolder per condition, then per-animal folders.
- The MATLAB pipeline requires every recording renamed `D1_habituation`, `D1_baseline`, `D1_stim`, `D1_post` etc.; the Python file selection should be smart enough to drop that requirement. NTG folder recordings were never renamed.
- New unprocessed data under `compiled\Treg` (animal/date/camera_index layout, not renamed). The DLC batch of the Treg videos ran in a separate session; the batch log and `paradigm.yaml` it wrote in `compiled\Treg` are the reference for Treg paradigms.
- Treg paradigm: 30 min ambient habituation, 30 min ambient test, 30 min dark habituation, 30 min dark recording; the first runs had only one 30 min dark block (really a habituation) — treated as the dark recording for now.
- Treg experiment days: NavTIGRE-2608-1 and NavTIGRETreg-2608-1 were the first run, with no baseline day, starting at D3; all other animals have a baseline day (1–2 days before D0) and their second time point is D3.
- Standard (opto) paradigm for all conditions except NTG and Treg: 30 min habituation, 30 min baseline, 1 h stimulation, 30 min post; each recording ~30 min, some a bit longer. A rare recording-PC performance bug drops the acquisition frame rate, so a 30 min recording plays shorter than 30 min.
- `compiled\High_led_test`: tests whether the set LED intensity is enough to trigger squinting — light-only pilots (1500 / 2000 lux, 2–3 h) and single-day NTG injection tests (30 min hab + 30 min baseline, injection, 30 min wait, 4 × 30 min post, then a 5000 lux + dark test or 2 × max-power ~6000–8000 lux test).
- Several squint paradigms run in parallel (Treg and others) and may be processed together, so the pipeline/GUI should recognise each recording's paradigm automatically.
- Goals: make squinting and wiping characterisation more convincing — normalisation rework, a visual atlas to add/remove detected events, and better detection.
- Concern: the eye points shift, mistrack and jump across recordings, so metrics built on individual point pairs may not be reliable.
- The 30 s minimum squint duration in the MATLAB pipeline was an arbitrary choice from observing test runs.

## Decisions
- Event review tool = desktop Python app (trace + synced video frame, add/remove events, writes a curated events CSV).
- Build order: file discovery + feature/normalisation layer first, validated against an already-analysed animal; then detection and the review tool.
- The Python pipeline is a GUI exposing all functions; it lives at `D:\Behavior\Chao\squint_pipeline` (tkinter app, `python -m squint gui`).
- Treg: the eye-openness reference uses both the ambient and the dark habituation blocks.
- The behavior DeepLabCut model (8-point eye tracking) lives in `D:\Behavior\Squint_8P_ispi\`; the separate pupil model is only for the 2-photon recordings.
- The squint GUI also includes the DeepLabCut processing step (`DLC_batch_gui_squint.py` in the pipeline folder), recognising the paradigm and processing videos from raw.

## CPP (conditioned place preference)
- Setting up the CPP paradigm. For now the chamber is recorded with a Logitech HD 720p USB webcam; the only device available to connect it to is a Samsung Galaxy tablet.
- The tablet's Chrome lists only its built-in cameras, not the Logitech, so recording uses the Logitech + a webcam app on the tablet. CPP recordings don't need audio.
- The chamber itself (Panlab-style place preference box) and its 3D-printed floor inserts are covered in the 3D-printer project.

## Added 10 Nov 2026 (from the project's chats and a check of the disk)
Chat-by-chat summaries with links: `docs\claude_notes\behavior\SESSIONS.md`.

### Squint DLC model
- `D:\Behavior\Squint_8P_ispi`: project Squint_8P_LargeSep11, iteration-1, trainset95 shuffle2, ResNet-50 (GN), 448×448 crops, 200 epochs. Run 28–30 Sep 2026 (two earlier attempts were interrupted); snapshots 050/100/150/200 + best-120.
- 1 Oct: `evaluate_network` on shuffle 2, best-120: test RMSE 7.73 px (6.64 px above pcutoff 0.6), mAP 97.3, mAR 98.0.
- Treg videos tracked with best-120 (early Oct); NTG (`compiled\NTG\1mg`, `10mg`) with the older best-60 (Jul 2026).

### Squint pipeline timeline (file dates + notes)
- By 1 Oct: v2 analysis and legacy mode. Legacy reproduces MATLAB on NavRearChr-2511-2-15ms D2–D5; MATLAB double-counted recordings that had two DLC snapshots (28 in Bone-400-2hz).
- 6 Oct: package, paradigm files and the squint DLC batch GUI on disk. 10 Oct: detectors and review code. 13 Oct: Treg NavTIGRE-2609-1 D-2/D3/D6 analysed and event review started.
- 2 Nov: GUI, discovery, run, plots and README updated. 5 Nov: `squint_gui.bat`, and the DLC step's python set to `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe` (`squint_dlc_settings.json`).
- 10 Nov: NavTIGRE-2609-1 days 261006 (D10) and 261014 (D18) not yet tracked. Possible stray output: `compiled\Treg\NavTIGRE-2609-1\261002\analysis_v2\` (`261002_D1_*`, 6 Oct).

### NTG naming (to confirm)
- The notes say NTG recordings were never renamed, but `compiled\NTG\1mg` and `10mg` hold `D1_habituation_…`, `D1_baseline_k_…` and `D1_post_k_…` videos with best-60 DLC outputs, MATLAB results and `Prism_Master_acuteD1.xlsx`. Their sidecars sit in `compiled\NTG\<animal>\<date>`.

### CPP recording (6–12 Oct 2026)
- Web recorder built in `D:\Behavior\CPP_recorder` (README there): timed auto-stop, files named `CPP_<date>_<time>_<animal>_<phase>[_D<day>][_Dots|_Stripes].mp4`, crash-safe chunks and a CSV session log; one-file copy in `github_pages\` via `build_single.py`. Tested only with a simulated camera.
- The Samsung tablet doesn't expose USB webcams to browsers, and Logitech has no Android webcam app. Decision 12 Oct: "USB Camera" by ShenYao at 1280×720 MJPG 30 fps, with audio, auto-delete and auto-segment off; camera overhead ~55–60 cm. Possible next step: a timer + session-log page for the tablet.

### Other behavior work
- 29 Sep: literature summary on the light levels behind migraine-like photophobia in mice (no code).
- 4 Nov: von Frey hindpaw figure, Control vs NTG 10 mg/kg (the script exists only in the chat).
- 10 Nov: method review of SUDO vs Chaplan for facial von Frey (no code).
- 11 Nov: facial von Frey after 0.1 mg/kg NTG, all three batches scored with SUDO as % change, baselines under 0.07 g left out. No effect in the new setup; −45 % in the previous one. Details in `docs\claude_notes\behavior\facial-von-frey-0.1NTG.md` (scripts and outputs only in the chat).
