# Behavior project: chats and what each one holds

Written 10 Nov 2026 (updated 11 Nov) in the old Claude project "Bahvior" (the behavior project), for the move to Claude Code.

Claude Code can't open these chats because they live inside a claude.ai project. Everything that matters is summarised below or already saved on disk. The links are only for opening a chat in the Claude app if a detail is missing.

Next to this file are the project's three notes, copied as they were: `squint-pipeline-status.md` (1 Oct), `cpp-recorder.md` (12 Oct) and `facial-von-frey-0.1NTG.md` (11 Nov). The project had no custom instructions, no uploaded files and no scheduled tasks.

## The chats

| # | Chat | When | What it was | Code / output |
|---|---|---|---|---|
| 1 | [Squinting behavior analysis pipeline](https://claude.ai/chat/940c54d7-c9f4-8dc8-b939-f556bca44110) | until 5 Nov 2026 | Main build of the Python squint pipeline | `D:\Behavior\Chao\squint_pipeline` |
| 2 | [Camera control GUI for tablet](https://claude.ai/chat/d4470f26-ae56-4f22-a049-c23986bfde91) | 6–12 Oct 2026 | CPP recording: web recorder + tablet webcam decision | `D:\Behavior\CPP_recorder` |
| 3 | [DLC model training status](https://claude.ai/chat/8491e557-92f6-42bd-aaf0-5a1782cde7fe) | 1 Oct 2026 | Checked the squint DLC retraining, then evaluation | `D:\Behavior\Squint_8P_ispi` |
| 4 | [DeepLabCut Qt bindings error](https://claude.ai/chat/60ea074f-2a72-49dd-a1a8-bbf35c746fc1) | 25 Sep 2026 | Fixed DLC GUI / labeling errors on Windows | none (install fixes below) |
| 5 | [Von Frey hindpaw plot](https://claude.ai/chat/4de76a70-0854-4791-b7c0-c822b1477675) | 4 and 11 Nov 2026 | Hindpaw figure (NTG 10 mg/kg); facial von Frey after 0.1 mg/kg NTG, SUDO-scored | only in the chat; rules in `facial-von-frey-0.1NTG.md` |
| 6 | [White light lux threshold photophobia](https://claude.ai/chat/d7d75c93-c2a0-470e-8960-adb66e97539e) | 29 Sep 2026 | Literature review (no code) | none |
| 7 | [Facial Von Frey method comparison](https://claude.ai/chat/da80fb54-caaf-462f-9cd8-202f746b4d58) | 10 Nov 2026 | SUDO vs Chaplan for facial von Frey (no code) | none |

Related work kept elsewhere:

- The DLC batch of the Treg squint videos ran in the deeplabcut project (its log and `paradigm.yaml` are in `compiled\Treg`). Pupil DLC history: `docs\dlc-pupil-history.md`.
- The CPP chamber and its 3D-printed floor inserts are in the 3D-printer project.

## 1. Squinting behavior analysis pipeline (main code chat, last active 5 Nov 2026)

This chat's messages can't be read from inside the project, so this summary comes from the code, `README.md`, the project notes and memory.

- Built the `squint` package and the tkinter app (`squint_gui.bat` or `python -m squint gui`, tabs 0–8). It is a port of the MATLAB `dlc_eyelid_area_auto_by_day` workflow, reworked as "v2". `README.md` is the manual.
- File dates: core package and paradigm files late Sep to 6 Oct; detectors and review code 10 Oct; GUI, discovery, run, plots and README last changed 2 Nov; `squint_gui.bat` and `squint_dlc_settings.json` on 5 Nov.
- Validation (recorded by 1 Oct): `--mode legacy` reproduces MATLAB on NavRearChr-2511-2-15ms D2–D5. It also found that the MATLAB glob double-counted recordings with two DLC snapshots (28 in Bone-400-2hz).
- Decisions are in `CLAUDE.md`, `squint-pipeline-status.md` (this folder) and the README section "Analysis".
- The DLC step (tab 0) runs `DLC_batch_gui_squint.py` with snapshot "best" of `D:\Behavior\Squint_8P_ispi\config.yaml` (= best-120). On 5 Nov its python was set to `C:\ProgramData\anaconda3\envs\DEEPLABCUT\python.exe`.
- On disk on 10 Nov (checked for NavTIGRE-2609-1 only):
    - 260924 (D-2), 260929 (D3) and 261002 (D6) are tracked with `snapshot_best-120` and analysed in `analysis_v2`. Event review started 13 Oct.
    - 261006 (D10) and 261014 (D18) have raw video only, with no DLC yet. Tab 0 "Scan ticked folders" lists what is pending for every animal.
- Probably a stray: `compiled\Treg\NavTIGRE-2609-1\261002\analysis_v2\` (`261002_D1_*`, 6 Oct) looks like a run on a single date folder. Ask before deleting.
- Next (from the README): a hand-scored validation clip to get precision/recall of the squint and wipe detectors; the squint-duration distribution across animals, to justify a report cut.

## 2. Camera control GUI for tablet (CPP, 6–12 Oct 2026)

Full notes: `cpp-recorder.md` (this folder).

- Built a browser recorder for CPP sessions in `D:\Behavior\CPP_recorder`: live view, timed auto-stop, files named from animal / phase / day / compartment, crash-safe saving and a CSV session log. A one-file version is in `github_pages\`, rebuilt with `python build_single.py`. It was tested only with a simulated camera, never with the real webcam.
- Found that the Samsung tablet doesn't pass USB webcams to Chrome or any browser, and Logitech has no Android app for its webcams.
- Decision (12 Oct): record with the Android app "USB Camera" by ShenYao.
    - Settings: 1280×720 MJPG at 30 fps (YUY2 manages only ~7–10 fps), audio off, auto-delete and auto-segment off.
    - The first test frame was 320×240 and filmed from the front. Advised: camera overhead, ~55–60 cm above the floor, centred over the corridor.
- The web recorder still works on a PC (double-click `index.html`). Offered but not built: a timer + session-log page to run beside the camera app.
- CPP analysis has not started.

## 3. DLC model training status (1 Oct 2026)

- Squint 8-point model in `D:\Behavior\Squint_8P_ispi`: project Squint_8P_LargeSep11, iteration-1, trainset95 shuffle2, ResNet-50. The 200-epoch run finished 30 Sep 02:31 despite a PC shutdown. Snapshots: 050/100/150/200 + best-120.
- Training-log test metrics at epoch 200: RMSE 7.74 px, mAP 97.1. Plateau from ~epoch 100.
- `evaluate_network` (shuffle 2) was run on 1 Oct for best-120: test RMSE 7.73 px (6.64 px above pcutoff 0.6), mAP 97.3, mAR 98.0. The other rows in `CombinedEvaluation-results.csv` come from earlier training runs.
- Treg videos were then tracked with best-120. NTG (`compiled\NTG\1mg`, `10mg`) had been tracked earlier with best-60.

## 4. DeepLabCut Qt bindings error (25 Sep 2026)

Fixes for DeepLabCut 3.0.0rc13 on Windows (conda env `deeplabcut`, working folder `D:\2photon`):

- `qtpy.QtBindingsNotFoundError`: DLC was installed without its GUI extras. Fix: `pip install "deeplabcut[gui]==3.0.0rc13"`.
- `OSError [WinError 1114]` in `c10.dll`: PyTorch loaded after Qt. Fix: import torch before starting the GUI, or pin `torch==2.8.0` with the matching CUDA/CPU wheel.
- napari `TypeError: object.__init__() takes exactly one argument` when labeling: napari 0.4.18 clashes with a newer PySide6. Fix: `pip install "pyside6==6.4.2"` (Python 3.11 or older).
- Harmless warnings: `QPixmap::scaled: Pixmap is a null pixmap` and the napari-console npe2 schema warning. Labels save as `CollectedData_<name>.h5`.

## 5. Von Frey hindpaw plot (4 Nov 2026, continued 11 Nov)

Hindpaw figure, 4 Nov:

- Input was an Excel workbook attached in the chat, with column D = 50 % threshold per mouse. The unlabeled top block was taken as pre and the "post" block as post.
- Groups: NavWT-L + WT-L = Control; NavWT-R + WT-R = NTG 10 mg/kg.
- Thresholds were recomputed from the raw o/x up-down marks and matched column D. One thing to check in the lab notes: WT-R pre and NavWT-L post have identical raw sequences (copy-paste or coincidence).
- Results: NTG 0.43 → 0.21 g (−51 %, both mice down); Control 0.73 → 0.56 g (−24 %, mostly one mouse).
- Figure: paired pre/post per group (mean bars + one line per mouse, marker = genotype), as a 300 dpi PNG and a PDF with editable text (`pdf.fonttype 42`, Arial-metric font). Colours were matched to the slide deck: Control `#111111`, NTG `#e3261b`, open outlined bars, genotype key in grey.
- Not on disk: the script (`plot.py`) and `hindpaw_von_frey_NTG.png` / `.pdf` exist only in that chat. Offered but not made: a percent-of-baseline panel.

Facial von Frey after 0.1 mg/kg NTG, 11 Nov. Scoring rules, findings and output names are in `facial-von-frey-0.1NTG.md` (this folder).

- Input: `Facial_von_Frey.xlsx`, attached in the chat (sheets 260729-0.1NTG-1stbatch, 260805-0.1NTG-2ndbatch and 261110-0.1NTG).
- Decisions: score every series with SUDO (Chaplan series too, from their first 5 pokes), report % change, and leave out mice with a baseline under 0.07 g. With the 0.07 g start, SUDO can't score above 1.18 g.
- Findings: no 0.1 mg/kg NTG effect in the new setup (261110: Vehicle −7 %, DT −4 %, p = 0.89). In the previous setup (260729 + 260805) all 8 mice dropped, mean −45 % (p = 0.008). The two setups also differ in baseline (median 0.14 vs 0.62 g). Setup is confounded with batch and with the vehicle/DT pretreatment.
- Not on disk: the scoring code, the figures (`facial_von_frey_0.1NTG_pct_change` and `facial_von_frey_0.1NTG_baseline_distribution`, PNG + PDF) and `facial_von_frey_0.1NTG_per_mouse.csv` exist only in that chat.

## 6. White light lux threshold photophobia (29 Sep 2026)

Literature review on the light levels behind the squint paradigms (no code).

- Unsensitized wild-type mice needed ~27,000 lux (measured at the floor) to show CGRP-induced light aversion. Sensitized preparations (central CGRP, hRAMP1 overexpression, NTG) bring that down to ~55–1,000 lux.
- Papers measure at the floor, not at the eye. Under overhead light the eye gets roughly half the floor value. A 2025 Kaiser/Aguirre preprint was the only source found with photopic and melanopic lux and a dose–response.
- For Treg: ambient room light (~100–1,000 lux) is below the aversion threshold for unsensitized wild-type mice, so squinting in the ambient blocks is more likely nociceptive than photophobic. A dark-to-light step (~1,000 lux) closes the eyelid even in naive mice, which matters at block boundaries.
- Offered but not made: a methods paragraph or a reference list.

## 7. Facial Von Frey method comparison (10 Nov 2026)

Method review, no code. Its conclusions shaped the scoring in section 5.

- SUDO (Bonin, Bories & De Koninck, Mol Pain 2014) uses five pokes and reads the threshold from the fifth, moved half a step. It matched Chaplan's up-down closely, but was validated only on the paw. On the face it is defensible but not proven; other groups have used it periorbitally.
- Risks on the face: responses are noisier (the filament is visible and can touch whiskers), and thresholds far above the starting filament hit SUDO's ceiling.
- To make it defensible: a pilot that scores the same full Chaplan sequences both ways and checks agreement, not just correlation; responses defined in advance; approach from outside the visual field; a blinded tester; analysis in log units. Mason et al.'s periorbital protocol (Dixon up-down from 0.07 g) is a template.
- One SUDO sequence per mouse per time point is enough. If two are averaged, use the same rest gap every time. Labs that repeat baselines spread them over separate days.
- Chaplan vs SUDO: Chaplan counts six responses from the first change in response (4–9 pokes in total) and scores with Dixon's k table; SUDO always stops at 5 pokes and needs no table.
