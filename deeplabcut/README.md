# deeplabcut

Pupil tracking with DeepLabCut 3.0 (PyTorch, CPU only on samson). `DLC_batch_gui.py` converts MJ2 to MP4, then runs DLC analysis, filtering, a labeled sanity video and likelihood QC.
It now handles the two-camera setup: `_eye` videos go to the frontal model `Pupil_9P_WB_260402` and `_face` videos to the angled model `Pupil_9P_angled_261016`, with per-folder `pupil_models.yaml` overrides.

| File | What |
|---|---|
| [`HANDOFF.md`](HANDOFF.md) | **Start here.** Layout, camera → model rule, how to run, QC, results of the 3–4 Nov batch, pending actions |
| [`../sessions/dlc-batch-gui/recovered/`](../sessions/dlc-batch-gui/recovered/) | Docs the session wrote that were rebuilt from its transcript: its own HANDOFF draft, the angled-model notes, the batch results, a `pupil_models.yaml` example |

**Code:** the session's export (44 files, including `DLC_batch_gui.py` identical to samson) is ready on branch `cowork-export/dlc-behavior` in that session's container.
It is waiting for the user to choose how to publish it (private artifact, or a push once the Claude GitHub App is installed on cwei6/GCAMP).
