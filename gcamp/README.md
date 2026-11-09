# gcamp

GCaMP imaging of visual cortex (V1 L2/3, jGCaMP8s) across optostimulation cohorts: the `GCAMP_pipeline_reworked` MATLAB pipeline (Suite2p, ROICaT tracking, steps 1–11) and the visual experiments in `D:\2photon\Chao\GCAMP_visual_exp` on samson.

| File | What |
|---|---|
| [`docs/artifacts/gcamp-pipeline-audit.md`](docs/artifacts/gcamp-pipeline-audit.md) | Code audit of `GCAMP_pipeline_reworked_20260829` (Aug 30). Blocking defect: PRE/POST ROIs are paired by row index instead of ROICaT track ID. Also covers missing RNG seeds, no FDR correction, and pupil-metric problems, with a prioritised fix list |
| [`docs/artifacts/v1-habituation-deficit.md`](docs/artifacts/v1-habituation-deficit.md) | Findings memo for run 260901_2114 (Sep 2): habituation deficit after 2 Hz optostimulation, animal-level audit, next analyses |

**Not transferred yet:** the Cowork session(s) where this work happened. Their links are needed (see `../sessions/README.md`). No HANDOFF exists until then.
Related: the pupil DLC outputs this pipeline reads come from [`../deeplabcut/`](../deeplabcut/).
