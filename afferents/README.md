# afferents

Two-photon GCaMP imaging of trigeminal afferents in awake mice. Events are found with AQuA2, units with CFU, and everything runs through the `afferent2p` pipeline (Python Tk GUI + MATLAB) on samson.
Two cohorts: the RearChr optostim animals (pre vs post optostimulation) and the Treg/DTx animals (chronic, across days).

| File | What |
|---|---|
| [`HANDOFF.md`](HANDOFF.md) | **Start here.** Data layout on samson, pipeline steps, decisions, status per animal, open questions, pitfalls. Checked by the Cowork session itself on 2026-11-10. |
| [`code/`](code/README.md) | The pipeline (`afferent2p` v45 as on samson), the Treg driver and MATLAB drivers, the label files, the analysis scratch scripts and the GUI render tests |
| [`docs/project-notes/`](docs/project-notes/) | The 22 notes of the claude.ai project, verbatim — [`v34-progress.md`](docs/project-notes/v34-progress.md) is the master continuation log (v34 … v45, Oct 19 – Nov 9), [`treg-processing-nov3.md`](docs/project-notes/treg-processing-nov3.md) the Treg state, [`cfu-rules-overlap-threshold-nov2.md`](docs/project-notes/cfu-rules-overlap-threshold-nov2.md) the CFU-rule decision, [`afferent-prepost-analysis-pipeline.md`](docs/project-notes/afferent-prepost-analysis-pipeline.md) the September history (v22 … v29) |
| [`docs/memory/`](docs/memory/) | The project memory files: the user's requests and decisions since Sep 2026, dated ([`gcamp-afferent-analysis.md`](docs/memory/gcamp-afferent-analysis.md)), and the Treg cohort facts ([`treg-project.md`](docs/memory/treg-project.md)) |
| [`docs/artifacts/afferent2p-how-to-run-aqua2-pipeline.md`](docs/artifacts/afferent2p-how-to-run-aqua2-pipeline.md) | How to run the AQuA2 pipeline (user guide, Sep 22; the package's own `code/afferent2p/GUIDE.md` is the current, much longer version) |
| [`docs/artifacts/navrearchrtigre-2608-1-pre-vs-post-optostim.md`](docs/artifacts/navrearchrtigre-2608-1-pre-vs-post-optostim.md) | NavRearChrTIGRE-2608-1 pre vs post optostim report (Oct 3; predates the re-tuning and the 0.3 CFUs) |
| [`docs/artifacts/mouse-trigeminal-fiber-planner.md`](docs/artifacts/mouse-trigeminal-fiber-planner.md) | Trigeminal ganglion optic-fiber planner: coordinates, approaches, light model (Oct 13) |

**Code on samson:** the working copy is `D:\2photon\Chao\GCAMP_affrerents\afferent2p\` (v45 + the user's `pipeline_gui_settings.json`), with versioned copies `afferent2p_vNN\` + `afferent2p_vNN.zip` next to it. `code/afferent2p/` here equals the v45 package; keep them in step when changing the pipeline (see HANDOFF §3 "Delivery convention").

Source session: [`../sessions/afferent-pipeline/`](../sessions/afferent-pipeline/). Imported on 2026-11-11 from the session's own export (private artifact "afferents export", 131 files, checksums verified against its `MANIFEST.sha256`). The earlier transcript-based draft is kept as [`docs/HANDOFF-transcript-draft.md`](docs/HANDOFF-transcript-draft.md).
