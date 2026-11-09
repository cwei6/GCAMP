# afferents

Two-photon GCaMP imaging of trigeminal afferents in awake mice. Events are found with AQuA2, units with CFU, and everything runs through the `afferent2p` pipeline (Python Tk GUI + MATLAB) on samson.
Two cohorts: the RearChr optostim animals (pre vs post optostimulation) and the Treg/DTx animals (chronic, across days).

| File | What |
|---|---|
| [`HANDOFF.md`](HANDOFF.md) | **Start here.** Data layout on samson, pipeline steps, decisions, status per animal, open questions, pitfalls |
| [`docs/artifacts/afferent2p-how-to-run-aqua2-pipeline.md`](docs/artifacts/afferent2p-how-to-run-aqua2-pipeline.md) | How to run the AQuA2 pipeline (user guide, Sep 22) |
| [`docs/artifacts/navrearchrtigre-2608-1-pre-vs-post-optostim.md`](docs/artifacts/navrearchrtigre-2608-1-pre-vs-post-optostim.md) | NavRearChrTIGRE-2608-1 pre vs post optostim report (Oct 3) |
| [`docs/artifacts/mouse-trigeminal-fiber-planner.md`](docs/artifacts/mouse-trigeminal-fiber-planner.md) | Trigeminal ganglion optic-fiber planner: coordinates, approaches, light model (Oct 13) |

**Code:** not in the repo yet. The current `afferent2p` (v45) lives on samson in `D:\2photon\Chao\GCAMP_affrerents\afferent2p\`.
The Cowork session "Afferent activity analysis pipeline" was asked to export it to branch `cowork-export/afferents`, and is waiting for the user to confirm.

Source session: [`../sessions/afferent-pipeline/`](../sessions/afferent-pipeline/).
