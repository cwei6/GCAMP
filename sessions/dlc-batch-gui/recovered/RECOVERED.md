# Recovered files: "DLC batch video analysis GUI"

**Caveat:** these files were rebuilt from the Write/Edit events in the transcript of Cowork session `session_01TXsqNvay4sMGNwGK4XyuVe`, not copied from disk. The rebuild was marked complete for each one: every edit applied and no later shell command touched it. They show the file as of its last change below.
**The session's own export (branch `cowork-export/dlc-behavior` of `cwei6/GCAMP`, or the private artifact "dlc-behavior export") is authoritative** wherever it differs.

| File here | Original path (in the Cowork session's container) | Last change (UTC) | What it is |
|---|---|---|---|
| `angled_pupil_model_notes.md` | `/home/claude/project_docs/angled_pupil_model_notes.md`; also the Claude project doc `claude/angled_pupil_model_notes.md` (doc `5219e46a-1556-43d9-a7c3-dda973ef433b`) | 2026-11-04 08:48:20 | Angled-camera pupil model notes: training, validation, point convention, training recipe, two-camera usage, and the 3–4 Nov camera check and results |
| `HANDOFF.md` | `/tmp/claude-0/-home-claude/4b7db539-9055-56fe-aefd-ca35393a5b5d/scratchpad/export/deeplabcut/HANDOFF.md` | 2026-11-09 22:30:32 | The session's own pupil-DLC handoff, written for its export (destined for `deeplabcut/HANDOFF.md`) |
| `batch_2026-11-03_results.md` | `.../scratchpad/export/deeplabcut/docs/batch_2026-11-03_results.md` | 2026-11-09 22:29:06 | Results of the two-camera pupil DLC batch of 3–4 Nov 2026, plus the 261014_001 archive failure |
| `pupil_models.example.yaml` | `.../scratchpad/export/deeplabcut/code/pupil_models.example.yaml` | 2026-11-09 22:28:29 | Documented template for the per-folder camera → model override used by `DLC_batch_gui.py` |

**Not recovered:**
- `/home/claude/DLC_batch_gui.py`: the rebuild is incomplete (edits only, no full write). The current version is on samson (`D:\Behavior\Chao\squint_pipeline\DLC_batch_gui.py`, 115,749 B, md5 `60cfe6672c6427edd68f38546c0607a9`) and in the session's export.
- `scratchpad/launch_gui.sh`: a mock-test helper.

**Checked against the transcript (the files above are left verbatim):**
- `HANDOFF.md` §4 says samson runs "Windows 10". The device info only reports `win32` / `x64`; the Windows version is not recorded.
- `HANDOFF.md` §2 says afferent recording folders hold `suite2p\`. The listing of `NavRearChrTIGRE-2608-1_261101_001` on 3 Nov shows only the `.mj2`, `.sbx`, `.mat` and `_quadrature.mat`.
- `HANDOFF.md` and `batch_2026-11-03_results.md` do not say that `docs\claude_notes\` (named in samson's `CLAUDE.md`) was missing on samson on 9 Nov; see `deeplabcut/HANDOFF.md` §2.
- Their versions, md5s, likelihood results, paths and commands match the transcript.

Other files the session copied into its export without rewriting them (code from samson, the angled-model scripts, the squint package) are not in the transcript's write events. They will arrive only with the export.
