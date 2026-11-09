# Working in this repo

This repo is the Claude Code home for lab analysis work that started in Cowork ("normal") Claude sessions.
It holds code, project notes and handoff documents. **It holds no data.**

## Start here
1. Read the `HANDOFF.md` of the project you are working on (`afferents/`, `gcamp/`, `deeplabcut/`, `behavior/`, `vasculature/`).
   It tells you the goal, where the data and code live on samson, how the user runs things, the current status and what is pending.
2. Then read that project's `docs/` (project notes and converted reports) and `code/`.
3. `sessions/<name>/` has the source Cowork session records: the session's own summary, the user's recent messages verbatim, and a digest.
   Use them when a HANDOFF is unclear.

## Where the data is
- All data is on the Windows lab PC **samson** (e.g. `D:\2photon\Chao\GCAMP_affrerents`, `D:\2photon\Chao\GCAMP_visual_exp`, `D:\Behavior`).
- A Claude Code cloud session cannot reach samson by default.
  If the user has shared a samson folder through Remote Control (Claude desktop app on samson → Settings → Claude Code → Remote Control),
  use `list_computer_folders` / `request_computer_folder` to attach it. Otherwise write or fix code here and have the user run it on samson.

## Rules
- The GitHub repo is **public**. Never commit data files, results, images, credentials, or personal information. `.gitignore` blocks the common data types.
- Keep each project's code under `<project>/code/` and notes under `<project>/docs/`.
- When you change a pipeline, update that project's `HANDOFF.md` (status, versions, decisions) in the same commit.

## Re-syncing a Cowork session that is still running
Each `sessions/<name>/README.md` records a sync cursor (the last transcript event already transferred).
To pull what happened since, page the session's events forward from that cursor with the `list_events` tool,
save each page to a folder, and run `tools/session-transfer/condense.py OUT pages/*` (then `split.py`, `replay.py`)
to get a readable transcript. Then update `sessions/<name>/digest.md`, the project's `HANDOFF.md`, and the cursor.
