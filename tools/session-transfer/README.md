# session-transfer tools

Scripts used to move Cowork session transcripts into this repo.

- `condense.py OUT_DIR PAGE_FILE...`: parses saved `list_events` result pages (JSON) into `OUT_DIR/condensed.md`, a readable chronological transcript.
  It also writes `human_messages.md` (the user's messages, verbatim), `compaction_summary.md` (the session's own summary at its last compaction),
  `file_writes.jsonl` (every Write/Edit with full content) and `stats.json` (includes `last_uuid`, the sync cursor).
- `split.py OUT_DIR [CHARS]`: splits `condensed.md` into chunks for reading.
- `replay.py OUT_DIR`: rebuilds files from Write and Edit events. It marks a file complete only when every edit applied cleanly and no later shell command touched it.
