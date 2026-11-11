# NavRearChrTIGRE-2609-1 — quick read of the three days (Oct 27, 2026)

Run remotely on samson while the user was away (computer use on the open v37 pipeline window + the open MATLAB R2021b; the pipeline code of v37 = v39 for every step used). No curation ("ignore curation" ticked, so every run is uncurated), no pupil (no DLC file for 261025; the older days still have theirs), bad frames excluded as usual.

## Days (same two FOVs on every day; FOV1, FOV2 pre → 1 h stimulation → FOV1, FOV2 post; 2 runs per FOV and condition)

- 261006 sham — settings of the Oct 25 v37 apply: FOV1 thr 2 / minSize 50 (matched), FOV2 thr 2 / 100 (common, TOO WEAK)
- 261014 "5 Hz max power + 40 ms" — FOV1 thr 2.5 / 50 (common, NOT COMPARABLE), FOV2 thr 1.5 / 50 (matched, level 0.5)
- 261025 "5 Hz max power" (the previous pulse width, same FOVs as 261014) — NO sweep: `aqua_opts\NavRearChrTIGRE-2609-1_261025_FOV1_Opt.csv` / `_FOV2_Opt.csv` are copies of the 261014 files (note in the thrARScl row says so), so 261014 and 261025 are detected at identical settings per FOV. `day_labels.csv` got `261025,5 Hz max power`.

## Numbers (experiments_summary.csv, whole-FOV rate, ≈31 min per condition)

| day | FOV | pre | post | ratio (95 % CI) | notes |
|---|---|---|---|---|---|
| sham | 1 | 19 ev, 0.611/min | 1 ev, 0.032/min | ×0.05 (0.00–0.33) | still-only ×0.05; running 12 % pre |
| sham | 2 | 11 ev, 0.350/min | 1 ev, 0.032/min | ×0.09 (0.00–0.63) | still-only ×0.07; running 30 % pre, 4 % post |
| 5 Hz + 40 ms | 1 | 0 | 1 ev, 0.032/min | – | silent both sides |
| 5 Hz + 40 ms | 2 | 0 | 11 ev, 0.350/min | 0 → 11, direction only | post decays: 8 → 3 events between the two post runs |
| 5 Hz | 1 | 40 ev, 1.281/min | 24 ev, 0.771/min | ×0.60 (0.35–1.02); GLM-adjusted ×0.68 (0.29–1.36) | 8 up / 30 down / 24 unchanged of 62 CFUs |
| 5 Hz | 2 | 4 ev, 0.127/min | 27 ev, 0.859/min | ×6.75 (2.35–26.6); still-only ×12.7 | 18 up / 0 down / 13 unchanged of 31 CFUs; flat over the hour (12 → 15; slope ×0.93 per 10 min, CI 0.77–1.18) |

Reading: FOV2 increases after both protocols and decreases on the sham day; at identical settings the post rate with the previous width is ×2.5 the 40-ms one (0.86 vs 0.35/min) and sustained rather than decaying. FOV1: nothing on the 40-ms day; a mild, non-significant decrease on the 5 Hz day, far smaller than the sham-day drop of the same field (×0.05) — not interpretable as a stimulation effect. Baselines swing between days (FOV1 pre 19 / 0 / 40 events; FOV2 pre 11 / 0 / 4). One animal → descriptive.

Caveats: (1) motion check flags 3 of 4 post runs of 261025 — 005: 9 of 14, 006: 6 of 12, 008: 6 of 15 event onsets on moving frames (2.4–3.5× chance; mouse still 100 %, so image motion/twitches) → look at those events before trusting the FOV2 increase; (2) FOV2 events at thr 1.5 are large and dim (median area 590–660 µm², peak dF/F 0.04) — check by eye (tab 3 Pick by eye / QC report); (3) sham day at different settings → only within-day ratios comparable with it; (4) no pupil; (5) the step-at-stimulation on 261025 (×1.34 mean over the two fields) only just exceeds the within-condition drift (×1.25).

The whole animal was re-analysed (the Oct 21 analysis used the first-sweep settings; the AQuA2 results on samson are from the Oct 25 re-detection), so the sham/261014 numbers changed from the Oct 21 report (e.g. 261014 FOV2 post 4 → 11 events). Outputs: `NavRearChrTIGRE-2609-1\aqua_analysis\` (summary.txt, csvs, pngs), `NavRearChrTIGRE-2609-1_report.html` rebuilt Oct 27; `_qc.html` NOT rebuilt (still Oct 21).

## Practical findings about running the pipeline remotely

- The Batch tab had "overwrite (re-queue tiffs that have results)" ticked from the user's last session — that is why the Oct 25 batch (built 22:50 UTC) had all 16 runs queued although their results existed. Unticked; the new batch queued only the 8 new movies.
- Tab 5 "Run both now" (matlab -batch started from the pipeline window) failed in 9 s: "AQuA2 not found" — the mapped drive H: (H:\MATLAB\AQuA2, and the user's startup.m paths) is not visible to a MATLAB process started by the GUI (`run_aqua_batch_log.txt`: every H:\MATLAB path "nonexistent"). Pasting `run('…\aqua_batch\run_aqua_and_cfu.m')` into the open MATLAB window works (AQuA2 ran 8 movies in ~9 min; CFU batch seconds). Possible fix for the pipeline: pass the AQuA2 folder as a UNC path, or `addpath` it explicitly in run_aqua_batch.m from a path MATLAB -batch can see.
- `run_gui.bat` double-clicked from Explorer fails ("python is not recognized") — python is only on PATH in the Anaconda Prompt. Use the open window or an Anaconda Prompt.
- Computer-use notes: MATLAB's command window accepts key events but not clipboard pastes from the tool — type the command in ≤ 14-character chunks. Timings: unpack 3.2 min/run (8 runs: 25.5 min), AQuA2 ≈ 1 min/movie here, Step 7 with the motion check on all 24 runs 76 min (the motion check is recomputed for every run after a re-detection, ≈3 min each), report 3 s.
