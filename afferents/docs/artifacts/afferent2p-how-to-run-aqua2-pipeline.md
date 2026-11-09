# afferent2p — how to run the AQuA2 pipeline

Source: https://claude.ai/artifact/HdCUxqgcFyeAn8uJSVBvUk (claude.ai artifact, last updated 2026-09-22)

Step-by-step guide to running the afferent2p pipeline on GCaMP afferent two-photon recordings: unpack and brightness-QC the `.sbx` files, tune AQuA2 per FOV, batch AQuA2 and CFU detection in MATLAB, curate CFUs by eye, and compare conditions with locomotion and pupil taken into account.

## Before you start (once)

Everything runs from an Anaconda prompt in the pipeline folder, plus two MATLAB scripts. Open the prompt and go there first, every session:

```
cd D:\2photon\Chao\GCAMP_affrerents\afferent2p
```

**One window for every step:** `python pipeline_gui.py` (from that folder) opens a window with one tab per step below — unpack, clean & check frames, AQuA2 GUI notes, batch, the MATLAB commands to copy, CFU curation, analysis — plus a table of every recording showing what the brightness step removed and the sanity-check verdict (double-click a row for its removed-frames sheet) and a log. It calls the same functions as the commands in this guide, so the two can be mixed.

1. Unzip the newest `afferent2p_vNN.zip` (currently v22, in `D:\2photon\Chao\GCAMP_affrerents\`) over that folder. Nothing in your data folders is touched by the unzip.
2. Install the Python packages once: `pip install tifffile scipy matplotlib pandas h5py sbxreader` (`h5py` reads AQuA2's `.mat` files).
3. Update AQuA2 from GitHub (yu-lab-vt/AQuA2). Your copy is older than 2025-12 and has the `curSz` save bug; the current version has it fixed. If you cannot update, add `curSz = []; dsSclXY = [];` above the line `if opts.sz(3)>1` in `src\+ui\+proj\saveExp.m`.

**Folder and name conventions the tools rely on**

| What | Convention |
| --- | --- |
| Animal folder | `D:\2photon\Chao\GCAMP_affrerents\<mouse>\` |
| Run folder | `<mouse>_<yymmdd>_<nnn>\` with the `.sbx`, `.mat`, `_quadrature.mat` and the DLC eye csv |
| Paradigm | no `RearChr` in the mouse name: FOV k = run 00k on every day, conditions = days. `RearChr` in the name: runs 1–4 = pre (FOV1, FOV2, FOV1, FOV2), runs 5–8 = post, one experiment per day and FOV |
| AQuA2 settings per FOV | `FOV1_Opt.csv`, `FOV2_Opt.csv`, … in the animal folder (any name containing `FOV1` works) |
| Binned tiff | `<run>_bin31f_1p98s.tif` next to the `.sbx` (older ones are `..._1.98s.tif`; both are accepted) |
| AQuA2 results | next to the tiff: `<run>_..._AQuA2.mat`, `_AQuA2_Ch1.csv`, `_AQuA2_res_cfu.mat` |
| Your decisions | `<run>_..._AQuA2_res_cfu_curated.json` next to the results; never inside AQuA2's files |
| Summaries | `<mouse>\aqua_batch\` (batch files) and `<mouse>\aqua_analysis\` (results) |

Folders the tools create (`aqua_batch`, `aqua_analysis`, `AQuA_input`, `afferent2p_out`) are never scanned for recordings, and a run that appears twice (a copy in `test_run`, say) is used once and reported.

## Step 1 — Unpack the recordings

One command per animal turns every `.sbx` into a 2 s/frame tiff, corrected for whole-field brightness changes, with the locomotion and pupil signals binned the same way:

```
python unpack_bin.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --sec 2
```

Runs that already have a binned tiff are skipped, so re-running after a new imaging day only unpacks the new folders. Add `--overwrite` to redo everything. `python unpack_bin.py` with no arguments opens a folder-picker GUI with the same options.

What appears next to each `.sbx`:

| File | Content |
| --- | --- |
| `<run>_bin31f_1p98s.tif` | the movie, mean of 31 frames per output frame (1.984 s), uint16, opens in ImageJ; whole-field brightness flattened; every frame kept, bad frames flagged |
| `<run>_..._ch2.tif` | the red tracer channel, same binning (only when recorded) |
| `<run>_..._meanproj.tif` | mean image, for a quick look |
| `<run>_..._json` | exact seconds per frame, frame count, size, what was corrected |
| `<run>_..._behavior.csv` | running speed, still fraction and pupil diameter per recorded frame — running from `<run>_quadrature.mat`, pupil from the DeepLabCut eye csv `<run>_eye_*DLC_*filtered.csv` (9 points, `pupil_9` = centre; diameter = 2 × median centre-to-edge distance over the points with likelihood ≥ 0.6, so a partly hidden outline still gives a radius); both files must sit in the run folder. Eye csv arrived later? `python unpack_bin.py "<mouse>" --behavior-only` (GUI tab 1: *Rebuild behaviour csv*) writes only this csv again, and Step 6 does the same by itself |
| `<run>_..._frames.csv` and `_globalfix.png` | the brightness QC: one row per frame with its level, class (artefact / state), stretch, decision and excluded flag; the png shows the level and histogram-shape traces (red = artefact, amber = state, black bar = excluded) with running speed underneath; \_flagged_frames.png/.csv is the check sheet described below |

**Nothing is removed from the movie.** The step flags two kinds of stretches and then checks every frame of them against the kept frames around it. *Artefact* (red): fast changes — a flash, a jump between two frames beyond 3 %, a spike shorter than 10 s, a content change. Only its **flagged frames** are excluded from the statistics by default: the frames that differ from the neighbours by ≥ 10 % of the frame and ≥ 3× the clean floor (the flash itself; `--core-min-pct` / `--core-ratio`), ±1 frame, their time leaving the denominator; the tails and margins around them and stretches without such a frame are kept. *State* (amber): slow changes — a level plateau, or the background/fluorescence ratio (p10/median) more than 5 % from its baseline; a laser or focus change, or widespread activity, the trace cannot tell, so they are kept (a note warns when the background ends up brighter than normal after the correction). In `260911_001` that means the flash frames 67–68, 76–77 and 99 (+25 to +46 % whole field) are excluded and the nine stretches' other 72 frames — the ±5–9 % neighbours, the mild stretches, the dim state 55–62 — are kept. `--no-global-fix` turns the step off, `--keep-raw` also keeps the uncorrected movie (`..._raw.tif`, twice the disk; lets you undo a drop later).

**Check frame by frame in ImageJ, then decide.** Per recording the step writes `_check_frames.tif`: a 16-bit ImageJ stack of every flagged stretch with 6 frames of context on each side, full resolution, each slice labelled (ImageJ title bar and a line burned into the top of the frame) with its recorded frame number, stretch and class, whole-field change, difference from the clean neighbours (`x` = relative to how much two clean frames differ), FLAGGED FRAME or not, EXCLUDED or kept. `_check_movie.tif` has the same slices as a composite: the frame (fixed grey scale, a flash shows as brighter), its % difference from the neighbours (red brighter, blue darker) and the level/shape traces with a cursor on the current frame, flagged frames ticked in red, excluded frames barred. `_flagged_frames.png/.csv` is the one-glance overview. In the GUI (tab 1b) set the ImageJ/Fiji executable once in the top row, select the recording, *Open frame stack in ImageJ*, scroll; the lower table lists each stretch with its flagged frames, verdict, decision and excluded frames. Decide per stretch: *Exclude flagged frames* (the default), *Exclude whole stretch*, *Keep all*, *Back to default*, or type frame numbers as in the labels (e.g. `68,76-77`) and *Exclude these frames*; or `aqua_clean.py "<run folder>" --decisions-only --keep 13-19 --exclude 63-87 --exclude-all 38-42 --exclude-frames "94-102:99"`. Decisions live in `<run>_..._frames_decisions.json`, take effect at Step 6 without re-running AQuA2, and survive re-runs of the step (*Make / refresh check stacks* rewrites the labels with the current decisions). Removing frames from the movie is an explicit choice (`--drop --force`, or the GUI button): then `_frames.csv` column `tiff_frame` maps movie frames to recorded frames and Step 6 converts event frames back.

**Tiffs unpacked before v22** (named `..._1.98s.tif`) get the same treatment without re-reading the 24 GB `.sbx`, about 15 s per recording; **tiffs treated by an older version** (v22 fillers, v23/v24 drops) get the current rules with `--force`. Frames those versions removed stay excluded — their content is gone — unless `--restore-raw` finds a `..._raw.tif` (then the run is restored and flagged again with every frame kept); otherwise re-unpack that run with `--overwrite`:

```
python aqua_clean.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --force
python aqua_clean.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --restore-raw
```

Without `--force` it never applies twice. Decisions never require an AQuA2 re-run; only a movie that was changed afterwards (`--drop`, `--fill`, `--restore-raw`) makes its AQuA2 result stale → Step 3 with `--overwrite`.

## Step 2 — Tune each FOV once in the AQuA2 GUI

The GUI is used only to decide the settings, on the first recording of each FOV; the batch does the rest. Tune on the corrected tiff from Step 1, so you are looking at what the batch will see.

1. `aqua_gui` in MATLAB, load `<run>_bin31f_1p98s.tif` of the FOV's first recording (`..._001` for FOV1, `..._002` for FOV2, …).
2. Project setup: temporal resolution **1.984** s per frame (the value is in the run's `.json`; 2.0 is 0.8 % off), **0.792** µm per pixel, border exclusion **5**. Leave `detectGlo` off: it did not catch the whole-field step, its one global event was a 310 s drift.
3. Registration on if you want AQuA2 to correct motion (the batch keeps whatever you choose); bleach correction can stay off, the level is already flat.
4. Work through Act → Temporal → Spatial until the events match what you see in the movie. The per-parameter reasoning is in the project note *AQuA2 GUI setup notes*.
5. Detect tab → **SaveOpts**. Save it **in the animal folder** with the FOV in the name:

```
D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\FOV1_Opt.csv
D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\FOV2_Opt.csv
```

The GUI's suggested name with `FOV1` added (`NavRearChrTIGRE-2608-1_260911_001_bin31f_1p98sFOV1_Opt.csv`) works as well. That file is then applied to FOV1 in every run and on every day. Frame time, pixel size and border are taken from the tiff and the batch command, not from this file, so a file saved before you set them is still safe.

You do not need to save the detection result from the GUI (the button that crashed with `curSz`); the batch writes all results itself. A GUI result folder left from earlier (`...\<run>_bin31f_1.98s\`) is ignored and can be deleted.

To give one day its own settings: `<mouse>_<yymmdd>_FOV1_Opt.csv` in the animal folder wins over `FOV1_Opt.csv` for that day.

## Step 3 — Build the batch

One command per animal pairs every tiff with its FOV's settings and writes the files MATLAB will run. Replace `D:\AQuA2` with the folder that contains `aqua_cmd_batch.m`:

```
python aqua_batch.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --aqua-dir "D:\AQuA2" --um-per-px 0.792 --border 5
```

Read the printout before moving on. It lists each experiment (FOV, and for RearChr mice the day) with its runs labelled pre/post or by day, which settings file it found, and the tiffs it queued. A FOV without a settings file is named and left out; a settings file saved on day 1 being applied to day 2 is noted (that is the intended use). Tiffs that already have an `_AQuA2.mat` are skipped unless you pass `--overwrite`.

It writes `<mouse>\aqua_batch\`:

| File | Purpose |
| --- | --- |
| `batch.csv` | one line per tiff and the preset (settings) number it gets |
| `parameters_for_batch.csv` | the presets, one column per FOV, in AQuA2's own format |
| `run_aqua_batch.m` | MATLAB: runs AQuA2 on `batch.csv` |
| `run_cfu_batch.m` | MATLAB: CFU detection on every result |
| `manifest.json` | which tiff belongs to which FOV, day and condition; read by Steps 5 and 6 |

Useful options: `--fov 2` limits the batch to one FOV (with `--overwrite`, redoes just that FOV after you changed its settings); `--sec-tol` if your bin length differs from 2 s by more than 0.25 s; `--stim-after-run N` if a RearChr day was not the standard 4 pre + 4 post; `python aqua_batch.py` alone opens the GUI version.

## Step 4 — Run AQuA2 and the CFU detection in MATLAB

In the MATLAB command window, two lines per animal, in this order:

```
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_aqua_batch.m')
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_cfu_batch.m')
```

The first script changes into the AQuA2 folder and calls `aqua_cmd_batch` with the two csv files. Each tiff takes about as long as a GUI run of the same movie; nothing needs clicking. Results are written **next to each tiff**, with fixed names:

| File | Content |
| --- | --- |
| `<run>_..._AQuA2.mat` | the full result (about 1 GB per 15-min run: it contains the movie); keep it, the CFU step reads it |
| `<run>_..._AQuA2_Ch1.csv` | one column per event: start frame, area, peak dF/F, duration, rise, decay, p value … (this is what Step 6 reads) |
| `<run>_..._AQuA2_Ch1_curves.xlsx` | the dF/F curve of every event |
| `<run>_..._AQuA2_risingMaps_*` | propagation onset maps (not used downstream) |

The second script takes the `_AQuA2.mat` of the batch's movies (not GUI test results elsewhere in the folder) and writes `<run>_..._AQuA2_res_cfu.mat` beside it: the CFUs (spatial units with at least 3 events), their maps, curves and event lists. It skips results whose CFU file is already newer than the AQuA2 file, so after re-running one FOV only that FOV's CFUs are redone. CFU settings (minimum events per CFU, overlap) are set when the script is generated: `aqua_batch.py --cfu-min-events 3 --cfu-overlap 0.5`.

Since v27 the same script also writes, beside each result, `<run>_..._AQuA2_events.csv` (one row per event: start and end frame, centre, footprint size and box, peak dF/F, area — the feature table has no location), `<run>_..._AQuA2_eventmap.tif` (how many events touched each pixel) and `<run>_..._AQuA2_meanimg.tif` (mean of the registered movie). Step 6's maps need them; results made with an older script get them by deleting their `_res_cfu.mat` and running the CFU script again.

Both scripts can be re-run at any time; nothing is lost, and your curation from Step 5 lives in separate files.

## Step 5 — Check the CFUs by eye

```
python cfu_curate.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1"
```

The window shows one run at a time: the mean image with every CFU outlined (red = real, grey = noise, cyan = added by you, yellow = selected), the selected CFU's dF/F curve with its event onsets, its numbers, and underneath a raster of all CFUs' events over time (an all-noise CFU is obvious there).

| Do | Key or button |
| --- | --- |
| Select a CFU | click its outline, or click again to cycle through overlapping ones |
| Mark it real / noise | `k` / `x` |
| Add a region AQuA2 missed | `a`, click the vertices, Enter to close (Esc cancels) |
| Next / previous run | `n` / `p` (saves first) |
| Save | `s` |
| Copy this run's decisions to the same FOV's other runs | "Propagate to this FOV" |

Decisions go into `<run>_..._AQuA2_res_cfu_curated.json`, so re-running Step 4 never loses them. Propagation matches CFUs across runs by overlap after aligning the mean images; a CFU marked noise in one run counts as noise in all of them. Hand-added regions are measured from the tiff with a simple threshold, not by AQuA2, and are labelled `manual` everywhere downstream.

Without the GUI: `python cfu_curate.py <mouse folder> --list` shows what exists per run; `python cfu_curate.py <run tiff> --reject 3,7` or `--keep-all` sets decisions from the prompt; `--propagate` copies them within each FOV.

## Step 6 — Compare conditions

```
python aqua_analyze.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --um-per-px 0.792
```

It reads every run's event table, the QC flags from Step 1 and the decisions from Step 5, and writes `<mouse>\aqua_analysis\`. Start with `summary.txt`, then the figures.

Events starting in an excluded stretch (±1 frame; artefacts and the state stretches you decided to exclude) are dropped and that time is not in `minutes_used`; `run_summary.csv` also says how many events start in state stretches you kept (`n_events_in_kept_state`) — the number to look at when deciding. Frame numbers in the event tables are recorded frames (`t0_frame`, aligned with `_behavior.csv`); when frames were dropped from a movie, `t0_frame_movie` is the frame AQuA2 saw.

| File | What it holds |
| --- | --- |
| `summary.txt` | per FOV and condition: events per minute, peak dF/F, area, duration; the change vs baseline; CFU counts and the paired CFU comparison; the pupil coupling per condition and the pre→post rate change at matched pupil size |
| `run_summary.csv` | one row per recording: minutes used, events, events/min (also per mm²), medians of every property, still minutes, pupil, plus how many frames were flagged and how many events the QC and your curation removed, and the `pupil_*` columns below. **This is the table to build statistics from.** |
| `condition_summary.csv` | one row per FOV × condition (pre/post or day), pooled and as median of runs |
| `events_<FOV>.csv` | every event with run, condition, time and the running/pupil state at its onset |
| `cfu_matched_<FOV>.csv` | CFUs matched across the runs of a FOV, with each CFU's events/min per run and per condition and the pre→post ratio; hand-added regions as `manual1`, `manual2` |
| `pupil_<FOV>.csv` | one row per run: pupil median and 10–90 % range, pupil–running r, pupil–activity r at its best lag with the shuffle p, partial r with running removed, still-only r, CFUs coupled up/down, events/min inside each pupil tercile |
| `pupil_cfu_<FOV>.csv` | one row per CFU and run: its dF/F–pupil r and `up` / `down` / `no` against the shuffle band |
| `<FOV>_overview.png` | event rate per run and the distributions of amplitude, area, duration, rise and decay per condition |
| `<FOV>_cfu_map.png` | the CFU maps of every run, one colour per matched identity |
| `<FOV>_pupil.png` | per run the pupil over the activity trace (grey = excluded frames) with the lag curve; then event rate per pupil tercile per condition, per-CFU r, pupil vs running, the pupil distributions with the tercile edges, activity vs pupil octiles |
| experiments_summary.csv | one row per experiment (FOV on one day, or FOV across days) and condition vs its baseline: raw rate ratio, still-only and running-only ratios, matched-pupil ratio, behaviour-adjusted ratio with 95 % CI, how much behaviour vs the condition explains, common-region ratio, map correlation, % of active tiles that moved, fibre segments up/down |
| \<mouse>\_effects.png | those ratios drawn, one row per experiment — the first figure to look at |
| \<FOV>\_adjusted.png | rate vs running per 20-s bin and condition, rate around running onsets, raw vs behaviour-adjusted ratio |
| \<FOV>\_maps.png, tiles\_\<FOV>.csv, segments\_\<FOV>.csv | runs aligned to the first, common region, fibre atlas, 25-um tile maps per condition, change map, per-fibre-segment rates |

What is applied by default: events starting in a flagged frame (± 1 frame) are dropped and that time is removed from the rate denominator; CFUs marked noise and their events are dropped; added regions are included. Options: `--min-p 2` keeps only events AQuA2 rates significant (−log10 p ≥ 2); `--gate still` restricts to quiet-wakefulness frames — not useful for `NavRearChrTIGRE-2608-1`, whose recordings are 97 % running; `--baseline d260917` changes the reference condition; `--no-curation`, `--keep-flagged-events`, `--event-source all` switch the filters off for comparison; --no-pupil skips the pupil block.

**Pupil and activity.** For every FOV whose runs have a pupil trace, Step 6 also asks two things. *Does arousal go with afferent activity?* Per run the mean CFU dF/F (the smoothed event rate when there is no `_res_cfu.mat`) is correlated with the pupil diameter at lags of −20…+20 s (Spearman; positive lag = activity follows the pupil): `pupil_dff_r_best` at `pupil_dff_lag_s`, a shuffle p from circular shifts of the pupil trace (`pupil_dff_p`; the traces are autocorrelated, a plain p would be far too small), the partial r with running taken out (`pupil_dff_partial_r` — pupil and running go together, r ≈ 0.45 in `260911_001`, so a raw pupil correlation may be a running correlation), r on still frames only (`pupil_dff_r_still`, needs 40+ still frames), and per CFU whether it is coupled (`pupil_cfu_<FOV>.csv`, counted in `n_cfu_pupil_up/down`). *Does the optostimulation change the pupil, and is the activity change just arousal?* `pupil_median_px` and `pupil_range_pct` pre vs post say whether the arousal state differs; the pupil terciles are then set once per FOV (pooled over pre and post) and the event rate is taken inside each (`rate_low/mid/high_pupil`), so `rate change at matched pupil: low x…, mid x…, high x…` compares pre with post at the same pupil size — a change that holds in every tercile is not explained by arousal, one that is only there because post had a larger pupil disappears within them. Excluded frames are left out; a tercile rate needs 15 s of data, the ratio 30 s in both conditions. Pupil pixels are comparable within one mouse and camera setting only.

**The comparison itself (v27; `ANALYSIS_PLAN.md` in the pipeline folder has the reasoning).** Three problems are built in. The field is never exactly the same: runs are aligned to the first run of the experiment and only the *common region* is compared (a run that does not match is reported `NOT ALIGNED`). Afferents do not respond at the same place twice, so CFU matching is not the unit: the field-level rate in the common region (also per mm² of fibre), 25-µm tiles on a *fibre atlas* segmented from the averaged image (map correlation, fraction of active tiles active in one condition only = activity that moved), and *fibre segments* (a fibre crossing another stays one segment; events assigned to the segment they lie on; segments up/down beyond the behaviour factor). Locomotion and arousal differ between pre and post, so every experiment gets five rate ratios in `experiments_summary.csv` and `<mouse>_effects.png`: raw, still frames only, running frames only, matched pupil size, and the **behaviour-adjusted** ratio from a Poisson regression on 20-s bins (condition + running speed + running state + pupil; 95 % CI from a 2-min block bootstrap) — the number to lead with; the matched-state ratios should agree with it and the tile/segment numbers are computed with that behaviour factor taken out. A real change: adjusted CI excludes 1, the matched ratios point the same way, and the change is not one hotspot. `--no-maps` skips the spatial part, `--tile-um` sets the tile size.

Read the statistics with the caveats in the summary: the event-level p values treat events as independent and are descriptive only. The unit that carries is one value per animal (average the runs of a condition within a FOV first, then across FOVs), compared between treated and control mice, or pre vs post within mice. Keep a FOV's AQuA2 settings identical across its recordings (Step 3 guarantees this) and prefer event rate and peak dF over peak dF/F for pre vs post, because mCitrine bleaches within a session.

## Adding a day, changing settings, re-running

Every step skips what is already done, so the normal loop after a new imaging day is the full sequence again and only the new recordings are processed:

1. `unpack_bin.py` — unpacks the new folders only.
2. `aqua_batch.py` — queues the new tiffs only (the old ones already have `_AQuA2.mat`), reusing each FOV's settings file.
3. `run_aqua_batch.m`, then `run_cfu_batch.m` — the CFU script only touches results without an up-to-date CFU file.
4. `cfu_curate.py` — curate the new runs; "Propagate to this FOV" carries the earlier decisions onto them.
5. `aqua_analyze.py` — recomputes everything with the new day as an extra condition.

A new FOV on a later day (a fourth FOV for a mouse that had three) needs its own `FOV4_Opt.csv`; its conditions are just the days it exists on.

**Changing one FOV's settings.** Re-tune in the GUI, overwrite that FOV's `_Opt.csv`, then redo only that FOV:

```
python aqua_batch.py "<mouse folder>" --aqua-dir "D:\AQuA2" --um-per-px 0.792 --border 5 --fov 2 --overwrite
```

then the two MATLAB scripts and `aqua_analyze.py`. Your curation of that FOV is kept but was made on the old CFUs — look through it again in `cfu_curate.py`.

**Redoing an animal from scratch:** `aqua_batch.py --overwrite` (all FOVs), the MATLAB scripts, curation, analysis. Deleting the `_AQuA2*` files next to a tiff has the same effect for that run.

**After `aqua_clean.py` on old tiffs:** all AQuA2 results of that animal are stale — redo with `--overwrite`.

## When something goes wrong

| You see | It means | Do |
| --- | --- | --- |
| `Unrecognized function or variable 'curSz'` when saving in the AQuA2 GUI | your AQuA2 predates the fix (2025-12) | update AQuA2, or add `curSz = []; dsSclXY = [];` above `if opts.sz(3)>1` in `saveExp.m`; the .mat and csv were already written when it crashed |
| AQuA2 files named `..._bin31f_1_Ch1.csv` | the tiff name had a dot (`1.98s`) and MATLAB cut it there | harmless for the GUI; the batch names its outputs correctly, and new tiffs are `1p98s` |
| `NO settings file` for a FOV in the `aqua_batch.py` printout | no `*_Opt.csv` naming that FOV in the animal folder | save one from the GUI (Step 2); the FOV is skipped until then |
| `was saved for day 260911; applying it to 260917 too` | the settings file names a day | intended: the day-1 settings are used for every day of that FOV |
| `found twice — using …, ignoring the copy in test_run` | a recording exists in two folders | fine; the copy in the run's own folder is used |
| `no binned tiff` for a run | Step 1 was not run, or the bin length is far from 2 s | run `unpack_bin.py`, or pass `--sec-tol` |
| `Missing Python package: h5py` (or pandas, tifffile) | the package is not in this Anaconda | `pip install h5py pandas tifffile scipy matplotlib` |
| `already corrected … skipping` from `aqua_clean.py` | the correction was applied before | nothing to do; `--force` re-does it with the current rules (needed once for tiffs treated by an older version) |
| `AQuA2 results exist … computed on the PREVIOUS movie` | the movie itself was changed after AQuA2 ran (only --drop, --fill or --restore-raw do that; decisions never do) | re-run Step 3 with `--overwrite` for that animal or FOV |
| `no _res_cfu.mat` in `cfu_curate.py` or `CFUs: no _res_cfu.mat` in the analysis | `run_cfu_batch.m` did not run for that result | run it; the event statistics are unaffected |
| `NOT aligned … matching CFUs on the unshifted maps` | the mean images of two runs of a FOV did not match well enough to trust a shift | check `<FOV>_cfu_map.png`; if the field really moved, the runs may be different planes |
| hundreds of events in a few frames of one run | a whole-field brightness step got through unflagged | open that run's `_globalfix.png`; if those frames are not flagged, lower `--min-rel` in `aqua_clean.py --force` (decisions apply at Step 6, no AQuA2 re-run) |
| the movie is shorter than the recording, or an AQuA2 frame number does not match the recorded frame | frames were dropped from that movie (explicit --drop, or v23/v24 defaults) | column tiff_frame of \_frames.csv maps movie frames to recorded frames, Step 6 converts event frames back and excludes the missing time; to get the frames back: --restore-raw if ...\_raw.tif exists, else re-unpack the run |
| frozen (still) stretches in the movie or the AQuA2 overlay, often with big slow events of 100 s or more | a v22 tiff: bad frames were replaced by a filler, and a long excursion was only partly flagged | --restore-raw or re-unpack that run, then the step again (every frame kept), then Step 3 with --overwrite |
| big slow events spanning half the field in one run | a state stretch that AQuA2 turned into events, or a slow change below the floors | in the GUI (tab 1b) mark that stretch Exclude whole stretch — no AQuA2 re-run needed; if it is not even flagged, aqua_clean.py --force --shape-rel 0.03 |
| CHECK: excluded but not visibly different (clean log, GUI table) | an excluded stretch looks like its kept neighbours | open the check sheet; if the frames are fine, select the stretch and Keep events (or raise --min-rel / --shape-rel and re-run with --force) |
| N frame(s) are no longer on disk | an older version dropped or replaced them; they stay excluded because their content is gone | aqua_clean.py --restore-raw if ...\_raw.tif exists, else re-unpack that run (unpack --overwrite); then decide per stretch |
| a long amber (state) stretch in the QC plot | a slow change the rules cannot classify (laser, focus, or widespread activity) | kept by default unless its background is brighter; look at the sheet and the running trace, decide (Exclude / Keep events) and note it; run_summary.csv says how many events start in kept state stretches |
| pupil: no pupil data (analysis) | no \*eye\*DLC\*.csv next to the .sbx of any run of that FOV | copy the DeepLabCut files into the run folders and run Step 6 again; it rebuilds the behaviour csv by itself (or unpack_bin.py --behavior-only) |
| pupil: N eye frames for M imaging frames — resampling | the eye camera and the 2P clock differ by a frame or so (14751 vs 14750 is normal) | nothing; a ratio far from 1 means the wrong DLC file is in the folder |
| still-only n/a in the pupil summary | fewer than 40 still frames with a pupil value in that run (a running animal) | use the partial r (running removed) instead |

`aqua_gui` fails at start with `Property assignment is not allowed when the object is empty` in `addCon_dat` line 123: your MATLAB creates the axes toolbar lazily, so AQuA2's `Toolbar.Visible = 'off'` lines fail. In MATLAB run `patch_aqua2_toolbar('D:\AQuA2')` (the script is in the afferent2p folder; it rewrites those seven lines to `Toolbar = []` and keeps a `.bak`), then restart MATLAB. Every fresh download of AQuA2 brings the lines back, so re-run it after each update.

If the GUI windows (`unpack_bin.py`, `aqua_batch.py`, `cfu_curate.py`) do not open, the command-line forms in this guide do the same thing.

## All commands in order

Anaconda prompt (replace the animal folder; `D:\AQuA2` = your AQuA2 folder):

```
cd D:\2photon\Chao\GCAMP_affrerents\afferent2p
set M=D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1

python unpack_bin.py "%M%" --sec 2
python unpack_bin.py "%M%" --behavior-only  :: only if the DLC eye csv arrived after unpacking (Step 6 also does this itself)
python aqua_clean.py "%M%" --force        :: once for tiffs treated by an older version (v29: per-frame check + ImageJ check stacks); then check the stacks and decide per stretch (GUI tab 1b)
:: --- AQuA2 GUI: tune each FOV once, SaveOpts -> %M%\FOV1_Opt.csv, FOV2_Opt.csv ...
python aqua_batch.py "%M%" --aqua-dir "D:\AQuA2" --um-per-px 0.792 --border 5
```

MATLAB:

```
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_aqua_batch.m')
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_cfu_batch.m')
```

Anaconda prompt again:

```
python cfu_curate.py "%M%"
python aqua_analyze.py "%M%" --um-per-px 0.792
```

Results: `%M%\aqua_analysis\summary.txt`, `run_summary.csv`, the `_overview.png` and `_cfu_map.png` figures, and \_pupil.png with pupil\_\<FOV>.csv when the runs have a pupil trace. The same sequence for the spontaneous mice (`NavTIGRETreg-2608-1`, three or four FOVs per day) — the tools tell the paradigms apart from the animal name.
