# afferent2p — pre/post optostim GCaMP afferent pipeline (Python)

Replaces the Excel-driven MATLAB front end (`UnpackSBX` + `ProcessScanboxData_concatFirst`)
with one config file and one command. Reads Scanbox `.sbx` directly, registers the
baseline and post runs **together to a single reference** (so the same afferents are
comparable pre vs post), and produces QC images, registered tiffs for AQuA2, Suite2p
ROIs/traces (including a MATLAB-readable `Fall.mat`), locomotion state, and a
first-pass pre/post event-rate comparison.

## Install (once, on the analysis PC)

Another PC that sees the data through a mapped drive (the server's `D:` as
`Z:`): see `REMOTE_PC.md` — copy this folder locally, `setup_remote.bat`,
then `remap_batch.py` (or rebuild the batch there) so the MATLAB scripts use
`Z:\...` and that PC's AQuA2 folder. Nothing in the data folders changes.

Python 3.9+ (3.10/3.11 recommended). Then:

```
pip install suite2p sbxreader tifffile pyyaml scipy matplotlib
```

Tested with suite2p 1.1.0. A GPU is not required (`torch_device: cpu` in the config);
if the PC has an NVIDIA GPU, set `torch_device: cuda` for much faster registration.

The AQuA2 route also needs `pandas` and `h5py` (reading AQuA2's v7.3 `.mat`
files): `pip install pandas h5py` — both are already in Anaconda except h5py in
some installs.

## Use — step 1 of the AQuA2 route (unpack + bin, no registration): `unpack_bin.py`  ← start here

Matches the colleague's pipeline: `UnpackSBX.m` (binT) → binned tif → AQuA2
(which does its own registration). Python replacement of `WriteSbxPlaneTif`:

```
python unpack_bin.py                       # GUI: add session/mouse folders, Run
python unpack_bin.py "D:\2photon\Chao\NavTIGRETreg-2608-1" --sec 2     # CLI
```

For every `.sbx` found (a mouse folder adds all its session folders) it writes,
**next to the .sbx** — no extra output folders:

| file | what |
|---|---|
| `<run>_bin31f_2.00s.tif` | mean of every 31 frames (2 s at 15.62 Hz), **uint16 ImageJ TIFF** — opens in plain ImageJ, Fiji, MATLAB; this is the AQuA2 input |
| `<run>_bin31f_2.00s_meanproj.tif` | mean image of the whole recording, for checking by eye |
| `<run>_bin31f_2.00s.json` | frame rate after binning etc., for AQuA2's parameters |

Options: `--sec` (seconds per output frame, default 2), `--frames N` (bin size in
frames, overrides `--sec`; e.g. `--frames 8` reproduces `binT=8`), `--chan`,
`--crop L,R,T,B`, `--overwrite`. Pixel values are the same as the MATLAB output
(65535 − raw, averaged, uint16), no scaling.

Each channel-1 tiff also gets `_frames.csv`, `_globalfix.png` and
`_flagged_frames.png/.csv`, `_check_frames.tif` and `_check_movie.tif` from the global brightness
step described in the next section (step 1b); every frame stays in the movie. An older dotted-name copy (`..._1.98s.tif`) next to the current `..._1p98s.tif` is ignored by every tool and can be deleted with its sidecars.

Two-channel recordings: by default every channel in the file is written —
channel 1 (green) keeps the plain name `<run>_bin31f_2.00s.tif`, channel 2 (red
tracer) goes next to it as `<run>_bin31f_2.00s_ch2.tif` with its own
`_meanproj.tif` and `.json`; the behaviour csv is written once. `--chan 1` /
`--chan 2` (or the GUI's channel box) limit it to one channel.

## AQuA2 route, end to end

```
unpack_bin.py (bins + flattens global brightness + flags bad frames)  →  aqua_clean.py (check & decide)
      →  aqua_tune.py --plan → MATLAB run_aqua_tune.m → aqua_tune.py --report --apply   [detection settings by measurement;
         optional: AQuA2 GUI once per FOV (SaveOpts); a FOV without any file gets default_Opt.csv]
      →  aqua_batch.py  →  MATLAB: run_aqua_and_cfu.m (AQuA2 batch, then CFU batch — one paste)
      →  cfu_curate.py  →  aqua_analyze.py  →  aqua_report.py (the report) / aqua_report.py --qc (the QC report)
```

The GUI's **▶ Run all (guided)** button runs everything after `unpack_bin.py`
in one go, for a Treg mouse as for a RearChr mouse: at the two MATLAB steps it
shows the command, copies it to the clipboard, checks the outputs every 30 s
and continues by itself (or starts `matlab -batch` from the window when that
box is ticked); it never stops on the AQuA2 folder — the generated scripts
find AQuA2 on MATLAB's path. It ends with both reports
(`aqua_analysis\<mouse>_report.html`: the general representation of the data
— sessions on the day's clock, activity by running / resting state, ECDFs,
acute / chronic change of every measure —, the CFU × time heatmaps with
up / down / unchanged CFUs, the locomotion-evoked activity (Yu's trial
pipeline), pupil, coupling GLMs, the behaviour-controlled pre/post or
across-day change with a verdict per experiment, the time course;
`<mouse>_qc.html`: the evidence — a QC table with warnings and, per recording,
every detected event on the real frames with clips, rasters and detection
maps); curation stays a manual look.

`python pipeline_gui.py` puts all of these in one window — tabs *1 Unpack ·
2 Clean & check · 3 Detection settings · 4 Batch · 5 MATLAB · 6 Curate ·
7 Analyze · 8 Report*, in the order they are used — with a per-recording table
of the flagged stretches, check verdicts and your exclude/keep decisions, the
MATLAB commands to copy, a log, and a status strip above the tabs that shows where the
animal is (one chip per step from the files on disk + the next action; also
`python pipeline_status.py <mouse folder>`); it calls the same functions.

Tiffs unpacked before v22 (named `..._bin31f_1.98s.tif`) get the brightness
step with `python aqua_clean.py <mouse folder>` — 15 s per recording, no
re-unpacking; tiffs treated by an older version are brought to the current
rules with `python aqua_clean.py <mouse folder> --force`.

This is the main route. Each FOV gets its detection settings decided **once** —
by measurement (tab 3) or by eye in the AQuA2 GUI; every other recording of that FOV is then processed with
exactly those settings, and the results are read back in Python for the
across-day / pre-post comparison. `aqua_tune.py` (GUIDE Step 3, GUI tab 3)
replaces the by-eye choice with measurements that run
overnight: for one movie per day and FOV it runs a grid of `thrARScl` × `minSize`
on the movie, on a null copy (per-pixel time shifts — every event found is
false) and on a copy with synthetic transients of known amplitude, size and
duration injected on the fibres, and recommends the setting with the highest
recall among those with a false-alarm share ≤ 10 % and most events on fibres,
then one common setting per FOV that meets the rule on every day; `--apply`
writes it as `aqua_opts\<experiment>_Opt.csv` (a Treg mouse's experiment is
the FOV across days, so its days always share the setting).

**1. Unpack** (2 s per frame, as the lab does):

```
python unpack_bin.py <mouse folder> --sec 2
```

Writes `<run>_bin31f_1p98s.tif` (uint16, ImageJ-readable) next to each `.sbx`
— the tag is spelled `1p98s`, not `1.98s`, because MATLAB's `fileparts()` cuts
a name at its last dot and the AQuA2 GUI then saved everything as
`..._bin31f_1_*`; every tool here accepts both spellings — plus `_meanproj.tif`,
a `.json` with the exact frame time, `_behavior.csv` (velocity / still fraction
/ pupil, binned the same way), and two-channel runs also get `_ch2.tif`.

**2 (Clean & check). Whole-field brightness changes are flattened, bad frames are flagged
and decided on — nothing is removed.** A step in global brightness (the brain
shifting when the animal brakes, a hemodynamic swing, a laser/PMT change,
stray light) is brighter than any afferent event and present everywhere at
once, so AQuA2 carves it into hundreds of "events". In a real recording
(`NavRearChrTIGRE-2608-1_260911_005`) **227 of 291 events started within 9
frames** of one such flash. No threshold can separate this from real events,
so `unpack_bin.py` (default on; `--no-global-fix` to disable, `--keep-raw` to
also keep the uncorrected movie as `..._raw.tif`) and `aqua_clean.py`:

* divide every frame by its own spatial median — the global level becomes
  flat (bleaching and steps included); dF/F and AQuA2's z-scores are
  scale-free, so events are untouched (verified: 0.69 → 0.685 on a synthetic
  event) while the median of ~400 000 pixels does not move for sparse events;
* flag two classes of stretches. **Artefact**: a jump between consecutive
  frames beyond max(4 robust SD, 3 %), a spike (≤ 10 s and returning), or a
  collapse of the correlation with the reference image; tails while the level
  is still moving, +1 frame each side, clean islands < 10 s between two
  stretches. **State**: a level plateau longer than a spike or not returning,
  or the p10/median ratio — which a pure gain change cannot move — more than
  5 % from its clean-frame 5-min moving median (iterated, edges extrapolated).
  A state change is a laser/focus change or widespread activity; the trace
  cannot tell which, so it is only annotated;
* check every frame of a flagged stretch against the kept frames around it
  (whole-field change; smoothed RMS % difference from the clean neighbours,
  relative to how much two clean frames differ). The **flagged frames** of an
  artefact stretch are the frames that differ by ≥ 10 % of the frame and ≥ 3×
  the clean floor — the flash itself, not the tails, margins and mild frames
  the flagging rules wrap around it (adjustable: `--core-min-pct`,
  `--core-ratio`);
* apply a **decision per stretch**: *exclude* (only its flagged frames:
  `aqua_analyze.py` drops events starting in them ± 1 frame and their time),
  *exclude-all* (the whole stretch), *frames:68,76-77* (exactly these), or
  *keep*. Defaults: artefact with flagged frames → exclude; artefact without →
  keep; state → keep (annotated; a note warns when its background ends up
  brighter than normal after the gain). Your own decisions (`pipeline_gui.py`
  tab 2, or `aqua_clean.py --decisions-only --keep 13-19 --exclude 63-87
  --exclude-all 38-42 --exclude-frames "94-102:99"`) live in
  `<run>_..._frames_decisions.json`, take effect at Step 7 (analyze) without re-running
  AQuA2, and survive re-runs of the step;
* keep every frame in the movie. `aqua_clean.py --drop --force` leaves the
  *excluded* frames out (the csv column `tiff_frame` maps movie frames to
  recorded frames; `_ch2.tif` loses the same frames) and `--fill` replaces
  them by a frozen filler — both explicit, never the default;
* write `<run>_..._frames.csv` (per recorded frame: level, gain, class,
  stretch, decision, flagged_frame, excluded, tiff_frame, reason, the
  per-frame check), `<run>_..._globalfix.png` (level and shape traces, red =
  artefact, amber = state, black bar = excluded, running speed underneath),
  and the check files for ImageJ: `<run>_..._check_frames.tif` (16-bit stack
  of every flagged stretch with 6 frames of context on each side, full
  resolution, each slice labelled with its recorded frame number, stretch,
  numbers, flagged / EXCLUDED — scroll through it as you would through the
  raw movie), `<run>_..._check_movie.tif` (the same slices as a composite:
  frame, % difference from the neighbours, traces with a cursor) and the
  overview sheet `<run>_..._flagged_frames.png/.csv`. `--report-only` makes
  all of this without touching the movie; `--stacks-only` rewrites the stacks
  with the current decisions in their labels.

What this gives on the real recording `..._001`: nine stretches flagged;
only the flash frames 67–68, 76–77 and 99 (+25 to +46 % whole field, 7–12×
the clean floor) are excluded — the ±5–9 % frames next to them, the mild
stretches (13–19, 38–42, 50–54, 389–394, 460–465) and the two state stretches
(55–62, the dim state; 88–93) are kept.

Why not drop by default: whether a plateau is an artefact or widespread
activity (locomotion) is not decidable from the movie, and a dropped frame is
gone; an excluded event is one decision away from being counted again.
Earlier versions dropped or filled the flagged frames: v22's filler ramped
from a clean frame into a dim state and AQuA2 found 30 "events" of 100–180 s
in the ramp; frozen fillers also dilute AQuA2's noise estimate (mean squared
frame difference over the whole movie). Tiffs treated by those versions get
the current rules with `--force`; frames they removed stay excluded (their
content is gone) unless `--restore-raw` finds a `..._raw.tif`, else re-unpack.

`aqua_analyze.py` converts AQuA2's frame numbers back to recorded frames when
frames were dropped (`t0_frame` = recorded frame, `t0_frame_movie` = frame in
the movie AQuA2 saw), drops events starting in excluded frames (± 1 frame,
`--flag-margin`; `--keep-flagged-events` keeps them), removes the excluded
time from the rate denominator, and reports per run `flagged_frames`
(artefact), `state_frames`, `excluded_frames`, `n_events_dropped_flagged` and
`n_events_in_kept_state` (events starting in state stretches you kept — the
number to look at when deciding). Step 7 (analyze) also warns when an AQuA2 result
predates a drop.

**AQuA2 itself:** after every download run `patch_aqua2_toolbar('D:\\AQuA2')` in
MATLAB (script in this folder) — otherwise `aqua_gui` stops at start with
*Property assignment is not allowed when the object is empty* (`addCon_dat`
line 123), because this MATLAB creates the axes toolbar lazily and AQuA2's
`Toolbar.Visible = 'off'` lines fail; the patch rewrites them to `Toolbar = []`.

**3. Detection settings — by measurement (`aqua_tune.py`, GUIDE Step 3) or by eye in the AQuA2 GUI.** By eye: load the first recording of that FOV
(after step 1b, so you tune on what the batch will see); set temporal resolution
from the `.json` (31-frame bins at 15.62 Hz = **1.98 s**, not 2.00), `0.792`
µm/px, border exclusion `5`; consider `detectGlo` off — it did not catch the
step above (its one global event was a 310-s drift); work through Pre → Act →
Temporal → Spatial until the detection looks right (see `aqua2-gui-setup-notes`
in the project docs for the parameter-by-parameter reasoning). Then **Detect tab
→ SaveOpts** and save the file **in the animal folder** with the FOV in its name:

```
D:\2photon\Chao\NavTIGRETreg-2608-1\FOV1_Opt.csv
D:\2photon\Chao\NavTIGRETreg-2608-1\FOV2_Opt.csv          ...
```

Any name containing `FOV1`, `FOV2`, … works, so the GUI's own suggestion is fine
as long as you add the FOV to it — e.g.
`NavRearChrTIGRE-2608-1_260911_001_bin31f_1.98sFOV1_Opt.csv`. The file is then
used for that FOV in **every** recording and on **every day**: FOV1 across days
−1/5/7/9 for the spontaneous mice, and FOV1's four runs before and after the
stim hour on each day for the RearChr mice. A file that names a day
(`..._260911_...FOV1_Opt.csv`) is still applied to the other days, and the log
says so. Other places that are searched, in order: `<mouse>\<experiment>_Opt.csv`
(per-day override, e.g. `NavRearChrTIGRE-2608-1_260917_FOV1_Opt.csv`),
`<mouse>\aqua_opts\`, a file naming one of the experiment's runs, then any
`*_Opt.csv` inside the run folders.

**4. Build the batch** (no MATLAB needed):

```
python aqua_batch.py <mouse folder> --aqua-dir "D:\AQuA2" --um-per-px 0.792 --border 5
python aqua_batch.py            # or the GUI: add mouse folders, press Build
```

It groups runs into FOV experiments with the same rules as the suite2p planner:

* **no "RearChr" in the animal name** — FOV *k* = run `00k` of each day, so
  `260911_001`, `260917_001`, `260919_001` are one experiment with the days as
  conditions and the earliest day as baseline. A day with a fourth FOV simply
  adds `FOV4`, whose conditions are only the days it exists on — nothing else
  shifts.
* **"RearChr"** — one experiment per (day, FOV): runs 1–4 of a day are the pre
  block (FOV1, FOV2, FOV1, FOV2) and 5–8 the post block, so FOV1 gets
  `001, 003` = pre and `005, 007` = post. A short day is treated as pre-only
  rather than split in half; `--stim-after-run` overrides.

Pipeline output folders (`afferent2p_out`, `AQuA_input`, `aqua_batch`,
`aqua_analysis`, `suite2p`, …) are never scanned for raw runs, and a recording
that appears twice (a copy in a `test_run` folder, say) is used once and the
duplicate reported. The tiff is matched by its actual bin length, so
`_bin31f_1.98s.tif` is found when you ask for 2 s (`--sec-tol`, default 0.25 s).

It then finds each FOV's `*_Opt.csv` and writes `<mouse>\aqua_batch\`:

| file | what it is |
|---|---|
| `batch.csv` | one row per tiff: `inputPath,preset,outputPath,cellMaskPath,landmarkMaskPath` (AQuA2's format; empty `outputPath` = results next to the input) |
| `parameters_for_batch.csv` | `Preset1..PresetN` columns, one preset per FOV, in AQuA2's own layout |
| `run_aqua_batch.m` | `cd` to AQuA2 and call `aqua_cmd_batch(batch.csv, parameters_for_batch.csv)` |
| `run_cfu_batch.m` | CFU detection on every `_AQuA2.mat` (mirror of `CFU_cmd_batch.m`, but skips files already done — with the same CFU rules — and does not rewrite the big result file). The CFU rule is AQuA2's published overlap clustering with the pipeline's overlap threshold 0.3 and min events 1 (`--cfu-overlap`, `--cfu-min-events`; AQuA2's own defaults 0.5 / 3) — GUIDE.md "CFU rules" says why, and why AQuA2's newer rules (28 Aug 2026 on: basin splitting + post-merge) are not used. Since v42 the script also writes `<run>_AQuA2_cfu_pre.mat` (AQuA2's pairwise event overlaps + per-event centre / size / frames) for `cfu_overlap_check.py` — the check that re-asks per animal whether the threshold fits (GUI tab 4) |
| `manifest.json` | tiff → experiment, FOV, condition, preset, frame time |

The frame time always comes from the tiff's `.json`, never from the Opt file (a
GUI file saved before you set the resolution would otherwise poison the batch);
`--um-per-px` / `--border` likewise override what the Opt file stored. A FOV with
no settings file is listed as such and left out of the batch. Tiffs that already
have an `_AQuA2.mat` are skipped unless `--overwrite`.

**5. Run it in MATLAB** (one paste; the two scripts can also be run one at a time):

```matlab
run('D:\...\<mouse>\aqua_batch\run_aqua_and_cfu.m')      % = run_aqua_batch.m, then run_cfu_batch.m
```

The scripts find AQuA2 on MATLAB's path (`which('aqua_cmd_batch')`) when the
folder given to `aqua_batch.py --aqua-dir` is not visible from that MATLAB.

**6. Check the CFUs by eye and mark noise: `cfu_curate.py`.**

```
python cfu_curate.py <mouse folder>                 # the GUI
python cfu_curate.py <mouse folder> --list          # what is there, per run
python cfu_curate.py <mouse folder> --propagate     # copy decisions within each FOV
```

The window shows one run at a time: the mean image with every CFU outlined (red
= real, grey = noise, cyan = added by hand, yellow = selected), the selected
CFU's dF/F curve with its event onsets, its numbers (events, events/min, size,
peak dF/F), and a raster of all CFUs' onsets underneath. Click an outline to
select it, then `k` for real or `x` for noise; `n` / `p` step through runs,
`s` saves. Decisions belong to one CFU result: after the CFU batch remakes a
run's `_res_cfu.mat` (other CFU rules, e.g. another overlap threshold) the
window ignores the earlier decisions, opens the run clean and keeps the old
file as `_res_cfu_curated_superseded_<date>.json` on the first save (v43);
`--list` ends with how many runs are curated under the current rules.

*Adding a CFU AQuA2 missed*: press `a`, click the vertices of a region, Enter to
close it. The trace is then measured from the binned tiff — mean over the region,
dF/F against a moving 10th-percentile baseline, events = crossings of z > 3 on
the MAD-scaled trace — and stored with `source: manual`, because those events did
not come out of AQuA2's detector and should not be pooled with its events
silently. (The region is placed correctly even though AQuA2 crops
`regMaskGap` pixels off each border: the offset between its frame and the tiff's
is measured and applied.)

*Propagating*: "Propagate to this FOV" copies the decisions of the current run to
the CFUs of the other runs of the same FOV that overlap them, using the same
matching as the analysis (mean images aligned, then IoU). A CFU identity marked
noise in one run is marked noise in all of them — noise does not become signal
because the next recording drew the region slightly differently.

Decisions live in a sidecar next to each result,
`<run>_..._AQuA2_res_cfu_curated.json`, holding `{"rejected": [...], "added":
[...]}`. AQuA2's own files are never touched, so re-running `run_cfu_batch.m`
does not lose the curation. Without the GUI:
`python cfu_curate.py <run tiff> --reject 3,7` / `--keep 3` / `--keep-all`.

(AQuA2's own `CFU_Viewer_Standalone.m` is a good cross-check for looking at a
single result in MATLAB, but it keeps no keep/reject record for a batch.)

**7. Read the results back, the general representation, locomotion trials, heatmaps and the comparison of conditions:**


```
python aqua_analyze.py <mouse folder> --um-per-px 0.792
python aqua_analyze.py <mouse folder> --gate still --min-p 2      # quiet wakefulness only, significant events only
python aqua_analyze.py <mouse folder> --no-glm                    # skip the Loco2Fluor / LocoPupil2Fluor GLMs (the slowest part)
```

Since v45 the check buttons of GUI tab 2 (*Show check sheet*, *Frame stack
(ImageJ)*, *Composite stack*, *QC plot*) and the per-stretch decision buttons
sit directly under their tables and can no longer be squeezed out of the tab
by a short window (from v36 on they were silently dropped at the default
window size). The check sheet window has a *Zoom / pan* button that opens the
sheet and the QC plot in the picture viewer.

Since v44 a run that breaks the naming rule (run k = FOV k of the day for a
Treg mouse; FOV1, FOV2, FOV1, FOV2 pre then post for a RearChr mouse) is
described in `<mouse>\run_labels.csv` — one row per such run: `fov`,
`condition` (`d<date>_pre` / `d<date>_post` for a second-injection day with
pre and post recordings), `use` = 0 to leave a run out (a lost FOV, a shifted
plane), `note`. Every step reads it through the planner (batch, tuning,
analysis, curation, reports); GUI tab 1 creates / completes the file and
prints the resulting plan; `day_labels.csv` labels the days. GUIDE.md "Which
run is which FOV".

Since v40 an output that cannot be overwritten no longer stops the batch or
the plan: a csv left open in Excel, or a file that belongs to another account
(a pipeline window started from an Administrator prompt writes files owned by
Administrators, which a window started the normal way cannot change) is
written under the next free name (`parameters_for_batch_2.csv`,
`run_aqua_batch_2.m`, …), `aqua_batch\manifest.json` records the names and
the GUI's MATLAB commands use them; only the manifest itself must stay
writable (the message then says what to close or whose files they are). The
sweep script does the same for its log (`tune_log_<date>.txt` when the old
`tune_log.txt` cannot be appended to).
`run_gui.bat` activates the base Anaconda by itself, so a double-click in
Explorer opens the window. GUIDE.md ("When something goes wrong") explains
the Administrator-prompt effects, including why `matlab -batch` started from
such a window cannot see a mapped H: drive.

Since v39 the tuning is presented in two numbers and one picture: CATCH (the
share of 72 planted transients a setting finds) and NOISE (the events it
invents in a scrambled copy of the movie, as a share of what it reports in
the real one); `<mouse>_tune.png` plots every setting as a dot at (noise,
catch) per day and FOV, pre and post together, the rule's pick starred and
the pair the day uses ringed, with plain sentences under each panel (the
same open `tune_recommendation.txt`; the curves moved to
`<mouse>_tune_detail.png`). Tab 3 *Pick by eye…* shows the few settings
worth a look side by side on the real movie — events as circles, the noise
events overlaid on request, a frame slider that plays the movie in every
panel — and writes the setting you pick (`aqua_tune.apply_pick`, rule
"picked" in the tuning json). The grace rule needs at least 4 real events.

Since v38 the GUI shows pictures in a viewer of its own (tab 3 *Figure*):
the tuning figure and every "why" picture at full resolution, zoom at the
cursor (Ctrl + wheel, + / −), pan by dragging, *Fit* / *1:1*, the file path
in a copyable box with *Copy path* / *Open folder* / *Open file*, and
*Redraw at 2× resolution* (`aqua_tune.py --figure-only --dpi 260`); the
figure itself got an explanatory header, the amplitude classes in shades
(light = dimmest) whatever their values, and the grid verdict in each row's
title.

Since v37 the tuning tells when it cannot reach a field, instead of choosing
among blind settings. The sweeps of Oct 21 looked insensitive on every
animal (recall ≤ 0.3 at every setting) because the injected transients of
v30–v36 (nominal 0.1 / 0.2 / 0.4 dF/F, reported at ~0.05–0.2) were dimmer
than nearly every real event; on NavRearChrTIGRE-2609-1 FOV 2 the grid also
stopped at minSize 50 while the field's events are 20–50 px, and the fibre
mask dropped the whole fibre network as a "round speck". v37: default
classes nominal 0.3 / 0.6 / 1.0 dF/F × 20 / 60 / 200 px, default grid
thrARScl 1.5 / 2 / 2.5 × minSize 20 / 50 / 100, the mask keeps networks
whole (old masks are made again), and the report gives a **grid verdict**
per movie (*GRID TOO STRICT* / *edge* / *TOO PERMISSIVE* / *WEAK*, with the
values to sweep next; `--plan --extend` or tab 3 *Extend the grid* adds
them and sweeps only the new combinations — finished results are matched by
their setting values and kept; *Run all* does one such round by itself), a
classes check (the brightest injected class as AQuA2 reports it against the
real events), a fibre-mask check, and a floor for the matched rule
(`--min-recall` 0.30: a pair matched below it is *TOO WEAK* and the day
not comparable). A result tuned by hand in the AQuA2 GUI, exported once with
`aqua_export_events.m` (tab 3 *Hand result → reference*), becomes the
reference of its run: per setting the report shows the share of the
hand-accepted events recovered (**hand%**) and the share of extra events
(**extra%**), names the setting closest to the hand result, and draws the
hand events in the pictures.

Since v36 the tuning of a RearChr day no longer forces pre and post onto
one setting: with `--day-rule matched` (the default; GUI tab 3 *RearChr day
(pre + post tuned)*) each condition gets the setting that detects the
injected transients with the **same measured recall** (within
`--tol-recall`, 0.05) **and the same false-alarm rate** on the null movie
as the other — the same setting when it serves both, one per condition when
the signal changed between the recordings (the manual "pre looks good but
post does not" problem), and **NOT COMPARABLE** when no pair of clean
settings matches, in which case the common setting is written and the
analysis / report flag the day instead of reporting a ratio. The choice is
drawn as pictures (`aqua_tune\<day>_FOV<k>_why_day.png`: pre above post at
the same thresholds, each movie's own setting framed green, at the top of
`tune_why.html`); `--apply` writes `aqua_opts\<exp>_pre_Opt.csv` /
`_post_Opt.csv` (or the single `<exp>_Opt.csv`) plus `<exp>_tuning.json`;
Step 4 picks the condition's file per run; Step 7 reports the settings,
recall and false-alarm rate per condition and the rates with the expected
false alarms subtracted (`ratio_fpcorr`), and the report's verdict carries
them. In the GUI the tabs and the log share a draggable divider (the log
never disappears) and **Stop** now ends any running step at its next log
line.

Since v35 `aqua_conditions.py` puts the stimulation protocols side by side
across animals — every RearChr animal's tab-7 tables under one parent folder,
each day × FOV experiment labelled by its day's protocol (sham / 2 Hz / 5 Hz
/ 5 Hz + 40 ms from `day_labels.csv`), the pre → post change of the event
rate (all / still / running frames), amplitude and CFU fate per protocol
with the behavioural state next to it, `<parent>\cross_animal\conditions.html`
(GUI tab 8, *Across animals*); the tuning report draws the sanity check
behind every recommended setting (`aqua_tune\<run>_why.png`,
`tune_why.html`: per setting what is found in the real, the null and the
injected movie, on the mean image); `reg_demo.py` turns the registration
comparison into a movie (as recorded / after AQuA2's own registration / the
legacy chain, side by side, zoomed, with a fixed outline of the reference
and the shift traces); the GUI shows what runs and how far it is (progress
bar, item i of n, elapsed time and the time left, from the `[i/n]` lines the
steps print and the result counts of the MATLAB steps; the running step's
chip is highlighted); `aqua_reset.py` / the GUI's *Fresh start…* parks
every output of an animal (stage by stage, into `_previous_<date>\`, or
deletes it) so that the whole pipeline runs again with the current version;
and an aborted recording (a `.sbx` without its info `.mat`, much smaller than
the others) is skipped by the planner and the unpack, its re-take (the next
free run number, e.g. 051 for 005) placed by recording time.

Since v34 the analysis compares the stimulation's effect on three triggers
side by side — running onsets (the locomotion trials), pupil dilations
during stillness (new; `afferent2p/triggers.py`) and resting activity — each
with the rate comparison's statistics (`afferent2p/stats.py`), so that a
session in which the animal never runs is read on its pupil dilations and
its rest instead of on nothing (`trigger_comparison.csv`,
`<exp>_triggers.png`, report section 3); and the tuning of an optostim mouse
is per day (the sweep runs on a pre and a post movie of every day and the
report says whether they agree; the day gets one setting).

Since v33 `reg_compare.py` compares up to five registrations of the same
recordings (AQuA2's own; the raw frames registered before binning in Python —
`afferent2p/fullreg.py`, also `unpack_bin.py --register`; both; the lab's
legacy MATLAB chain itself, CorrectExpt + RegisterSBX, run per recording by
`legacy_register_run.m` with `--add-legacy`; and that with AQuA2's on top) on
what AQuA2 detects, and the analysis
also runs a motion / registration check per run
(`afferent2p/motion.py`: how far the raw frames moved, what AQuA2's own
registration left behind, whether the event onsets fall on moving frames more
often than chance — `motion_summary.csv`, the QC report), and the curation
window and the reports show AQuA2's own registered movie (read out of the
`_AQuA2.mat`) instead of the raw binned tiff; the curation window's zoomed
view can be moved (right-drag, "move view" + left-drag, Ctrl+arrows).

Since v32 the analysis also writes the lab's presentation measurements —
sessions on the day's clock, activity by running / resting state with event
overlays, ECDFs, the CFU × time heatmaps (identities followed across runs,
rows read from the movies; up / down / unchanged CFUs), the locomotion-evoked
activity (strict and relaxed running-onset trials, three responsive rules,
per-CFU responsiveness, rest-only activity, behaviour check, GLMs), pupil
measures and the acute / chronic change of every measure — see GUIDE Step 7.
**8.** `python aqua_report.py <mouse folder>` writes the report, `--qc` the QC
report, `--all` both (GUIDE Step 8).

It reads each run's `_AQuA2_Ch1.csv` event feature table (rows = features,
columns = events — that is how `writetable` saves it), the `.json`, the
behaviour csv and, when present, `_AQuA2_res_cfu.mat`, and writes
`<mouse>\aqua_analysis\`:

| file | contents |
|---|---|
| `events_<experiment>.csv` | every event with run, condition, start time, and the locomotion/pupil state at its onset |
| `run_summary.csv` | one row per recording: events/min (and per mm²), median peak dF/F, dF, area, duration, rise, decay, AUC, plus still minutes and median pupil — **this is the unit to average for group statistics** |
| `condition_summary.csv` | one row per FOV × condition, pooled and as median-of-runs |
| `cfu_matched_<experiment>.csv` | CFUs matched across the runs of one FOV, with each CFU's events/min per run and per condition and the pre→post ratio |
| `<experiment>_overview.png` | event rate per run + ECDFs of amplitude, area, duration, rise, decay per condition |
| `<experiment>_cfu_map.png` | the CFU maps of every run, coloured by matched identity |
| `pupil_<experiment>.csv` | one row per run: pupil median / range, pupil–running r, pupil–activity r at its best lag with shuffle p, partial r vs running, still-only r, CFUs coupled up/down, events/min inside each pupil tercile (terciles pooled over the FOV) |
| `pupil_cfu_<experiment>.csv` | one row per CFU and run: Spearman r of its dF/F against the pupil and `up`/`down`/`no` against the shuffle band |
| `<experiment>_pupil.png` | per run pupil over activity with the lag curve; event rate per pupil tercile per condition; per-CFU r; pupil vs running; pupil distributions; activity vs pupil octiles |
| `experiments_summary.csv` | one row per experiment × condition vs its baseline: raw ratio, still-only and running-only ratios, matched-pupil ratio, behaviour-adjusted ratio with 95 % CI, deviance explained by behaviour / condition, common-region ratio, map r, % tiles switched, fibre segments up/down |
| `<mouse>_effects.png` | those ratios drawn, one row per experiment — the first figure to look at |
| `<experiment>_adjusted.png` | rate vs running per 20-s bin and condition, rate around running onsets, raw vs adjusted ratio |
| `<experiment>_maps.png`, `tiles_<experiment>.csv`, `segments_<experiment>.csv` | runs aligned, common region, fibre atlas, 25-µm tile maps per condition, change map, per-fibre-segment rates |
| `summary.txt` | the numbers in prose, with the change vs baseline, every behaviour control, where the activity is, the pupil coupling, and the statistics caveats |

`ANALYSIS_PLAN.md` explains the comparison: why CFU matching is not the unit
(afferents respond at different places), how the common region, the fibre
atlas, the tiles and the fibre segments replace it, and how the four
behaviour controls (still-only, running-only, matched pupil, Poisson
regression with a block-bootstrap CI) are meant to be read together.

The pupil block (`--no-pupil` skips it) needs `pupil_diam_px_mean` in the
behaviour csv, i.e. the DeepLabCut eye file `<run>_eye_*DLC_*filtered.csv` in
the run folder. If it was not there at unpack time, `unpack_bin.py <mouse>
--behavior-only` (or the GUI button on tab 1) rewrites the behaviour csv from
it; `aqua_analyze.py` also does this by itself when it finds the DLC file.
Positive lag = activity follows the pupil; the shuffle p comes from circular
shifts of the pupil trace (autocorrelated traces make a plain p far too
small); the partial r removes running, which the pupil follows (r ≈ 0.45 on
the real recording); the tercile rates compare pre with post at the same
pupil size, which is the test of whether a post change is only arousal.

**Curation is honoured automatically.** CFUs marked noise are dropped from the
CFU tables and figures, and their events are dropped from the event statistics
too (`--event-source drop-rejected-cfu`, the default). `--event-source
kept-cfu-only` keeps only events that belong to a kept CFU — stricter, but note
that AQuA2 only forms a CFU from at least `minNumEvt` events, so isolated real
events belong to no CFU and would be lost. `--event-source all` ignores the
curation for events (it still applies to the CFU tables), and `--no-curation`
ignores the sidecars entirely. `run_summary.csv` records `curated`,
`n_cfu_rejected`, `n_cfu_added` and `n_events_dropped_by_curation` per run, so
what the curation did to the numbers stays visible. Hand-added regions appear in
`cfu_matched_<experiment>.csv` as `manual1`, `manual2`, … with `source =
manual`; a region drawn in the same order in several runs of a FOV is treated as
the same region across them.

Event rates use the **analysed** frame (AQuA2's border crop removed), not the
full tiff, so `events_per_min_per_mm2` is comparable across settings.

Conditions come from the planner: `pre`/`post` for RearChr mice, `d260911`,
`d260917`, … (earliest = baseline) for the spontaneous ones; `--baseline` picks
another reference.

**CFU matching across runs.** AQuA2 finds CFUs per recording, so a CFU has no
identity across runs. `aqua_analyze.py` aligns each run's mean image to the
first run's (rotation ± 4° + translation, the same routine as the suite2p
route's `prealign`), warps the CFU maps, and matches them greedily by
intersection-over-union (`--cfu-iou`, default 0.3) — the Python equivalent of
AQuA2's `tools/multi_session_FOV_registration.m` + `multi_session_CFU_comparison.m`,
done for all runs of a FOV at once. An alignment is only applied when it clearly
improves the match of the two mean images, otherwise the maps are matched
unshifted and the log says so. CFUs present in *every* run give a paired
pre→post comparison (Wilcoxon over CFUs, exploratory) — that is the main thing
CFU matching buys you over pooled event statistics.

**What to compare, and what not to.** Event rate is events per minute of the
time actually used (all frames, or still frames with `--gate still`).
Amplitudes are medians over events *and* median-of-run-medians, because the
event-level p values in `summary.txt` treat events as independent and are
descriptive only — the statistic that carries is one value per animal (or per
FOV, declared as such). Keep the settings identical within a FOV (which is the
point of the preset mechanism) and report them; event counts move with
`thrARScl`, `minDur` and `minSize`. And because mCitrine bleaches within a
session, prefer event rate and peak dF over peak dF/F for pre-vs-post amplitude.

## Use — AQuA route with suite2p registration first: `aqua_prep.py`

```
python aqua_prep.py
```

The quickest way to get AQuA2 input. Click **Add folder** and pick either a
*session* folder (one recording, e.g.
`D:\2photon\Chao\NavTIGRETreg-2608-1\NavTIGRETreg-2608-1_260911_001`) or a
*mouse* folder (`...\NavTIGRETreg-2608-1`) — in which case every session folder
inside it is added at once. Add as many as you like, then **Run**. Each session
is unpacked, registered (no suite2p ROI detection — not needed for AQuA, saves
most of the time), binned to the AQuA rate, cropped, and written as uint16.
Every AQuA tiff (+ its `.json` with the frame rate) is also copied into one
collection folder (default `<mouse folder>\AQuA_input`) so AQuA2's batch mode can
be pointed at a single directory.

Options: AQuA rate (default 1 Hz), crop L,R,T,B (default 80,80,40,40 = MATLAB),
alignment temporal smoothing (default 2 frames — needed for dim axon data),
*rigid only* (closest to the MATLAB `method='rigid'`), GPU, and **register all
together** — tick this when the selected sessions are the *same field of view*
(pre / post optostim): they are then registered to one common reference, the
first labeled pre and the rest post. Leave it unticked for independent
recordings (each session gets its own registration).

## Use — suite2p route by paradigm: `plan.py` / GUI "Add mouse folder (auto-plan)"

Recordings are named `<mouse>_<yymmdd>_<nnn>.sbx`. The planner groups them into
FOV-level experiments (nothing is moved; one output folder per experiment under
`<mouse>\afferent2p_out\<experiment>\`):

| paradigm (from the animal name) | grouping | conditions |
|---|---|---|
| **spontaneous** (no `RearChr`): n FOVs per day, several days | one experiment per FOV = run `00k` of every day | the dates (`d260911`, `d260917`, …); earliest = baseline |
| **optostim** (`RearChr`): FOV1, FOV2, FOV1, FOV2 → 1 h stim → FOV1, FOV2, FOV1, FOV2 | one experiment per (day, FOV): first half of the day's runs = pre, second half = post, FOV = position in the block mod n_fov (001,003,005,007 → FOV1) | pre / post |

```
python plan.py "D:\2photon\Chao\NavTIGRETreg-2608-1"                 # show the plan
python plan.py "D:\2photon\Chao\NavTIGRETreg-2608-1" --run           # write configs + run all
python plan.py <mouse> --paradigm optostim --n-fov 2 --stim-after-run 4 --run
```

In the GUI, **Add mouse folder (auto-plan by paradigm)** does the same and lists
each planned experiment; double-click a run to exclude it (skip) or restore it.
All runs of an experiment are registered together to one reference (see
"Tracking afferents across sessions" below), so the same afferent has the same
ROI in every condition; the analysis compares every condition against the
baseline (`analysis.baseline_condition` to override).

## Tracking afferents across sessions (pre/post, different days)

The pipeline does not "match" ROIs between sessions — it makes matching
unnecessary: all runs of an experiment are aligned into one coordinate frame
and ROIs are detected once on the combined movie, so an afferent's ROI *is* the
same object in every condition. Two stages make that work when the FOV was
re-found by hand:

1. **Pre-alignment** (`prealign:`): the structural image of each run (red
   vessels when every run has them, else the green 10th-percentile image =
   mCitrine + resting GCaMP) is aligned onto the first run by a search over
   rotation (±6°) and translation. The transform is applied to that run's frames
   before suite2p (which handles translation and small deformation but not
   rotation). `qc/prealign.png` shows before/after overlays (magenta =
   reference, green = run; white = overlap) and the log/`prealign.json` report
   shift, rotation and the structural correlation achieved.
2. **suite2p** then registers residual motion frame by frame (on the red
   channel when available) to a single reference.

Check `qc/qc_report.txt`: the mean-image correlation of each run vs the first
run is the same-field score (≥ 0.8 OK). The analysis summary also counts
afferents with events in only one condition — if many, alignment probably left
session-specific ROIs. ROICaT-style matching of separately detected ROIs is not
used: its ROI network is trained on somata, and axon fragmentation differs
between sessions, so footprint matching is unreliable for afferents; joint
detection on aligned data is the approach used in the published afferent work.

## Use — GUI (full pipeline incl. suite2p ROIs)

```
python gui.py
```

Click **Add experiment folder…** and pick a session folder — the one that contains the
.sbx runs, e.g. `...\GCAMP_afferents\NavRearChrTIGRE-2608-1\NavRearChrTIGRE-2608-1_260831_001`.
Add as many folders as you want (each folder = one experiment); they queue up and run
in order. Runs inside each folder are auto-discovered and labeled pre/post
(first half pre, second half post when there are more than two) — **double-click a run
to flip its label**. Adjust the settings row if needed, hit **Run all**, and watch the
log. The GUI writes an `afferent2p_config.yaml` into each folder, so any experiment can
be re-run later from the command line without the GUI.

## Use — command line

1. Copy `config_example.yaml` → e.g. `mouse1_0616.yaml`, edit the `experiment:` block
   (folder containing the .sbx/.mat/.quad files, and which run is pre / post).
2. Run:

```
python run_pipeline.py mouse1_0616.yaml
```

That runs all steps: `convert, register, export, loco, qc, analyze`.
Rerun any subset with `--steps`, e.g. after curating ROIs in the Suite2p GUI:

```
python run_pipeline.py mouse1_0616.yaml --steps analyze
```

## Outputs (`<data_dir>/afferent2p_out/`)

| Path | What it is |
|---|---|
| `manifest.json`, `config_used.yaml` | record of exactly what was run |
| `raw_tifs/*_raw.tif` | unpacked raw movie per run (viewable in Fiji; replaces UnpackSBX output) |
| `suite2p/plane0/` | registration + ROIs, traces (`F.npy`, `Fneu.npy`, `stat.npy`, `iscell.npy`) and **`Fall.mat`** (works with the lab's existing MATLAB post-processing, e.g. `GCAMP_concat_Fall`) |
| `aqua_input/*_reg.tif` | **registered** movie per run — the input for AQuA2, already in a common reference frame. uint16 by default (MATLAB-safe); optional temporal binning / edge crop via `export.aqua_bin_t` / `export.aqua_crop` |
| `qc/qc_projections.png` | mean + structural (10th-percentile) image per run: check the *same* axons are visible pre & post |
| `qc/qc_report.txt` | pre-vs-post mean-image correlation — a z-drift score with OK/CHECK/BAD flags |
| `qc/qc_shifts.png`, `qc/qc_locomotion.png` | registration shifts; velocity with still bouts shaded |
| `suite2p/plane0/runs.json`, `frame_labels.csv`, `afferent2p_runs.mat` | which frames of the concatenated suite2p output belong to which run / condition (see *Joint registration* below) |
| `analysis/roi_curation.png/.csv/.txt` | the afferent filter: features, thresholds, example traces, pass/fail per ROI |
| `analysis/prepost_per_roi.csv` | events/min, mean peak dF/F, and baseline F0 per ROI, pre and post |
| `analysis/prepost_summary.txt`, `analysis/prepost_overview.png` | FOV-level summary + figures |

## Recommended workflow per experiment

1. Run the pipeline (all steps).
2. Open `qc/qc_report.txt` + `qc_projections.png`. If the pre-vs-post correlation is
   flagged BAD, the imaging plane likely drifted in z during the stim hour — treat the
   per-afferent comparison for that session with caution or exclude it.
3. Check the afferent filter: `analysis/roi_curation.png` (kept vs rejected traces).
   Tune `analysis.curate_*` + `--steps curate` if needed, or open `python -m suite2p`
   → load `afferent2p_out/suite2p` and use the probability slider (0.5 = the filter)
   / flip individual ROIs.
4. Rerun `--steps analyze`.
5. If using AQuA2: run it in MATLAB on `aqua_input/*_reg.tif`, with the SAME AQuA
   parameters for pre, post, and every animal.
6. Group statistics: one value per animal/FOV, compare the pre→post change between
   ChR2 and control mice, within matched locomotion state (`loco_gate: still`).

## Opening the tiffs / feeding AQuA2

If ImageJ says *format not supported*, the usual cause is **signed 16-bit** data:
ImageJ reads 8-bit, 16-bit *unsigned* and 32-bit float only, and suite2p's
registered binary is signed int16 (SampleFormat=2). The exports are therefore
**uint16** by default (`export.aqua_dtype`) and open normally. Size alone is not
the problem — ImageJ handles large/BigTIFF stacks fine.

**Binning for AQuA (lab practice: 1 Hz).** The AQuA tiff is temporally binned
exactly the way the MATLAB pipeline does it (`GenerateExptProjections` →
`WriteSbxPlaneTif`): bin factor = `round(fs / rate)`, each output frame is the
**mean** of that many consecutive registered frames, written as uint16 with the
edges cropped by `[80,80,40,40]` (MATLAB `projParam.edge`). Defaults:
`aqua_rate_hz: 1.0` → at 15.62 Hz that is 16 frames per output frame (0.976 Hz).
Change the rate in the config or the GUI's *AQuA rate (Hz)* box; use the *same*
value for every recording you compare. Each tiff gets a `.json` sidecar with the
resulting frame rate, bin factor, and crop for entering into AQuA2.

## How this maps onto the MATLAB pipeline (and how to check)

| MATLAB (`ProcessScanboxData_concatFirst`) | afferent2p |
|---|---|
| `ConcatenateRuns` – join runs, one coordinate frame | all runs registered together to one reference (suite2p) |
| `CorrectExpt` – per-15-frame-chunk DFT rigid pass (2× spatial bin, σ=1 blur, edge crop) | suite2p rigid phase-correlation (same math as `dftregistration`) with spatial smoothing + temporal smoothing for alignment |
| `AlignPlanes`/`RegisterSBX` – per-frame `imregtform` **rigid** to a still-period reference | suite2p reference from high-correlation frames; nonrigid block refinement (set `nonrigid: false` for rigid-only, closest to MATLAB `method='rigid'`) |
| `GenerateExptProjections` – bin to `rate_target` (1 Hz), crop `[80,80,40,40]`, uint16 tif | `export` step: identical rule (`aqua_rate_hz`, `aqua_crop`, uint16) |

The registrations are not bit-identical (different implementations), but they do
the same job. To verify on your own data run:

```
python run_pipeline.py <config> --steps compare --ref "<MATLAB projection or registered tif>"
```

`--ref` can be the mean projection UnpackSBX/WriteSbxProjection wrote, or a
MATLAB-registered movie tif (its mean is used). You get `qc/compare_matlab.png`
(side by side + overlay) and `qc/compare_matlab.txt` with the correlation after
removing any global offset, and a sharpness score for each — the registered
Python mean should be at least as sharp as MATLAB's raw projection, and r ≳ 0.9
means both pipelines saw the same field the same way.

## Joint registration of pre and post (and of days) — why, and how it is labelled

suite2p runs on **all runs of an experiment concatenated** (pre + post of one FOV;
or the same FOV across days). That is deliberate: one ROI set, detected on the
whole recording, gives every afferent one identity across conditions, so a
pre-vs-post comparison is a paired comparison of the same pixels — no ROI
matching step, and an afferent silent pre and active post is still found (its
ROI comes from the post frames). The cost is that it assumes the runs show the
same plane. The pipeline guards that in three places:

1. **`prealign`** — before conversion, each run's structural image (red mean, or
   green 10th percentile) is aligned to the first run's by rotation (±6°) and
   shift, which absorbs the offset from lifting/re-lowering the objective for the
   stim hour; `qc/prealign.png` / `prealign.json` show what was applied.
2. **`qc/qc_report.txt`** — correlation of each run's registered mean image with
   the first run's: OK ≥ 0.8, CHECK, BAD. A BAD post run means a different plane
   (z-drift): the joint ROIs are then not the same afferents, and that session
   should be analysed per run or excluded.
3. **Per-frame labels**, so nothing downstream has to guess where pre ends:
   `suite2p/plane0/runs.json` (runs, conditions, start/stop frame, sbx path,
   prealign transform), `frame_labels.csv` (frame → run, condition, time in run)
   and `afferent2p_runs.mat` (the same for MATLAB, 1-based, plus a
   `frame_condition` vector) — next to `Fall.mat`; `analysis/afferent_traces.npz`
   also carries `bounds`, `run_conditions` and `frame_condition`; suite2p's own
   `ops['frames_per_folder']` has the counts. The AQuA exports stay one tiff per
   run, named `<run>_<condition>_reg.tif`.

If you prefer to detect ROIs per run (e.g. after a BAD drift flag), make one
experiment per run in the GUI ("Add folder as ONE experiment…" on the run folder)
— then afferents are matched across runs only by position, which is what the joint
run avoids.

## Manual ROIs (recommended for sparse, low-SNR FOVs)

Hand-drawn ROIs along visible fiber segments often beat automatic detection when
the FOV is one dim bundle. Two-phase workflow (needs `pip install roifile`):

1. `python run_pipeline.py <config> --steps manualroi`
   → writes `manual_roi/mean_for_rois.tif` (and `struct_for_rois.tif`).
2. Open `mean_for_rois.tif` in Fiji. Draw an ROI along each fiber segment
   (polygon / freehand / oval — not line tools), press `t` to add each to the
   ROI Manager, then ROI Manager → More → Save → `RoiSet.zip` into the
   `manual_roi` folder. (Alternative: a label-mask tiff named `roi_mask.tif`.)
3. Run the same command again → traces are extracted from the registered movie
   under each ROI (ring-neuropil subtracted), `manual_roi/roi_overlay.png` shows
   the ROIs on the mean image for verification, and the standard event analysis
   is written to `analysis_manual/`.

This runs alongside (not instead of) suite2p's automatic ROIs, so you can compare
`analysis/` and `analysis_manual/` for the same recording.

## Two-channel recordings (green GCaMP + red vascular tracer)

Runs may have 1 or 2 channels, mixed within a mouse folder — the analysis always
uses channel 1 (green). When **every** run of an experiment has the red channel,
`imaging.register_on_red: auto` (default) computes the motion correction on the
bright, static vessels (suite2p `align_by_chan2`) and applies it to the GCaMP
channel — far more reliable than aligning dim axons; the log says which channel
was used. If any run lacks red, registration falls back to green. Force with
`register_on_red: true/false`. suite2p then also writes `data_chan2.bin`,
`F_chan2.npy` (red traces per ROI) and `redcell.npy`.

## Behaviour: locomotion + pupil (both routes)

Both signals live next to the `.sbx` and are found automatically:

| file | content | alignment |
|---|---|---|
| `<run>_quadrature.mat` | `quad_data`, int32 cumulative wheel counts | one sample per 2P frame |
| `<run>_eye_*DLC_*filtered.csv` | DeepLabCut: `pupil_1..8` = edge points (45° apart), `pupil_9` = centre | one eye frame per 2P frame (camera on the frame clock) |

Pupil diameter = 2 × median distance from the centre point to the edge points
with likelihood ≥ `pupil.likelihood_min` — the model's design: occluded or
mistracked edge points just drop out (≥ `min_edge_points` needed; otherwise the
frame is a blink → NaN, short gaps interpolated). An ellipse fit is stored as a
supplementary measure. Validated on a paired recording: pupil vs running r ≈ 0.65,
pupil lagging running by ~0.8 s; `qc/qc_behavior.txt` prints this check per run
(a low r flags an eye/quad misalignment).

Where it goes:
- **suite2p route**: `behavior/<run>_<cond>_behavior.csv` (frame, time, velocity,
  still, pupil_diam_px, pupil_valid); `qc_locomotion.png` shows pupil over velocity;
  `prepost_per_afferent.csv` gains `r_pupil`, `r_loco` (coupling of each afferent's
  dF/F to pupil size / running speed, per condition) and `pupil_median`; the summary
  reports pre vs post pupil and warns when arousal differs by >15%.
  `analysis.pupil_gate: match` compares pre/post only in the overlapping pupil-size
  range (matched arousal), on top of the still-gating.
- **AQuA route** (`unpack_bin.py`): `<run>_bin31f_1p98s_behavior.csv` — velocity,
  still fraction, pupil diameter and valid fraction averaged over the *same*
  bins, so row i = recorded frame i (`afferent2p/behavior.py`; `--behavior-only`
  rewrites it for existing movies). `aqua_analyze.py` uses it for the
  locomotion/pupil state at each event onset, still-gating, and the pupil
  block (`afferent2p/pupil_coupling.py`): pupil–activity lag correlation with
  shuffle p and partial r vs running, per-CFU coupling, and pre vs post event
  rate at matched pupil size (terciles pooled per FOV).

## Per-afferent analysis (suite2p route) — what the `analyze` step does

AQuA2 reports events and CFUs; suite2p reports ROIs. Neither is an *afferent*: a
single axon is usually split into several ROIs (fragments), and a silent axon
produces no ROI at all. The analysis step turns ROIs into afferents and compares
each one pre vs post:

1. **dF/F** — neuropil-corrected F, moving 10th-percentile baseline, MAD z-score.
2. **Curation (the afferent filter)** — suite2p's built-in classifier is trained
   on somata and is useless for axons, and on a real afferent FOV detection
   returns thousands of ROIs, most of them background speckle whose "signal" is
   photon noise. ROIs are therefore classified by their **response profile**, the
   way AQuA decides an event is real (noise-scaled significance + duration):

   | criterion | default | what it rejects |
   |---|---|---|
   | `curate_min_autocorr`: lag-1 autocorrelation of the high-passed F − 0.7·Fneu (the same quantity as suite2p's stat `snr`) — the fraction of variance that is temporally smooth | `auto` = max(0.05, median + 4 robust SD over all ROIs) | photon-noise ROIs. On a real afferent FOV (4015 ROIs) the junk bulk sits at 0.00–0.06, the auto threshold lands at ≈ 0.07, and the afferents at 0.1–0.7; `roi_curation.png` shows the histogram with the cut |
   | `curate_min_sig_events` significant transients: 0.5-s smoothed z > `curate_sig_z_on` (3) for ≥ `curate_sig_min_dur_s` (0.5 s), peak z ≥ `curate_sig_z_peak` (4) | 1 | smooth-but-eventless traces |
   | `curate_max_red_ratio`: ROI intensity in the red (tracer) mean image / surrounding annulus | 2.0 | ROIs sitting on a vessel (vessel walls make fake transients). Needs channel 2; suite2p's own `redcell.npy` is ignored when it flags > 50 % of ROIs (it flagged 100 % on real data) |
   | `curate_min_npix` | 6 | single-pixel speckle |
   | `curate_min_asymmetry`, `curate_skew_min`, `curate_snr_min`, `curate_min_aspect`, `curate_max_neuropil_r` | 0 = off | optional; all values are still reported. On real data the asymmetry test (fast rise / slow decay) rejected the afferents themselves — their responses are plateau-like — and skewness cannot separate anything because Poisson noise is skewed too (4001/4015 junk ROIs had skew > 0.5) |

   Thresholds live in `afferent2p/curate.py` (`DEFAULTS`); a config lists a
   `curate_*` key only to override it, so pipeline updates take effect. Configs
   written by v17/v18 carry old values — run with `--default-curate` once, or
   delete the `curate_*` lines. Any value can also be overridden on the command
   line: `--set analysis.curate_min_autocorr=0.1`.

   Every criterion gives a soft score (0.5 at its threshold); the ROI score is the
   minimum, and it is written into `iscell.npy[:, 1]`, so in the **suite2p GUI the
   probability slider at 0.5 reproduces this selection and higher = stricter**.
   `analysis/roi_curation.png` shows the feature scatter with thresholds and example
   traces (kept, borderline-rejected, clearly rejected); `roi_curation.csv` has every
   feature and pass/fail per ROI; `roi_curation.txt` the counts. On simulated traces
   (Poisson noise, 1/f background, motion steps, jGCaMP8s transients of various
   rates/amplitudes) the defaults reject ≥ 93 % of noise ROIs and keep ≥ 90 % of
   axons except very dim / very sparse ones, which cannot be told from noise at
   the single-ROI level (they can still merge into afferents from stronger
   fragments). To tune on real data run only the filter — it takes seconds:

   ```
   python run_pipeline.py <config> --steps curate      # rewrites iscell.npy + roi_curation.*
   python run_pipeline.py <config> --steps analyze
   ```

   `python debug_dump.py <experiment output folder>` writes `curation_debug.npz`
   (~40 MB: dF/F of the 300 smoothest ROIs + kept + 100 random, features, outlines,
   mean image) — a shareable file for tuning the filter offline.

   If you then curate in the suite2p GUI (`python -m suite2p`, load
   `afferent2p_out/suite2p`), the change is detected on the next `--steps analyze`
   and your curation wins (`iscell_auto.npy` keeps the automatic verdict,
   `iscell_classifier.npy` suite2p's original one). After curating a few FOVs you can
   also save a suite2p classifier trained on your own afferents and point
   `roi_detection.classifier_path` at it.
3. **ROIs → afferents** — hierarchical clustering by activity correlation of the
   smoothed dF/F traces during still frames (`cluster_r: auto` picks the threshold
   from a shuffle null; `roi_correlation.png` shows observed vs null pairwise r,
   so you can see whether fragments separate from chance at your SNR).
   `afferent_map.png` shows the result on the mean image.
4. **Per-afferent trace** = npix-weighted average of member ROIs (higher SNR).
5. **Per condition, still-gated**: events/min, mean peak dF/F, mean peak ΔF
   (absolute, robust to F0 changes), fraction of time active, deconvolved
   `spks`/min, baseline F0 → `prepost_per_afferent.csv`, paired plot
   `prepost_overview.png`, and `prepost_summary.txt` with a within-FOV Wilcoxon
   (exploratory; real statistics are one value per animal, ChR2 vs control).

Detection settings are axon-oriented (`spatial_scale: 1`, `threshold_scaling: 0.8`,
`soma_crop: false`, PCA `denoise`, ~1 s temporal binning *for detection only*).

**Limits to know**: (a) only *active* axons become ROIs — for a "fraction of
afferents that changed" denominator, add manual ROIs on the structural image
(`--steps manualroi`; they go through the same analysis); (b) fragment merging is
SNR-limited — when observed and null correlations overlap, afferents will be
over-counted (fragments counted separately), which inflates n but does not bias
the pre/post direction; (c) mCitrine bleaching lowers F0 post — check
`F0_ratio_post_pre` and prefer events/min, spks/min and peak ΔF over peak dF/F.

## Event definition (analysis step)

Follows the published meningeal-afferent approach (Blaeser et al., eLife 2024):
neuropil-corrected F, moving 10th-percentile baseline (32 s window) → dF/F → z-score;
an event = z > 1 for ≥ 1 s with peak dF/F ≥ 5%. All thresholds are in the config.
Event *frequency* is the primary metric — with mCitrine (ReaChR tag) sharing the green
channel, dF/F amplitudes are compressed and amplitude increases can be mimicked by
bleaching of the static background; the per-ROI `F0` columns in the CSV are there to
check for exactly that.

## Things to verify on first real use

- **.quad format**: the loco loader handles a .mat with a counts vector (cumulative
  encoder position by default; set `quad_kind: rate` if the file already stores speed)
  and falls back to raw binary. Verify `qc_locomotion.png` against what you saw the
  animal do; if it looks wrong, the format needs pinning to one real example file.
- **counts_per_cm**: set the wheel calibration to get cm/s; otherwise thresholds are in
  counts/s.
- Single plane, single functional channel for now. Optotune volumes: convert.py is the
  place to extend (plane loop), suite2p handles multiplane natively.

## What this replaces / keeps from the MATLAB pipeline

- No Excel sheet, no fixed folder-naming convention — point `data_dir` at wherever the
  files are; runs are discovered or listed explicitly.
- Unpack settings (`binT` etc.) are gone: raw tiffs are always full-rate; nothing
  downstream depended on those settings anyway (the MATLAB registration read raw .sbx).
- Registration is Suite2p rigid + nonrigid to one reference across all runs — same
  intent as `ProcessScanboxData_concatFirst`, standard implementation.
- AQuA2 stays available via `aqua_input/`; Suite2p outputs feed the ROI-based route,
  and `Fall.mat` keeps MATLAB post-processing usable.
