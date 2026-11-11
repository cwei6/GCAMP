# afferent2p — how to run the AQuA2 pipeline

(Offline copy of the guide; the same text lives as a Claude doc you can edit.)

## Before you start (once)

Everything runs from an Anaconda prompt in the pipeline folder, plus two MATLAB
scripts. Open the prompt and go there first, every session:

```
cd D:\2photon\Chao\GCAMP_affrerents\afferent2p
```

**On a different PC** (data reached as `Z:\...` instead of `D:\...`): `REMOTE_PC.md`
— copy the folder locally, `setup_remote.bat`, and rewrite the animal's batch
paths once with `remap_batch.py --from D:\2photon --to Z:\2photon --aqua-dir <AQuA2 there>`.

**One window for all steps:** `python pipeline_gui.py` opens a window with one
tab per step below, in the order they are used — *1 Unpack · 2 Clean & check ·
3 Detection settings · 4 Batch · 5 MATLAB · 6 Curate · 7 Analyze · 8 Report* —
a table of every recording with its flagged stretches, check verdicts and your
decisions (double-click a row to open its frame stack in ImageJ), and a log.
The strip under it shows **where the animal is** — one chip per step, read from the
files on disk (green done, amber pending or partial, red stale, grey not yet), the
next step framed, and the next action written under it ("next: 5 MATLAB: AQuA2 6/6,
CFU 0/6 — paste run_aqua_and_cfu.m"); it refreshes after every step, *Refresh* reads
the folder again, *Details → log* prints every step, and clicking a chip opens its
tab. The same text comes from `python pipeline_status.py <mouse folder>`. So the
settings → batch → MATLAB loop can be dropped and picked up again without
remembering what was done. Under the chips a progress row says what runs
and how far it is (item i of n, elapsed and remaining time, from the steps'
own `[i/n]` lines). The tabs and the log share the window with a draggable
divider between them (v36: the log never disappears on a small screen; drag
the divider to give a tall tab or the log more room). **Stop** (next to the
log, and in the *Run all* row) ends the running step at its next log line —
a recording being unpacked or analysed finishes first, a MATLAB started from
the window is closed, a MATLAB script pasted into a MATLAB window keeps
running there (close it yourself); the log then says `stopped: by the Stop
button` and the chips are refreshed (since v36 — before, Stop only closed a
MATLAB started from the window).

**▶ Run all (guided)** above the tabs runs 4 → 3 → 4 → 5 → 7 → 8 in one go: at
each MATLAB step it shows the command to paste, copies it to the clipboard,
checks the outputs every 30 s and continues by itself (tick *start MATLAB from
here* to have it run `matlab -batch` instead). It calls exactly the same
functions as the commands in this guide, so both can be mixed.

1. Unzip the newest `afferent2p_vNN.zip` over that folder. Nothing in your data
   folders is touched by the unzip.
2. Install the Python packages once:
   `pip install tifffile scipy matplotlib pandas h5py sbxreader`.
3. Update AQuA2 from GitHub (yu-lab-vt/AQuA2). A copy older than 2025-12 has the
   `curSz` save bug; if you cannot update, add `curSz = []; dsSclXY = [];` above
   the line `if opts.sz(3)>1` in `src\+ui\+proj\saveExp.m`.
4. After EVERY AQuA2 download, in MATLAB: `patch_aqua2_toolbar('D:\AQuA2')`
   (script in this folder). Otherwise `aqua_gui` stops at start with
   `Property assignment is not allowed when the object is empty` (addCon_dat
   line 123): this MATLAB creates the axes toolbar lazily, so AQuA2's
   `Toolbar.Visible = 'off'` lines fail; the patch turns them into
   `Toolbar = []`, keeping a `.bak`.

Conventions the tools rely on:

| What | Convention |
|---|---|
| Animal folder | `D:\2photon\Chao\GCAMP_affrerents\<mouse>\` |
| Run folder | `<mouse>_<yymmdd>_<nnn>\` with the `.sbx`, `.mat`, `_quadrature.mat`, DLC eye csv |
| Paradigm | no `RearChr` in the name: FOV k = run 00k every day, conditions = days. `RearChr`: runs 1–4 pre (FOV1, FOV2, FOV1, FOV2), 5–8 post; one experiment per day and FOV |
| AQuA2 settings per FOV | `FOV1_Opt.csv`, `FOV2_Opt.csv`, … in the animal folder (any name containing `FOV1` works) |
| Binned tiff | `<run>_bin31f_1p98s.tif` next to the `.sbx` (older `..._1.98s.tif` also accepted) |
| AQuA2 results | next to the tiff: `_AQuA2.mat`, `_AQuA2_Ch1.csv`, `_AQuA2_res_cfu.mat` |
| Your decisions | `<run>_..._AQuA2_res_cfu_curated.json` next to the results |
| Summaries | `<mouse>\aqua_batch\` (batch files), `<mouse>\aqua_analysis\` (results) |

Folders the tools create (`aqua_batch`, `aqua_analysis`, `AQuA_input`,
`afferent2p_out`) are never scanned for recordings; a run present twice (e.g. a
copy in `test_run`) is used once and reported.

## Step 1 — Unpack

```
python unpack_bin.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --sec 2
```

Skips runs that already have a tiff (`--overwrite` redoes them); no arguments =
GUI. Per `.sbx` it writes the 1.984 s/frame uint16 tiff (whole-field
brightness flattened; every frame kept), `_ch2.tif` when red was recorded,
`_meanproj.tif`, `.json`, `_behavior.csv` (speed, still fraction, pupil per
frame), and the brightness QC: `_frames.csv`, `_globalfix.png`,
`_flagged_frames.png/.csv`, `_check_frames.tif`, `_check_movie.tif`. `--keep-raw` also keeps the uncorrected movie
(`..._raw.tif`, twice the disk; lets you undo a drop later), `--no-global-fix`
skips the step.

The behaviour csv takes running from `<run>_quadrature.mat` and the pupil from
the DeepLabCut eye file `<run>_eye_*DLC_*filtered.csv` (9 points: `pupil_1..8`
on the edge, `pupil_9` the centre; diameter = 2 × median centre-to-edge
distance over the points with likelihood ≥ 0.6, so a partly occluded outline
still gives a radius). Both files must sit in the run folder next to the
`.sbx`. If the DLC file arrives after unpacking, nothing needs re-unpacking:

```
python unpack_bin.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --behavior-only
```

(GUI: *Rebuild behaviour csv (running + pupil) for existing movies*) writes
only the csv again — and Step 7 does the same by itself when it finds a DLC
file next to a movie whose csv has no pupil columns.

## Step 2 — Clean & check: the flagged frames and your decisions

Nothing is removed from the movie. The step flags two kinds of stretches and
then checks every frame of them against the kept frames around it:

| class | what | default |
|---|---|---|
| **artefact** (red) | fast changes: a flash, a jump between two frames beyond 3 %, a spike shorter than 10 s, a content change (shutter, brief shift) | its **flagged frames** — the frames that differ from the neighbours by ≥ 10 % of the frame and ≥ 3× the clean floor (the flash itself) — are **excluded**: events starting in them (±1 frame) leave the statistics, their time leaves the denominator. The tails and margins around them, and stretches without any such frame, are **kept** |
| **state** (amber) | slow changes: a level plateau, or the background/fluorescence ratio (p10/median) more than 5 % away from its baseline — a laser or focus change, or widespread activity: the trace cannot tell | **kept** (counted); the verdict says how different the frames are, and the stretch note warns when its background ends up brighter than normal after the correction (check that AQuA2 does not turn it into events) |

Check them frame by frame in ImageJ, like you would scroll through the raw
movie — but only the flagged parts, labelled. Per recording the step writes:

| file | what |
|---|---|
| `_check_frames.tif` | 16-bit ImageJ stack: every flagged stretch with 6 frames of context on each side, at recorded brightness and full resolution. The slice label (ImageJ title bar) and a text line burned into the top of each frame say: recorded frame number, stretch and class, whole-field change, difference from the clean neighbours (`x` = relative to how much two clean frames differ), **FLAGGED FRAME** or not, **EXCLUDED** or kept |
| `_check_movie.tif` | the same slices as a composite: the frame (one grey scale for the whole stack, so a flash shows as brighter), its smoothed % difference from the neighbours (red brighter, blue darker), and the level/shape traces of that segment with a cursor on the current frame — the stretch shaded, flagged frames ticked in red, excluded frames barred in black |
| `_flagged_frames.png` / `.csv` | the one-glance overview sheet per stretch, and the per-frame table (field %, diff rms %, floor, ratio, flagged, decision, excluded) |

GUI tab 2: select the recording, then the row of buttons directly under the
recordings table — *Show check sheet* (the `_flagged_frames.png` pictures of
every flagged stretch with the per-frame table, inside the GUI; *Zoom / pan*
opens it in the picture viewer together with the QC plot), *Frame stack
(ImageJ)* (set the ImageJ / Fiji executable once in the top row; a
double-click on the recording does the same), *Composite stack*, *Make check
stacks*, *QC plot* (`_globalfix.png`: the brightness trace with the flagged
stretches), *Run folder*. The lower table lists every stretch with its
flagged frames, verdict, decision and excluded frames; the row under it
decides per stretch: *Exclude flagged frames* (the default for artefacts),
*Exclude whole stretch*, *Keep all*, *Back to default*, or type frame numbers
(as in the stack labels, e.g. `68,76-77`) and *Exclude these frames*.

> v45: these two button rows had become invisible — they were the last
> things in the tab, and once the window gained the *Run all* row, the status
> strip and the log pane (v35/v36) the tab no longer had room for them, so Tk
> dropped them without a word. They now sit directly under their tables and
> are laid out first, so a short window shrinks the tables instead; the
> window also opens as tall as the screen allows. If a row is ever missing
> again, drag the divider above the log down.

Or on the command line:

```
python aqua_clean.py "<run folder or tiff>" --decisions-only --keep 13-19 --exclude 63-87
python aqua_clean.py "<run folder or tiff>" --decisions-only --exclude-all 63-87
python aqua_clean.py "<run folder or tiff>" --decisions-only --exclude-frames "63-87:68,76-77;94-102:99"
```

Decisions live in `<run>_..._frames_decisions.json` (`keep`, `exclude` = the
flagged frames only, `exclude-all`, `frames:68,76-77`), take effect at Step 7
(no AQuA2 re-run) and survive re-runs of the step. `_frames.csv` has one row
per frame with class, stretch, decision, `flagged_frame` and `excluded`. The
stack labels show the decision at the time they were written — *Make /
refresh check stacks* (or `aqua_clean.py --stacks-only`) rewrites them.

The rule for a flagged frame is adjustable (`--core-min-pct 10 --core-ratio
3`, or the two boxes in the GUI); on the real animal the defaults pick the
flash frames (+25–45 % whole field) and leave the +5–9 % neighbours in.

Tiffs corrected by an older version get the current rules and the stacks with
`--force` (frames an older version dropped or replaced stay excluded — their
content is gone; `--restore-raw` puts `..._raw.tif` back first, otherwise
re-unpack that run). `--report-only` previews everything without touching the
movie; `--no-stacks` skips the two stacks (they take 20–200 MB per recording,
`--stack-ds 2` halves the 16-bit one). Removing frames from the movie is an
explicit choice: `aqua_clean.py --drop --force` (or the GUI button) leaves the
*excluded* frames out; the csv column `tiff_frame` then maps movie frames to
recorded frames and Step 7 converts event frames back. Any AQuA2 result made
from a movie that was changed afterwards is stale → Step 4 with `--overwrite`.

```
python aqua_clean.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --force
```

## Step 3 (optional, by eye) — Tune a FOV in the AQuA2 GUI

1. `aqua_gui`, load the FOV's first recording (`..._001` for FOV1, `..._002` for FOV2).
2. Project: temporal resolution **1.984** s/frame (see the `.json`), **0.792**
   µm/px, border **5**; `detectGlo` off.
3. Registration as you prefer (kept by the batch); bleach correction can stay off.
4. Act → Temporal → Spatial until the events look right (see the project note
   *AQuA2 GUI setup notes*).
5. Detect tab → **SaveOpts** → save in the animal folder as `FOV1_Opt.csv`
   (or the suggested name with `FOV1` in it). Applied to that FOV on every day.
   Frame time, pixel size and border come from the tiff/batch command, not this file.

No need to save the detection result from the GUI; the batch writes everything.
A per-day override: `<mouse>_<yymmdd>_FOV1_Opt.csv`.

**Since v30 this step is optional.** A FOV that has no settings file of its
own gets `default_Opt.csv` from the pipeline folder when the batch is built
(always, since v31 — tab 4's "settings for a FOV that has no *_Opt.csv" field
and `aqua_batch.py --opts` name another file); that file holds the settings
this animal was analysed with (registration 2, bleach correction 2, median
0.5, smoothing 0.5, minDur 1, sigThr 3.5, source sensitivity 8, thrARScl 2,
minSize 50). The sweep of Step 3 then chooses thrARScl and minSize by
measurement, so for a new animal the route is **Step 1 → ▶ Run all (guided)**
(= Step 4 → Step 3 plan → sweep → apply → Step 4 → Step 5 → Step 7 → Step 8,
with two MATLAB pastes) **→ Step 6 → Steps 7 and 8 again**, for a Treg mouse
exactly as for a RearChr mouse. Open the AQuA2 GUI only to look at a field, or
to change one of the other parameters by hand (save it as `FOV<k>_Opt.csv` as
before — a FOV's
own file always wins over the default; a tuned `aqua_opts\<experiment>_Opt.csv`
wins over both).

## Step 3 — Detection settings by measurement: `aqua_tune.py`

"Change the threshold until activity shows on the afferent" has no criterion,
every change costs a re-run with you present, and different days may need
different settings. `aqua_tune.py` replaces this with three measurements per
setting that run unattended (GUI tab *3*; the FOV's current settings — its
Opt file or `default_Opt.csv` — are the starting point, so the batch of Step
3 is built first when it is missing). What is swept are AQuA2's detection
parameters, by default `thrARScl` × `minSize`, each setting a complete AQuA2
run with the other 45 parameters as in the base file; the fibre mask is not
given to AQuA2 at all — it only places the injected transients and decides
which detections count as on-fibre in the scoring:

* **false alarms** — the setting is run on a *null* copy of the movie: every
  pixel's time series is shifted in time by its own random offset, so each
  pixel keeps its brightness, noise and slow drift while nothing is left that
  is coherent across pixels. Every event found there is false. Counted on the
  fibre mask: null events on the fibres / real events on the fibres (the
  *false-alarm share*), plus null events/min.
* **recall** — the setting is run on an *injected* copy: the real movie plus
  72 synthetic transients on the fibres (4 per class: 0.1 / 0.2 / 0.4 dF/F of
  the pixel's own baseline; 80 / 200 / 400 px footprints grown along the
  fibre; 4 or 12 frames above half amplitude; `--amplitudes 0.3,0.6,1.0
  --sizes 200,400,800` changes the classes, and a changed class set makes the
  injected movie again). The fraction recovered, overall and per class,
  corrected for chance coincidence with real events. Mind the units: AQuA2
  grows an event's footprint 2–3× and reports the footprint's average dF/F,
  so a transient injected at 0.4 appears in its tables at about 0.2 — the
  real events it detects in 2608-1 (reported 0.3–0.46 dF/F, 140–1100 px)
  correspond to injected amplitudes of roughly 0.6–0.9.
* **specificity** — the share of the real movie's events whose centre lies on
  the fibres of the mean image, the event rate on and off the fibres, and how
  the rate changes with the threshold (a plateau is a robust setting, a cliff
  is not). Everything the choice depends on is counted on the fibres, so
  activity elsewhere — autofluorescent cells or debris, vessel walls, dura —
  is shown as the off-fibre rate but neither helps nor hurts a setting; what
  sits on the fibres themselves (an autofluorescent speck on a bundle) is
  indistinguishable from an afferent event here and is left to the CFU
  curation of Step 6.

```
python aqua_tune.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --plan --as-root "Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1"
    MATLAB:  run('Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_tune\run_aqua_tune.m')
python aqua_tune.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --report            # then --report --apply
```

`--plan` takes the first run of every condition of each experiment — for a
RearChr mouse the first **pre and the first post run** of each day and FOV
(since v34; so that the report can say whether pre and post agree, below),
for a Treg mouse (whose experiment is the FOV across days) the first run of
every day — and lists them first; a test day you do not want tuned is left
out with `--days 260911,260917,260921` (`--runs 260921_002`, `--one-per-fov`,
`--baseline-only` = pre only, half the sweep, also narrow it; the GUI's
*movies* box has the same choices). It writes the null and injected movies (~400 MB each, 30 s per movie) and the
truth table into `<mouse>\aqua_tune\<run>\`, the grid as AQuA2 presets
(`parameters_for_tune.csv`: the FOV's own Opt file with the swept values
replaced; default grid `thrARScl` 1.5 / 2 / 2.5 / 3 / 4 × `minSize` 50 / 100 /
200, plus the current setting; `--grid sigThr=3,3.5,4` adds a variable, a
single value such as `--grid minSize=50` fixes one; `--grid maxSize=inf,1500`
caps the footprint — note that a finite maxSize also switches AQuA2's
active-region step to its multi-threshold mode, which splits large regions),
one
batch file per job and `run_aqua_tune.m`. `--as-root` writes the MATLAB paths
for the PC that runs MATLAB (it sees D: as Z:). The MATLAB script runs the jobs
(`tune_one_job.m` does one: AQuA2, export of the events and the mean image,
deletion of the 1.2 GB `.mat`), logs to `tune_log.txt`, and skips finished
jobs, so it can be stopped and restarted. **What a re-plan keeps and what it
redoes (v35 rule):** a job is "done" when its `sNN\<variant>\…_events.csv`
exists; the planner records each movie's settings list in
`aqua_tune\<run>\settings.json` and, on the next plan, keeps every `sNN`
folder whose setting VALUES are unchanged (their jobs are skipped) and moves
the others to `stale_<date>\` to be redone — so after an *Apply* only the
"current" setting (s00) is redone, not the grid; results without a
`settings.json` (a sweep from a version before v30) are all moved aside and
redone; the null movie is reused; the injected movie and its results are
remade only when the fibre mask or the injected classes changed (the recall
must be measured against the same transients). Movies planned earlier and
not selected now stay in the plan with their finished jobs. Nothing of
Steps 4–7 is touched by the sweep: the run's own `_AQuA2.mat` / CFU files
change only when *Apply* writes a newer settings file and Step 4 re-queues
that run (see below). About 1.4 min per job: 6 movies × 16
settings × 3 variants ≈ 7 h one at a time — start it in the evening. AQuA2
itself uses one core most of the time, so the script runs several jobs at
once when it can: with the Parallel Computing Toolbox it opens a pool
(`--workers N` / the GUI field "MATLAB jobs at once"; 0 = chosen from the
cores and the free memory, up to 6 — each job needs roughly 6–8 GB) and
finishes in about 1/N of the time; without the toolbox,
`run_aqua_tune_parallel.bat` starts N MATLAB processes that each take every
Nth job (`run_aqua_tune_<k>ofN.m`; `matlab` must be on the PATH, else open N
MATLAB windows and run one script in each; each writes its own
`tune_log_<k>ofN.txt`). Watch the first minutes once: if the network drive or
the memory is the limit, use fewer.

`--report` scores every setting and prints one table per movie (also
`tune_recommendation.txt`, `tune_report.csv` and `<mouse>_tune.png`: recall,
recall by amplitude and duration, false-alarm share on the fibres, events/min
on and off the fibres with the on-fibre fraction, per movie). **Since v35 it
also draws the sanity check behind the choice, `aqua_tune\<run>_why.png`
per tuned movie, collected in `aqua_tune\tune_why.html` (GUI tab 3: *Why
this setting? (pictures)*)**: one column per setting of the sweep (the
thresholds at the chosen minimum size, plus the current setting), three
rows — the real movie with every detected event drawn as a circle of its
size on the mean image (green on the fibres, orange off them), the null
movie where every detection is a false alarm (red), and the injected movie
with the 72 synthetic transients found (filled cyan) or missed (hollow
magenta) — with the numbers the rule uses under each panel. The chosen
column is framed green; a setting that fails the false-alarm rule is tagged
red ("too permissive: false alarms 40 % > 10 %"), the others say what they
lose ("too strict: recall 0.61 < 0.82"). Going left to right, red
disappears faster than green appears, then cyan turns to magenta: the
chosen column is the last one before the magenta starts. The recommended
setting is the **highest recall among the settings whose false-alarm share on
the fibres is ≤ 10 % (`--max-false`) and whose events are ≥ 80 % on fibres
(`--min-on-fibre`)**; `--target-amp 0.2` ranks by
the recall at one amplitude. If no setting keeps the false alarms that low, the
report shows the best trade-off instead, marked *rule not met* — accept it with
a higher `--max-false` knowingly, or conclude that this field is not reliably
detectable at these settings. For a FOV recorded on several days the report
also shows every day's recall and false-alarm share at each day's chosen
setting, and then **one common setting for the FOV**: the highest mean recall
among the settings that keep the false-alarm share ≤ `--max-false` on *every*
day (the best common trade-off, marked *rule not met*, when none does). A
Treg mouse compares its days with each other, so its days must share one
setting — `--apply` always writes the common choice for an experiment that
spans days (`<mouse>_FOV<k>_days`). A RearChr mouse (one experiment per day
and FOV, pre vs post within the day) gets **each day its own setting** (since
v34 the default: the field is not the same on another day, and on a "new
FOV" day it is a different field altogether; `--common-settings` / GUI tick
*one common setting per FOV for all days* forces one setting per FOV across
the days as before). Within a day, pre and post are both swept, and the
report's block *"pre and post tuned separately — do they agree?"* shows each
one's own choice and the other's recall / false-alarm share at that choice
(**AGREE** = the same optimum; **DIFFER** = the signal changed between pre and
post more than the settings tolerate — SNR, focus, the field after the lens
came back). What the day then gets is decided by the **day rule (v36)**,
`--day-rule matched`, the default (GUI tab 3, row *RearChr day (pre + post
tuned)*): **matched sensitivity — one setting per condition, chosen so that
pre and post are detected with the same measured sensitivity.** Tuning pre
and post apart by their own optima would build the pre-vs-post difference
into the detection (a more permissive post threshold finds more events
whether or not anything changed); one common setting avoids that, but on a
DIFFER day it detects one recording worse than the other (recall 0.82 for pre
and 0.54 for post at the common setting means a third of the real post events
are missing before the comparison starts). The matched rule takes the middle
road with a measurement instead of a guess: among the settings that pass the
false-alarm rule on their own movie, it takes the pair — one per condition —
whose recalls of the injected transients are within `--tol-recall` (0.05) of
each other AND whose false-alarm rates on the null movies are within 10 % of
the lower real event rate (at least 0.05 / min) of each other, at the highest
common sensitivity; when pre and post agree, that pair is simply the same
setting twice. So the two conditions are detected with equal sensitivity and
equal false-alarm rate — a difference in the counts cannot come from the
detector — and Step 7 additionally subtracts each condition's own expected
false alarms. When no pair of clean settings matches (one movie's best clean
recall is far below the other's), the day is marked **NOT COMPARABLE**: the
common setting is written so the batch can run, and Step 7 / the report say
so in the verdict ("DETECTION NOT COMPARABLE") instead of reporting a ratio
as if nothing were wrong. `--day-rule common` (GUI: *one common setting for
the day*) keeps the v34/v35 behaviour.

**Does the sweep reach this field at all? (v37)** The sweep can only choose
among the settings it swept and can only measure recall on the transients it
injected — three things went wrong at once on NavRearChrTIGRE-2609-1 FOV 2
(261014), which a hand-tuned AQuA2 run (thrARScl 2, minSize 20) showed to be
active (144 events, 9 per min, 20–50 px, 0.18–0.41 dF/F) while the sweep's
choice found nothing: (1) the injected classes of v30–v36 (nominal 0.1 / 0.2 /
0.4 dF/F, which AQuA2 reports at ~0.05–0.2 because it grows the footprint and
averages over it) were dimmer than nearly every real event, so EVERY setting
looked blind (recall ≤ 0.3 on every animal swept on Oct 21) and the choice
among blind settings meant nothing; (2) the grid never went below minSize
50 px while that field's events are 20–50 px; (3) the automatic fibre mask
dropped the whole fibre network as a "round speck" (crossing fibres form one
big connected piece that is as wide as it is long), so the transients were
injected on specks and the real fibre events counted as "off the fibres".
v37 changes all three and adds the checks that catch them: the default
classes are nominal **0.3 / 0.6 / 1.0 dF/F** (reported ≈ 0.15 / 0.3 / 0.5 —
below, at and above the real events) and **20 / 60 / 200 px** (a small class
tells whether minSize lets small events through); the default grid is
thrARScl 1.5 / 2 / 2.5 × minSize **20** / 50 / 100 (3 and 4 were never
chosen on any animal); the fibre mask keeps every piece whose skeleton is ≥
10 µm long, whole networks included (a saved mask of the old version is made
again by the next *1 Plan*, and the injected movie with it). The report now
gives, per movie, a **grid verdict**: *GRID TOO STRICT* when the best clean
setting is the most permissive value swept of some variable and still finds
< 50 % of the injected transients (the optimum lies beyond the grid — it says
which values to sweep next, one step further: thrARScl − 0.25/0.5, minSize
÷ 2.5), *edge* for the same with recall 50–70 %, *GRID TOO PERMISSIVE* when
no setting keeps the false alarms low, *WEAK* when more permissive settings
were swept but fail the rules; the suggestions are written to
`aqua_tune\grid_suggestion.json` and **`--plan --extend`** (GUI: *Extend the
grid (as suggested)*) adds them to the grid — finished jobs are kept whatever
their number (a result folder is matched to the new list by its setting
VALUES and renumbered), only the new combinations are swept; *Run all*
does one such round by itself. It also says when the **injected classes
probe below the real signal** (the brightest class, as AQuA2 reports it,
dimmer than the real events — the hand result's when there is one — with the
`--amplitudes` to plan with), when **the fibre mask may be incomplete** (most
events at the chosen setting off it — look at `fibre_mask_raw.png`, draw the
mask yourself if so), and when a setting **finds no event at all**. The
matched-sensitivity rule refuses a pair matched below `--min-recall` (0.30):
"TOO WEAK — the movies can be matched only at a level where most injected
transients are missed" → the day is NOT COMPARABLE and gets the common
setting until the grid or the classes are fixed.

**A result tuned by hand as the reference (v37).** The best yardstick is what
you accept yourself: save a result from the AQuA2 GUI (its default save
folder is a subfolder named after the movie, `<run folder>\<movie>\`), export
its events once in MATLAB — `aqua_export_events('<that folder>\<movie>_AQuA2.mat')`
(tab 3 *Hand result → reference* copies the command; AQuA2 and the afferent2p
folder on the MATLAB path) — and *3 Report* then shows, for that movie, two
more columns per setting: **hand%** = the share of the hand-accepted events
the setting recovers (same onset within 3 frames, same place) and **extra%**
= the share of the setting's events the hand result does not have, a line
"hand-made result (144 events): this setting recovers 85 % of them and adds
12 %; closest to the hand result: setting 7", and the hand events as white
dashed circles in the pictures. The hand result also serves as the real
signal for the classes check. The run must be one of the tuned movies (the
plan takes the first pre and the first post run of a day; a hand result on
another run is listed at the end of the report — put it in the *runs* field,
*1 Plan*, sweep it). The hand result does not decide the setting — the
injected / null measurements do, because the hand result has no measured
false-alarm rate — but a chosen setting that recovers < 50 % of it is flagged.

**Three refinements of the rules (v38), after the first v37 sweep of
2609-1.** (1) *One null event of grace*: the false-alarm share is a ratio of
small counts (4 null events among 36 real = 11 %, 1 among 4 = 25 %), so a
setting also passes when the share without one null event would — a single
null event never decides (with at least 4 real events to lean on; a setting
that finds 1 or 2 events is judged on the raw share). (2) *The on-fibre rule gives way when it would
leave a day unmatchable*: when the ≥ 80 % on-fibre requirement allows only a
match below `--min-recall`, it is relaxed for the movie it limits (the
report says so: "the on-fibre rule was relaxed for 261006_002 (61 % on
fibres)") — a setting that finds the fibre events and some off-fibre ones is
better than one that finds nothing, and the off-fibre events are what the
curation is for. (3) *The same setting for pre and post whenever it matches
within the tolerance*: equal recall of the injected transients makes two
settings equivalent by our measure, and one setting cannot favour either
condition on the real events — on the sham day 261006 FOV 2, pre at thr
1.5 / minSize 50 and post at thr 2 / 20 (both recall ≈ 0.55) gave 1.1 vs
2.3 events/min, a "doubling"; thr 2 / 20 for both gives 2.29 vs 2.29.
Different settings for pre and post are kept for the case the rule was made
for — the same setting cannot match them at all.

**The day as pictures (v36)** — `aqua_tune\<day>_FOV<k>_why_day.png`, at the
top of `tune_why.html`: pre above post at the same settings, one column per
threshold, rows = real / null / injected of each movie. Each movie's own
matched setting is framed green (pre may get one column and post another),
the common setting is framed dashed, and a setting that fails the false-alarm
rule on that movie is tagged red — so it is visible at a glance why pre keeps
one threshold and post another, and what either would find at the other's.

`--report --apply` writes the settings: with the matched rule
`aqua_opts\<experiment>_pre_Opt.csv` and `..._post_Opt.csv` (one file per
condition; when the same setting serves both, the single
`aqua_opts\<experiment>_Opt.csv` as before) plus
`aqua_opts\<experiment>_tuning.json` — the rule, the settings and the
measured recall / false-alarm rate per condition, read by Step 7. Step 4
looks for a run's settings in this order: `aqua_opts\<experiment>_<condition>_Opt.csv`,
`<mouse>\<experiment>_<condition>_Opt.csv`, then the experiment file
`aqua_opts\<experiment>_Opt.csv`, the FOV file, the default; the batch log
says which ("[pre] settings of this condition: …") and the manifest records
the file per run. Rebuild the batch (tab 4; a run whose settings file is
newer than its result is queued again by itself, no `--overwrite` needed) and
re-run the MATLAB scripts of Step 5 to apply them (`--apply-all` also writes
the trade-off settings). `--progress` prints the sweep's done / failed /
pending jobs.

Step 7 then carries the detection into the results: `summary.txt` gets one
line per comparison ("detection (matched sensitivity): pre … (recall 0.82,
false alarms 0.06/min) / post … ; with the expected false alarms subtracted:
1.31 -> 4.56/min = x3.48"), `experiments_summary.csv` the columns
`detection_rule`, `detection_comparable`, `setting_baseline` /
`setting_condition`, `recall_*`, `fp_rate_*` (false alarms per minute of the
whole frame, from the null movie), `sensitivity_gap`, `sensitivity_matched`,
`rate_*_fpcorr` and `ratio_fpcorr` (the two rates with the expected false
alarms subtracted, and their ratio — compare it with `ratio_raw`: when they
differ much, the false alarms were a large share of one condition's events).
The Step 8 report's verdict line ends with "detection (matched sensitivity):
… false-alarm-subtracted ratio ×r", or "DETECTION NOT COMPARABLE" when the
day had no matched pair.

**The tuning at a glance (`<mouse>_tune.png`, v39).** Two numbers per setting
and movie, nothing else: CATCH — 72 fake transients are planted on the fibres
(three brightnesses × three sizes × short/long); catch = the share the
setting finds. NOISE — in a copy of the movie with every pixel's time course
scrambled (no real events left) everything found is invented; noise = those
events as a share of what the setting reports in the real movie. The figure
is one panel per day and FOV (pre and post together; one per movie for a
Treg animal): every setting a dot at (noise, catch) — up = finds more of
the planted events, right = invents more — labelled with its values
("2/20" = thrARScl 2 / minSize 20), sized by the events per minute it
reports on the fibres, hollow when fewer than 80 % of its events lie on the
fibres; the rule's pick is the black star = the highest dot left of the red
line (noise ≤ 10 %); the green rings are the pair the day uses (one setting
for pre and post whenever one matches both); the red dotted ring is the
setting in the batch now. Under each panel, plain sentences: what the day
uses and why, what each movie's pick finds, and the more sensitive setting
the rule rejected with the reason ("35 % of its events would be noise; to
take it anyway: max false-alarm share 0.35"). The same sentences open
`tune_recommendation.txt` (IN PLAIN WORDS, before the tables). The curves
are still there as `<mouse>_tune_detail.png`.

**Pick by eye (tab 3 *Pick by eye…*, v39).** The measurement narrows the
choice to a few settings; the last word can be yours. The window shows, for
a tuned movie, up to four settings worth a look — the rule's pick, the more
sensitive one it rejected, the setting in the batch now, the next cleaner
one, and the day's matched choice — side by side on the real movie: every
event the setting found is a circle (green on the fibres, orange off), the
events the same setting found in the scrambled movie can be overlaid in red
(noise), a hand result as white dashed circles, and the other condition's
numbers at each setting are written under the title. The frame slider and
*Play* run the movie in every panel at once, an active event filled, so a
detection can be seen happening; the toolbar's magnifier and cross zoom and
pan all panels together. *Pick …* writes that setting for the whole day (pre
and post) or for this condition only (`aqua_opts\<exp>_Opt.csv` or
`<exp>_<cond>_Opt.csv`, note "picked by eye") and records the measured
catch and noise of every movie at it in the tuning json, so Step 7 still
reports the sensitivity and the false-alarm rate; tab 4 then queues the
runs again.

**Reading the detail figure (`<mouse>_tune_detail.png`, v38 viewer).** One
row per tuned movie, four panels, thrARScl along the x axis, one line per
minSize. Left: the share of the injected transients the setting recovers
(recall) — rising towards the top-left corner of the grid means "more
permissive finds more". Second: the same split by injected amplitude class
(light = dimmest) and by duration — the dim class is the sensitivity to weak
events. Third: the false-alarm share on the fibres, null-movie events as a
share of the real events — 0 is clean, 1 means the setting invents as many
events in a movie without signal as it reports in the real one; this is what
stops the rule from taking the most permissive setting. Right: the events
per minute the setting reports in the real movie on the fibres (circles) and
off them (triangles), with the on-fibre fraction (squares, right axis). The
red dotted line is the setting in the batch now, the black line and star the
chosen one; the title carries the grid verdict. The viewer (v38) shows the
file at full resolution: Ctrl + wheel or + / − zoom at the cursor, drag to
pan, *Fit* / *1:1*, ◀ ▶ flip through the tuning figure and every "why"
picture, the path is in a box you can select (double-click copies it),
*Copy path* / *Open folder* / *Open file*, and *Redraw at 2× resolution*
draws the figure again at 260 dpi from `tune_report.csv` (seconds;
`python aqua_tune.py <mouse> --figure-only --dpi 260` does the same; the
report's own figure is 150 dpi, `--dpi` sets it).

In the GUI (tab *3 Detection settings*) the same flow is the numbered
buttons: *1 Plan* (movies / days / runs / grid fields above it; builds tab 4
first when the batch is missing), *2 Run the sweep in MATLAB now* (starts
MATLAB without its desktop from the executable in the top row — found by
itself under `C:\Program Files\MATLAB` — and streams its output into the
log and `aqua_tune\run_aqua_tune_log.txt`; *Stop* ends it and a later start
continues; *or copy the command* for a MATLAB window, also listed in tab 5),
*Progress* (done / failed / pending jobs, measured minutes per job and the
hours left, the last lines of the MATLAB log — for the morning after), *3
Report* (the recommendation text appears in the tab, the figure opens in its
own window; tick *apply* to write the Opt files, *one common setting per FOV
for all days* is off by default since v34 — a Treg mouse gets one anyway, a
RearChr mouse per-day settings unless ticked —, *trade-off if the rule is not
met* writes the flagged choice too; the row *RearChr day (pre + post tuned)*
chooses the day rule — *matched sensitivity (one setting per condition)*, the
default since v36, or *one common setting for the day* — and its *recall gap
≤* field is `--tol-recall`) and *4 Apply → next: 4 Batch*, which switches to tab 4
and builds the batch — then tab 5 (one MATLAB paste) and Step 7 as usual.

**▶ Run all (guided)** above the tabs does all of it in one go — tab 4
(settings file or the default) → plan → the sweep in MATLAB → report + apply
(per-day settings for a RearChr mouse, one per FOV for a Treg mouse; when no setting meets the rule the best trade-off,
flagged in the log and the Opt file's note) → tab 4 again → the AQuA2 batch +
the CFU batch in MATLAB → Step 7 → Step 8 (both reports). At the two MATLAB
steps a small window shows the command, copies it to the clipboard and waits:
paste it into any MATLAB window; the outputs are checked every 30 s and the
chain goes on by itself when they are complete (*Continue now* goes on at
once, *Stop* ends it). Tick *start MATLAB from here (matlab -batch)* to have
the window start MATLAB itself (needs the MATLAB field in the top row and the
"mouse folder as the MATLAB PC sees it" field empty). It never stops on the
AQuA2 folder: the generated scripts find AQuA2 on MATLAB's path. It takes a
night for a typical animal (the sweep is the long part) and can be stopped
and started again: finished sweep jobs, results and tuned settings are kept
(untick *skip the sweep if every FOV is tuned* to sweep again). Curation
(Step 6) stays a manual look afterwards; run Steps 7 and 8 again after it.

Where the fibres are is found automatically, no input needed: a ridge filter
on the run's mean projection picks out thin bright line-like structures
(0.8–2.4 µm wide) that stand clearly above the local background, specks under
20 µm² and round pieces whose skeleton is shorter than 10 µm are dropped
(v37: a connected network of fibres is kept whole — before, it could be
dropped as one big "round" piece), the frame's 6-px border is never a fibre,
the injection sites go on the long pieces, and `aqua_tune\<run>\fibre_mask_raw.png`
shows the result over the mean image — look at it once per field. On a dim or noisy field the
filter finds little, so the mask is relaxed and, failing that, the brightest
5 % of the flattened mean image is used — the log says which, and how many
sites fit (fewer than 72 on a sparse field). To draw the fibres yourself,
save a binary image of the movie frame's size (e.g. draw on
`<run>_meanproj.tif` in ImageJ, Edit → Selection → Create Mask, save) as
`<run>_bin31f_1p98s_fibre_mask.tif` next to the movie, or as
`FOV<k>_fibre_mask.tif` in the animal folder for every day of that FOV; it
replaces the automatic mask, and a changed mask makes the injected movie
again. The same mask gives the on-fibre counts in the report.

What the numbers mean and do not mean: the null movie measures the false
alarms that come from pixel noise — motion artefacts are coherent across
pixels and are not in it, so the false-alarm share is a lower bound and the
on-fibre fraction and the registration carry that part; recall is measured on
synthetic transients, so it is comparable between settings and days rather
than an absolute detection rate of real afferent events.

## Step 4 — Build the batch: `aqua_batch.py`

```
python aqua_batch.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --aqua-dir "D:\AQuA2" --um-per-px 0.792 --border 5
```

Read the printout: experiments, run labels, which settings file was used, what
was queued. A FOV without a settings file gets `--opts <file>`, else the
pipeline's `default_Opt.csv`; tiffs with an `_AQuA2.mat` are skipped unless
`--overwrite`, or unless their FOV's settings file is newer than the result
(a tuned file from Step 3) — then they are queued again, which the printout
says. Writes `<mouse>\aqua_batch\`: `batch.csv`,
`parameters_for_batch.csv`, `run_aqua_batch.m`, `run_cfu_batch.m`,
`manifest.json`. Options: `--fov 2` (one FOV; with `--overwrite` = redo it),
`--sec-tol`, `--stim-after-run N`, `--cfu-min-events 1` (events a place needs before it becomes a CFU; 1 = every event gets a unit, so single events are not lost), `--cfu-overlap 0.3` (see *CFU rules* below).

### CFU rules (v41) — what a unit is, and why the overlap threshold is 0.3

A CFU (consensus functional unit) is AQuA2's "one place that produces
repeated events". The pipeline applies AQuA2's **published rule**
(`cfu.CFU_tmp_function` + `cfu.hierarchicalClusteringMinMeasure`): every event
gets a duration-weighted footprint (pixel weight = the fraction of the event's
frames in which the pixel was active); two events overlap by the weighted
intersection over union of these footprints; events whose overlap is at least
the **overlap threshold** join one unit (single linkage, so chains join too);
a unit needs **min events** events. AQuA2's own defaults are 0.5 and 3; the
pipeline uses **0.3 and 1** (tab 4, `--cfu-overlap`, `--cfu-min-events`).

Why 0.3: on the three RearChr animals (96 runs, 3,860 events) the repeated
events of one bouton have footprints of different extent, so their overlaps
are only 0.25–0.45 — at 0.5 a spot that fired seven times came out as seven
one-event units (97 % of all units were single events). The pooled
distribution of all 10,429 overlapping pairs falls off smoothly between 0 and
0.5 with no gap, so the value had to be chosen by a criterion, and it was
optimised on all 96 runs (2 Nov 2026): every overlapping pair was labelled
without using the overlap — *same place* when the two events' duration-
weighted centres lie within half the smaller event's radius of each other
(AQuA2's "consistent spatial foundation"), *two places* when each centre lies
outside the other event's extent — and for every threshold each recording was
clustered exactly as AQuA2 does (single linkage, chains included). The
agreement of the resulting units with the labels (F1) peaks at 0.3 pooled
(0.275 in 2608-1 and 2608-2, 0.3 in 2609-1): at 0.3, 74 % of same-place pairs
end up in one unit and 8 % of two-place pairs do; at 0.2 the chains pull 44 %
of two-place pairs together, at 0.5 only 0.4 % of the repeats are joined.
At 0.3 the merged groups are compact (centres a median 3.5 µm apart, 7 % of
groups wider than 12 µm) and 70 % of the 0.5 units remain. The figures are in
the project notes (`overlap_threshold_optimisation.png`,
`overlap_threshold_all_runs.png`, `260911_001_cfu_rules_compared.png`).
In a methods section: *"CFUs were defined with AQuA2's overlap clustering;
the overlap threshold was set to 0.3 (AQuA2 default 0.5), the value at which
the clustering best reproduced same-location repeats across all recordings
of the three animals, because bouton events of varying extent overlapped only
0.25–0.45."* Run the final pre/post statistics once at 0.5 as well (uncurated
is enough) to show the conclusion does not hinge on it.

**Does 0.3 fit this animal? — the overlap-threshold check (v42).** The value
is one number for the dataset, but the question can be re-asked on any animal
from its own results: tab 4 *Check the overlap threshold (from the CFU
results)* or `python cfu_overlap_check.py <mouse folder> [--thr 0.3]`. The CFU
batch (v42 scripts) writes next to each result `<run>_AQuA2_cfu_pre.mat` —
AQuA2's own pairwise event overlap table (`cfu.CFU_tmp_function`'s `s_t0`)
plus, per event, the duration-weighted centre, the footprint size and the
first/last frame. The check labels every overlapping pair without using the
overlap (*same place* = the two centres within half the smaller event's
equivalent radius of each other; *two places* = each centre outside the other
event's extent; the rest unscored), clusters the events of each recording as
AQuA2 does at every threshold 0.10 … 0.50 (single linkage, chains included)
and measures the agreement of the units with the labels (F1 of "same unit"),
per experiment (day × FOV, pre and post together) and pooled. Each experiment
gets **OK** (the chosen threshold reaches ≥ 90 % of that experiment's own best
F1), **CHECK** (it does not — the experiment's optimum is printed) or *too few
pairs* (fewer than 30 same-place pairs: nothing can be said; most single
recordings are in this class, which is why the threshold is not tuned per
recording). Outputs: `aqua_batch\cfu_overlap_check.txt` (plain words; the
report's section 9 quotes its verdict), `.csv` (one row per experiment + pooled,
with the whole F1 curve) and `.png`. The check changes nothing — change the
field in tab 4 only if several experiments say CHECK, knowing that every CFU
file and the curation are then redone. Results made before v42 have no
`_cfu_pre.mat`: the check lists them and asks for a CFU re-run.

Changing the overlap or min events makes the CFU batch **redo every result's
CFU file** (the script compares the `cfuOpts` saved in each `_res_cfu.mat`
with the current values) — and the **curation has to be redone**: the
decisions in `_res_cfu_curated.json` refer to the old units. Since v43 the
curation window recognises them (file dates + the rules recorded in the
file), ignores them, opens those runs clean and keeps the old file as
`_res_cfu_curated_superseded_<date>.json` on the first save (Step 6).

**AQuA2's newer CFU rules (not used).** AQuA2 copies from 28 Aug 2026 on
(GitHub `cfu/algo` branch, merged 14 Sep 2026) add two rules after the
overlap clustering: every cluster is cut into watershed basins of its weight
map (hard-coded parameters), and touching units whose movie-derived dF/F
curves correlate at r ≥ 0.85 + 0.15 × the GUI's *Post-merge correlation*
control are merged (default control 0.5 → r ≥ 0.925). On this data the basin
splitting puts a single event with an uneven footprint into up to 13 units
(260911_001: 53 events → 154 units, 28 events in 2–13 units each), so a unit
no longer equals a set of events and one event is counted several times. The
pipeline therefore runs the CFU step with the same (mid-August 2026) AQuA2
copy that detects the events (`H:\MATLAB\AQuA2` on samson). A copy of the
newer release is installed next to the animal folders
(`GCAMP_affrerents\AQuA2_2026-10`, event detection identical) and can be
named in tab 4 ("Other AQuA2 copy for the CFU step only") to compare; the CFU
script prints which rules it is running, and a result made by the other rule
set is redone automatically when the batch runs again with the pipeline's.
The curation window's title says which rule made the result it shows.

## Step 5 — MATLAB (one paste)

```
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_aqua_and_cfu.m')
```

One script runs the AQuA2 batch and then the CFU batch (`run_aqua_batch.m`
and `run_cfu_batch.m`, which can still be run one at a time). Paste it into any
MATLAB window that can see the mouse folder: the scripts find AQuA2 by
themselves (`which('aqua_cmd_batch')`) when the folder set in the GUI is not
visible from MATLAB — the "AQuA2 folder cannot be located" stop of v31 is gone.
In the GUI, tab 5 shows and copies the command, *Check what MATLAB has
produced* counts the results (every queued tiff needs an `_AQuA2.mat` newer
than the batch build, every result a `_res_cfu.mat` and an `_AQuA2_events.csv`
newer than its `.mat`), and *Run both now* starts `matlab -batch` from the
window (the MATLAB field in the top row; output in the log and
`aqua_batch\run_*_log.txt`; *Stop* ends it). *Run all (guided)* does the same
check every 30 s and goes on by itself when everything is there (or when the
outputs stop changing for 20 min after some progress — MATLAB finished with
something missing; the log says what).

Results next to each tiff: `_AQuA2.mat` (~1 GB, keep it), `_AQuA2_Ch1.csv`
(events; read by Step 7), `_Ch1_curves.xlsx`, rising maps; then, from the
CFU script, `_AQuA2_res_cfu.mat` (CFUs), `_AQuA2_events.csv` (every event's
frames, centre, footprint size and box — the feature table has no
location), `_AQuA2_eventmap.tif` (events per pixel) and `_AQuA2_meanimg.tif`
(mean of the registered movie). The last three feed the maps of Step 7 and
the QC report. The CFU script skips results whose CFU file is newer than the
AQuA2 file, so only re-run FOVs are redone (delete a run's `_res_cfu.mat` to
redo it).

## Step 6 — Curate the CFUs: `cfu_curate.py`

*Several at once, and a first opinion (v31):* drag a box with the left mouse
button over either image and every CFU whose centre (or half of whose
footprint) lies in it is selected — the title says how many — then `k` / `x`
act on all of them (shift-drag adds to the selection, ctrl-click toggles one,
Esc clears; ctrl/shift-click works in the list too). The list's **looks**
column rates every CFU against the fibre mask of the registered mean image
(`afferent2p/cfu_shape.py`): *fibre 85 %* = that share of the footprint lies
on fibres; *off-fibre 5 %* = it sits beside them (vessel wall, dura, debris);
*blob* = a compact bright spot off the fibres (cell body, autofluorescent
speck). Rows that look non-fibre are tinted; **Select suggested** (`u`) picks
them all so that one `x` marks them noise, and `f` shows the mask in orange.
It is a suggestion to look at, not a decision: a bright varicosity ON a fibre
stays "fibre", a faint fibre the mask misses can be flagged wrongly, and a
CFU on the edge of a bundle may read 30 %. `python cfu_curate.py "<mouse>"
--list --suggest` prints the suggestions per run without the window. On
2608-1 the rule flags exactly the vertical vessel-wall band of 260921_004/006
(the tuning sweep's hotspot) — 21 of the 39 CFUs of run 004.

*Tuning overlay (v30):* `python cfu_curate.py "<mouse>" --tune 1` (GUI tab 6,
field "overlay the events of tuning setting #") draws the events AQuA2 found
with setting 1 of the tuning sweep (`aqua_tune\<run>\s01\real`) as magenta
boxes T1… on both panels, each with its onset, a trace of its box read from
the movie, and `[` `]` jumping to it — so a candidate setting can be judged
on the movie before it is applied. They are read-only: `s` saves only the CFU
decisions. A run without a CFU result but with tuning results opens on the
sweep's mean image with the boxes alone. The setting numbers are the
"setting" column of `tune_recommendation.txt` (0 = the current setting).


```
python cfu_curate.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1"
```

The window shows the mean image with the CFU outlines (left) and the binned
movie (right) with the same outlines, the selected CFU's dF/F curve, its
numbers, and the event raster. Every CFU has its own colour (overlapping or
touching CFUs never share one); noise = thin grey, hand-drawn = cyan, the
selected one is thick white over black. The list on the left ("CFUs of this
run") is the easy way to pick one of several overlapping CFUs; clicking the
same spot again cycles through them. `Tab` / `Shift-Tab` step through the CFUs.

Zoom: scroll wheel over either image (around the cursor, both images move
together), `+`/`-`, `0` = whole field. Move the zoomed view (v33): drag with
the right (or middle) mouse button, or tick **move view** (`h`) and drag with
the left button, or `Ctrl` + arrow keys (a quarter of the view per press).
The view is kept while you step through CFUs and runs of the same field.

Movie: `space` play / pause, `←`/`→` one frame (`Shift`: 10), `Home`/`End`,
the slider, the fps box, `raw` / `dF/F` view (`v`), a contrast slider,
"outlines on movie" on/off, and `[` / `]` = jump to the previous / next event
of the selected CFU (3 frames before its onset — then `space` shows it).
Click in the trace or in the raster to move the movie to that time.

*Which movie (v33):* by default the panel shows **AQuA2's own registered
movie**, read straight out of the run's `_AQuA2.mat` (`res.datOrg1`: the
binned frames after AQuA2's rigid registration, its border crop and its
bleach correction — the frames the events were detected on, in the same
coordinates as the CFU maps, so the outlines sit exactly on it and the field
does not wander). "show: as recorded" switches to the binned tiff as it was
recorded (the crop offset applied, the motion not) — useful to see how much
the field moved. Needs `h5py` (installed with the pipeline) and a `-v7.3`
result, which `aqua_cmd_batch` always writes; an older result without the
movie falls back to the tiff. Regions you add (`a`) are measured on the
movie shown. The dF/F view is (F − F0) / F0 per pixel with F0 the median of
~60 frames of the run.

Decisions: click or pick a CFU → `k` real / `x` noise; `a` + vertices + Enter
adds a region; `n`/`p` runs; `s` saves; "Propagate to this FOV" copies
decisions to matched CFUs of the other runs. Decisions live in
`_res_cfu_curated.json`. Headless: `--list`, `--reject 3,7`, `--keep-all`,
`--propagate`; `--fps 15` sets the playback speed.

*After a CFU re-run (v43):* the decisions name CFU numbers of one CFU result.
When the CFU batch remakes a run's `_res_cfu.mat` (other CFU rules — e.g. the
overlap threshold 0.5 → 0.3 of Nov 2026 — or a redone AQuA2 result), the
units are renumbered and the old decisions would land on the wrong units.
The window therefore checks every run's curation file against its CFU file
(file dates and the CFU rules recorded in the file) and **ignores decisions
made for an earlier result**: those runs open clean, the title and the run
list say "earlier decisions ignored", a message at start lists the runs whose
old file held decisions, and the old file is kept beside the new one as
`_res_cfu_curated_superseded_<date>.json` the first time the run is saved
(nothing is deleted; delete the superseded files yourself when they are no
longer wanted). `--list` (GUI tab 6 "List CFU results and decisions") prints
per run which rules made the CFU file and ends with "N of M runs … are curated
under the current CFU rules" — that line is the curation status. Before v43
only a redone `_AQuA2.mat` was recognised, not a redone CFU file.

## Step 7 — Analyze: `aqua_analyze.py`

```
python aqua_analyze.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --um-per-px 0.792
```

Writes `<mouse>\aqua_analysis\`: `summary.txt`, `run_summary.csv` (one row per
recording — the table for statistics), `condition_summary.csv`,
`events_<FOV>.csv`, `cfu_matched_<FOV>.csv`, `<FOV>_overview.png`,
`<FOV>_cfu_map.png`. Defaults: events starting in an excluded frame (±1
frame; the flagged frames of artefact stretches, or what you decided in tab 2) dropped and
that time taken out of the denominator, with a note of how many events start
in state stretches you kept; AQuA2 frame numbers converted back to recorded
frames when frames were dropped (`t0_frame`; `t0_frame_movie` is the frame in
the movie AQuA2 saw); CFUs marked noise and their events dropped; added
regions included. Options: `--min-p 2`, `--gate still` (not for the
97 %-running animal), `--baseline d260917`, `--no-curation`,
`--keep-flagged-events`, `--event-source all`, `--no-pupil`. Event-level p
values are descriptive; the unit for statistics is one value per animal.

**Pupil and activity.** For every FOV whose runs have a pupil trace, Step 7
also writes `pupil_<FOV>.csv` (one row per run), `pupil_cfu_<FOV>.csv` (one
row per CFU) and `<FOV>_pupil.png`, and adds `pupil_*` columns to
`run_summary.csv`. It answers two questions:

*Does the pupil (arousal) go with afferent activity?* Per run, the mean CFU
dF/F (or the smoothed event rate when there is no `_res_cfu.mat`) is
correlated with the pupil diameter at lags of −20…+20 s (Spearman; positive
lag = activity follows the pupil). Reported: the strongest r and its lag
(`pupil_dff_r_best`, `pupil_dff_lag_s`), a shuffle p from circular shifts of
the pupil trace (`pupil_dff_p`; the traces are autocorrelated, so a plain p
value would be far too small), the partial r with running speed taken out
(`pupil_dff_partial_r` — pupil and running go together, r ≈ 0.45 in the real
recording, so a raw correlation with the pupil may be a correlation with
running), and r on still frames only (`pupil_dff_r_still`, needs 40+ still
frames). Per CFU, `pupil_cfu_<FOV>.csv` gives r against the pupil and whether
it is `up`/`down`/`no` (outside the 95 % band of shuffled traces);
`n_cfu_pupil_up/down` count them per run.

*Does the optostimulation change the pupil, and is the activity change just
an arousal change?* `pupil_median_px` and `pupil_range_pct` pre vs post say
whether the animal's arousal state differs between the conditions. The
pupil-size terciles are then defined once per FOV (pooled over pre and post)
and the event rate is computed inside each tercile
(`rate_low/mid/high_pupil`), so `post vs pre … rate change at matched pupil:
low x…, mid x…, high x…` compares pre with post at the same pupil size. A rate
change that holds in every tercile is not explained by arousal; one that is
only there because post had a larger pupil disappears within the terciles.
Excluded frames are left out of all of this; a tercile rate
needs at least 15 s of data, the pre-vs-post ratio at least 30 s in both
conditions.

`<FOV>_pupil.png`: per run the pupil trace above the activity trace (grey =
excluded frames) with the lag curve on the right; then the event rate per
pupil tercile per condition (bars = pooled, dots = runs), the per-CFU
correlations, pupil vs running, the pupil distribution per condition with
the tercile edges, and mean activity against pupil-size octiles per
condition. Pupil size in pixels is comparable within a mouse and camera
setting only.

**The comparison itself (v27).** `ANALYSIS_PLAN.md` explains the reasoning;
in short, every experiment (one FOV on one day for RearChr, one FOV across
days for Treg) is compared with its baseline in five ways, all in
`experiments_summary.csv` and drawn in `<mouse>_effects.png`: the raw rate
ratio; the ratio inside still frames only and inside running frames only;
the ratio at matched pupil size; and the **behaviour-adjusted ratio** from a
Poisson regression on 20-s bins (condition + running speed + running state
+ pupil; 95 % CI from a 2-min block bootstrap) — the number to lead with,
with how much of the rate variation behaviour explains versus the condition.
`<FOV>_adjusted.png` shows rate against running per bin, the rate around
running onsets, and raw vs adjusted. Then *where* the activity is: the runs
are aligned to the first (rotation + shift of the mean images), the common
region is the part of the field every run shows, the fibres are segmented
from the averaged image (fibre atlas), and events are counted on 25-µm
fibre tiles (`--tile-um`) and per fibre segment (a fibre crossing another
stays one segment). `<FOV>_maps.png`, `tiles_<FOV>.csv`,
`segments_<FOV>.csv`; the summary gives the rate in the common region (also
per mm² of fibre), the map correlation, the fraction of active tiles active
in one condition only (activity that moved), and segments up/down beyond the
behaviour factor. `--no-maps` skips this part; it needs the
`_AQuA2_events.csv` / `_meanimg.tif` exports of Step 5.

**At matched running speed, and whether matching is possible at all (v30).**
The regression above can only take the behaviour out where the two
conditions overlap in behaviour; when the mouse ran through all of post and
sat through all of pre, no model can. Step 7 therefore also sorts the 20-s
bins of each experiment into four speed bands — still, and slow / medium /
fast running (terciles of the running bins' speed, set once per experiment)
— and reports per band the minutes and the rate per condition
(`speed_bands.csv`, `<mouse>_speed_tuning.png`). The **matched-speed ratio**
(`ratio_speed_matched` in `experiments_summary.csv`) averages the band rates
with one fixed set of weights (the pooled minutes per band) over the bands
both conditions visited for ≥ 1 min — a direct standardisation. Next to it
the **behaviour overlap** says what fraction of each condition's minutes lies
in such shared bands; below 50 % the line reads `NOT COMPARABLE` and no ratio
is given, because the two conditions were recorded in different states and
only the design can fix that (record pre and post in the same state, or long
enough that both states occur in both).

**Time since the stimulation (v30, optostim days).** The two fields alternate,
so the post runs of a day sit at four different times after the stimulation
and the pre runs at four times before it. Step 7 puts every run on the day's
clock — the end of a recording is the modification time of its `.sbx`
(`_quadrature.mat` / `.mat` as fallbacks) next to the movie, the start that
minus the recorded duration; without such files it reads `run_times.csv`
(`run`, `end` [, `start`] as ISO times or epoch seconds) from the mouse
folder, else it places the runs in order with assumed gaps and says
"ESTIMATED". Outputs: `timecourse_runs.csv` (per run: start, minutes since
the first pre / first post recording, rate, ratio to its field's pre rate,
behaviour, event properties), `timecourse_bins.csv` (2-min bins, `--bin-min`),
`timecourse_steps.csv` (consecutive runs of one field: pre→pre, pre→post
across the stimulation, post→post), `timecourse_summary.csv` (one row per
day) and `<mouse>_timecourse.png` (top: rates on the session clock with the
stimulation gap; bottom: ratio to the field's pre rate against minutes since
the first post recording, with the fitted slope). The summary gives per day:
the gap between the last pre and the first post run; the **slope of the post
rate** per 10 min (Poisson regression with a per-field intercept; 95 % CI
from a block bootstrap; also with speed / running / pupil taken out) — does
the effect decay, build up or stay flat over the hour; the second post run
vs the first; and the **step at the stimulation** (last pre run → first post
run of the same field) against the largest change between consecutive runs
of one field *within* a condition — if the step is not larger than that
drift, a pre/post difference of the same size can be the session's trend.
The stimulation itself is not timed; the axis is minutes since the first post
recording. Name the days' conditions with `<mouse>\day_labels.csv` (`day,label`,
e.g. `260917,sham`) or `--day-label 260917=sham`; the labels appear in the
summary and the figures. `--no-timecourse` skips this part.

**Scale (v31).** A ratio of event rates means little when one side has a
handful of events: 1 → 10 events reads ×10, 51 → 144 reads ×2.8, and 0 → 8 has
no ratio at all, although in events per minute the changes are +0.3, +3.0 and
+0.26. Step 7 therefore writes, per comparison, the **rate difference**
(events/min, the one measure on the same scale for every experiment), the
ratio's **exact 95 % interval** (two Poisson rates; infinite when the baseline
has no events), a **log2 ratio with a pseudo-count** of 0.5 events (defined
for every experiment), the **modulation index** (condition − baseline) /
(condition + baseline) in −1 … +1, and a **few-events flag** (fewer than 10 in
total or fewer than 3 in a condition — the direction may be clear, the size
is not; `low_counts` in `experiments_summary.csv`). `summary.txt` has a
"scale:" line per comparison that says which of these applies. The figure
`<mouse>_rates.png` puts every condition of every experiment on ONE rate axis
with its Poisson interval and the single runs, and the change underneath on
one log axis with its interval; the per-experiment overview figures share one
y axis across the animal; the effects figure draws the raw ratio with its
interval and keeps experiments without baseline events (hollow, pseudo-count).
When comparing days or FOVs, use the difference and the intervals, not the
bare ratios.

### v34 — the three triggers: running onsets, pupil dilations, rest

The locomotion trials answer "does the stimulation change the running-evoked
activity?" — and nothing when the animal does not run (no trials), or runs
all the time (no still 10 s before a bout, no rest). v34 compares the same
question on three triggers side by side, so that every session has a
comparison it can carry:

* **running onsets** — the locomotion trials below (10 s still, 30 s
  response window);
* **pupil dilations during stillness** — the same windows around a rise of
  the pupil (DLC diameter per 2-s frame, scaled 0–1 between the run's 1st
  and 99th percentile) of ≥ 20 % of that range within 10 s after a stable
  (range ≤ 0.15), still, fully valid 10-s baseline, one onset per 30-s
  window; a dilation with no running in the response window either is a
  *still* dilation — arousal without movement, the set the comparison uses
  (dilations during running are kept as *moving* in `pupil_trials_<exp>.csv`
  and the figure). Whether a dilation triggers afferent activity is an open
  question this answers per session; the pupil remains a control elsewhere
  (section 4's matched-pupil terciles, the LocoPupil2Fluor GLM, the
  behaviour-adjusted ratios);
* **resting activity** — kept events in still frames ≥ 20 s away from any
  running (the whole run when the animal never runs), per minute of rest.

Per trigger, condition vs baseline: the **evoked event rate** (kept events
starting in the response windows per minute of window; rest: per minute of
rest) with the rate difference, the exact Poisson interval of the ratio, the
pseudo-count log ratio and the modulation index — the statistics of the rate
comparison in section 1/6 (`afferent2p/stats.py`) —, the **responsive-trial
fraction** (z rule) with a risk ratio and Fisher's exact p, and the
**population response amplitude** (peak − baseline dF/F per trial) with a
bootstrap interval and the Mann–Whitney p. A comparison is **testable** when
both conditions have ≥ 3 trials (rest: ≥ 1 min each); otherwise the row says
why — *the animal never ran*, *ran 95 % of the time: no still 10 s before a
bout*, *no still pupil dilation*, *no pupil data*, *less than 1 min of rest*
— and the other triggers carry the comparison. The verdict reads the rate
first (an increase / decrease the interval supports; direction only on few
events; not distinguishable) and the amplitude when the rate cannot speak.
Files: `trigger_comparison.csv` (one row per experiment × condition ×
trigger), `pupil_trials_<exp>.csv`, `<exp>_triggers.png` (onset-aligned
population dF/F for running and still-pupil onsets, resting and evoked rates
with exact intervals, the verdicts); the report's section 3 opens with the
table of all comparisons and that figure, the locomotion trial details follow.
The thresholds (`PUPIL_RISE` 0.2, `PUPIL_STABLE` 0.15, `MIN_TRIALS` 3) are at
the top of `afferent2p/triggers.py`.

### v33 — the motion / registration check

AQuA2 does its own registration before detecting: with `registrateCorrect =
2` every binned frame is shifted by one whole-pixel translation that best
matches the mean of the first 10 frames (`reg.regCrossCorrelation`), then
`regMaskGap` plus the motion margin is cropped off the border; the detection,
the CFUs and the event boxes all live in that registered frame. What that
cannot undo: movement *within* a 2-s bin (31 raw frames are averaged before
AQuA2 sees them), deformation that is not one translation of the whole field,
and z-motion. The analysis now measures all of this per run
(`afferent2p/motion.py`; "motion / registration check" on the Analyze tab,
`--no-motion` to skip):

* from the binned tiff: the rigid shift of every frame against AQuA2's
  reference (`dy`, `dx`, `|shift|`), the jump between consecutive frames, the
  correlation to the reference, a sharpness index;
* from the registered movie in the `_AQuA2.mat`: the rigid residual (should be
  ≤ 0.5 px) and a non-rigid index — the field is cut into 3 × 4 blocks, each
  block matched on its own; the spread of the block shifts around the
  whole-field shift is the deformation one translation could not remove, and
  how far the blocks moved since the previous frame is the deformation that
  happened at that frame;
* against the behaviour csv: everything split into still / running frames
  and the correlation of |shift| with running speed;
* against the events: how many event onsets fall within ±1 frame of a jump
  (≥ 2 px, or more when the run is restless) or of a deformation jump, versus
  the number expected if events and motion were unrelated. A ratio ≥ 2 with
  several events = **warn** ("events follow the motion"); the QC report lists
  and flags those runs, and its gallery shows every event's frames so you can
  see whether the "event" is a smear at a fibre edge.

Per run: `<run>_motion.csv` (one row per frame) and `_motion.json` (the
summary) next to the tiff, rebuilt when the tiff or the result is newer;
per mouse `aqua_analysis\motion_summary.csv`; the QC report has the figure
(shift traces, jumps, residual after AQuA2, running speed, event onsets) and
the numbers for every recording, and the main report's QC paragraph sums
them up. Noise floor of the indices on these movies: ~0.2 px (white noise
at the real level); a 0.5 % stretch of the field reads 1.1 px.

### Is AQuA2's registration enough? `reg_compare.py` (v33)

Up to five ways of registering the same recordings, the same AQuA2 settings,
compared on what they detect — the experiment that answers "do we need the
legacy registration as well, or both?" for this data:

| variant | movie given to AQuA2 | AQuA2's own registration |
|---|---|---|
| **A_aqua2** | the binned tiff as the pipeline makes it (the current route) | on (`registrateCorrect = 2`) |
| **B_prereg** | the raw 15.6-Hz frames registered first — one rigid translation per raw frame to a still-period reference (`afferent2p/fullreg.py`, the legacy operation in Python, also `unpack_bin.py --register`) — then binned | off (`= 1`) |
| **C_both** | the same pre-registered tiff | on again |
| **D_legacy** (`--add-legacy`) | the lab's legacy MATLAB registration itself, run by MATLAB on the raw frames (`legacy_register_run.m`, below) and binned the same way | off |
| **E_legacy_both** (with D) | the same lab-registered tiff | on again |

```
python reg_compare.py pick    "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1"
python reg_compare.py prepare "D:\...\NavRearChrTIGRE-2608-1" --runs 260921_004,260921_006,260921_002,260917_006 --add-legacy
   -> in MATLAB:  run('D:\...\NavRearChrTIGRE-2608-1\reg_compare\run_reg_compare.m')
python reg_compare.py report  "D:\...\NavRearChrTIGRE-2608-1"
```

`pick` ranks the recordings by running time, locomotion onsets, events and
(where the motion check has run) movement; `prepare` registers the raw
frames of the chosen runs in Python (a few minutes each; the `.sbx` must be
next to the tiff), builds the variant folders under
`<mouse>\reg_compare\<run>\` with the run's own tuned settings (from
`aqua_batch\manifest.json`) and one MATLAB paste: with `--add-legacy` the
legacy registration first (`run_reg_compare_legacy.m`), then the AQuA2 + CFU
batch on every variant (results stay in the variant folders, the run folders
are untouched, and every other tool skips `reg_compare\`; `--variants A,D,E`
runs a subset). `report` writes `reg_compare_summary.csv`,
`reg_compare_events.csv` and `reg_compare.html` — per recording and variant:
how much the movie AQuA2 was given still moved, what is left after AQuA2
(rigid residual, deformation), events per minute, the share of events on the
fibre mask, event onsets on moving frames versus chance, CFUs, and which
events every variant finds versus one variant alone (centres aligned across
the variants' frames and matched within 15 px with overlapping times), plus
a figure with the mean images, the event boxes (red = found by that variant
only) and the motion traces. Reading it: if A's events are also found by the
pre-registered variants, A has no excess of onsets on moving frames and no
lower share of events on the fibres, AQuA2's own registration is enough; if
B/C (or D/E) find the same events with fewer motion-locked extras, switch
the unpack to **register the raw frames before binning** (tab 1 checkbox /
`unpack_bin.py --register`) for every run of the animal and redo Steps 4–8
(same treatment for every run) — or, if only the legacy chain does it, run
`legacy_register_run.m` on every run and point the batch at its tiffs.
Measured so far on 2608-2 and the running runs of 2608-1 (Oct 14): the 2-s
frames sit 1–3 px (max 4–8) from AQuA2's reference, register well, blur
within the bins is ~1 % of contrast, and event onsets do not sit on moving
frames more often than chance — the prior is that A is enough; the comparison
is the proof for the lab.

**The legacy chain, as `legacy_register_run.m` runs it.** The lab's
`ProcessScanboxData_concatFirst.m` (dataSet `Afferents`, `regParam.method =
'rigid'`) registers a day's concatenated runs in two stages and bins
afterwards: **CorrectExpt** takes every raw frame against the mean of its
own 15-frame chunk (edges `[80 80 40 40]` px cropped, 2× `imresize`, central
95 %, Gaussian blur 1 px, `dftregistrationAlex` to 1/100 px), subtracts the
global median and applies the sub-pixel shift with `imtranslate` — the
within-second jitter goes, chunks are not aligned to each other; then
**AlignPlanes / RegisterSBX** builds the reference from the longest
zero-speed epoch (≥ 30 s, 5 s trimmed off each end, `DetermineReference`)
and aligns every frame to it with MATLAB's intensity-based `imregtform
'rigid'` (rotation + translation; RegularStepGradientDescent, 1000
iterations, maximum step 0.001, MeanSquares) on edge-cropped, 2×-binned
frames, applied with `imwarp`; **GenerateExptProjections** bins to 1 Hz
(15 frames; the pipeline's 31 here) and crops the projection edges. The
wrapper calls the same legacy functions with the same settings on one run
and one channel, chunk by chunk in memory (no `.sbxdft` / `.sbxreg` files),
keeps the first bin, crops the output by the shifts (99.5th percentile)
instead of the fixed edges so the frames line up with the other variants,
and normalises the brightness level of every binned frame as the pipeline's
unpack does (the frame flagging is not ported). It needs the legacy code
(found automatically at `lagacy_matlab\Matlab` next to the mouse folders —
`MovieProcessing\`, `Misc\`, `Locomotion\` — or `--legacy-matlab DIR`), the
Image Processing Toolbox, and time: stage 2 is `imregtform` on every raw
frame (~15 000 per run), 0.2–0.5 s each without a parallel pool (an hour or
two per run; `parpool` first and it is minutes). Standalone:
`legacy_register_run('…\<run>.sbx', '…\D_legacy', 'legacyDir', '…\lagacy_matlab\Matlab', 'binFrames', 31, 'tag', 'bin31f_1p98s', 'behaviorCsv', '…\<run>_bin31f_1p98s_behavior.csv')`
— `'maxScans', 3000` for a quick test, `'method', 'translation'` / `'dft',
false` to switch stages off, `'outputCrop', [80 80 40 40]` for the legacy
projection edges. Output next to the tiff: `_legacy_shifts.csv` (per raw
frame: the DFT shift, the rigid translation and rotation) and `_legacy.png`.
Differences between the Python stand-in (B) and the real chain (D): whole-
pixel instead of sub-pixel shifts (a blur below 0.3 px RMS after 31-frame
averaging), no rotation, a Hann window instead of the edge crop for the
estimate, and the reference from this run instead of the whole day; if B and
D disagree on a run, D is the lab's answer.

**What the comparison found (2608-1, four heavy-running runs, Oct 2026).**
The field drifts by a few px over a 16-min run and jitters by less than half
a pixel within a second; running and still frames move the same. AQuA2's own
registration (A) takes the drift out to the nearest pixel and gives the
sharpest registered mean image; the Python pre-registration adds nothing
(its per-raw-frame estimates are noisy at single-frame SNR); the legacy chain
cannot see the slow drift (stage 1 is intra-chunk, stage 2's optimiser moves
the frames by 0.02–0.05 px) and its two bilinear resamplings smooth the pixel
noise, so AQuA2 reports 2–4× more, mostly weak, events on its movie at the
same settings. Decision: the pipeline stays as it is (no `--register`, no
legacy chain). Details: `reg_compare\reg_compare_verdict.md` + `.png`.

**See it rather than read it: `reg_demo.py` (v35).** A movie of one run with
the registrations side by side, zoomed on the region with the most
structure, with a fixed cyan outline of each movie's reference (the mean of
its first 10 frames, AQuA2's reference) on every frame — whatever slides off
the outline is motion the route left in — plus the per-frame shift traces,
the running speed and a cursor:

```
python reg_demo.py "%M%" --run 260921_006                 -> %M%\reg_compare\<run>_reg_demo.mp4 (+ _excerpt.gif, .png)
python reg_demo.py "%M%" --run 260921_006 --variants raw,aqua,legacy --scale 2.5 --fps 20 --region 200,300,160,240
```
Panels: *as recorded* (the binned tiff), *after AQuA2's own registration*
(`res.datOrg1` of the run's `_AQuA2.mat` — the frames the events were found
on; when the .mat is not readable the same step is emulated in numpy and the
panel says "emulated"), *legacy MATLAB chain* and *Python pre-registration*
when the run went through `reg_compare.py prepare` (their tiffs in
`reg_compare\<run>\D_legacy\` / `B_prereg\`). The mp4 needs `imageio` +
`imageio-ffmpeg` (`pip install imageio imageio-ffmpeg`); the GIF excerpt (the
5 s in which the as-recorded frames move most) and the PNG need only pillow.
Any run works, not only the compared ones: with the run's own `_AQuA2.mat`
the movie shows exactly what the detection saw.

### v32 — the general representation, the heatmaps and the locomotion-evoked activity

The analysis now also writes, per experiment, what the lab's presentation
shows (figures in `aqua_analysis\`, the numbers in csvs the report reads):

* **Sessions on the day's clock** — `<mouse>_session_<day>.png` (optostim,
  both FOVs of the day) / `<exp>_session.png` (Treg, one block per day):
  events per minute (baseline grey, condition coloured, the stimulation gap
  hatched, excluded frames orange), the running speed with the running state
  shaded, the bouts and the locomotion trials (▼), the pupil.
* **Activity by behavioural state** — `<exp>_state.png`, `state_rates.csv`:
  events/min and median event max dF/F in running vs resting frames (dot =
  run), % running, and the event footprints that started while running vs
  resting on the mean image; `<exp>_ecdf.png`: ECDFs of amplitude (÷ the
  FOV's 99th percentile), duration and area by condition × state.
* **Heatmaps** — `<exp>_heatmap.png`, `cfu_traces_<exp>.npz`: CFU identities
  (matched across the runs) × time, every row read from every run's binned
  movie under the identity's footprint (dF/F, F0 = 20th percentile), scaled
  0–1 per CFU, k-means (k = 3) groups ordered by their change, the event
  onsets and the speed underneath; `<exp>_updown.png`, `cfu_updown_<exp>.csv`:
  CFUs up / down / unchanged (20-s bin means, Mann–Whitney p < 0.05) and
  their event rates. Without readable movies the rows are AQuA2's CFU dF/F
  where the CFU was detected (`--no-movies` forces this).
* **Locomotion-evoked activity** (Yu's trial pipeline) —
  `loco_trials_<exp>.csv` (one row per trial), `loco_cfu_<exp>.csv`,
  `loco_glm_<exp>.csv` (+ kernels), `loco_summary.csv` (per run),
  `loco_conditions.csv` (per condition), `<exp>_loco.png`,
  `<exp>_loco_cfu.png`. A trial is 10 s still → 30 s running → 20 s still
  (strict) or ≥ 10 s still before a bout ≥ 4 s (relaxed); responses = the
  population dF/F (mean of the kept CFUs), the event-count trace and every
  CFU; a trial is responsive by three rules reported side by side — *yu*
  (peak ≥ 1.5 × baseline on F/F0, the slide's rule), *z* (peak − baseline
  ≥ 2 × the still-period noise; the default), *ev* (≥ 1 kept event in the
  running window) — tell the pipeline which one the lab uses. Also: per-CFU
  locomotion responsiveness against circularly shifted onsets, rest-only
  activity (events ≥ 20 s away from any running), the behaviour check
  (% running, bouts, pupil dilation), Loco2Fluor and LocoPupil2Fluor ridge
  GLMs with held-out deviance explained and leave-one-predictor-out
  (`--no-glm` skips them; they are the slowest part; `--no-loco` skips the
  whole block).
* **Pupil measures** — `pupil_measures.csv`: pupil 0–1 per FOV and day,
  variation (SD), resting pupil (valid with ≥ 10 % resting frames),
  locomotion-related dilation per trial, spontaneous dilations at rest.
* **Acute and chronic change** — `acute_chronic.csv`,
  `<mouse>_acute_chronic.png`: for every measure above, post / pre per day
  (acute, RearChr) and the pre value of each day / the first day's (chronic,
  both paradigms), per FOV, on a log axis with 1 = no change.

These animals run most of the time, so strict trials and rest-only minutes
are few — read every number with its N.

## Step 8 — The reports: `aqua_report.py`

```
python aqua_report.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1"          -> <mouse>\aqua_analysis\<mouse>_report.html
python aqua_report.py "D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --qc     -> <mouse>\aqua_analysis\<mouse>_qc.html
python aqua_report.py "..." --all                                                           both
```

Two self-contained HTML files per animal (GUI tab *8 Report*: *Build the
report*, *Build the QC report*; the last step of *Run all*). Open them in a
browser; click any figure to zoom and pan; print to PDF from the browser if a
PDF is needed. Both are built from `aqua_analysis\` (and, for the QC
report, the movies) — nothing is recomputed.

**Across animals — the stimulation protocols side by side (v35):**
`python aqua_conditions.py "D:\2photon\Chao\GCAMP_affrerents"` (GUI tab 8,
*Across animals* block: parent folder = the folder that holds the mouse
folders, *Compare the protocols*). It reads every RearChr animal's tab-7
tables (`aqua_analysis\experiments_summary.csv`, `condition_summary.csv`,
`heatmap_summary.csv`) and `day_labels.csv`, labels each day × FOV
experiment with the protocol of its day — *sham*, *2 Hz*, *5 Hz*, *5 Hz +
40 ms* are read from the label text ("2 Hz max power" → 2 Hz, max power;
"5 Hz max power + 40 ms (new FOV)" → 5 Hz + 40 ms) — and compares the
pre → post change per protocol: the event rate on all / still / running
frames (log2 post/pre with pseudo-count, exact intervals, few-events flag,
the pooled rate ratio per protocol), the event amplitude and duration, the
CFU fate (up / down / unchanged) and the behavioural state (time still pre
and post, the confound to read the raw ratio against). Mean ± SEM per
protocol, Wilcoxon vs 0 and Mann-Whitney vs the sham days when there are
enough experiments — descriptive with 2–6 experiments per protocol. No pupil
step. Writes `<parent>\cross_animal\conditions.html` (figure + tables),
`conditions.png`, `conditions_experiments.csv`, `conditions_groups.csv`,
`conditions.txt`, and lists the days that are not analysed yet (raw only /
unpacked but no tab 7). Run tab 7 for every animal first; curation is not
required for a first look (tab 7 with *no curation* ticked).

**The report** (`<mouse>_report.html`, seconds, ~30 MB) reads in this order:
a QC paragraph with a link to the QC report; **1 General representation** —
the sessions on the day's clock, activity by running / resting state with
the event overlays, the ECDFs, the acute / chronic change of every measure
(figure + table); **2 Heatmaps** — CFU × time, up / down / unchanged, the
tile maps; **3 Triggered activity** (v34) — the three triggers (running
onsets, still pupil dilations, rest) compared per experiment on the evoked
rate, the responsive fraction and the response amplitude, with why a
comparison is not testable; then the locomotion trial definitions, the three
responsive rules, one table per experiment × condition (trials, responsive %,
response, events per trial, responsive CFUs, rest-only rate, behaviour check,
GLM) and the figures; **4 Pupil**; **5 Coupling GLMs** (table + kernels);
**6 Behaviour-controlled change** — the table with the rate difference, the
raw ratio with its exact interval (*few events* where it cannot be sized), the
modulation index and every behaviour control, the verdicts in words, the "on
one scale" and effects figures; **7 Time course** on the day's clock and the
matched-speed bands; **8 Per experiment** (run table, overview, CFU map,
maps, pupil, the experiment's part of `summary.txt`); **9 Settings and
tuning**; **10 Methods**.

**The QC report** (`<mouse>_qc.html`, minutes; one page per experiment plus
an index when it would exceed ~80 MB) is the evidence that the detection is
right: a summary table with one row per recording (minutes, excluded frames,
events detected / kept / dropped, CFUs kept / noise, how many look like fibre /
off-fibre / blob, curation status, crop match, movie and export present,
locomotion trials) with warning chips (stale curation, > 30 % non-fibre CFUs,
no event export, > 10 % frames excluded, 0 events, < 3 trials), sortable by
clicking a column, each row linking to its recording; a navigation panel on
the left; and per recording, collapsed until opened: the mean image with
**every detected event's box** (blue kept, red dashed in a CFU marked noise,
grey dashed in an excluded stretch) and the CFU outlines coloured by what
they look like (yellow fibre, red off the fibres, pink blob); the events on
the recording's clock against running speed (▼ = locomotion trials) and pupil;
**clips** of the movie around a few kept events; a **gallery with a montage of
every detected event** (12 at a time; filters *all / kept / dropped / in
non-fibre CFU / in a locomotion trial / top amplitude / long*, sort by onset,
amplitude or duration); and the CFU table (events, looks, on-fibre %,
brightness, locomotion-responsive, status). A real afferent transient rises
over one or two frames on a fibre-shaped footprint and decays; a flicker, a
motion smear or a bright non-fibre structure is a noise detection — mark its
CFU in Step 6, run Step 7 and build both reports again. `--max-events N`
limits the montages per recording (0 = all), `--clips N` / `--no-clips` the
clips. The montages need the binned movies next to the results; without them
the sections have the maps and rasters only.

## Adding a day / changing settings

Run the sequence again: each step processes only what is new. A new FOV gets
the default settings until it is tuned (or its own `FOV4_Opt.csv`). To change
one FOV's settings: re-save its `_Opt.csv` (or tune again), `aqua_batch.py`
(a run whose settings file is newer than its result is queued by itself;
`--fov 2 --overwrite` forces it), the MATLAB scripts, re-check its curation,
`aqua_analyze.py`. Whole animal from scratch: `--overwrite` without `--fov`.
A new day of a RearChr mouse (e.g. 2608-1's 261002): Step 1 for the new
runs, then *Run all (guided)*: the plan adds that day's first pre run per
FOV, the sweep runs only its jobs (the old ones are done), the report
recomputes the common setting per FOV over all days and rewrites only the
settings files whose values changed — so only the new day's runs go
through AQuA2, unless the setting itself moved (then the log says which
files changed and every run of that FOV is queued again). Add the day's
condition to `day_labels.csv` first.
A new day of a Treg mouse: Step 1, then *Run all (guided)* again — with *skip the sweep
if every FOV is tuned* ticked it keeps the FOV's common setting and only runs
the new day through AQuA2 (untick it to check the setting on the new day too:
the plan adds that day's movie, the sweep runs only its jobs, and the report
chooses again over all days). If the new day's run numbers do not follow
the rule (a lost FOV, a recapture, pre + post of a second injection), say so
in `run_labels.csv` first (next section).

## Which run is which FOV: `run_labels.csv` and `day_labels.csv` (v44)

The planner reads only file names. Its rule: a spontaneous (Treg) mouse
records FOV1, FOV2, FOV3 … as runs 001, 002, 003 … of every day, condition =
the date (`d260925`), the earliest day = baseline; an optostim (RearChr) mouse
records FOV1, FOV2, FOV1, FOV2 before the stimulation and the same after it.
Real days break the rule: a FOV lost to bone growth, a plane that shifted and
was recaptured under the next run number, a second injection with pre- and
post-injection recordings on one day. `<mouse>\run_labels.csv` says what the
rule cannot — one row per run that deviates (the others follow the rule):

```
run,fov,condition,use,note
NavTIGRE-2609-1_261010_001,1,d261010_pre,1,
NavTIGRE-2609-1_261010_002,2,d261010_pre,1,
NavTIGRE-2609-1_261010_003,,,0,lost FOV (bone growth)
NavTIGRE-2609-1_261010_004,,,0,shifted plane - replaced by 006
NavTIGRE-2609-1_261010_005,2,d261010_post,1,
NavTIGRE-2609-1_261010_006,1,d261010_post,1,recapture of FOV1
```

`run` = the run name (or just `<yymmdd>_<nnn>`); `fov` = the FOV the run
belongs to (empty = unchanged); `condition` = its label (empty = the rule's:
`d<date>`, or pre / post) — a same-day pre / post pair is written
`d<date>_pre` / `d<date>_post`, date first, so the day's label still applies
and the two sort in recording order; `use` = 0 leaves the run out of every
step (batch, tuning, analysis, curation list); `note` is for you. GUI tab 1,
block *Which run is which FOV*: **Create / complete run_labels.csv** writes
the file with every run and the rule's values (an existing file keeps its
rows, new runs are appended — do it again after a new day), **Open** edits
it in Excel / Notepad, **Show the plan** prints what every step will use
(`python plan.py <mouse>` does the same). Rebuild the batch (tab 4)
afterwards — the status strip says "run_labels.csv newer than the batch"
until you do. Letters, digits, `_`, `+`, `-` only in a condition name.

`day_labels.csv` (`day,label`) names the days for the analysis and the
reports: `260925,Baseline`, `261001,D5`, `261010,D14 (second injection)`; a
row `261010_post,D14 after the injection` names one same-day condition
exactly (otherwise the day's label plus the suffix is used: "D14 (second
injection) post"). Both files survive a fresh start (`aqua_reset.py` keeps
them).

## Starting over: a fresh start (v35)

When results of several versions are mixed up in an animal and you want one
defined state — every run processed by the current version — press **Fresh
start…** (top row, next to *Run all*), or on the command line

```
python aqua_reset.py "%M%" --dry-run                   :: list what would go
python aqua_reset.py "%M%"                             :: park everything in %M%\_previous_<date>\
python aqua_reset.py "%M%" --stages aqua,analysis      :: only the AQuA2 + CFU results and the analysis
python aqua_reset.py "%M%" --discard-decisions         :: also the frame decisions and the CFU curation
python aqua_reset.py "%M%" --delete                    :: delete instead of parking (no way back)
```

The dialog lists, per stage, what would go and how big it is — *1 unpacked
movies + sidecars*, *3 tuning sweep + its settings* (`aqua_tune\`,
`aqua_opts\`), *4-5 AQuA2 + CFU results + batch* (`aqua_batch\` and every
`_AQuA2*` file of every run), *7-8 analysis + reports* (`aqua_analysis\`) —
and what the re-run costs (unpack ~3 min per run; the sweep ~4.5 h of MATLAB
per day of recordings with one worker; AQuA2 + CFU ~6 min per run). The old
files are **moved** into `<mouse>\_previous_<yymmdd_HHMM>\` with their
folder structure — a rename on the same drive, instant and reversible; the
pipeline never looks into that folder — and you delete it when the new run
is through (tick *DELETE* to skip the parking when the disk is short).
Never touched: the recordings (`.sbx`, info `.mat`, quadrature, eye videos,
DLC files), `day_labels.csv`, settings files you saved by hand (`*_Opt.csv`
in the mouse folder or a run folder), hand-drawn `*_fibre_mask.tif`,
`reg_compare\`, `afferent2p_out\`. Your decisions — the tab-2 frame
decisions (`*_frames_decisions.json`) and the tab-6 CFU curation
(`*_res_cfu_curated.json`) — stay unless you tick *also discard my
decisions*; they apply again to the re-made results (the curation by CFU
number only if the detection gives the same CFUs — re-check tab 6 after a
re-run with other settings). *Reset and run everything (guided)* then runs
tab 1 (unpack) and the whole *Run all* chain — tuning plan, the sweep in
MATLAB, report + apply, batch, AQuA2 + CFU in MATLAB, analysis, reports —
with the tab-3 grid and movie choice as set; *Reset only* just parks, and
you run the tabs yourself. *Preview* writes the list to the log without
moving anything.

## When something goes wrong

| You see | Do |
|---|---|
| `Property assignment is not allowed when the object is empty` at `aqua_gui` start | run `patch_aqua2_toolbar('D:\AQuA2')`, restart MATLAB (needed after every AQuA2 update) |
| `curSz` error when saving in the GUI | update AQuA2 (or the two-line patch); the .mat/csv were already written |
| files named `..._bin31f_1_Ch1.csv` | old dotted tiff name; harmless in the GUI, correct in the batch |
| `NO settings file` for a FOV / `NOTHING QUEUED` | the pipeline's `default_Opt.csv` is missing — unzip the package again, or save `FOVk_Opt.csv` from the GUI |
| `run_aqua_tune.m not found` in MATLAB | the plan was not written: tab 4 queued nothing (see above) or *1 Plan* was not pressed — press it; it builds the batch itself |
| *Auto* said `The AQuA2 folder is not on this PC` (v31) | gone in v32: *Run all (guided)* never stops on it; the MATLAB scripts find AQuA2 on MATLAB's path when the folder in the top row is not visible from Python |
| *Run all* waits at a MATLAB step although MATLAB has finished | something expected is missing (tab 5 *Check what MATLAB has produced* says what, e.g. a run whose CFU step failed): press *Continue now* in the small window to go on with what exists, or fix it and let the 30-s check pick it up; after 20 min without change it continues by itself |
| `AQuA2 not found: set the AQuA2 folder in the GUI (top row) or addpath it in MATLAB` | MATLAB cannot see the folder in the top row and has no AQuA2 on its path: `addpath(genpath('<AQuA2 folder>'))` in that MATLAB, or set the folder as MATLAB sees it |
| *Run … in MATLAB now* says no MATLAB executable | top row: MATLAB → Browse… (`C:\Program Files\MATLAB\R20xx\bin\matlab.exe`); needs R2019a+ for `-batch` |
| MATLAB finished with exit code ≠ 0 | the script hit an error — the last lines in the log / `run_*_log.txt` say where; a sweep job that failed leaves `FAILED.txt` in its folder and the others go on |
| `plan: skipping <run>.sbx — no <run>.mat next to it` | an aborted recording (Scanbox writes the info .mat only when a recording ends; e.g. 261016_005 of 2608-1, 266 frames) — not a run; its re-take has the next free number (051) and is placed by recording time: `plan: run numbers [...] are not 1..8 — ordered by recording time` |
| tab 3 *Why this setting? (pictures)* says `tune_why.html does not exist yet` | the sweep report was made before v35 — press *3 Report* again (seconds; the MATLAB sweep is not repeated) |
| `matched sensitivity: NOT POSSIBLE` in the tuning report / `DETECTION NOT COMPARABLE` in the verdict (v36) | no pair of clean settings detects pre and post with the same recall (within `--tol-recall`) and the same false-alarm rate: one movie's best clean recall is far below the other's (look at `<day>_FOV<k>_why_day.png` — usually a dim or unfocused recording, or the field moved). The day gets the common setting and its ratio is flagged; a larger `--tol-recall` (0.10) or `--max-false` accepts the mismatch knowingly; otherwise treat that day's change as not interpretable |
| `[pre] settings of this condition: …` in the batch log, and `aqua_opts\` holds `<exp>_pre_Opt.csv` and `<exp>_post_Opt.csv` | the matched-sensitivity rule wrote one file per condition (v36): pre and post are detected at different thresholds ON PURPOSE, with the measured sensitivity equal — not a mistake. `--day-rule common` (GUI: *one common setting for the day*) + *4 Apply* writes one file again and removes the pair |
| the buttons at the bottom of a tab are hidden behind the log | drag the divider between the tabs and the log (v36), or make the window taller; the window opens at 1180 × 900 |
| the sweep's choice finds no (or almost no) events on a field that has activity by eye — recall ≤ 0.3 at every setting (v37) | the injected transients are dimmer or smaller than the field's events, the grid stops short of the setting that works, or the fibre mask missed the fibres: read the report's grid verdict, classes check and mask note; set the injected dF/F to 0.3, 0.6, 1.0 and sizes to 20, 60, 200 (the v37 defaults), add minSize 20, press *1 Plan* (only the changed jobs are swept) — and save a hand-tuned result of that run in the AQuA2 GUI, export it (*Hand result → reference*) so the report scores every setting against it |
| `GRID TOO STRICT … sweep minSize=20` / `thrARScl=1.25` in the report | the best clean setting is the most permissive one swept and still misses most injected transients: *Extend the grid (as suggested)* (tab 3) or `--plan --extend`, run the sweep (new combinations only), *3 Report* again |
| `matched sensitivity: TOO WEAK` / `DETECTION NOT COMPARABLE` with recall < 0.3 | both movies can only be matched at a level where most transients are missed — the day is not comparable as it stands; fix the grid / the classes (the lines above) before reading the day |
| `the injected classes probe BELOW the real signal` | the brightest injected class, as AQuA2 reports it, is dimmer than the real events: plan again with the `--amplitudes` the line gives (tab 3 *injected dF/F*); only the inject jobs are run again |
| `the fibre mask may be incomplete` | most events at the chosen setting lie off the mask: open `aqua_tune\<run>\fibre_mask_raw.png`; if the fibres are not outlined, draw `FOV<k>_fibre_mask.tif` (above) and plan again |
| a run folder has files of several pipeline versions and you no longer know which result is current | *Fresh start…* (top row) parks everything the pipeline made in `_previous_<date>\` and runs it all again with this version; `--dry-run` / *Preview* lists it first |
| the progress row says `no item count for this step` | that step prints no `[i/n]` lines (batch build, CFU list, matlab -batch from here …) — the log shows what it does; the MATLAB paste steps count their results instead |
| `mp4 writer not available` from reg_demo.py | `pip install imageio imageio-ffmpeg`; the GIF excerpt and the PNG are written anyway |
| `was saved for day 260911; applying it to 260917 too` | intended |
| `found twice — using …` | a copied recording; the run folder's copy is used |
| `no binned tiff` | run Step 1, or `--sec-tol` |
| `Missing Python package` | `pip install h5py pandas tifffile scipy matplotlib` |
| `PermissionError … aqua_batch\manifest.json is locked by another program or owned by another user` (v40), or a `WARNING: … written as parameters_for_batch_2.csv instead` | an existing output cannot be overwritten. Two causes seen on samson: a csv open in Excel (close it), or the file was written by a pipeline window started from an **Administrator** Anaconda Prompt — such files belong to BUILTIN\Administrators and a window started the normal way (double-click `run_gui.bat`, a plain prompt) cannot change or delete them, nor can MATLAB. v40 writes the batch/plan files under the next free name (`…_2.csv`, `run_aqua_batch_2.m`; the manifest records the names and the GUI's MATLAB commands use them), but the manifest itself and the analysis outputs must be overwritable: either keep using an Administrator window for that animal, or once, in an Administrator prompt, `icacls "D:\2photon\Chao\GCAMP_affrerents" /grant <user>:(OI)(CI)F /T` and from then on start the window the normal way |
| *Run both now* / *run the sweep from here* (matlab -batch) ends in seconds with `AQuA2 not found`, and `run_*_log.txt` says every `H:\MATLAB\…` folder is `nonexistent` | the window runs as Administrator: a mapped network drive (H:) belongs to the normal logon and is invisible to elevated processes (Windows "linked connections"). Paste the command into a MATLAB window that sees H: (the normal route), start the pipeline window from a plain prompt, or use the share's UNC path / a local AQuA2 copy in the top row |
| the sweep stops at once with `Error using diary … tune_log.txt: Cannot open file: permission denied` | the old log was written by another account (see the PermissionError row): a sweep script written by v40 or later falls back to `tune_log_<date>.txt` by itself; for an older script, add the try/catch around `diary(logFile)` (or delete the old log from an Administrator prompt) and paste the command again — finished jobs are skipped |
| `run_gui.bat` says `'python' is not recognized` | python is only on the PATH inside the Anaconda Prompt; v40's `run_gui.bat` activates the base Anaconda itself when it is in `C:\ProgramData\Anaconda3` or `%USERPROFILE%\anaconda3` — otherwise edit `PY=` in the file |
| `already corrected … skipping` | fine; `--force` re-does it with the current rules (needed once for tiffs corrected by an older version) |
| `results … computed on the PREVIOUS movie` | only after `--drop`/`--fill`: Step 4 with `--overwrite` |
| `N frame(s) are no longer on disk` | an older version dropped or replaced them; `--restore-raw` if `..._raw.tif` exists, else re-unpack that run (`--overwrite`) |
| frames in the check stack look smoothed / flattened, label says `filler` | the movie holds fillers there (v22 replaced flagged frames by a blend of clean frames — no noise, so they look smooth); they count as flagged and stay excluded. Typically the old `..._1.98s.tif` copy of a day-1 run: the current `..._1p98s.tif` next to it has the real frames, and the old copy with its sidecars can be deleted (the tools skip it when both exist) |
| `CHECK (old version)` in tab 2 | the recording was corrected by an older version that excluded whole stretches; run the step with *force* once to get the per-frame check and the stacks |
| a long amber (state) stretch in the QC plot | a slow change the rules cannot classify; kept by default — scroll through it in the frame stack and look at the running trace, decide (*Exclude whole stretch* if it is a laser/focus change), and note it |
| the movie is shorter than the recording / AQuA2 frame ≠ raw frame | only after `--drop`: `_frames.csv` column `tiff_frame` maps them, Step 7 converts event frames back |
| the movie has frozen stretches | filled frames (`--fill`, or a v22 tiff) — `--restore-raw` or re-unpack, then the step again |
| big slow "events" of 100 s or more spanning half the field | a state stretch that AQuA2 turned into events: in the GUI mark that stretch *Exclude whole stretch* (no AQuA2 re-run needed); if it is not even flagged, `aqua_clean.py --force --shape-rel 0.03` |
| `no _res_cfu.mat` | run `run_cfu_batch.m` |
| `NOT aligned … unshifted maps` | check `_cfu_map.png`; the field may have moved |
| hundreds of events in a few frames | open `_globalfix.png`; if those frames are not flagged, `aqua_clean.py --force --min-rel 0.02` (decisions apply at Step 7, no AQuA2 re-run) |
| `pupil: no pupil data` | no `*eye*DLC*.csv` next to the `.sbx` of any run of that FOV; copy the DLC files into the run folders and run Step 7 again (it rebuilds the behaviour csv itself) |
| `pupil: N eye frames for M imaging frames … resampling` | the eye camera and the 2P clock differ by a frame or so (14751 vs 14750 is normal); a ratio far from 1 means the wrong DLC file is in the folder |
| `still-only n/a` in the pupil summary | fewer than 40 still frames with a pupil value in that run (a running animal) — use the partial r instead |
| `where: no per-event locations` | the AQuA2 results predate v27's CFU script; delete the runs' `_res_cfu.mat` and run `run_cfu_batch.m` again (it rewrites CFUs and the exports) |
| `NOT ALIGNED` in the maps line / figure | that run's mean image does not match the first run's (different plane or field): its events are still counted for the field-level numbers but the maps of that experiment are unreliable |
| `not enough bins with behaviour data` | the behaviour csvs lack running or pupil for most frames — `unpack_bin.py --behavior-only` |
| `still: not enough still time in both conditions` | the matched-state ratio cannot be formed; use the adjusted ratio (it uses every bin) |
| QC chip *events follow the motion* / motion **warn** (v33) | the run's event onsets sit on moving frames ≥ 2× more often than chance: open the run's gallery in the QC report and look at those events (a smear along one side of a fibre that comes and goes with a jump = motion, mark its CFU noise); if most runs of a FOV warn, the field moves too much for the 2-s bins — register at full rate before binning (the suite2p route of the old pipeline exports registered, binned tiffs) and redo Step 4 for that FOV |
| the curation window shows *as recorded* although "AQuA2-registered" is ticked | the run's `_AQuA2.mat` has no readable `res.datOrg1` (older AQuA2 output, or `h5py` missing — `pip install h5py`); the log line says which. Re-running the AQuA2 batch for that run writes a movie the window can read |
| a FOV whose CFU result has no or almost no CFUs although the fibres are there | look at the sweep report (`tune_recommendation.txt`) for that FOV: a chosen setting with recall 0.00 and 0.0 events/min detects nothing by construction (this happened to 2608-2 FOV 2: the one-setting-for-all-days rule found no sensitive setting that kept the false-alarm share ≤ 10 % on every day). For an optostim mouse untick *one common setting per FOV for all days* before *4 Apply* (each day is compared within itself) and let Step 4 → 5 redo the re-queued runs |

## All commands in order

```
cd D:\2photon\Chao\GCAMP_affrerents\afferent2p
set M=D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1
python pipeline_status.py "%M%"                                      :: where is this animal? (one line per step + the next action)
python aqua_reset.py "%M%" --dry-run                                 :: v35: a fresh start — park (or --delete) every output and run the steps again with this version
python unpack_bin.py "%M%" --sec 2                                   :: 1 Unpack
python unpack_bin.py "%M%" --behavior-only                           :: only if the DLC eye csv arrived after unpacking (Step 7 also does this itself)
python aqua_clean.py "%M%" --force                                   :: 2 Clean & check — once for tiffs corrected by an older version; then decide per stretch (GUI tab 2)
python aqua_batch.py "%M%" --aqua-dir "D:\AQuA2" --um-per-px 0.792 --border 5    :: 4 Batch (a FOV without a settings file gets default_Opt.csv)
:: 3 Detection settings by measurement (overnight; "Run all (guided)" in the GUI does this block and the rest):
python aqua_tune.py "%M%" --plan                                     :: then MATLAB: run('%M%\aqua_tune\run_aqua_tune.m')
python aqua_tune.py "%M%" --progress                                 :: the morning after
python aqua_tune.py "%M%" --report --apply                           :: per-day settings (RearChr) / one per FOV (Treg); also aqua_tune\tune_why.html (v35: the pictures behind the choice); then 4 Batch again
:: v36: a RearChr day gets ONE SETTING PER CONDITION with matched sensitivity by default (aqua_opts\<exp>_pre_Opt.csv / _post_Opt.csv + <exp>_tuning.json);
python aqua_tune.py "%M%" --report --apply --day-rule common         ::   the v34/v35 rule instead (one common setting per day); --tol-recall 0.05 = the allowed recall gap
python aqua_tune.py "%M%" --plan --extend                            :: v37: add the values the last report suggested (grid verdict) and sweep only the new combinations
python aqua_tune.py "%M%" --figure-only --dpi 260                    :: v38: the tuning figure again at higher resolution (no rescoring)
python -c "from afferent2p.pick import main; main()" "%M%"          :: v39: pick a setting by eye among the measured candidates, on the real movie (GUI: tab 3 "Pick by eye…")
python cfu_overlap_check.py "%M%"                                   :: v42: does the CFU overlap threshold fit this animal? (GUI: tab 4 "Check the overlap threshold")
:: v37: a hand-tuned AQuA2 result as the reference — in MATLAB: aqua_export_events('%M%\<run>\<movie>\<movie>_AQuA2.mat'), then --report again
```
MATLAB (5): `run('…\aqua_batch\run_aqua_and_cfu.m')` — the AQuA2 batch, then the CFU batch
```
python cfu_curate.py "%M%"                                           :: 6 Curate
python aqua_analyze.py "%M%" --um-per-px 0.792                       :: 7 Analyze (--no-glm skips the GLMs, --no-loco / --no-heatmaps the blocks)
python aqua_report.py "%M%" --all                                    :: 8 -> %M%\aqua_analysis\<mouse>_report.html + <mouse>_qc.html
:: optional (v33): is AQuA2's registration enough? up to five registrations of the same runs, compared
python reg_compare.py pick "%M%"                                     :: which recordings move most
python reg_compare.py prepare "%M%" --runs 260921_004,260921_006 --add-legacy  :: registers the raw frames, builds the variants + one MATLAB paste (legacy chain + AQuA2 + CFU)
python reg_compare.py report "%M%"                                   :: after MATLAB -> %M%\reg_compare\reg_compare.html
python unpack_bin.py "%M%" --register --overwrite                    :: only if the comparison says so (it did not): raw frames registered before binning, whole animal
python reg_demo.py "%M%" --run 260921_006                            :: v35: the registrations as a movie, side by side -> %M%\reg_compare\<run>_reg_demo.mp4
:: v35: the stimulation protocols across animals (after 7 Analyze on each animal)
python aqua_conditions.py "D:\2photon\Chao\GCAMP_affrerents"          :: -> ...\cross_animal\conditions.html (+ png, csv, txt)
```
A Treg mouse (`NavTIGRE-2608-1`, runs 001–003 = FOV1–3 on every day) takes the
same commands: its experiments are `<mouse>_FOV<k>_days` with the dates as
conditions (baseline = the earliest day; `day_labels.csv` names them, e.g.
`260911,day -1`), the plan tunes one movie per day and FOV, and `--apply`
writes one setting per FOV for all its days.
