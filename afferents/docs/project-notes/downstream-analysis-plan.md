# Downstream analysis — did the afferents change, and is it real?

What Step 6 (`aqua_analyze.py`, GUI tab 6) computes, why, and how to read it.
Two paradigms share it: **RearChr** (per day: runs 1–4 *pre*, 1 h optostim,
runs 5–8 *post*; FOV1 = odd runs, FOV2 = even runs; each day is its own
condition) and **Treg** (baseline day, DTx, then D5/D7/D9; FOV k = run k of
every day; each day is compared with the baseline day). One *experiment* =
one FOV on one day (RearChr) or one FOV across days (Treg); the baseline is
`pre` or the earliest day.

## The three problems and what is done about them

**1. The field of view is never exactly the same.** Every run's mean image
(the registered movie's average, exported by `run_cfu_batch.m`) is aligned to
the first run of the experiment (rotation ±4°, shift). The **common region**
is the part of the field present in every run; every spatial number is
computed inside it only, and the summary says how much of the field that is
and whether a run could not be matched (`NOT ALIGNED` — usually a different
plane; treat that comparison with suspicion). Rates are also given per mm² of
fibre inside the common region, so a run that shows less tissue is not
penalised.

**2. An afferent does not respond at the same place twice.** A CFU (AQuA2's
"functional unit") is a hotspot; the same fibre can have a different hotspot
active in the next recording, so CFU-to-CFU matching across runs mostly
fails and is kept only as an exploratory extra. Three levels replace it:

* *field level* — events/min of the whole common region (and per mm² of
  fibre): the primary number; needs no matching at all;
* *tile level* — the aligned mean images of all runs are averaged into a
  **fibre atlas** (bright ridges = the mCitrine-filled fibres), and events
  are counted on 25-µm tiles of that atlas (only tiles on fibres, inside the
  common region). 25 µm is large enough that a hotspot drifting a few µm along
  a fibre stays in its tile, small enough to show redistribution. Between
  conditions: the map correlation (how much of the pattern is kept), the
  fraction of active tiles that are active in one condition only (activity
  moved), and a tile-wise log2 ratio with a spatial block-bootstrap CI;
* *fibre-segment level* — the atlas skeleton is cut into segments at branch
  points, with a fibre that merely crosses another kept as one segment
  (straight continuation through the crossing). Every event is assigned to
  the segment it lies on. A segment is the afferent-level unit: whichever
  part of it is active, it is the same fibre. Per segment: events/min per
  condition, up/down beyond what behaviour explains (binomial test on the
  split of its events between the conditions), and a Wilcoxon over segments.
  Segments of one axon are not independent, so this is descriptive.

**3. Locomotion and arousal.** Afferent events follow running (the
running-onset-triggered rate in `<exp>_adjusted.png` shows how much), and
the animal is not in the same state pre and post — on the real animal the
still fraction went 86 % → 10 % on one day and 3 % → 38 % on the next. Four
numbers, all in `experiments_summary.csv` and `<mouse>_effects.png`:

* *raw ratio* — events/min post / pre; confounded whenever the still
  fraction differs;
* *still frames only* and *running frames only* — the same ratio inside one
  behavioural state; honest but thin when one condition hardly has that
  state (the minutes are printed);
* *matched pupil size* — the ratio inside pupil-size terciles (defined once
  per experiment, pooled over conditions), combined over the terciles both
  conditions populate;
* *behaviour-adjusted* — a Poisson regression on 20-s bins: events ~
  minutes × exp(condition + running speed + running state + pupil size).
  exp(condition) is the rate ratio at the same behaviour, using every bin;
  its 95 % CI is a block bootstrap (2-min blocks within each run, because
  neighbouring bins are not independent). The summary also says how much of
  the rate variation behaviour explains versus the condition. **This is the
  number to lead with**; the matched-state ratios should agree with it, and
  the tile/segment comparisons are computed with this behaviour factor
  taken out of the post map.

Frames in excluded stretches (flashes, state changes you excluded in tab 1b)
are left out of every count and every denominator.

## Reading the outputs

| file | what it answers |
|---|---|
| `<mouse>_effects.png` | one row per experiment: the ratio by every control, adjusted with CI — the picture to look at first |
| `experiments_summary.csv` | the same as numbers, plus map r, % tiles switched, segments up/down, common-region fraction, minutes per state |
| `<exp>_adjusted.png` | rate vs running per bin and condition; rate around running onsets; raw vs adjusted ratio |
| `<exp>_maps.png` | aligned images with the events, the fibre atlas, tile maps per condition, the change map, per-segment scatter |
| `tiles_<exp>.csv`, `segments_<exp>.csv` | the maps and segments as tables (positions in µm) |
| `<exp>_pupil.png`, `pupil_<exp>.csv` | pupil–activity coupling (lag, partial r vs running, per CFU) and the pupil-tercile rates |
| `run_summary.csv`, `condition_summary.csv` | per recording / per condition: rates, medians of the event properties, still and running minutes and rates |
| `events_<exp>.csv` | every kept event with its recorded frame and the behaviour at its onset |
| `summary.txt` | all of the above in sentences |

How to read a result: the adjusted ratio's CI excludes 1 **and** the
matched-state ratios point the same way **and** the tile/segment level shows
the change is not one hotspot → a real change of afferent activity. Adjusted
ratio near 1 while the raw ratio is large → the raw change was behaviour.
Map r low with many switched tiles but adjusted ratio ≈ 1 → the activity
moved (different functional units) without changing in amount — report that
as its own finding.

## What the design can and cannot support

* One animal, two runs per condition per FOV: everything here is
  descriptive for that animal. The unit for a group statistic is one value
  per animal (or per FOV, stated as such): take the adjusted log ratio of each
  experiment, average within animal, then test across animals (paired for
  pre/post, unpaired for treated vs control). Days are different conditions
  in RearChr, so they are not replicates of each other — report each.
* Event rates depend on the AQuA2 threshold; keep the settings fixed within
  a FOV (the batch does) and quote them. Peak dF/F drifts with mCitrine
  bleaching; prefer rate and peak dF for pre/post.
* Pupil size is in pixels and comparable within one mouse and camera
  setting only; the pupil follows running (r ≈ 0.5), which is why it enters
  the regression together with running rather than alone.
* The regression assumes the behaviour effect is the same in both
  conditions (no interaction) and log-linear in speed; when the two
  conditions hardly overlap in behaviour (day 1 of the real animal) the
  adjustment relies on that assumption — the wide CI says so.
* The fibre atlas needs the mean image to show the fibres; a very dim FOV
  gives few segments and the tile/segment numbers become noisy (printed
  counts tell).
* `NOT ALIGNED` runs, a common region below ~60 %, or fewer than ~20 fibre
  tiles: report the field-level numbers only.
