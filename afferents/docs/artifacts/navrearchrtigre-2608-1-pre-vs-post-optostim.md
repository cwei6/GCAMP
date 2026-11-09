# NavRearChrTIGRE-2608-1 — pre vs post optostim report

Source: https://claude.ai/artifact/BKpBo2S4arAicBipQpubrU (claude.ai artifact, last updated 2026-10-03)

AQuA2 pipeline report (Steps 1–6) for mouse NavRearChrTIGRE-2608-1 comparing afferent event rates before and after 2 Hz optostimulation on a sham, a 3.3 mW and a max-power day, which finds a roughly two-fold behaviour-adjusted increase only in the one interpretable field (260921 FOV2) and shows that this increase is as much a session-long upward trend as a step at the stimulation.

Report heading: *NavRearChrTIGRE-2608-1 — afferent activity before and after optostimulation* (dated 2026-10-03).

Report of the AQuA2 pipeline results (Steps 1–6) for the three RearChr days — 260917 (sham), 260911 (2 Hz, 3.3 mW) and 260921 (2 Hz, max power): each effect with its figure and what it means for this animal.

**In short.** Three stimulation days, each with two fields of view and two 15.7-min runs before and after the stimulation: **sham** (260917), **2 Hz at 3.3 mW** (260911) and **2 Hz at maximum power** (260921). Only one experiment has both enough events and the same behaviour before and after — FOV2 on the max-power day (51 pre, 144 post events). There the afferent event rate after the stimulation is about twice its baseline at matched behaviour (× 2.1, 95 % CI 1.1–4.7; raw × 2.8), the increase is spread over the fibres of the field rather than confined to a few sources, the events become longer and lower in amplitude, and the activity switches from positively to negatively coupled to pupil size. FOV1 of the same day agrees in direction (× 2.7 on running frames, but only 3 → 8 events). Put on the day's clock, however, the max-power FOV2 rate rose through the whole session — × 3.7 between its two pre runs, × 0.86 across the stimulation itself, × 3.2 between its two post runs — so the pooled × 2.8 is a day-long upward trend as much as a stimulation effect (Effect 8). The sham day (0 → 0 and 0 → 8 events) and the 3.3 mW day (5 → 1 and 1 → 10 events, with the mouse resting in one condition and running in the other) have too few events to measure a ratio; they neither support nor contradict the max-power result, and the sham day shows that a handful of events appearing after a procedure is not by itself evidence of an effect. An HTML version with zoomable figures is in the animal's `aqua_analysis` folder (`NavRearChrTIGRE-2608-1_report.html`).

## What was analysed

Three imaging days of NavRearChrTIGRE-2608-1 (jGCaMP8s in the afferents, ReaChR-mCitrine), each a different stimulation condition, plus a test run on 260831 without stimulation that has no pre/post pair and was left out:

| Day | Condition | Runs | Events pre → post, FOV1 | Events pre → post, FOV2 |
| --- | --- | --- | --- | --- |
| 260831 | test run, no stimulation | 001 only | not analysed | not analysed |
| 260917 | sham stimulation | 8 (4 pre, 4 post) | 0 → 0 | 0 → 8 |
| 260911 | 2 Hz, 3.3 mW | 8 | 5 → 1 | 1 → 10 |
| 260921 | 2 Hz, max power | 8 | 3 → 8 | 51 → 144 |

Each day has eight 15.7-min runs of 475 binned frames (31 frames → 1.984 s per frame, 512 × 796 px at 0.792 µm/px). Runs 1–4 were recorded before the stimulation (**pre**) and runs 5–8 after it (**post**); odd runs are field of view 1, even runs field of view 2. That gives six experiments (3 days × 2 FOVs), each with 2 pre and 2 post runs, about 31 min of recording per condition.

Events were detected by AQuA2 in batch with one preset per FOV, identical for pre and post: registration and bleach correction on, smoXY 0.5, active-voxel threshold (thrARScl) 2, significance threshold 3.5, source sensitivity 8, no global-signal detection; minimum event size 100 px for FOV1 and 50 px for FOV2. CFUs (recurring sources) were grouped with a minimum of 1 event. Before the statistics, the Step 1b brightness QC removed the frames with a flash or a shutter artefact (0–9 excluded frames per run, out of 475); events starting in or next to an excluded frame were dropped (at most 1 per run) and their time was taken out of the denominator. No run had been hand-curated, so no CFUs were rejected.

The table lists every run with the events that entered the statistics (kept / detected by AQuA2), the resulting rate, the excluded frames (and the number of frames the QC had flagged in total), and the behavioural state of the run: time spent still, median running speed (encoder units) and median pupil diameter.

| Day | Run | FOV | Condition | Events kept / detected | Events / min | Excluded frames (flagged) | Time still | Median speed | Pupil (px) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 260911 | 001 | 1 | pre | 5 / 5 | 0.33 | 5 (63) | 79 % | 0 | 98.5 |
| 260911 | 002 | 2 | pre | 0 / 0 | 0.00 | 2 (9) | 75 % | 0 | 98.4 |
| 260911 | 003 | 1 | pre | 0 / 0 | 0.00 | 2 (85) | 92 % | 0 | 95.8 |
| 260911 | 004 | 2 | pre | 1 / 1 | 0.07 | 5 (19) | 94 % | 0 | 91.2 |
| 260911 | 005 | 1 | post | 1 / 1 | 0.06 | 3 (29) | 3 % | 30 | 101.9 |
| 260911 | 006 | 2 | post | 9 / 9 | 0.59 | 7 (32) | 3 % | 35 | 103.2 |
| 260911 | 007 | 1 | post | 0 / 0 | 0.00 | 5 (26) | 18 % | 30 | 101.9 |
| 260911 | 008 | 2 | post | 1 / 2 | 0.07 | 9 (28) | 36 % | 16 | 100.7 |
| 260917 | 001 | 1 | pre | 0 / 0 | 0.00 | 3 (46) | 0 % | 46 | 102.1 |
| 260917 | 002 | 2 | pre | 0 / 1 | 0.00 | 1 (20) | 1 % | 35 | 103.3 |
| 260917 | 003 | 1 | pre | 0 / 0 | 0.00 | 5 (21) | 5 % | 34 | 101.5 |
| 260917 | 004 | 2 | pre | 0 / 0 | 0.00 | 1 (12) | 9 % | 31 | 100.7 |
| 260917 | 005 | 1 | post | 0 / 0 | 0.00 | 4 (13) | 18 % | 33 | 98.4 |
| 260917 | 006 | 2 | post | 6 / 7 | 0.38 | 1 (17) | 40 % | 18 | 95.6 |
| 260917 | 007 | 1 | post | 0 / 0 | 0.00 | 4 (11) | 35 % | 19 | 96.1 |
| 260917 | 008 | 2 | post | 2 / 3 | 0.13 | 3 (13) | 60 % | 0 | 90.8 |
| 260921 | 001 | 1 | pre | 3 / 4 | 0.19 | 2 (14) | 4 % | 61 | 94.7 |
| 260921 | 002 | 2 | pre | 11 / 11 | 0.70 | 0 (5) | 1 % | 62 | 96.4 |
| 260921 | 003 | 1 | pre | 0 / 0 | 0.00 | 1 (8) | 20 % | 43 | 99.8 |
| 260921 | 004 | 2 | pre | 40 / 40 | 2.56 | 1 (35) | 12 % | 57 | 96.3 |
| 260921 | 005 | 1 | post | 6 / 6 | 0.39 | 2 (13) | 9 % | 74 | 94.2 |
| 260921 | 006 | 2 | post | 34 / 35 | 2.19 | 2 (28) | 10 % | 72 | 93.1 |
| 260921 | 007 | 1 | post | 2 / 2 | 0.13 | 4 (87) | 14 % | 70 | 94.8 |
| 260921 | 008 | 2 | post | 110 / 110 | 7.05 | 1 (23) | 18 % | 67 | 92.2 |

Two things stand out before any effect is read. Detection is sparse: 20 of the 24 runs have fewer than 10 events, 10 have none at all, and only the max-power day's FOV2 (195 events) has enough events for the behaviour-controlled comparisons below. And the days differ in what the mouse was doing: on the 3.3 mW day it was still for most of the pre runs and ran through the post runs, on the sham day it ran before and was half still after, and only on the max-power day was it running in both conditions. These two facts shape everything that follows.

## Effect 1 — event rate before and after the stimulation (raw)

> **Figure — event rate per condition, day and FOV.** Pooled pre and post bars with the single runs as dots (left, square-root axis so the sparse experiments stay visible), and the raw post/pre ratio with a descriptive 95 % interval (right). Key numbers: only 260921 FOV2 (1.63 → 4.63 / min, × 2.8) has an interval well clear of 1; 260911 FOV2 (× 10.2) and FOV1 (× 0.20) rest on 1–5 events, and sham FOV2 goes 0 → 0.26 / min.

The raw event rate (events per minute of usable recording, both runs of a condition pooled) is the first number out of the pipeline, before any behavioural control. It goes up after the stimulation in four of the six experiments and down in one; the sixth (sham, FOV1) has no events in either condition.

| Condition — experiment | Pre | Post | Raw ratio post / pre |
| --- | --- | --- | --- |
| sham — 260917 FOV1 | 0 events | 0 events | — |
| sham — 260917 FOV2 | 0 events | 8 events, 0.26 / min | — (0 → 8) |
| 2 Hz, 3.3 mW — 260911 FOV1 | 5 events, 0.16 / min | 1 event, 0.03 / min | × 0.20 |
| 2 Hz, 3.3 mW — 260911 FOV2 | 1 event, 0.03 / min | 10 events, 0.33 / min | × 10.2 |
| 2 Hz, max power — 260921 FOV1 | 3 events, 0.10 / min | 8 events, 0.26 / min | × 2.7 |
| 2 Hz, max power — 260921 FOV2 | 51 events, 1.63 / min | 144 events, 4.63 / min | × 2.8 |

How much weight each ratio carries depends on the counts behind it. With fewer than ten events on one side, a ratio is mostly noise: the × 10 of the 3.3 mW day's FOV2 rests on a single pre event, the × 0.20 of its FOV1 on five pre events and one post event, and the descriptive intervals on the right of the figure span more than an order of magnitude for every experiment except the max-power day's FOV2. The sham day's FOV2 went from no events to 8 (0.26 / min) with no stimulation at all — the same size of post-only appearance as on the 3.3 mW day (0.33 / min) — which is the floor against which "events appeared after" has to be read. The only experiment with enough events to measure a change is **260921 FOV2 (2 Hz, max power)**: 51 events before, 144 after, a raw increase of about 2.8-fold, with the interval well clear of 1.

The single-run dots add a caution even there. The two runs of each condition are not alike: the second pre run (004) has 3.7 times the rate of the first (002), and the second post run (008) 3.2 times the first (006), so the rate also rises *within* each condition over the session, and the lower post run (2.2 / min) is below the higher pre run (2.6 / min). Part of the pre-to-post difference could therefore be a trend over time rather than a step at the stimulation. The sham day would be the natural reference for such a drift, but with 0–8 events it cannot show one either way (its FOV2 post events fall 6 in run 006 and 2 in run 008, no rising trend). Note also that the baseline rate of "FOV2" differs enormously between days (0 events in 31 min on the sham day, 51 on the max-power day, with the mouse running on both), so days are not comparable with each other in absolute rate; only the within-day pre/post change is.

## Effect 2 — the behavioural state was not the same before and after

> **Figure — behaviour of the pre and post runs per experiment** (ordered sham → 3.3 mW → max power): time spent still, median running speed and median pupil diameter (bars = mean of the two runs, dots = single runs). Key numbers: time still 3 % → 27 % and 5 % → 50 % (sham), 85 % → 10 % and 84 % → 20 % (3.3 mW), 12 % → 12 % and 6 % → 14 % (max power).

The afferent event rate in this animal depends on what the mouse is doing: in every experiment with enough events, running frames carry more events per minute than still frames (max-power FOV2 pre: 1.7 / min running vs 0.5 / min still; 3.3 mW FOV1 pre: 0.44 vs 0.11), and the pupil is larger while running (pupil–running correlation 0.4–0.6 in the 3.3 mW day's pre runs). A pre/post comparison is therefore only as good as the match in behaviour between the two conditions, and that match differs sharply between the three days.

| Condition — experiment | Time still, pre → post | Median speed, pre → post | Pupil, pre → post |
| --- | --- | --- | --- |
| sham — 260917 FOV1 | 3 % → 27 % | 40 → 26 | 101.8 → 97.2 px (−4 %) |
| sham — 260917 FOV2 | 5 % → 50 % | 33 → 9 | 102.0 → 93.2 px (−9 %) |
| 2 Hz, 3.3 mW — 260911 FOV1 | 85 % → 10 % | 0 → 30 | 97.2 → 101.9 px (+5 %) |
| 2 Hz, 3.3 mW — 260911 FOV2 | 84 % → 20 % | 0 → 25 | 94.8 → 102.0 px (+8 %) |
| 2 Hz, max power — 260921 FOV1 | 12 % → 12 % | 52 → 72 | 97.2 → 94.5 px (−3 %) |
| 2 Hz, max power — 260921 FOV2 | 6 % → 14 % | 59 → 69 | 96.4 → 92.7 px (−4 %) |

On the **3.3 mW day** the mouse sat still through most of the pre runs (median speed 0) and ran through the post runs with a larger pupil, so the post condition is confounded with running and arousal in the direction that raises the rate. On the **sham day** it is the other way round: the mouse ran almost continuously before and spent 27–50 % of the post runs still with a smaller pupil, a confound in the direction that lowers the rate. Only on the **max-power day** are the two conditions in the same state — running in both, with similar time still (12 % → 12 % and 6 % → 14 %) and a pupil that shrank by only 3–4 %. This is the main reason the max-power FOV2 result can be interpreted and the other days' results cannot, independent of their event counts: on the sham and 3.3 mW days the raw ratios compare a resting mouse with a running one.

Pupil diameters are camera pixels and were not calibrated; they are comparable within a day, and across days only loosely (the camera framing may differ).

## Effect 3 — the change at matched behaviour (still frames, running frames, pupil terciles)

> **Figure — post/pre event-rate ratio per experiment** (log scale, 1 = no change; ordered sham → 3.3 mW → max power): raw, still frames only, running frames only, at matched pupil size, and behaviour-adjusted by regression with its 95 % interval; a triangle at the left edge is a subset with no post events, a missing dot a ratio that could not be formed (no pre events). Key numbers for 260921 FOV2: raw × 2.8, still × 18.3, running × 2.3, pupil terciles × 5.8 / × 2.4 / × 3.7, adjusted × 2.14 (95 % CI 1.13–4.66).

The first behaviour control recomputes the ratio on frames of the same kind only: pre still frames against post still frames, pre running frames against post running frames, and pre against post within each pupil-size tercile (terciles set on the pooled frames of the FOV). A change that survives in every subset is not a change of state; one that vanishes in all of them was.

| Condition — experiment | Still frames only | Running frames only | Pupil terciles low / mid / high |
| --- | --- | --- | --- |
| sham — 260917 FOV2 | 0.19 / min vs 0 | 0.32 / min vs 0 | — / 0.21 vs 0 / 0 vs 0 |
| 2 Hz, 3.3 mW — 260911 FOV1 | 0 vs 0.11 / min (3 vs 26 min) | × 0.08 (0.04 vs 0.44 / min; 28 vs 4.5 min) | × 3.9 / × 0 / — |
| 2 Hz, 3.3 mW — 260911 FOV2 | 0 vs 0.04 / min (6 vs 26 min) | 0.41 / min vs 0 (25 vs 4.9 min) | × 0 / 0.25 vs 0 / 0.42 vs 0 |
| 2 Hz, max power — 260921 FOV1 | 0 vs 0 | × 2.7 (0.29 vs 0.11 / min; 27 vs 27 min) | × 0 / 0.38 vs 0 / × 11.9 |
| 2 Hz, max power — 260921 FOV2 | × 18.3 (9.2 vs 0.5 / min; 4.3 vs 2.0 min) | × 2.3 (3.9 vs 1.7 / min; 27 vs 29 min) | × 5.8 / × 2.4 / × 3.7 |

For **260921 FOV2 (max power)** the increase holds in every subset. The running frames are the bulk of the data (27–29 min on each side) and give × 2.3; the still frames give a much larger ratio (× 18) but rest on 2–4 min and about 1 pre event against 40 post events; and the rate is higher after the stimulation in each pupil tercile (× 2.4 to × 5.8). So neither locomotion nor arousal explains the change in this field. One detail is worth keeping: before the stimulation the rate was higher while running than while still (1.7 vs 0.5 / min), after it the still frames had the highest rate (9.2 / min; 37 of the 110 events of run 008 occurred while the mouse was still). The activity became less tied to running, which Effect 7 sees again from the pupil side.

For **260921 FOV1 (max power)** the one well-matched number is the running-only ratio (× 2.7, 27 min on each side), the same direction and size as FOV2 on the same day, but it rests on 3 against 8 events and the pupil terciles scatter (× 0, 0.38 vs 0, × 11.9), as expected from so few events.

For the **3.3 mW day** the matched comparisons collapse: the pre condition has only 4.5–4.9 min of running and the post condition only 3–6 min of stillness, so each subset compares a few events on one side with a few minutes on the other (the × 0.08 of FOV1 is 2 pre events in 4.5 min of running). The **sham day's FOV2** has no pre events at all, so no ratio can be formed in any subset; its 8 post events are spread over still and running frames alike (0.19 and 0.32 / min). These experiments cannot tell whether anything changed.

### At matched running speed — and whether matching is possible at all

> **Figure — event rate per running-speed band** (still, and slow / medium / fast running = terciles of the running bins' speed, set per experiment), pre vs post, with the minutes spent in each band above the bars; each panel title gives the matched-speed ratio (band rates averaged with one fixed set of weights over the bands both conditions visited) and the behaviour overlap (share of each condition's minutes in such shared bands). Key numbers: 260921 FOV2 × 2.9 at 97 % / 89 % overlap (slow × 3.4, medium × 3.7, fast × 1.6), 260921 FOV1 × 4.8; 3.3 mW day post overlap only 27 % (FOV1) and 33 % (FOV2), marked NOT COMPARABLE.

Still / running is a coarse split; this finer version (added 5 October) sorts every 20-s bin into four speed bands and compares pre with post inside the bands both conditions visited. It also makes explicit *how much* of each condition can be matched at all. On the max-power day the overlap is near-complete (FOV2: 97 % of the pre minutes and 89 % of the post minutes lie in shared bands) and the matched-speed ratio is × 2.9, with the increase present in every band (slow × 3.4, medium × 3.7, fast × 1.6; 4.0 vs 9.4 / min in the few still minutes); FOV1 of the same day gives × 4.8 on full overlap but 3 vs 8 events. On the 3.3 mW day only 27 % (FOV1) and 33 % (FOV2) of the post minutes fall in bands the pre runs also visited — the mouse sat still for 26 of the 31 pre minutes and ran at medium or fast speed for 20–22 of the 31 post minutes — so the comparison is marked NOT COMPARABLE and no ratio is given. The sham day's FOV2 has the overlap (97 % / 54 %) but no pre events in the shared bands. The rule this gives for future sessions: when the overlap is below about half, no regression or matching can take the behaviour out, because the two conditions never visited the same state; only the recording design can.

## Effect 4 — the behaviour-adjusted rate ratio (260921 FOV2, 2 Hz max power)

> **Figure — 260921 FOV2, activity and running.** Left: event rate of every 20-s bin against running speed, pre (blue) and post (orange); middle: event rate around the 7 post running onsets that followed at least 20 s of stillness; right: the raw post/pre ratio and the ratio at the same behaviour from the Poisson regression, with its 95 % block-bootstrap interval. Key numbers: post bins above pre at every speed including 0; 13 / min in the 40 s before onset vs 8 / min in the first 20 s of running; adjusted × 2.14 (1.13–4.66) vs raw × 2.83.

The second control fits one model to all four runs instead of slicing the frames: the event count in each 20-s bin (192 bins) is regressed on running speed, a running indicator, pupil size and the condition. The condition coefficient is then the pre-to-post change at the same behaviour, and a block bootstrap over the bins gives its interval. This needs at least 5 events in each condition, so it was computed for the max-power day's FOV2 only.

For 260921 FOV2 the **adjusted ratio is × 2.14 (95 % CI 1.13–4.66)** against × 2.83 raw on the same bins. Behaviour accounts for 15 % of the bin-to-bin variation in rate and the condition for 11 %, so about a quarter of the raw increase is attributable to the slightly different behaviour of the post runs and the rest is not. The model's rates tell the same story as Effect 3: at a given behaviour the rate roughly doubles (still 4.6 → 10.0 / min, running 1.6 → 3.4 / min). The interval is wide and its lower end is close to 1, which is what 51 and 144 events across four runs can support; it excludes no change, but not by much.

The left panel shows why the adjustment moves the ratio only modestly: post bins (orange) sit above pre bins (blue) at every running speed, including the bins at speed 0, so the post increase is not confined to a particular behaviour. The middle panel checks the opposite confound, that running itself triggers the events: around the 7 post running onsets the rate does not rise — it is 13 / min in the 40 s before the onset and 8 / min in the first 20 s of running — so the events are not locked to the start of locomotion. (With 7 onsets this is a sanity check, not a measurement.)

For the other five experiments the adjusted ratio is not available (fewer than 5 events in a condition), and the behaviour mismatch of the sham and 3.3 mW days described in Effect 2 could not be corrected by a model anyway: a regression cannot estimate the effect of the condition when the condition and the behaviour change together.

## Effect 5 — where the activity is (260921 FOV2, 2 Hz max power)

> **Figure — 260921 FOV2, where the activity is.** Top: mean image with the fibre atlas (yellow) and the common region (dashed), then each run aligned to the first with its events at their centres (blue = pre, orange = post); bottom: events per minute per 25-µm fibre tile before and after, the log2 post/pre change per tile (red = more after), and the rate per fibre segment before against after with the behaviour factor taken out. Key numbers: active tiles 23 → 55, tile-wise log2 ratio +1.08 (CI +0.66 to +1.46), active segments 14 → 36, median segment change log2 +2.47 (× 5.5).

The four runs of the field were aligned to the first (shifts of 1–2 px and up to 0.25°) and the comparison was restricted to the region all four runs cover (94 % of the field). Fibres occupy 15 % of that region (34,000 µm², 145 segments), and the fibre area was cut into 221 tiles of 25 µm. Inside the common region the rate rises **× 2.76** (1.63 → 4.50 / min; 48 → 132 events per minute per mm² of fibre), of which the behaviour factor explains × 1.33, leaving **× 2.1** beyond it — the same answer as Effects 3 and 4 by a different route.

The maps say the increase is spread over the field rather than concentrated in one hotspot. Active tiles go from 23 to 55; 82 % of the tiles active in either condition are active in only one of them (65 % newly active after, 17 % lost), and the pre and post tile maps are only weakly correlated (r = 0.25). The tile-wise log2 ratio is +1.08 (spatial bootstrap CI +0.66 to +1.46, i.e. × 2.1). At the level of fibre segments — the closest thing to an afferent-level unit — 14 segments were active before and 36 after; of the 19 segments with at least 3 events, 4 are up and 1 is down beyond the behaviour factor (binomial p < 0.05), and the median segment change is log2 +2.47 (× 5.5; Wilcoxon p = 0.002 across the 19 segments). In the scatter, most segments sit above the no-change line, with a cluster of segments that were silent before and active after along the left edge.

> **Figure — 260921 FOV2, CFU maps per run** (colour = matched identity across runs) and the area covered by CFUs in all four runs (empty: no CFU recurs in every run). Key numbers: 11 and 39 CFUs in the pre runs, 35 and 106 in the post runs, 189 identities, none recurring in all four runs at IoU ≥ 0.3.

The CFU maps (recurring sources from AQuA2's CFU step) show the same thing from the source side: 11 and 39 CFUs in the two pre runs, 35 and 106 in the two post runs, lying along the fibre bundles. No CFU recurs in all four runs at IoU ≥ 0.3 (189 identities in total), so the sources are not fixed hotspots that fire more after the stimulation; they are mostly different stretches of fibre becoming active. A paired CFU-level test is therefore not possible, and the segment-level comparison above is the per-afferent statistic to use.

For the other experiments the maps exist but carry too few events to read: active tiles 4 → 1 on the 3.3 mW day's FOV1, 1 → 7 on its FOV2, 0 → 4 on the sham day's FOV2, none on its FOV1. The max-power day's FOV1 has no map yet because two of its runs (001 and 003) lack the per-event location files (`_AQuA2_events.csv` and `_meanimg.tif`) that the current `run_cfu_batch.m` writes; re-running that script for them and then Step 6 produces the map.

## Effect 6 — the events themselves: smaller in amplitude, twice as long (260921 FOV2, 2 Hz max power)

> **Figure — 260921 FOV2, event properties by condition:** event rate per run, and cumulative distributions of peak dF/F, area, duration, rise time (10–90 %) and decay time (90–10 %) for the 51 pre and 144 post events. Key numbers: median peak dF/F 0.34 → 0.25, duration 42 → 91 s, rise 8 → 22 s, decay 16 → 24 s, area 198 vs 223 µm².

Besides being more frequent, the events after the stimulation look different. Their peak amplitude is lower — median dF/F 0.34 before, 0.25 after (× 0.74; peak dF 796 → 642), and the post curve in the figure sits to the left of the pre curve over its whole range — while they last about twice as long: median duration 42 → 91 s, with both the rise (8 → 22 s) and the decay (16 → 24 s) slower and a tail of events longer than 100 s that does not exist before the stimulation. The area is unchanged (198 vs 223 µm²) and the integrated signal per event (dF/F × time) nearly doubles (2.2 → 4.2), so the per-event calcium load went up through duration, not amplitude. The event-level p values printed by the pipeline (amplitude 0.002, duration 5 × 10⁻⁵, area 0.5) treat events as independent and are descriptive only.

Three readings are compatible with this, and the data do not yet separate them. The first is a genuine change in the afferents: longer, lower-amplitude calcium transients, as expected if the fibres fire in longer, lower-rate bouts after the stimulation. The second is a detection effect of the higher event density: when more sources are active and overlap in time, AQuA2 merges neighbouring transients into one longer event with a lower peak relative to its baseline, which would lengthen the duration and lower the dF/F without any change in the single transients. The third is drift over the session: the second post run (008), which holds 110 of the 144 post events, has by far the longest events (median 100 s against 52 s in run 006) and the lowest amplitude (0.24), so the properties change *within* the post condition as much as between the conditions, and mCitrine/GCaMP bleaching lowers F0 over a session, which biases dF/F upward rather than downward and so does not explain the drop, but reminds that peak dF is the safer amplitude measure across time.

For the other experiments the event numbers (0–10 per condition) do not support a comparison of properties; their medians are in `run_summary.csv` for completeness.

## Effect 7 — the coupling between pupil and activity reversed (260921 FOV2, 2 Hz max power)

> **Figure — 260921 FOV2, pupil and afferent activity.** Per run: pupil diameter and mean CFU dF/F over time with the cross-correlation by lag (+ = activity follows the pupil); below: event rate at matched pupil size (terciles), per-CFU correlation with the pupil, pupil against running speed, the pupil distribution per condition, and mean activity against pupil size (octiles). Key numbers: cross-correlation +0.41 / +0.24 in the pre runs vs −0.27 / −0.31 in the post runs; CFUs coupled 20 positive / 15 negative of 50 before, 3 / 30 of 141 after.

Pupil size is the arousal read-out that does not depend on the wheel, and the pipeline asks two things of it: whether the rate change survives at matched pupil size (Effect 3: yes, × 2.4–5.8 in every tercile), and how the activity is coupled to the pupil in each condition. The second shows a qualitative change.

Before the stimulation, activity and pupil were **positively** coupled: cross-correlation +0.41 (shuffle p < 0.01) and +0.24 (p = 0.07) in the two pre runs, peaking with the activity 12–14 s ahead of the pupil, 20 CFUs positively and 15 negatively coupled of 50, and the mean CFU dF/F rising with pupil size (bottom-right panel, blue). After the stimulation the coupling is **negative**: −0.27 (p = 0.02) and −0.31 (p = 0.05) in the two post runs, 3 CFUs positively and 30 negatively coupled of 141, and the mean activity highest at the smallest pupil sizes and falling as the pupil dilates (orange). The partial correlations with running taken out go the same way (+0.21 / +0.05 before, −0.24 / −0.17 after), so this is not just the running–pupil relation showing through. Together with the still-frame result of Effect 3 (the post rate is highest while the mouse is still), the picture is that the activity after the stimulation no longer rides on arousal; if anything it is strongest when the animal is quiet.

Two cautions. The pupil range within a run is narrow (10–90 % range 4–8 % of the median), so these correlations are over small fluctuations and the per-run values are modest. And the post runs had a smaller pupil overall (96.4 → 92.7 px) with a wider spread; the terciles used for the matched comparison were set on the pooled frames, so the low tercile is mostly post frames and the high tercile mostly pre frames (visible in the pupil-distribution panel), which weakens the matched comparison at the extremes but not the sign change in the coupling, which is a within-run quantity.

For the other experiments the pupil–activity correlation is undefined or rests on fewer than 10 events and is not reported here.

## Effect 8 — time since the stimulation: step or trend, decay or build-up

> **Figure — every run of a day on the day's clock** (from the recording files' time stamps). Top: event rate per 2-min bin (dots) and per run (circles, number = events) for both fields, pre period shaded blue, post period orange, the stimulation gap between them, grey line = running fraction; bottom: post rate divided by the same field's pre rate against minutes since the first post recording (log scale; empty bins at the floor), with the fitted slope. Key numbers: max-power FOV2 0.70 → 2.56 → 2.19 → 7.05 / min per run; post slope × 1.33 per 10 min (0.98–1.85) on the max-power day vs × 0.48 (0.11–0.87) on the 3.3 mW day; pre → post gap 73–74 min (stimulation days) and 49 min (sham).

Because the two fields alternate, the post runs of a day sit at four different times after the stimulation (their centres about 8, 24, 41 and 57 min after the first post recording started) and the pre runs at four times before it. The recording files carry their own clock: on the two stimulation days the first post run started 73–74 min after the last pre run ended (the hour of stimulation plus setup), on the sham day 49 min. Three readings come out of placing the runs on that clock.

| Day | Gap pre → post | Post slope (× per 10 min, 95 % CI; behaviour-adjusted) | Second post run vs first | Same field, run by run |
| --- | --- | --- | --- | --- |
| sham — 260917 | 49 min | not fitted (8 post events) | × 0.34 (FOV2: 6 → 2 events) | FOV2: 0 → 0 → 0.38 → 0.13 / min |
| 2 Hz, 3.3 mW — 260911 | 73 min | × 0.48 (0.11–0.87); adjusted × 0.44 — 11 events, weak | × 0.10 (FOV2: 9 → 1) | FOV2: 0 → 0.07 → 0.59 → 0.07 / min |
| 2 Hz, max power — 260921 | 74 min | × 1.33 (0.98–1.85); adjusted × 1.19 — 152 events | × 2.80 (FOV1 × 0.34: 6 → 2; FOV2 × 3.21: 34 → 110) | FOV2: 0.70 → 2.56 → 2.19 → 7.05 / min (× 3.7, × 0.86, × 3.2) |

**Decay or build-up.** On the max-power day the post rate did not decay over the hour; it rose, × 1.33 per 10 min since the first post recording (the interval just touches 1), × 1.19 once speed, running and pupil are taken out, and × 2.8 between the second and the first post run of a field. If the stimulation has an after-effect in this field, it is one that builds over the hour rather than one that fades. On the 3.3 mW day the few post events (11) are concentrated in the first post run of FOV2 and gone by the second, a decaying pattern, but on 11 events the slope is a weak estimate and most of the drop coincides with the mouse slowing down. The sham day has too few events for a slope.

**Step or trend.** The run-by-run sequence of the max-power FOV2 is the important one: 0.70 → 2.56 / min between the two pre runs (× 3.7, 11 → 40 events), 2.56 → 2.19 across the stimulation itself (× 0.86, 40 → 34 events, with 90 min between those two runs), 2.19 → 7.05 between the two post runs (× 3.2, 34 → 110). The step at the stimulation is smaller than the changes between consecutive runs within a condition, so the pooled post/pre ratio of × 2.8 (Effect 1) is produced by a rate that climbs through the whole session, with the last run contributing most of it. That is compatible with an effect that builds up after the stimulation, but equally with a session-long trend that had already started before it — bleaching changing the detection, the animal's state, the tissue settling under the objective — and the pre/post design cannot tell these apart. The within-condition drift was flagged in Effect 1; the clock makes it concrete.

**Reading for the behaviour-controlled result.** Effects 3–5 established that the post increase is not a change of behavioural state. Effect 8 adds that it is not a step either. The defensible statement is therefore: in this field the event rate at matched behaviour was about twice as high in the hour after the stimulation as in the hour before it, rising through both hours, and a sham day with events on the same schedule (or a longer post period) is needed to say how much of that rise the stimulation caused. The stimulation itself is not time-stamped in the recording files, so the axis is minutes since the first post recording; the delay between the end of the stimulation and that recording should be added from the lab notes. Outputs: `timecourse_runs` / `bins` / `steps` / `summary.csv` and `NavRearChrTIGRE-2608-1_timecourse.png` in `aqua_analysis`.

## Caveats, and what to do next

**What the numbers can say.** For one field on the max-power day (260921 FOV2) the afferent event rate after 2 Hz stimulation at maximum power is about twice what it was before at the same behaviour (× 2.1, 95 % CI 1.1–4.7), the increase is spread over the fibres of the field rather than confined to a few sources, the events are longer and lower in amplitude, and the activity is no longer coupled to arousal. That is a consistent picture, and it is the template for what an effect looks like in this preparation. It is also n = 1: one animal, one field, one day, with the lower end of the interval near 1. FOV1 of the same day points the same way (× 2.7 on running frames) on 3 against 8 events.

**What they cannot say.** The sham and 3.3 mW days do not test the hypothesis either way: they compare a resting mouse with a running one, or the reverse, and have 0–10 events per condition. Their ratios should not be averaged with the max-power result, and a dose relation (sham < 3.3 mW < max power) cannot be claimed from them — the 3.3 mW day may have had no effect or an effect that is invisible at 1–10 events. The sham day cannot serve as the drift reference either: the pre/post design cannot separate a step at the stimulation from a trend over the session (on the max-power day the rate rose within each condition too, run 004 > 002 and 008 > 006, and so did the event duration; Effect 8 puts the step across the stimulation itself at × 0.86 against within-condition changes of × 3.2–3.7), and a sham day can only show such a drift if its baseline carries events, which 260917 did not (0 pre events in both fields). The baseline rate of the same FOV also differs by more than an order of magnitude between days (FOV2: 0, 1 and 51 pre events on the three days, with the mouse running on two of them), so absolute rates are not comparable across days; the field, depth, expression or imaging settings may differ. The event-level p values in the pipeline output treat events as independent and are descriptive; the planned group statistic is one value per animal (or per FOV, stated as such), compared pre→post across animals or against sham-stimulated animals.

**Detection.** FOV1 found 17 events in its 12 runs against 214 in FOV2's, with the same threshold (thrARScl 2) but a larger minimum size (100 vs 50 px). Whether FOV1 is quieter or just under-detected can be checked by re-running its four 260921 runs with the FOV2 preset (one preset per FOV must then be kept for pre and post alike, as the batch enforces). If the FOV1 fields are simply dimmer or have fewer fibres, they may not be usable for this comparison. The changed event properties after stimulation (longer, lower) can be partly a merging artefact at higher event density; scrolling through the long post events of run 008 in `cfu_curate`, or in the AQuA2 GUI, will show whether they are single transients or chains.

**Housekeeping.** The max-power day's FOV1 maps are missing because runs 001 and 003 have no per-event location files; the current `run_cfu_batch.m` writes them for every run, so a re-run of that script followed by Step 6 (with the updated `spatial.py`) completes the set. All numbers in this report are from `aqua_analysis\summary.txt`, `experiments_summary.csv` and `run_summary.csv` in the animal's folder; the figures are the pipeline's own outputs (effects, overview, adjusted, maps, cfu_map, pupil, speed_tuning and timecourse figures) plus three summary figures built from `run_summary.csv` and `experiments_summary.csv`. The HTML version (`NavRearChrTIGRE-2608-1_report.html` in the same folder) has every figure at full resolution with zoom, including all six experiments' pipeline figures in an appendix.

**Next.** The design changes that would make the next animal interpretable are cheap: record pre and post in the same behavioural state (either habituate until the mouse runs throughout, or wait for stillness in both), and make sure the sham day's fields have a baseline rate comparable to the stimulation days' (same fields, same imaging settings, enough events) so that it can show the within-session drift. With two or three animals recorded that way, the per-animal pre→post ratios for sham, low and high power can be compared, which is the statistic this pipeline was built to feed. For the time course, note the clock time at which the stimulation ends so that the post runs can be placed in minutes after it, and consider a longer post period (or a second post block an hour later) if the build-up seen on the max-power day is to be followed to its end.
