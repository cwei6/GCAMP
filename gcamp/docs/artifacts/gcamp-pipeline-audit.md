# GCaMP Pipeline Audit

Source: https://claude.ai/artifact/5kyEZDibv8pRUuMjXmptUS (claude.ai artifact, last updated 2026-08-30)

A code audit of the two-photon GCaMP analysis pipeline (GCAMP_pipeline_reworked_20260829) finding that Steps 5, 7 and 10 pair PRE and POST cells by matrix row index instead of ROICaT track ID (171 of 177 tracked neurons mispaired on the sample map), plus 16 further statistical, pupil and data-contract findings and an ordered fix plan.

*Code audit · Two-photon GCaMP analysis pipeline*

## The row index is not the neuron

One defect explains most of what you are seeing. Steps 5, 7 and 10 pair PRE and POST cells by their position in the matrix instead of by their ROICaT track ID. On your sample map, that mispairs 171 of 177 tracked neurons.

- **Reviewed:** GCAMP_pipeline_reworked_20260829 · 34 files · 32,262 lines
- **Data:** Navcontrol-2504-1 D1 250602 _001/_002
- **Method:** static analysis + direct verification against your data files

### Verdict at a glance

| Value | Meaning |
|---|---|
| **96.6%** | of tracked ROI pairs mispaired in Step 5 |
| **0** | multiple-comparison corrections in 32k lines |
| **0.5%** | of frames yield a pupil value at the default threshold |
| **177** | real tracked pairs the ROICaT map does give you |

> Note on attachments: the audit refers to three attached MATLAB files (`GCAMP_pair_tracked_rois.m`, `GCAMP_preflight_contract_check.m`, `GCAMP_plot_tracked_heatmaps.m`). They were attachments to the original session, not part of the artifact page, and are not reproduced here.

## 01 · The blocking defect

Your ROICaT run is fine. The tracking map is well formed, the labels align to the full Suite2p ROI list, the `-1` sentinel is correct, and the curated-baseline policy did exactly what you intended. The problem is what Step 4 hands downstream and how Step 5 reads it.

Step 4's tracked-variant writer keeps ROIs with `keep = ismember(lbl, trackedGlobal); source_roi_index = find(keep);`. `find` returns rows in ascending *original Suite2p row order of that session*. Each session therefore emits its tracked cells in its own row order, and the ROICaT label sequences end up permuted relative to one another. Both files have the same number of rows, so nothing errors.

Step 5 then does `n = min(numel(preInfo.roiType), numel(postInfo.roiType))` and compares `preX(1:n)` against `postX(1:n)` — row *i* to row *i*.

**Figure (diagram): row-index pairing vs. label pairing.** Two columns show the first five rows of the supplied map. PRE variant `_001` (701 ROI → 270 curated → 177 tracked): rows 1–5 carry labels 0, 1, 2, 3, 4. POST variant `_002` (596 ROI → all eligible → 177 tracked): rows 1–5 carry labels 1, 0, 10, 8, 18. Dashed arrows show the current pairing (row *i* to row *i*); solid green arrows show the correct pairing by equal ROICaT label (PRE row 1/label 0 → POST row 2; PRE row 2/label 1 → POST row 1). Labels 2, 3 and 4 match POST rows 30, 11 and 7, outside this window.

*Caption:* The first five rows of the supplied map, verbatim. Row-index pairing joins PRE label 0 to POST label 1 — two different neurons. Across the whole map, 171 of 177 tracked pairs are joined to the wrong cell. Only 6 land correctly, by coincidence.

The PRE sequence happens to be monotonic (`0,1,2,3,4,…`) because ROICaT's sequential-Hungarian labels the reference session in order. The POST sequence is not (`1,0,10,8,18,23,…`). That asymmetry is exactly what makes the bug invisible: the PRE side looks perfectly sensible on inspection.

> **What this means for your results**
>
> Every PRE→POST transition matrix, gain/loss/switch fraction, ΔP value and chronic trajectory in Step 5 is computed on shuffled pairs. That is statistically close to a permutation null — which is why effects appear, disappear and reverse between reruns, and why they would have looked obviously wrong as a heat map. Your PI's instinct was right, and for a more specific reason than they knew.

Affected sites:

```text
GCAMP_step5_core.m:803      computeOnePairMetrics_main       — acute + all four baseline pairs
GCAMP_step5_core.m:983      computeOnePairMetrics_suppEver   — supplementary ever-responsive
GCAMP_step5_core.m:1063     buildOneAnimalChronicTrajectory  — pMat(roi, day) across days
GCAMP_step7_..._core.m:2358 local_infer_track_key(meta, scenario, r) — r is the row index, so TrackKey
                             groups unrelated cells under (1|TrackKey)
GCAMP_step10_roi_qc_heatmap.m:213 deltaC = profPostC - profPreC — row-wise; roiTrackID is
                             stored for annotation but never used to reorder
```

Step 3 already computes and saves `roiTrackID` for every ROI. Step 5 simply never loads it — `loadStep3_ROIinfo` reads ten fields and `roiTrackID` is not among them. Steps 8 and 11, and Step 7's change-class branch, *do* use it, so the fix is to make the broken sites behave like the working ones.

## 02 · Where the damage is

Reachability matters here — you mentioned not every file is in play. This is the live call graph from the entry points, with each stage's status.

**Figure (pipeline strip): status of each live stage.**

| Step | Stage | Status |
|---|---|---|
| Step 1 | Quadrature → loco | warning: unchecked clock |
| Step 2_1 | Extract + dF/F | warning: silent truncation |
| Step 2_2 | Photobleach QC | ok: QC-only, fine |
| Step 3_1 | Responsiveness | bad: test degenerate |
| Step 4 | ROICaT + compile | ok: reader correct |
| Step 5 | Pair metrics | bad: mispaired |
| Step 6 | Signal properties | warning: unseeded RNG |
| Step 7 | Mixed models | bad: TrackKey wrong |
| Step 8 | Unbinned curves | ok: uses track IDs |
| Step 10 | QC heat maps | bad: rows mispaired |
| Step 11 | Balanced rest | ok: uses track IDs |

### Files nothing calls

Legitimate user-invoked entry points aside (`step4_run_multiCondition_withQC`, `step5_compile`, `step6…`, `step7_mixed_model_analysis`, `step8_final_figure_run`, `step10_unified_qc`, `step11…`, `GCAMP_trace_tracked_roi`), these five are referenced by nothing and appear superseded. I did not audit them:

- `GCAMP_step5_replot_refined.m` — 1,116 lines, duplicates Step 5 plotting
- `GCAMP_step5_transition_acute.m` — 584 lines
- `GCAMP_step7_plot_focus_easyplot.m` — 643 lines
- `GCAMP_step4_backfill_dlc_filtered_csv.m` — superseded by Step 4's CSV auto-discovery
- `DEBUG_WHY_LOCO_BOUTS_NOT_THERE.m` — a debugging script; keep, but it is not pipeline

Deleting them removes 2,400 lines of code that will otherwise drift out of sync with the fixes below.

## 03 · Does the ROICaT output fit the pipeline?

Mostly yes. I loaded `Navcontrol25041_D1_250602_001_002_tracking_map.mat` and checked it field by field against Step 4's reader.

| Contract | Expected by Step 4 | Actual in your map | Verdict |
|---|---|---|---|
| Top-level variable | `clusters` | `clusters` | pass |
| Session names from path | folder above `suite2p` | `…_250602_001 / _002` | pass |
| `labels_bySession` shape | 1×N cell, one vector per session | 1×2 cell, 701 and 596 | pass |
| Labels span full ROI list | = rows of `Fall.mat` | 701 / 596, matches `stat.npy` | pass |
| Untracked sentinel | `-1` | `-1` | pass |
| Label uniqueness per session | exactly once | 177 labels, no duplicates | pass |
| Curated-PRE policy | `iscell==1` only in `_001` | `curated_001_iscell_only` | pass |
| Labels are zero-based | label 0 is a real cluster | label 0 present in both | see F9 (high) |
| `iscell` source of truth | Step 4 reads `Fall.mat` | ROICaT read `iscell.npy` | see F8 (medium) |
| `labels_bool_tracking_input` | logical array | a 899-char text blob | see F16 (low) |

> **Yield**
>
> 270 curated cells in the D1 baseline, 177 of them tracked into the post recording — a 65.6% match rate. That is a reasonable yield for a same-day pair and gives you a usable *n*. The tracking is not your problem.

## 04 · Findings, ranked

### F1 · BLOCKER · ROIs paired by row index, not track ID

Detailed in §1. Affects Step 5 (all pair families and chronic trajectories), Step 7's main modelling path, and the Step 10 heat map.

**Fix:** Use `GCAMP_pair_tracked_rois.m` (attached) at all five sites. It intersects on `roiTrackID`, treats label 0 as valid, drops repeated labels, and refuses to fall back to index pairing.

### F2 · HIGH · The trial-level permutation test cannot go below p = 0.048

`trialPermTestOneVsPool` compares one stim value against `trialNullSubsetN = 20` null values by permuting which of the 21 is "the stim". Because `allVals(i) − mean(rest)` is monotonic in `allVals(i)`, the test is exactly "is the stim window the largest of 21?" — the floor is p = 1/21 = 0.0476, against a threshold of 0.05.

Two consequences. It is a binary max-rank test with almost no power gradient, so effect size is invisible to it. And with `B = 1000` Monte-Carlo draws the p-value at that floor has a standard deviation of about 0.007, so a genuine "max" trial fails the 0.05 cut a substantial fraction of the time, at random, on every rerun. Set `trialNullSubsetN` to 18 and *no trial in the dataset could ever be significant*.

**Fix:** Compute the rank analytically instead of by Monte Carlo (`p = (1 + #{null ≥ stim}) / (K+1)`), raise `trialNullSubsetN` to at least 99 so the floor is 0.01, and stop treating single-trial significance as the unit of evidence — the ROI-level test is the one that should carry the claim.

### F3 · HIGH · No RNG seed, so results change on every run

`rng(...)` is called exactly once in the entire pipeline, in Step 8. Steps 3, 6, 7 and 10 all call `randperm` — for the permutation nulls, for the random 20-window null subset, and for trial subsampling. Rerunning Step 3 on identical inputs produces a different set of responsive cells.

**Fix:** Add `rng(opts.seed, 'twister')` at the top of every step that randomises, default the seed to a fixed value, and record it in the output struct. This alone will make "not robust" measurable rather than anecdotal.

### F4 · HIGH · No multiple-comparison correction anywhere

A search for `fdr`, `bonferroni`, `holm`, `benjamini` and `mafdr` across all 32,262 lines returns nothing. `roiIsResp` is `p < 0.05` per ROI, applied to several hundred ROIs per session. Roughly 5% of your "responsive" cells in every session are noise, and that floor propagates into every fraction, transition matrix and GLMM outcome downstream.

**Fix:** Apply Benjamini–Hochberg within session across ROIs to `pPermExc_roi` and `pPermInh_roi`, store both the raw and adjusted p, and drive `roiIsResp` from the adjusted value. Report the q threshold in the figure legends.

### F5 · HIGH · The ROI-level null pool double-counts the same windows

In the deep-still branch of the ROI-level test, `baseVals` is built by looping over trials and appending every null window whose pupil Z falls inside the caliper. Trials with similar pupil size match overlapping window sets, so the same null window is appended many times.

```text
step3_1_trial_responsiveness_boutBaseline.m:702–718
  for ev = stimIdx        mask = |pz − trialPupilPreZ(ev)| <= caliper
                          baseVals = [baseVals, pv(r, mask)];   ← duplicates
```

`permTestMeanDiff` then treats those duplicates as independent observations. The permuted null distribution of the mean difference narrows roughly as 1/√N, so with ~60 stim values against a few thousand duplicated base values the p-values are strongly anti-conservative. Combined with F4, this is the most likely source of "many positive results" that do not replicate.

**Fix:** Deduplicate by null-window global index before testing (`unique(matchedIdx)`), and weight or cluster-resample by window rather than pooling raw values.

### F6 · HIGH · Runaway dF/F from non-positive baselines is flagged but never excluded

Step 2_1 computes `Fcorr_F0_bad` for ROIs whose maximin baseline is non-finite or ≤ 0 — which happens whenever `F − 0.7·Fneu` goes negative, and is more likely for the dim sourcery detections you merged in. It saves the flag and nothing ever reads it. Meanwhile `Fcorr_dff_maximin_t = (Fc − F0) ./ max(F0, 1e-6)` gives those ROIs dF/F values on the order of 10⁶.

A handful of such ROIs will dominate any population mean, any permutation statistic, and the colour scale of any heat map you draw.

**Fix:** Propagate `Fcorr_F0_bad` into Step 3 and exclude those ROIs from `roiIsTestable`. Print the count. The attached heat-map function already applies this guard.

### F7 · HIGH · The pupil metric is confounded with tracking confidence

Covered in detail in §5. The short version: at the shipped threshold the metric is undefined on 99.5% of frames, and where it is defined its value scales with how many body parts DLC happened to be confident about.

### F8 · MEDIUM · Step 4 reads `iscell` from a file ROICaT never saw

ROICaT loaded `suite2p/plane0/iscell.npy`. Step 4 loads `iscell` out of `Fall.mat`. The Suite2p GUI updates `iscell.npy` on every save but only rewrites `Fall.mat` when you explicitly export. Since you have just finished a large manual curation pass, a stale `Fall.mat` is a live risk — and it would silently change which cells the "curated PRE" filter keeps.

**Fix:** The map already records the answer: `clusters.tracking_n_roi_bySession(1)` is 270 for this pair, so `nnz(Fall.iscell(:,1)==1)` must also be 270. The attached preflight script asserts this for every session and diffs `Fall.mat` against `iscell.npy` directly.

### F9 · MEDIUM · `id > 0` silently discards ROICaT label 0

ROICaT labels are zero-based; `-1` is the not-tracked sentinel. Three places filter with `> 0` instead of `>= 0`, dropping one genuine tracked cell from every pair.

```text
GCAMP_step8_unbinned_curve_analysis.m:872   preGood = isfinite(preTrack) & (preTrack > 0)
GCAMP_step11_balanced_rest_magnitude.m:1146 okPre   = isfinite(idPre)  & idPre  > 0
GCAMP_step7_..._core.m:2647                 goodA   = isfinite(idA)   & idA    > 0
```

**Fix:** Change to `>= 0`. Also remove the silent `(1:n)` index fallbacks in Step 7's change branch and Step 11 — an unmatched pair should be skipped loudly, not paired wrongly.

### F10 · MEDIUM · The F0 window has the same period as the stimulus

The maximin baseline uses `win_baseline_sec = 60`, and your trial block is 2 s ON + 58 s OFF = 60 s. A rolling minimum whose window exactly matches the stimulus period produces an F0 that is itself phase-locked to the stimulus, so dF/F carries a systematic, stimulus-locked offset that is not a response. The null pool is drawn from bins 31–60 and the stim statistic from bins 1–2, so the comparison sits at two fixed, different phases of that artifact.

**Fix:** Set the baseline window to something incommensurate with the trial period — 90 s or 150 s — and verify by recomputing responsiveness on circularly shifted stimulus times, which should return a null result.

### F11 · MEDIUM · The 60 "one-second" bins are proportional, not fixed

`make_edges = @(s,e,N) round(linspace(s, e+1, N+1))` always yields exactly 60 bins over whatever the ON+OFF block happens to span. Bin *k* is one second only if the interval is exactly 60 s. For the last event, `off_end` is clipped to `nFrames`, so its 58 OFF bins are compressed into however many frames remain — and if the tail is short, individual bins come back `NaN` or duplicate frames. Averaging bin *k* across events then mixes different post-stimulus times.

**Fix:** Bin on fixed frame counts (`round(1*fps)` per bin) and mark trailing bins as missing rather than stretching. Drop the final event from block averages unless its full window was recorded.

### F12 · MEDIUM · Modality lengths are reconciled by silent truncation

Step 1 assumes one quadrature sample per imaging frame — `Vinst = diff(quad) · step · frameRate`, `Tinst = (0:N−1)/frameRate` — and never checks the length against `Fall.mat`. Step 2_1 then does `minLen = min(nFrames, numel(locostate))` and crops whichever is longer, with no warning printed. If the encoder is on a different clock, this produces a plausible-looking but time-warped locomotion trace, and every phase classification built on it is wrong.

Your sample quad file has 58,500 samples, which at 15.49 Hz is 62.9 minutes — consistent with ~63 trials of 60 s, so the assumption is probably right. But nothing in the code verifies it.

**Fix:** Assert that the lengths agree to within 1% and error otherwise. The attached preflight does this per session.

### F13 · MEDIUM · Day labels can disagree between ROICaT and Step 4

Both derive D-labels as calendar days from the first available session, but they disagree on what counts as available. The tracker requires `stat.npy` and `ops.npy`; Step 4's `findSessions` additionally requires the event file `<name>.mat`. A session with imaging but no event file shifts the anchor in one tool and not the other, and Step 4's `sameDaySet` check will then silently skip the pair.

**Fix:** Put the D1 anchor date in an explicit per-animal metadata file that both tools read, rather than inferring it from whichever folder sorts first.

### F14 · MEDIUM · A curated PRE cell can be matched to a duplicate POST mask

By design the post recording is uncurated, so every sparsery *and* sourcery detection is eligible for matching — including the near-duplicate pairs your merge produced. Hungarian matching is 1:1, so it picks one of the two essentially arbitrarily. The chosen mask may have different pixels and a different neuropil ring than its twin, which changes `F − 0.7·Fneu` for that cell in POST but not in PRE.

This does not invalidate the design — you deliberately ignore POST `iscell`, which is reasonable — but it adds variance to exactly the comparison you care about, and it can bias it if sourcery masks are systematically larger.

**Fix:** After tracking, compute mask IoU and area ratio between each matched PRE/POST pair from `stat`, and drop or flag pairs below an IoU threshold. Report the distribution — it is also a good sanity figure for the PI.

### F15 · MEDIUM · Pupil frame alignment falls back to a silent 1:1 map

`resolveFrameAlignment` accepts DLC frame IDs only if they fit `[0, nFrames−1]` or `[1, nFrames]` exactly. On any other relationship it returns `mode = 'sequential'` and maps CSV row *k* to imaging frame *k* — with no error and no warning. If the eye camera is not frame-locked to the microscope, the pupil trace is stretched across the recording and every pupil-matched null is wrong.

Separately, `copyDlcToVariant` takes the *newest* matching CSV. With an 8-point and a 9-point export sitting in the same folder, which model gets used depends on file timestamps.

**Fix:** Require `abs(nRows − nFrames) ≤ 1` and error otherwise; if the camera is genuinely on its own clock, resample by timestamp instead. Pin the DLC model by an explicit scorer string, not by recency.

### F16 · LOW · The tracker serialises two fields as text

In `ROICaT_Batch_Tracker.py`, `_safe_to_numpy` falls through to `np.array(str(x))` for anything it does not recognise. `labels_bool_tracking_input` is a scipy sparse matrix, so it lands in the `.mat` as the literal string `"<Compressed Sparse Row sparse matrix of dtype 'bool' with 866 stored elements…"`. `labels_tracking_input` is stored as an object array rather than `int64`. Nothing downstream reads either field, so this is cosmetic — but the file currently misrepresents its own contents.

**Fix:** Add a sparse branch (`if scipy.sparse.issparse(x): return x.toarray()`) and cast label arrays with `np.asarray(..., dtype=np.int64)` before saving.

### F17 · LOW · Two-sided p is taken as `min(pExc, pInh)`

`pVals(:,ev) = min(pExc, pInh)` without doubling. Where this feeds a two-sided claim it is anti-conservative by a factor of two. It is used for reporting rather than for the responsiveness call, so the impact is small.

## 05 · The pupil data

> **First, a mismatch**
>
> The CSV you sent as the 9-point sample has eight body parts — `pupil_1` through `pupil_8` — and the scorer string is the same `shuffle4_snapshot_best-130` model as before. There is no centre point in it. Either this is an old export, or the 9-point model writes under a different filename. `GCAMP_pupil_9point_metrics` handles this gracefully, but it means the 9-point upgrade is not being exercised.

I ran the reader's exact geometry over the file. With the shipped default `likelihoodThresh = 0.80`:

| Likelihood threshold | Mean visible points | Frames with ≥5 points | Frames with all 8 |
|---|---|---|---|
| 0.30 | 5.94 | 73.1% | 55.6% |
| 0.50 | 3.29 | 36.5% | 15.7% |
| 0.70 | 1.01 | 6.8% | 0.7% |
| **0.80** ← shipped default | 0.33 | **0.5%** | 0.0% |

The polygon fallback needs five points. At the default threshold it fires on 0.5% of frames, so `pupilAreaZAligned` is `NaN` almost everywhere, every deep-still trial fails the `pupilMinValidFrac = 0.50` gate, and pupil matching quietly does nothing. The median per-part likelihood in this session ranges from 0.32 to 0.70 — every body part sits below the threshold most of the time.

Lowering the threshold does not fix it, because the metric itself is not scale-invariant. Grouping frames by how many points were visible:

| Visible points (thr 0.50) | Frames | Median polygon area | Implied radius | Circle-fit radius |
|---|---|---|---|---|
| 5 | 2,176 | 32.3 | 3.21 | 3.94 |
| 6 | 3,405 | 55.1 | 4.19 | 4.86 |
| 7 | 6,584 | 64.1 | 4.52 | 4.96 |
| 8 | 9,206 | 150.0 | 6.91 | 7.45 |

Measured area is 4.6× larger when eight points are visible than when five are. Some of that is real — a small pupil genuinely is harder to track — but a large part is pure geometry: an inscribed pentagon covers 2.38*r*² where an octagon covers 2.83*r*², so dropping vertices shrinks the number with no change in the pupil. The correlation between visible-point count and measured area is **r = 0.75**.

Two further problems in the same function. `geometryModeRaw` can switch between centre-radius (area = π*r*² = 3.14*r*²) and polygon (≈2.83*r*² for an octagon) from frame to frame within one session, which puts a ~10% step into the trace at every switch and more when fewer vertices are visible. And 2.0% of frames contain a negative coordinate — DLC extrapolating off the image — which nothing rejects.

> **Recommendation**
>
> Replace both geometry modes with a single algebraic least-squares circle fit over whatever boundary points are visible. It needs only three points, so coverage rises from 36.5% to 43.5% at the same threshold, it is invariant to how many points survive, and the DLC centre becomes a cross-check rather than a dependency: if the fitted centre and the predicted centre disagree by more than a fraction of *r*, reject the frame. Then report per-session the fraction of frames retained, and treat pupil as unavailable for a session below some floor rather than letting it silently degrade.

## 06 · Why the figures look the way they do

Your PI's objection is quantifiable. Across the pipeline there are roughly 91 `bar`/`errorbar` calls against 19 `imagesc` calls. `GCAMP_step5_plot.m`, the main figure generator, has 37 bars and 2 heat maps. Of the 19 `imagesc` calls, 8 are in `GCAMP_step10_roi_qc_heatmap.m` — which is the right tool, but it is only reachable through `GCAMP_step10_unified_qc`, and it draws its PRE/POST difference on row-index pairs (F1), so it has been showing you a scrambled matrix.

There is also a display trap worth naming before you build the replacement. Sorting ROIs by their PRE response and then plotting POST for the same trials produces a peak in PRE that decays in POST by regression to the mean alone, whether or not anything changed. It looks exactly like the effect you are hoping to find.

The attached `GCAMP_plot_tracked_heatmaps.m` sorts on odd PRE trials and displays even ones, uses one row order for all three panels, pairs on track ID, and drops the runaway-dF/F ROIs from F6. Its panels are PRE, POST, POST−PRE, population mean ± SEM, per-ROI amplitude scatter against unity, and the sorted per-cell change.

## 07 · What to do, in order

1. **Run the preflight on one animal.** Drop `GCAMP_preflight_contract_check.m` on the path and run it against a single animal folder. It will tell you within a minute whether `Fall.mat` is stale relative to `iscell.npy` (F8), whether the frame clocks agree (F12), and it will print the exact mispairing percentage for each of your maps (F1). Do this before regenerating anything.
2. **Fix the pairing at all five sites.** Add `GCAMP_pair_tracked_rois.m` and route Step 5's three pairing functions, Step 7's `local_infer_track_key`, and Step 10's heat map through it. Have `loadStep3_ROIinfo` return `roiTrackID`. Nothing downstream of Step 4 is trustworthy until this lands.
3. **Make the run reproducible.** Seed every RNG (F3), fix the trial test's p-floor (F2), deduplicate the ROI-level null pool (F5), apply BH-FDR (F4), and exclude bad-F0 ROIs (F6). These five together determine whether a rerun gives you the same answer.
4. **Look at the heat maps before rerunning the statistics.** Run `GCAMP_plot_tracked_heatmaps.m` on two or three pairs. If the effect is real it will be visible here, and this is the figure to take to your PI. If it is not visible here, no model will rescue it — and that is worth knowing before you spend a week on Steps 6–11.
5. **Settle the pupil pipeline separately.** Find the actual 9-point CSV, confirm the ninth part's name, switch to the circle fit, and set the likelihood threshold from the model's own calibration rather than from 0.80. Until then, run the neural analysis with pupil matching explicitly disabled rather than silently inert.
6. **Then revisit the baseline window and binning.** F10 and F11 change the numbers but not the direction of the analysis. Fix them once the pipeline is reproducible, so you can measure what they change.

## 08 · What I could not verify

There is no MATLAB in this environment, so everything above is static analysis plus direct inspection of the four data files you sent. These are the questions whose answers would change what I wrote.

### Is the eye camera frame-locked to the microscope?

Your DLC CSV has 58,501 rows and your quadrature file 58,500 samples — consistent with one sample per imaging frame at 15.49 Hz over a 62.9-minute session. If that is right, F12 and the alignment half of F15 are latent rather than active. If the camera runs on its own clock, both are live and the pupil trace is currently stretched.

### Where is the 9-point DLC export, and what is the ninth part called?

The file I received is 8-point. If the ninth part is not named `pupil_center`, set `pupilCenterBodypartName` explicitly — the geometric inference fallback requires nine parts and a clear separation in radial score, and it will silently give up rather than guess wrong.

### Is the inter-stimulus interval fixed at 60 s?

The whole 60-bin representation assumes it. If the ISI is jittered, bin *k* means a different post-stimulus time on every trial and F11 becomes serious rather than an edge case.

### Did you re-export `Fall.mat` after the manual curation pass?

Saving in the Suite2p GUI updates `iscell.npy` but not `Fall.mat`. If you only saved, ROICaT tracked your curated set and Step 4 is filtering on the pre-curation one. `nnz(Fall.iscell(:,1)==1)` should be 270 for the D1 baseline.

### How many frames does `Fall.mat` actually contain?

I inferred ~58,500 from the quadrature file, but I never saw a `Fall.mat`. If it is materially different, Step 2_1 has been cropping one modality against the other on every session.

---

*Audit of GCAMP_pipeline_reworked_20260829 · findings verified against Navcontrol25041 D1 250602 tracking map, quadrature and DLC exports · 30 August 2026*
