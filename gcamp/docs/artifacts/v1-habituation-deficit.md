# V1 Habituation Deficit

Source: https://claude.ai/artifact/TRGTrkvShrt3DKcFN7rjJs (claude.ai artifact, last updated 2026-09-02)

A findings memo for the optostimulation cohort (run 260901_2114) arguing that the 2 Hz optogenetic protocol drives V1 L2/3 into a non-habituating, migraine-like regime (d ≈ 2 with animals as the unit), while the 5 Hz arm's cell-level effects do not survive animal-level permutation and the spontaneous-activity suppression still awaits its animal-level test.

*Findings memo · optostimulation cohort*

## A habituation deficit in V1

Across 14 mice, optogenetic stimulation left visual responses intact while control animals' responses weakened — the signature of a failure to habituate, which is the cardinal interictal finding in migraine. An animals-as-unit audit added after the first draft sharpens this considerably: the effect is carried by the 2 Hz arm, and two claims that looked strong at cell level do not survive.

- **Run:** 260901_2114 (deduplicated)
- **Animals:** 4 CTRL · 4 EXP-2 Hz · 7 EXP-5 Hz
- **Cells:** 8,695 tracked
- **Area:** V1 L2/3, jGCaMP8s

## The result in one picture

Two quantities move independently after the stimulation hour: the size of the visual response, and the spontaneous activity the cell shows before each stimulus. Plotting one against the other separates the three groups cleanly, and the axis a group moves along is more informative than how far it moves.

**Figure (scatter, one point per group):** x-axis = change in spontaneous activity (log rate, −0.10 to +0.10); y-axis = change in evoked response (per-cell SD, about −0.05 to +0.05); upper-left quadrant labelled "deficient habituation", lower-right "normal habituation". CTRL (no opsin, 4 animals): evoked −0.065 SD, spontaneous +0.073 log rate (down and right). EXP 2 Hz (4 animals): evoked +0.036 SD, spontaneous +0.003 log rate (up, no spontaneous change). EXP 5 Hz (7 animals): evoked −0.011 SD, spontaneous −0.085 log rate (left, little evoked change).

*Caption:* Each point is one group's state-adjusted PRE→POST change, from the mixed models (17,376 cell×session observations). Controls sit in the habituation quadrant: responses down, spontaneous activity up. The 2 Hz group sits opposite on the response axis; the 5 Hz group sits opposite on the spontaneous axis. All values are tabulated in Claims 1 and 2 below. Note that only the 2 Hz vertical displacement is confirmed with animals as the unit — the 5 Hz position on the horizontal axis rests on the one test still outstanding.

## Claim 1 — Stimulated animals fail to habituate

**Status: Well supported**

> Over the stimulation hour, control responses weakened on every measure. In both stimulated families they did not — and at 2 Hz they grew.

| Measure | CTRL | 2 Hz vs CTRL | 5 Hz vs CTRL |
|---|---|---|---|
| Evoked ΔF/F (per-cell SD) | −0.065 | +0.101 [+0.057, +0.146] | +0.054 [+0.020, +0.087] |
| Evoked rate (deconvolved) | −0.079 | +0.126 [+0.078, +0.175] | +0.057 [+0.021, +0.094] |
| Trial-to-trial reliability | −3.05 pp | +3.26 pp [+2.5, +4.0] | +1.79 pp [+1.2, +2.4] |
| Responsive cells retained | 56% | 49% | 65% |
| Silent cells recruited | 12% | 16% | 20% |

*State-adjusted mixed models (running and pupil regressed out); 95% CI in brackets. Turnover rows are observed counts, no model.*

The observed data says the same thing without any model: session by session, the median response strength fell in controls (`−0.64 z`) and rose at 2 Hz (`+0.39 z`). Treating each animal as one observation — the conservative standard — an exact permutation test over all 70 possible 4-vs-4 splits gives `p = 0.057` for amplitude (`d = 2.03`) and `p = 0.057` for reliability (`d = 1.93`). Both sit one swap away from the smallest p that 4 versus 4 animals can produce at all (`0.029`). The effect size is large and the separation is nearly complete; what limits the p-value is the number of animals, not the strength of the effect.

**Why this is the interesting framing**

Control behaviour here is textbook. In mouse V1 L2/3, seven days of repeated visual stimulation produces habituation carried by a *rise* in pre-stimulus activity (+35–40%) while evoked responses hold, collapsing the signal-to-noise distinction between baseline and stimulus; optogenetically raising pre-stimulus activity was sufficient to induce habituation on its own ([Cell Reports, 2022](https://www.cell.com/cell-reports/fulltext/S2211-1247(22)00515-0)). Your controls reproduce that: spontaneous `+8%`, evoked down, reliability down.

Deficient habituation of evoked responses — "potentiation instead of habituation" — is described as a neurophysiological hallmark of migraine between attacks, and it normalises just before and during an attack ([J Headache Pain, 2009](https://link.springer.com/article/10.1007/s10194-008-0095-x)). Your stimulated animals show the deficit; your controls show intact habituation. The mouse mechanism and the human phenotype point at the same variable.

**What would break this:** the 2 Hz family starts `0.29 SD` below controls before any light (`p ≈ 0.02`), so a floor effect could inflate a "preserved" response. The fix is a re-run restricted to cells already responsive at PRE. Note this objection does *not* touch the reliability or turnover rows, and does not touch 5 Hz.

## Claim 2 — Spontaneous activity is suppressed, ordered by light dose

**Status: Strong at cell level · animal test outstanding**

> Pre-stimulus activity rose in controls and fell in stimulated animals, with the 5 Hz family showing roughly twice the 2 Hz effect.

| Group | Duty cycle | Light on per hour | Spontaneous activity change |
|---|---|---|---|
| CTRL | — | — | +8% |
| EXP 2 Hz | 2–3% | 72–108 s | −7% vs CTRL |
| EXP 5 Hz | 5% | 180 s | −15% vs CTRL |

*Duty cycle is arithmetic from the protocol (2 Hz × 10–15 ms; 5 Hz × 10 ms). Estimates are state-adjusted mixed models on 17,406 cell×session means.*

Two things argue for this being real. It survives state adjustment untouched, so it is not running or arousal. And unlike the response measures, the groups have *no* baseline offset in spontaneous activity before stimulation (`−0.02, p = 0.93`), so the floor argument that shadows Claim 1 cannot apply here. The ordering also lines up with dose: at matched pulse width the 5 Hz protocol delivers 2.5× the photostimulation per unit time and shows roughly twice the suppression.

**Why the status is amber, not green.** This is the one headline measure that has *not* yet been tested with animals as the unit — the pipeline does not currently export a per-animal spontaneous-activity change. Every other 5 Hz claim that looked overwhelming at cell level evaporated under that test (see the audit below), so the cell-level `p` here should not be quoted until the animal-level version exists. That single export is the highest-priority analysis in the project.

A second caveat: deconvolved rate is in arbitrary units and sensitive to the fluorescence baseline, so a systematic change in F0 could in principle masquerade as reduced spontaneous activity. Controls received the same light and moved the other way, which argues against it, but a ΔF/F-based spontaneous measure plus an F0 audit would close the door.

**Mechanistic anchor**

Sustained suppression of spontaneous firing with a shift toward inhibition is the documented aftermath of cortical spreading depolarisation: spontaneous potentials and action potentials are abolished during the event and remain significantly reduced beyond 60 minutes, through reduced release probability presynaptically and increased inhibitory drive postsynaptically ([Cerebral Cortex, 2017](https://academic.oup.com/cercor/article/27/2/bhv328/3056266)).

Independently, Schoenen's preactivation model explains the migraine habituation deficit by *reduced* baseline activation of sensory cortex between attacks — a low starting point leaves a wide range before the response ceiling, so amplitudes do not decline with repetition ([J Headache Pain, 2009](https://link.springer.com/article/10.1007/s10194-008-0095-x)). Claims 1 and 2 are therefore not two findings but one: suppressed preactivation is the proposed *cause* of the habituation deficit, and you measured both in the same cells.

## The audit that reorganised this memo — Which claims survive with animals as the unit

**Status: Method note**

> A mixed model with four animals per arm estimates between-animal variance from four numbers. Where that is checked directly by permuting animal labels, the 2 Hz results hold and the 5 Hz results do not.

| Claim | Cell-level model | Animals as unit (exact permutation) | Verdict |
|---|---|---|---|
| 2 Hz — evoked amplitude | +0.101, p ≈ 2e−7 | p = 0.057, d = 2.03 | holds (at the 4v4 floor) |
| 2 Hz — reliability | +3.26 pp, p < 1e−4 | p = 0.057, d = 1.93 | holds |
| 5 Hz — evoked amplitude | +0.054, p = 0.0016 | p = 0.99, d = 0.01 | does not survive |
| 5 Hz — reliability | +1.79 pp, p < 1e−4 | p = 0.86, d = 0.18 | does not survive |
| 5 Hz — chronic baseline drift | −0.035 SD/day, p < 1e−4 | p = 0.81, d = −0.17 | does not survive |
| Spontaneous activity (both arms) | −7% / −15%, p < 1e−4 | not yet computed | outstanding |

*Permutation is exact over all splits (70 for 4v4, 330 for 7v4) on the per-animal mean change. Chronic drift compares per-animal slopes of mean amplitude across baseline days. Every group's baseline declines at about the same rate: CTRL `−0.0091`/day versus EXP-5 Hz `−0.0105`/day.*

This is not a power problem. The 5 Hz arm has *seven* animals against the 2 Hz arm's four, so it has more power, not less — and it shows nothing. The cell-level significance came from counting thousands of cells within a handful of animals, which measures how consistent an effect is inside an animal, not whether animals differ. The 5 Hz acute and chronic claims should be withdrawn from the headline and reported as null pending more animals.

One animal illustrates the fragility. Dropping NavRearChr-2506-3 — whose late-day POST sessions are already flagged as near-zero-responsive — moves the 5 Hz acute result from `p = 0.99, d = 0.01` to `p = 0.09, d = 1.21`. A single animal spanning that range means the arm is not yet measuring anything stable.

## Claim 3 — Responses return to baseline more slowly after stimulation

**Status: Supported, with a caveat**

> The tail of the response — 4 to 10 s after onset, i.e. 2 to 8 s after the screen goes blank — shrinks in controls and is relatively preserved in both stimulated families.

| Measure | CTRL | 2 Hz vs CTRL | 5 Hz vs CTRL |
|---|---|---|---|
| Decay index (late ÷ peak) | −0.100 | +0.055, p = 0.017 | +0.053, p = 0.0009 |
| Late-window ΔF/F (4–10 s) | −0.013 | +0.013, p < 1e−4 | +0.010, p < 1e−4 |

*Decay index computed on responding cells only (5,390 observations); at PRE it sits near 0.20, the expected jGCaMP8s tail.*

This was the analysis added for the aura question, and it returns a real but modest signal in the predicted direction. Unlike Claims 1 and 2 it does not separate the frequencies (`p = 0.71` for the contrast), which is consistent with it being a consequence of the shared change rather than an independent axis.

**What would break this:** the decay index is conditioned on cells that responded in that session, and which cells qualify differs by group and by PRE/POST — a selection asymmetry could generate the effect. Treat it as supporting, not load-bearing, until it is recomputed on a fixed cell set defined at PRE.

## Claim 4 — The acute effect does not build across days

**Status: Not supported**

> In the previous run the 2 Hz effect appeared to grow with repeated stimulation. Removing one duplicated week erased it.

The role × group × day term for 2 Hz went from `+0.022/day, p = 0.017` to `+0.008/day, p = 0.38`. For 5 Hz it is unchanged and marginally negative (`−0.011/day, p = 0.055`), i.e. if anything converging. There is no evidence of hypersensitisation in this cohort on this measure.

What does survive is that the 2 Hz effect is largest on the naive first day (`+0.26 SD [+0.09, +0.43]` vs controls) — the only day whose confidence interval excludes zero. Read that as a design instruction for the next cohort, not as a result: choosing the best-looking day after seeing the data and then testing it would be circular.

## Nulls worth reporting

- **Wheel lock does not change the acute effect.** `role × locked = +0.032, p = 0.31` (pupil-adjusted `p = 0.35`), 4 animals. Locomotion is not driving the response result — though with this N the test is underpowered, and locking removes running without removing arousal.
- **Population turnover does not separate groups at the animal level.** Retention `51% vs 58%, p = 0.41`; recruitment `14% vs 18%, p = 0.66`. The pooled counts point the right way but four control animals span `25–71%` retention, and that spread swamps the group difference.
- **The 5 Hz arm is null on every measure that has been tested with animals as the unit** — acute amplitude `p = 0.99`, reliability `p = 0.86`, chronic drift `p = 0.81` (see the audit above). Every per-day contrast straddles zero. Whether anything survives for 5 Hz now depends entirely on the outstanding spontaneous-activity test.

## What changed when the duplicate came out

One recorded week of NavRearChr-2511-1 had been sitting in two condition folders and was counted twice. Removing it took the 2 Hz family from 9,247 to 5,041 testable cell-sessions. The audit is short:

| Quantity | Before | After | Verdict |
|---|---|---|---|
| 2 Hz acute ΔF/F | +0.112 | +0.101 | holds, CI wider |
| 2 Hz reliability | +3.60 pp | +3.26 pp | holds |
| 2 Hz spontaneous | −5% | −7% | holds, slightly stronger |
| 2 Hz effect × day | p = 0.017 | p = 0.38 | did not survive |
| All 5 Hz and CTRL terms | — | — | unchanged to 3 decimals |

Only one of five conclusions moved, and it moved to null. That is the behaviour you want from a duplication: it inflated precision without inventing the main effects.

> **One housekeeping item.** The copy that was kept is the one in `2hz test`, which has no `gcamp_notes.json`, so the D5 field-of-view shift you confirmed by eye is no longer excluded — that day now enters the 2 Hz chronic maps. It cuts against you rather than for you (a shifted field can only add apparent loss, and 2 Hz still shows no chronic decline), so the conclusion is safe, but copying the notes file across and re-running would make it exact. The condition label also records 10 ms rather than the true 15 ms pulse width for that animal; it affects the families table only, not any model.

## Where this leaves the project

The premise was that optogenetic stimulation of V1 produces a migraine-like state. The data support a narrower version of that than the first draft of this memo claimed: **the 2 Hz protocol drives V1 into a non-habituating regime** — responses and reliability preserved where controls lose both, with an effect size around `d = 2` and near-complete separation of animals. Whether that sits on a background of suppressed preactivation, as the migraine literature would predict, hangs on one test that has not been run yet.

The 5 Hz arm, despite having more animals, currently supports nothing at the animal level. That inverts the dose reading: it is not obviously that more light does more. It may be that the two protocols do different things, that 5 Hz effects live only in spontaneous activity, or that the arm is too heterogeneous to read. Frequency and duty cycle are also fully confounded in the present design (2 Hz × 10 ms = 2%, 5 Hz × 10 ms = 5%), so even a clean difference could not be attributed to frequency — the 15 ms and 50 ms arms were the right idea for breaking that confound but have one animal between them.

Two honest limits on the migraine reading. Published optogenetic CSD induction from V1 uses brief high-power continuous light — around 4 mW for 10 s, with a threshold near 4 s at that power — and produces migraine-like pain and photophobia lasting 24–48 h after a *single* event ([J Headache Pain, 2025](https://link.springer.com/article/10.1186/s10194-025-01983-8)). Your protocol is different in kind: an hour of low-duty-cycle pulses, repeated daily for two weeks. Whether that is repeated sub-threshold activation or repeated CSD is currently unmeasured, and it is the single most valuable thing you could add. Second, imaging is off during the stimulation hour by design, so any event during the hour is invisible; everything here is an aftermath measurement.

### Next, in order

1. **Export per-animal spontaneous activity and permute it.** One column in the pipeline. It decides whether the 5 Hz arm has any finding at all, and whether the preactivation half of the migraine story stands.
2. **Make animal-level permutation the headline statistic** for every claim, with the cell-level model reported beside it as the precision estimate. State the achievable floor (`p = 0.029` at 4v4) so no reader over-interprets a marginal value.
3. **Close the floor objection.** Model POST with PRE as a covariate rather than differencing — the standard correction when baselines differ between arms — and repeat on cells responsive at PRE only.
4. **Establish what the stimulation actually does.** Record during the stimulation hour in one or two animals, or add a wide-field channel, and look for propagating events. It converts "migraine-like" from analogy to measurement, and it is the only item here that adds new physics rather than new statistics.
5. **Add 2 Hz controls, and break the frequency/dose confound.** A 2 Hz × 25 ms arm matches the 5 Hz duty cycle at a different frequency: if it behaves like 5 Hz the variable is total light, if like 2 Hz it is frequency. One condition, decisive.
6. **Power the naive day.** The largest, cleanest effect is on D1 (`+0.26 SD [+0.09, +0.43]`); a cohort designed around first exposure needs far fewer animals than one designed around the two-week mean.

## Sources

1. An increase in spontaneous activity mediates visual habituation. *Cell Reports*, 2022 — [cell.com](https://www.cell.com/cell-reports/fulltext/S2211-1247(22)00515-0) ([preprint full text](https://www.biorxiv.org/content/10.1101/507814v1.full)). Mouse V1 L2/3; habituation over 7 days of daily stimulation; pre-stimulus activity +35–40%; evoked activity stable; raising pre-stimulus activity optogenetically is sufficient to induce habituation.
2. Cortical inhibition and habituation to evoked potentials: relevance for pathophysiology of migraine. *The Journal of Headache and Pain*, 2009 — [springer.com](https://link.springer.com/article/10.1007/s10194-008-0095-x). Interictal habituation deficit as a hallmark across sensory modalities; Schoenen's reduced-preactivation model; normalisation just before and during attacks.
3. Mechanisms of neuronal silencing after cortical spreading depression. *Cerebral Cortex*, 2017 — [academic.oup.com](https://academic.oup.com/cercor/article/27/2/bhv328/3056266). Spontaneous activity abolished then persistently reduced >60 min; presynaptic release probability down, inhibitory drive up.
4. Optogenetic cortical spreading depression originating from the primary visual cortex induces migraine-like pain and anxiety behaviors in freely moving C57BL/6J mice. *The Journal of Headache and Pain*, 2025 — [springer.com](https://link.springer.com/article/10.1186/s10194-025-01983-8). 465 nm, 4 mW, 10 s continuous; CSD threshold 4.13 ± 0.58 s at 4 mW; single events; photophobia and anxiety persisting 24–48 h.
5. Representational drift in the mouse visual cortex. *Current Biology*, 2021 — [sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0960982221010526). Context for the substantial baseline cell turnover seen in all groups, including controls.

---

Prepared from run 260901_2114 (`stats_report.html`, `models.json`, and the per-animal CSV bundle). Effect sizes are state-adjusted mixed-model estimates in per-cell SD units unless noted; brackets are 95% confidence intervals. Animal-level tests treat each animal as one observation.
