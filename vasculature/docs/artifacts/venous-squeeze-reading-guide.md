# Venous Squeeze Reading Guide

Source: https://claude.ai/artifact/Y9VygbGVbm6f1R3s5XQRdn (claude.ai artifact, last updated 2026-09-24)

A figure-by-figure reading guide to the Zhang ... Drew (2026) bioRxiv preprint showing that abdominal pressure squeezes the bridging veins and superior sagittal sinus within about 100 ms, producing a fake GCaMP-like fluorescence transient at movement onset, with notes on what this means for our locomotion/GCaMP pipeline.

---

*Reading guide · bioRxiv preprint, May 2026*

## Ultrafast venous and sagittal sinus constrictions driven by abdominal pressure

Zhang Q, Garborg CS, Frank N, Salehi Shahrbarbaki F, Turner KL, Drew PJ (Penn State / Michigan State). [doi:10.1101/2026.05.02.722426](https://doi.org/10.1101/2026.05.02.722426) (v1). In the original artifact the figures were reproduced from the PDF for private reading; they are not reproduced here (each figure is described in text instead).

> When a mouse tenses its abdominal muscles (to run, to startle at an air puff, or on waking), pressure travels up the spinal canal and raises intracranial pressure. The largest, lowest-pressure veins, the bridging veins and the superior sagittal sinus, get squeezed within about 100 ms. Less blood means less light absorbed, so fluorescence goes up. The result is a fake, GCaMP-like transient that coincides with movement onset.

### Contents

Main figures (figures marked "priority" bear on our locomotion problem):

1. [Figure 1 - Locomotion](#figure-1--veins-constrict-at-run-onset-and-the-biggest-constrict-most) (priority)
2. [Figure 2 - Whisker puff](#figure-2--a-whisker-puff-squeezes-the-sinus-too-without-any-running)
3. [Figure 3 - Fluorescence artifact](#figure-3--the-constriction-makes-a-fake-calcium-like-signal) (priority)
4. [Figure 4 - ICP & breathing](#figure-4--intracranial-pressure-jumps-with-running-and-ripples-with-each-breath) (priority)
5. [Figure 5 - Sleep -> wake](#figure-5--veins-dilate-in-sleep-and-snap-shut-on-waking)
6. [Figure 6 - Abdominal cuff](#figure-6--squeezing-the-abdomen-reproduces-it-in-an-anaesthetised-mouse)
7. [Figure 7 - The model](#figure-7--the-model-pressure-from-the-abdomen-overtakes-the-lowest-pressure-veins) (priority)
- [Critical read](#critical-read--what-to-question-in-journal-club)
- [For your pipeline](#take-home--what-to-take-back-to-the-gcamp-analysis)

Supplements sit under the figure they support.

### If you have 20 minutes

Follow the causal chain, then the part that touches your data:

Fig 1 -> Supp 4 -> Fig 4C -> Fig 3 -> Fig 7

(In the original page, green dots in the index marked the figures that bear on the locomotion problem; these are the "priority" figures above.)

### Colour code used in every figure

| Label | Vessel | Colour in figures | Direction during running |
|---|---|---|---|
| SSS | superior sagittal sinus | pink | down (constricts) |
| BV | bridging vein | blue | down (constricts) |
| FL/HL | pial vein over limb cortex | green | up (dilates) |
| Wh | pial vein over whisker cortex | orange | up (dilates) |

Ordered largest to smallest. Pink and blue go *down*; green and orange go *up*.

### Measures

| Measure | Meaning |
|---|---|
| ΔD/D₀ | vessel diameter change, FWHM across the vessel |
| ΔHbT | total haemoglobin, a proxy for blood volume |
| ΔHbO−HbR | oxygenation index |
| IOS | intrinsic optical signal, 530/470 nm reflectance |
| EMG_Abd | external oblique abdominal muscle |
| ICP | intracranial pressure (mmHg) |
| LDF / CBF | laser Doppler flow, parenchymal blood flow |

### Preparation

Awake head-fixed C57BL/6 on a spherical treadmill, polished thinned-skull window over both hemispheres, widefield camera at 60 Hz alternating green/blue. Plus 2-photon (FITC-dextran), one miniscope mouse, CAG-EGFP mice for the artifact test, and ICP recordings re-analysed from earlier papers.

24 mice total · 4–7 per experiment · mostly male

---

## Figure 1 — Veins constrict at run onset, and the biggest constrict most

**Read first.**

*What do the different veins do when a mouse starts running?*

n = 7 mice · widefield IOS · abdominal + neck EMG

**Figure (not reproduced):** Figure 1 of the paper, venous diameter changes during voluntary locomotion: anatomy schematic, window image with measurement sites, a single run bout, onset/offset-triggered averages, group bar plots, and the diameter-vs-constriction scatter (slope −0.017 %/µm, R² 0.63).

- **A** — Anatomy. Cortical veins (blue) drain into bridging veins, which cross the arachnoid into the SSS, a sinus inside the dura. Almost all dorsal-cortex blood leaves through this route, which makes it a bottleneck.
- **B** — Setup and the measurement sites. The coloured lines on the window image are where each diameter is taken; this map recurs in every later figure.
- **C** — One run bout (grey). The SSS and BV diameters drop right at onset, while the FL/HL vein rises slowly. ΔHbT falls over SSS/BV and rises elsewhere. Both EMG channels burst with running.
- **D** — Onset- and offset-triggered averages. The pink SSS trace is the fastest mover: down within about 100 ms, before any cortical dilation. Everything returns to baseline after the bout ends.
- **E** — At 0–1 s: SSS −1.18 ± 0.42 %, BV −0.76 ± 0.84 %; pial veins not significant. At 2–5 s the SSS stays down (−2.13 %), the BV recovers.
- **F** — The key plot. The larger the vein's baseline diameter, the more it constricts (slope −0.017 %/µm, R² 0.63). Vessel size predicts the direction of the response, which is what the pressure model in Fig 7 explains.

> **Take-away:** A size-graded constriction arrives within about 100 ms of run onset. That is too fast for smooth muscle, so the authors look for a mechanical cause.

### Supplementary Figure 1 — Same result with 2-photon and in a freely moving mouse

**Figure (not reproduced):** Supplementary Figure 1, two-photon linescans and miniscope imaging (arteriole +2.0 % / +3.8 %, bridging vein −0.76 %, SSS −3.0 % sustained; LDF +~19 %).

- **A–C** — 2-photon linescans of FITC-dextran (n = 5). A pial arteriole dilates (+2.0 %, p = 0.09 at onset; +3.8 % sustained) while the bridging vein constricts (−0.76 %). With sustained running the SSS reaches −3.0 %.
- **D** — Laser Doppler flow rises about 19 % for both sustained running and a whisker puff.
- **E** — Miniscope in an open field: the bridging vein narrows when the mouse runs, so head fixation is not the cause.

### Supplementary Figure 2 — Blood-volume maps show the midline drop within 100 ms

**Figure (not reproduced):** Supplementary Figure 2, haemoglobin responses to locomotion: ΔHbT frame sequence after run onset and group ΔHbT/oxygenation plots (SSS −11 µM, BV −15 µM at 0–1 s).

- **B** — ΔHbT frames every 67 ms after run onset (54 events, one mouse). Blue (less blood) appears along the midline sinus first; red (more blood) spreads over cortex a second or more later.
- **C–F** — ΔHbT at 0–1 s: SSS −11 µM, BV −15 µM. Oxygenation rises everywhere during sustained running, even where volume fell.

### Supplementary Figure 4 — Abdominal muscle activity comes first; the wheel comes last

**Figure (not reproduced):** Supplementary Figure 4, abdominal EMG leads venous constriction: EMG power around locomotion onset, puff-evoked EMG, and cross-correlations of vessel signals with EMG and with wheel speed.

- **B** — Abdominal EMG power starts rising before wheel-detected locomotion onset.
- **C** — Contralateral and ipsilateral whisker puffs both evoke a sharp EMG burst. The auditory control does not.
- **D** — Cross-correlation with EMG: the SSS trace is negatively correlated (constricts) with the peak just after zero, so EMG leads. Pial veins correlate positively (dilate) and more slowly.
- **E** — Against wheel speed, the SSS dip is roughly simultaneous with speed. The wheel is a late marker of an event that started in the trunk.

> **For your pipeline:** Your locomotion state comes from the wheel, which on this evidence is the last signal to move. The body event (EMG, then pressure) starts earlier. A trial whose wheel onset falls just after your 4 s response window can still carry the onset transient inside it.

---

## Figure 2 — A whisker puff squeezes the sinus too, without any running

*Does a brief sensory stimulus trigger the same constriction?*

n = 4 mice · 0.1 s air puff at 10 psi · left, right or auditory control, 1:1:1, every 30 s

**Figure (not reproduced):** Figure 2 of the paper, venous constriction after whisker stimulation: single trials, stimulus-triggered averages with a 100 ms zoom, and group bar plots (SSS −1.24 % at 0–0.5 s; FL/HL +1.27 % at 0.5–3 s).

- **B** — Single trials. Contralateral and ipsilateral puffs both dip the SSS and both evoke an abdominal EMG burst; the auditory control does neither. The ipsilateral result matters: the constriction is a body reflex, not a barrel-cortex response.
- **C** — Averages, with a 100 ms zoom on the right. The SSS drops within 0.2 s. Pial veins dilate about a second later (ordinary neurovascular coupling).
- **D** — At 0–0.5 s: SSS −1.24 %, BV −0.25 %, FL/HL −0.56 % (all p = 0.014). At 0.5–3 s the FL/HL vein dilates +1.27 %.

> **Take-away:** Any stimulus that startles the animal can recruit the abdominal reflex and produce a fast artifact locked to stimulus onset.

### Supplementary Figure 3 — Blood volume and oxygenation after the puff

**Figure (not reproduced):** Supplementary Figure 3, haemoglobin responses to whisker stimulation: ΔHbT frames every 33 ms and group ΔHbT/oxygenation plots (SSS ΔHbT −12.8 µM at 0–0.5 s).

- **B** — ΔHbT frames every 33 ms: a midline decrease within the first 100–300 ms, then a broad increase.
- **C–F** — SSS ΔHbT −12.8 µM at 0–0.5 s. A brief drop in oxygenation follows. The authors link it to the "initial dip" seen in fMRI and suggest it happens because an unexpected puff is not accompanied by the faster breathing that comes with running.

> **For your pipeline:** A visual stimulus resembles their auditory control more than their puff: it produces cortical activity without a trunk reflex. Your real-vs-sham check (the stimulus does not make a still animal move) points the same way. Any tactile or startling stimulus would deserve its own check.

---

## Figure 3 — The constriction makes a fake calcium-like signal

**Read first.**

*How big is the fluorescence artifact, and can standard corrections remove it?*

n = 6 CAG-EGFP mice · head-fixed in a tube · 470 nm fluorescence + 530/630 nm reflectance, 10 Hz per channel

**Figure (not reproduced):** Figure 3 of the paper, fluorescence artifact from venous constriction in GFP mice: ΔHbT and raw ΔF/F frame sequences, three correction methods, example and group traces, and bar plots (uncorrected BV +1.86 % at 0–0.8 s; −1.3 to −2.3 % at 1–3 s).

- **A** — The trick: GFP does not report activity, so every change in its fluorescence here comes from haemoglobin absorbing light.
- **B** — Frames 0–1 s after a puff (47 events, one mouse). Row 1 is ΔHbT, row 2 the raw ΔF/F. Where blood volume falls along the midline (blue in row 1), the raw fluorescence rises (yellow in row 2). Rows 3–5 are the three corrections.
- **C–D** — Traces for one mouse and the group. The uncorrected ΔF/F mirrors ΔHbT: a fast positive transient, then a dip at 1–3 s as arterial dilation brings more blood in.
- **E** — At 0–0.8 s, uncorrected: SSS +0.85 %, BV +1.86 %, FL/HL +0.75 %. The spatial model and pixel-wise regression reduce these to about 0.05–0.15 %. The Ex-Em (Beer-Lambert) method overcorrects to −0.84 and −1.35 %.
- **F** — At 1–3 s the uncorrected signal is −1.3 to −2.3 % at every site. All three corrections remove most of it.

> **Take-away:** In widefield imaging the artifact is 1–2 % ΔF/F with a sub-second onset, largest over big veins. It looks like an indicator response and needs a reflectance channel to correct.

### Supplementary Figure 5 — The same animals, in the tube

**Figure (not reproduced):** Supplementary Figure 5, haemodynamic responses in the tube: puff-evoked abdominal EMG and 530 nm reflectance, oxygenation and diameter traces for the Fig 3 events.

- **B** — Puff-evoked abdominal EMG in the tube is smaller than on the treadmill, but still present for both puffs and absent for the tone.
- **C–D** — 530 nm reflectance (ΔR/R₀), oxygenation and diameter for the same events: the inputs the corrections in Fig 3 work from.

> **For your pipeline:** These are one-photon widefield numbers. In 2-photon the excitation path and depth differ, so the size will not carry over directly; the paper's reference 37 (Yogesh et al., *eLife* 2025) measures haemodynamic occlusion in 2-photon imaging and is the one to read next.
>
> You cannot apply their corrections: you have a single green channel and no reflectance. What you can do is look for the artifact's signature. It is spatially broad and strongest near large vessels, so it should appear in neuropil rings and non-cell background as much as in cells.

---

## Figure 4 — Intracranial pressure jumps with running and ripples with each breath

**Read panel C.**

*Is ICP a plausible driver, and does it have the right timing?*

ICP: n = 10 mice, catheter in cortex, 1 kHz, re-analysed from earlier papers · respiration: n = 4, thermocouple at the nose

**Figure (not reproduced):** Figure 4 of the paper, intracranial pressure during locomotion and respiration: ICP trace with run bouts, resting power spectrum with 1/f fit, locomotion-triggered ICP (+~8 mmHg, rising 0.56 ± 0.29 s before wheel onset), and breathing-locked diameter oscillations and cross-correlations.

- **A** — ICP rises by more than 10 mmHg during each run bout (black ticks). The inset shows breathing-driven oscillations of well under 1 mmHg.
- **B** — Resting power spectrum with a 1/f fit. The residual shows peaks at breathing (~3 Hz, 0.41 ± 0.11 mmHg) and heartbeat (~10 Hz, an order of magnitude smaller).
- **C** — Locomotion-triggered ICP: up about 8 mmHg, and it **starts rising 0.56 ± 0.29 s before wheel onset**. Falls back after the bout.
- **D** — At rest, vein diameters oscillate with breathing: narrower on inhalation, wider on exhalation.
- **E** — Cross-correlation of breathing with diameter and with ΔHbT, strongest for the SSS, within ±0.3 s.

> **Take-away:** ICP has the right size and the right timing. Locomotion produces large steps, breathing produces small ripples, and the veins follow both.

> **For your pipeline:** This half-second lead is the concrete number to use. Pad the pre-onset side of your state-change exclusion by about 0.5 s, or test whether the "starts moving later in the window" trials still differ from still trials when the window ends 0.5 s before wheel onset.
>
> Breathing at ~3 Hz is below your 7.8 Hz Nyquist (15.62 Hz frames). Because it is not locked to your stimulus, it averages out as noise rather than bias.

---

## Figure 5 — Veins dilate in sleep and snap shut on waking

*Does the effect generalise to arousal transitions without running?*

n = 4 (NREM) and 2 (REM) mice · head-fixed in a tube · states scored from whisking, body pressure sensor and EMG

**Figure (not reproduced):** Figure 5 of the paper, venous changes across sleep and waking: a long awake/NREM/REM recording, transition-triggered averages of EMG, ΔHbT, oxygenation and diameter, and group comparisons (waking constricts the SSS by several percent).

- **B** — A long recording across awake, NREM and REM (top bar). Veins and blood volume rise during sleep, then drop abruptly at waking, when the abdominal EMG surges.
- **C** — Transition-triggered averages: EMG power, ΔHbT, oxygenation and diameter around sleep onset and waking. Note that the y-scales differ between panels.
- **D** — Waking constricts the SSS by several percent, more than locomotion does, and unlike locomotion lowers oxygenation around these veins.

> **Take-away:** Any shift in trunk muscle tone moves the large veins, so vascular state follows arousal as well as locomotion.

> **For your pipeline:** Sessions where the animal is drowsy (small pupil, no running) differ in vascular state as well as neural state. This supports treating your "field does not respond" sessions as a state difference rather than a loss of responsive cells.

---

## Figure 6 — Squeezing the abdomen reproduces it in an anaesthetised mouse

*Is abdominal pressure enough on its own, without muscle activity or behaviour?*

n = 4 mice · 1 % isoflurane · pneumatic belt, 10 psi for 1 s, every 45 s · laser Doppler flow

**Figure (not reproduced):** Figure 6 of the paper, abdominal compression causes SSS constriction: example inflation trial with flat EMG, group diameter (SSS −3.18 %), ΔHbT (−19 / −12 µM) and oxygenation plots, and LDF (+~3 % at 0.2–1.2 s).

- **B** — Example: the SSS dips during inflation (grey) while the abdominal EMG stays flat. The muscles are not involved; pressure alone does it.
- **C** — SSS −3.18 %, BV −0.53 %; pial veins unchanged. The same size ordering as during running.
- **D–E** — Blood volume falls over SSS/BV (−19 and −12 µM); oxygenation rises over BV and pial veins.
- **F** — Parenchymal blood flow rises about 3 % at 0.2–1.2 s: squeezing the outflow briefly pushes more blood through cortex.

> **Take-away:** This is the causal test. External pressure produces the same size-graded constriction with no behaviour, so the mechanism is mechanical.

---

## Figure 7 — The model: pressure from the abdomen overtakes the lowest-pressure veins

**Read first.**

*Why do the largest veins constrict while small ones dilate?*

**Figure (not reproduced):** Figure 7 of the paper, schematic model of abdominal pressure transmission: abdomen -> vertebral venous plexus -> dural sac/CSF -> ICP, and pressure-vs-vessel diagrams at rest, on inhalation (+~1 mmHg) and on abdominal contraction (+~10 mmHg).

- **A** — Abdominal muscles contract, intra-abdominal pressure pushes blood into the valveless vertebral venous plexus inside the spinal canal, which presses on the dural sac and drives CSF toward the skull. ICP rises.
- **B** — At rest, pressure inside veins falls along the drainage tree (venules > BV > SSS) while ICP is the same everywhere and stays below all of them.
- **C** — Inhalation adds about 1 mmHg of ICP, enough to exceed the SSS's pressure only, so the SSS narrows a little with each breath.
- **D** — Abdominal contraction adds about 10 mmHg, which exceeds both the BV and SSS pressures, so both compress. Small veins stay above ICP and remain free to dilate with local activity.

> **Take-away:** A single pressure threshold explains the size ordering in Fig 1F, the breathing ripple in Fig 4, and the cuff result in Fig 6.

> **Why your PI might care beyond artifacts:** The discussion notes that the SSS is innervated by trigeminal and cervical nerves and is pain-sensitive, and cites work showing trigeminal afferents respond to locomotion-related meningeal deformation. If your project has a migraine or meningeal angle, running and arousal changes are a mechanical load on those structures, which makes the PRE/POST shifts in running and pupil worth reporting as findings.

---

## Critical read — What to question in journal club

- **Small samples.** Four to seven mice per experiment, mostly male, and two mice for REM. p = 0.0143 recurs across unrelated measures, a sign that it is the floor of the test at n = 4 rather than a measure of effect strength.
- **The chain is assembled across cohorts.** ICP comes from different animals in earlier papers. The 0.56 s ICP lead and the imaging latencies were never measured together.
- **The test choice.** The methods use a Wilcoxon *rank-sum* "against zero". Testing one sample against zero normally calls for the signed-rank test. It probably does not change the conclusions, but it is worth asking about.
- **Diameter from FWHM through thinned skull.** The authors note the SSS is triangular in cross-section, so percent-diameter is approximate. The same abdominal pressure also moves the brain (their ref 23). The widefield data come without a motion control, while the 2-photon data excluded z-motion.
- **The artifact test is one-photon.** How much of this reaches a 2-photon GCaMP measurement at cellular resolution is not addressed here.
- **The cuff runs under isoflurane,** which dilates vessels and changes baseline pressures. The match with the awake data is qualitative.

---

## Take-home — What to take back to the GCaMP analysis

**Your pipeline.**

The paper does not offer a new way to analyse locomotion; its event definitions (≥ 3 s rest before, ≥ 5 s running after, rest starting 4 s after a bout) are simpler than your state-change flag. It is useful as a mechanism for the movement-alone signal you measured, and as three checks:

1. **Pad the wheel onset.** Physiology leads the wheel by about 0.5 s (Fig 4C, Supp 4B). Extend the state-change exclusion to wheel onsets up to ~0.5 s after the response window ends, and require about 4 s of stillness after a bout before a baseline counts as still.
2. **Test whether the movement-alone signal is vascular.** Re-run the movement-onset-triggered average (panel D of your state check) on neuropil rings and non-cell background pixels, and split ROIs by distance to large vessels. An onset-locked bump in the background, largest near vessels, would point to haemodynamics rather than afferent activity.
3. **Check registration around onsets.** The same pressure moves the brain. Look at the frame-to-reference correlation or z-motion proxy in the second after each movement onset.
4. **Read Yogesh et al. 2025 (*eLife*)** if you want a correction rather than an exclusion for single-channel 2-photon data.

Your headline comparison should be largely protected: the stable view already drops trials with a state change, the visual stimulus does not trigger movement, and any residual artifact adds to CTRL and EXP alike. Absolute responder counts and time-to-peak are more exposed than the family deltas.

---

*Figures © the authors, taken from the bioRxiv PDF your PI shared, for private reading (not reproduced in this Markdown copy). Numbers are mean ± SD as reported in the figure legends.*
