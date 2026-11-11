# Oct 12 task list — DONE in v32 (Oct 13 2026)

**Status: implemented and delivered as afferent2p v32** — see `claude/pipeline-v32-report-redesign.md` for what was built, the decisions still open (which responsive-trial rule the lab uses; the optional cohort summary / non-fibre event source / fibre-area normalisation were NOT done, ask first) and what to run next on samson. The original spec is kept below for reference.

---

# DEFERRED TASK LIST (written Oct 12–13 2026, run Oct 13 after the weekly usage reset)

User's words:
1. (Oct 12) "the pipeline in general looks ok, but might need some changes in the report: lets first give a general representation of the data, can you look at the attached presentation and use the same or similar measurements? add more if you think it might help, the current time course and loco-matching can be in the later panel." → v32 main report sections 1 (general representation), 7 (time course / matched speed later).
2. (Oct 12) "can you give a heatmap activity comparison?" → v32 section 2 (`afferent2p/heatmaps.py`).
3. (Oct 12) "the pipeline 2b auto tune always give the aqua folder cannot be located error. let make it with manual copying commands to matlab" → v32 "Run all (guided)" with manual MATLAB pastes; scripts self-locate AQuA2.
4. (Oct 12) "the tuning step looks a bit not intuitive, so it is step 3 then 2b then back to step 3?" → v32 tabs 1 Unpack · 2 Clean & check · 3 Detection settings · 4 Batch · 5 MATLAB · 6 Curate · 7 Analyze · 8 Report.
5. (Oct 12) "the monitor session in the report is good but maybe make it into a separate and better structured report? currently it is too much to read as a whole" → v32 QC report (`aqua_report.py --qc`, `afferent2p/qc_report.py`).
6. (Oct 13) "YU's report shows analyze on the locomotion-induced activity, i think we could use the similar pipeline as we are testing spontaneous activiy and the locomotion is the major factor of triggering it" → v32 section 3 (`afferent2p/locotrials.py`).

## The presentation (source for items 1, 2 and 6)
`LabMeeting_261012_Yu.pdf` (29 slides, "Interactions Between Cortical Astrocytes and Meningeal Afferents in Migraine", Yu Yong, 10/12/2026). Mostly astrocyte GCaMP6f data, same rig/paradigm (2 s frames, pre 30 min → 1 h optostim → post 30 min; Nav1.8Cre control vs Nav1.8:ReaChR; weekly days), plus meningeal-afferent heatmaps. Measurements, slide by slide:
- s4: AQuA event overlays running vs resting; event rate/min running vs resting (bar mean ± sem + dot per recording/animal); ECDF of normalized dF/F, running vs resting.
- s5–s7 "locomotion trials": trial = 30 frames (1 min) = pre-bout window 5 frames (10 s) no running + running window 15 frames (30 s) running + post-bout window 10 frames (20 s) no running. Responsive trial: peak in the running window ≥ 1.5 × baseline (pre-bout window). Per day × pre/post: onset-aligned mean dF/F trace ± sem, pie of responsive-trial %, onset-aligned "#Events trace"; number of responsive loco trials per animal. Acute effect = Peak_Post / Peak_Pre (per day); chronic effect = Peak_Pre(Day#) / Peak_Pre(Day1); same for responsive-trial calcium event number. Stats: mean ± sem, dot = responsive trial, two-way ANOVA + Sidak, Mann–Whitney.
- s8 locomotion itself: % of running and number of bouts, pre vs post across days (the behaviour check).
- s9–s10 "Locotrial_Loco2Fluor GLM": predictor speed, response CFU × trial; trainFrac 0.75, Ncycle 50, CV fold 5, window [−5, 10] s, minDev 0.2; deviance explained per well-fitted CFU pre vs post, peak coefficient.
- s11–s12 "LocoPupil2Fluor GLM": predictors locomotion + pupil, response CFU dF/F trace; trainFrac 0.75, Ncycle 20, CV fold 10, window [−20, 20] s, minDev 0.2; deviance explained; LOPO relative explanatory value of speed and pupil; kernels.
- s13–s15 arousal: pupil variation (std of the 0–1 normalized diameter); resting pupil diameter (valid recording ≥ 10 % resting frames); locomotion-related pupil dilation; spontaneous pupil dilation; each acute (post/pre) and chronic (Day# pre / Day1 pre).
- s19, s20, s28: event-count histogram over the whole session; CFU × frame heatmap (0–1) with the event raster under it; running % pre vs post; event rate/min and event max dF/F in running and in resting, pre vs post; speed vs time with locomotion-state shading and bouts.
- s21, s27 (meningeal afferents): CFU × frame heatmaps of Pre | Post 0–15 min | Post 30–45 min, "three groups via k-means on mean levels", rows clustered & sorted, 0–1; pie of Up / Down / unchanged % of CFUs with "Total = N CFUs".

## Spec summary (as implemented; details in the v32 doc)
- Item 1: main report order 0 QC paragraph · 1 general representation (G1 sessions on the clock, G2 rates by state + overlays, G3 ECDFs, G6 pupil, G8 acute/chronic) · 2 heatmaps · 3 locomotion-evoked · 4 pupil · 5 GLMs · 6 behaviour-controlled change · 7 time course + matched speed · 8 per experiment · 9 settings · 10 methods. Module `afferent2p/general.py`, csvs read by the report.
- Item 2: H1 CFU × time from the movies under matched footprints (fallback AQuA2 CFU dF/F), k-means k=3 ordering, raster + speed; H2 up/down/unchanged (20-s bins, Mann–Whitney) pie + csv + rate scatter; H3 tile maps shown there; H4 peri-onset heatmap (`<exp>_loco_cfu.png`).
- Item 3: (a) no AQuA2-folder block, (b) `which('aqua_cmd_batch')` in every generated script, (c) guided chain with manual pastes, 30-s polling (`aqua_tune.sweep_progress`, `aqua_batch.batch_done`), Continue now / Stop, checkbox for `matlab -batch`, (d) `run_aqua_and_cfu.m`.
- Item 4: tabs renumbered; "Run all (guided)" above the tabs; references in code/GUIDE/README updated; CLI unchanged.
- Item 5: `aqua_report.py --qc` → `<mouse>_qc.html` (+ per-experiment pages when > 80 MB); summary table with warning chips, sticky nav, collapsed recordings with templates, gallery filters/sort/paging, CFU looks table.
- Item 6: `afferent2p/locotrials.py` L1–L8 (strict + relaxed trials, three responsive rules, per-CFU responsiveness vs shifted onsets, rest-only, GLMs, behaviour check, Treg handled by the same code); two-way ANOVA not implemented (cohort-level).

## Optional additions (not done — ask the user before doing)
- Cohort summary across animals (`aqua_cohort.py`; animals on samson: NavRearChrTIGRE-2608-1, -2608-2, -2609-1; NavTIGRE-2608-1, -2609-1; NavTIGRETreg-2608-1, -2609-1, -2609-2 — groups needed).
- Analysis event source "drop suggested non-fibre CFUs" (`cfu_shape`); main-report CFU outlines by looks.
- Fibre-area normalisation (events/min per mm² of fibre mask).
