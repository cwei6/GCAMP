> Copied on 11 Nov 2026 from the Claude project doc `claude/facial-von-frey-0.1NTG.md` (written 11 Nov 2026 in the chat "Von Frey hindpaw plot"; see `SESSIONS.md`).

# Facial von Frey, 0.1 mg/kg NTG: scoring rules and setup comparison

Last updated: 2026-11-11. Source workbook: Facial_von_Frey.xlsx (sheets 260729-0.1NTG-1stbatch, 260805-0.1NTG-2ndbatch, 261110-0.1NTG).

## Scoring rules in use
- Filaments: 0.008, 0.02, 0.04, 0.07, 0.16, 0.4, 0.6, 1.0, 1.4, 2.0 g. Every series starts at 0.07 g. Pokes at 1.4 and 2.0 g are written in the 1.0 g row.
- 261110 sheet: measurement 1 = columns C-H, measurement 2 = J onward. A multi-letter cell is consecutive pokes; each poke's filament follows the up-down rule from the previous poke.
- Primary method = SUDO: 5th poke moved half a filament step (geometric midpoint with the neighbouring filament). Chaplan series are SUDO-scored from their first 5 pokes. With a 0.07 g start, SUDO can't score above 1.18 g.
- Each mouse per timepoint = mean of its SUDO series. % change = (post - pre) / pre.
- Mice with baseline < 0.07 g are excluded from % change and from the baseline comparison.
- Earlier batches: series with pokes recorded but no score typed in the sheet are not counted (5 series: 3-L pre, 3-R pre, 4-L pre and post, 4-N pre). 260805 post blocks labelled cage 1/2 = baseline cages 3/4.
- 261110 groups: Vehicle = 833070 L/R/B + 837175 N/2L; DT = 833070 N/2L + 837175 L/R/B. All got NTG before the post test. The pre block labelled cage "5" is 837175-N.
- Chaplan alternative (used for the 261110 Chaplan-only plots): final filament x 10^(k*delta), k from Dixon's table, delta = 0.266 (mean log10 step of the filament set).

## Findings (SUDO, baseline >= 0.07 g)
- New setup (261110), n = 4 Vehicle + 4 DT: no 0.1 mg/kg NTG effect. Mean change Vehicle -7%, DT -4% (medians 0% and +8%). Vehicle vs DT: Mann-Whitney p = 0.89.
- Previous setup (260729 + 260805), n = 8: all 8 mice dropped, mean -45% (Wilcoxon vs 0, p = 0.008).
- All 16 mice combined: mean -26% (Wilcoxon vs 0, p = 0.017), driven by the previous setup.
- Previous vs new, % change: Mann-Whitney p = 0.037. Baseline: median 0.14 vs 0.62 g, Mann-Whitney p = 0.014.
- All mice tested at baseline: previous median 0.08 g (5 of 13 under 0.07 g); new median 0.44 g (2 of 10 under 0.07 g).
- Counting the 5 unscored earlier-batch series gives previous mean -52%, p = 0.028 (% change) and p = 0.026 (baseline). Same conclusions.
- Caveat: setup is confounded with batch/date and with the pretreatment injection (every new-setup mouse got vehicle or DT i.p. beforehand).

## Output files
- facial_von_frey_0.1NTG_pct_change (.png/.pdf): new setup Vehicle vs DT, all mice combined, previous vs new % change, previous vs new baseline.
- facial_von_frey_0.1NTG_baseline_distribution (.png/.pdf): baseline histogram per setup, including mice under 0.07 g.
- facial_von_frey_0.1NTG_per_mouse.csv: every mouse, inclusion status and notes.
