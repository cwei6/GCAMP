# First cross-animal protocol comparison (Oct 20, 2026) — 2 Hz vs 5 Hz vs 5 Hz + 40 ms vs sham

Made from the tab-7 tables of NavRearChrTIGRE-2608-1 and 2608-2 (the user ran both to step 7 on Oct 20; no curation, no per-day re-tuning of the new days, no pupil step) with `aqua_conditions.py` (v35). Outputs on samson: `D:\2photon\Chao\GCAMP_affrerents\cross_animal\conditions.html` (+ png, csv, txt). 2609-1 not analysed yet.

## Numbers (events pre → post per day × FOV, 2 pre runs vs 2 post runs of ~15.5 min each)

| protocol | day / FOV | events pre → post | rate ratio (95 % CI) | per-run events pre | per-run events post |
|---|---|---|---|---|---|
| sham | 2608-1 260917 FOV1 | 3 → 0 | 0 | 2, 1 | 0, 0 |
| sham | 2608-1 260917 FOV2 | 0 → 8 | – | 0, 0 | 6, 2 |
| sham | 2608-2 260922 FOV1 | 73 → 35 | 0.49 (0.31–0.74) | 51, 22 | 33, 2 |
| sham | 2608-2 260922 FOV2 | 0 → 1 | – | 0, 0 | 0, 1 |
| 2 Hz 3.3 mW | 2608-1 260911 FOV1 | 2 → 1 | – | 2, 0 | 0, 1 |
| 2 Hz 3.3 mW | 2608-1 260911 FOV2 | 2 → 20 | 10 (2.5–90) few | 1, 1 | 17, 3 |
| 2 Hz max | 2608-1 260921 FOV1 | 16 → 19 | 1.19 (0.58–2.48) | 5, 11 | 5, 14 |
| 2 Hz max | 2608-1 260921 FOV2 | 51 → 144 | 2.84 (2.05–3.99) | 11, 40 | 34, 110 |
| 5 Hz 3.3 mW | 2608-2 260915 FOV1 | 28 → 0 | 0 (0–0.14) | 14, 14 | 0, 0 |
| 5 Hz 3.3 mW | 2608-2 260915 FOV2 | 0 → 0 | – | 0, 0 | 0, 0 |
| 5 Hz max | 2608-1 261002 FOV1 | 30 → 13 | 0.43 (0.21–0.85) | 13, 17 | 12, 1 |
| 5 Hz max | 2608-1 261002 FOV2 | 24 → 89 | 3.74 (2.36–6.15) | 1, 23 | 40, 49 |
| 5 Hz max | 2608-2 261001 FOV1 | 23 → 182 | 7.94 (5.13–12.8) | 13, 10 | 89, 93 |
| 5 Hz max | 2608-2 261001 FOV2 | 33 → 9 | 0.27 (0.12–0.59) | 22, 11 | 0, 9 |
| 5 Hz max + 40 ms (new FOV) | 2608-1 261016 FOV1 | 2 → 5 | 2.5 (0.4–26) few | 2, 0 | 1, 4 |
| 5 Hz max + 40 ms (new FOV) | 2608-1 261016 FOV2 | 25 → 13 | 0.52 (0.24–1.05) | 8, 17 | 9, 4 |
| 5 Hz max + 40 ms (new FOV) | 2608-2 261018 FOV1 | 51 → 32 | 0.63 (0.39–1.00) | 19, 32 | 15, 17 |
| 5 Hz max + 40 ms (new FOV) | 2608-2 261018 FOV2 | 3 → 0 | 0 | 1, 2 | 0, 0 |

Pooled rate post/pre: sham ×0.58 (0.39–0.86); 2 Hz ×2.61 (1.97–3.48); 5 Hz ×2.14 (1.74–2.64); 5 Hz + 40 ms ×0.62 (0.42–0.89). Still-frames-only tells the same story (2608-2 is still ~99 % of the time; 2608-1 runs a lot — its 2 Hz still-only ratios are even larger). Amplitude (median peak dF/F): no systematic change; 40 ms days slightly down (log2 −0.44 ± 0.30, n = 3).

## Reading (as told to the user)

* The 40 ms pulse-width days show NO increase in either animal or FOV: both FOVs flat or down, like the sham days (whose own drift is downward — 2608-2 sham FOV1 73 → 35). So, on these data, 5 Hz + 40 ms does not work better; it does not visibly work at all.
* The increases are confined to the max-power days of the standard protocols, and always to ONE of the two FOVs of the day, with the other FOV flat or down: 2 Hz max 260921 FOV2 ×2.8 (the second post run carries it: 34, 110); 5 Hz max 261002 FOV2 ×3.7 (both post runs up); 5 Hz max 261001 FOV1 ×7.9 (both post runs ~90 events). 3.3 mW days: nothing (2608-2 260915 FOV1 even 28 → 0).
* Caveats: the 40 ms days are new fields (different afferents) — protocol and field changed together, one day per animal; those days ran with the old FOV settings, not re-tuned (2608-1 261016 FOV1 barely detects anything: 2, 0, 1, 4 events); no curation (the big post counts, e.g. 110 and 93 events in one run, include whatever vessel-band / blob CFUs there are — the 260921_004/006 vessel-band CFUs were never curated); n = 2 days per protocol per power.
* Next: 2609-1 (sham 261006 + 5 Hz + 40 ms 261014) through step 7 adds a third animal to the 40 ms group; re-tune the new-FOV days (tab 3, per day); curate the max-power post runs; then tab 8 → Compare the protocols again.
