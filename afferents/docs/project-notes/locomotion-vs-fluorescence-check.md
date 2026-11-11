# Locomotion vs whole-FOV fluorescence — NavTIGRETreg-2608-1 (Sep 24 2026)

Prompted by Zhang et al., bioRxiv 2026.05.02.722426. The paper reports ultrafast (~100 ms) constriction of the SSS and bridging veins at locomotion onset. The resulting drop in blood volume gives a spurious fluorescence increase.

Data: 14 spontaneous runs (Baseline 260911, D5 260917, D7 260919, D9 260921). Inputs were the *_bin31f_1.98s.tif stacks (ch1) and *_behavior.csv files. Report PDF is in the animal folder: NavTIGRETreg-2608-1_locomotion_vs_fluorescence.pdf

Findings at 2-s resolution:
- No consistent locomotion-locked whole-FOV ΔF/F. After drift removal, r ranged from -0.31 to +0.45 (4 runs significantly positive, 4 significantly negative). Pooled over 159 onsets, the change was about 0% (p≈0.96). Running minus still was under 1% of F in 12/14 runs.
- Single-bin flashes (>5 MAD): 47 in total, 32 of them during running. Running covered 44% of bins, so the flash rate is 2.7× higher when running (binomial p=0.001). 31 of the 47 are from the baseline day.
- The flashes are additive: dark and bright pixels gain similar counts even though bright pixels are ≥13× brighter at baseline. 33 of 47 fit a uniform offset better than a gain. This points to stray light or a PMT baseline shift, not hemoglobin absorption.

Open checks:
- 15.6 Hz frame means from the .sbx to look for the fast artifact
- ch2 vascular channel at running onsets
- a laser-shuttered recording while the mouse runs
- exclude or regress flash bins before AQuA2
- per-CFU locomotion modulation, with running as a covariate for day comparisons (running % differs a lot between days)
