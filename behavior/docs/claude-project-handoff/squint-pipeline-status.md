> Copied on 10 Nov 2026 from the Claude project doc `claude/squint-pipeline-status.md` (last written 1 Oct 2026). Newer than this note: `evaluate_network` ran on 1 Oct (best-120: test RMSE 7.73 px, mAP 97.3), Treg was tracked with best-120, and the review tool, clips and atlas now exist (see `README.md` and `SESSIONS.md`).

# Squint pipeline (Python port) — status and decisions

Code: `D:\Behavior\Chao\squint_pipeline` (run with `python -m squint ...`, see README.md there).
Data: `D:\Behavior\Squinting\compiled\<condition>\<animal>\`.

## Decisions
- Primary openness metric = 8-point polygon area / (canthal width, 5 s rolling median)². Area is the most mistracking-tolerant metric (drop-2-worst-points corr 0.97 vs 0.80 for EAR); width from eyelid_1/eyelid_5 (best-tracked points) removes most of the head-distance confound (raw area r = 0.88 with apparent eye size). EAR exported as secondary channel only.
- Per-point cleaning before any metric: likelihood < 0.25 or jump > 6 MAD from 1 s rolling median → reject; gaps ≤ 0.5 s interpolated; up to 2 missing points imputed from block median shape; frame valid if ≥ 6 tracked points.
- Reference = p95 of scaled area in habituation+baseline of the same day. Per-day references vary 0.47–0.56 across days (camera angle), so cross-day pooling is off by default.
- Squint candidates kept from 1 s; duration stored, 30 s cut only in `*_ge30s` columns (30 s was an arbitrary observation-based choice). Threshold-free metrics (`frac_below_*`) added.
- Wipes: legacy hysteresis detector kept + occlusion-based candidate detector (≥ 4 points lost and trace dips/blanks). Both go to the review tool.
- File discovery handles renamed `D#_phase_k` and raw `animal_date_cam_block` layouts; `paradigm.yaml` per condition maps block order → phase (Treg: hab, ambient rec, dark hab, dark rec; 260915 lacks the dark recording).

## Validation
- `--mode legacy` reproduces MATLAB on NavRearChr-2511-2-15ms D2–D5 (areas 1e-11 px, means 1e-6, identical event counts).
- Bug found in MATLAB workflow: `D*_filtered.csv` glob concatenated every snapshot; 28 recordings in Bone-400-2hz had two snapshots → phases double-counted (D1 of the validation animal: 167,724 vs 139,770 frames).
- Legacy likelihood gates: normalization gate passed 38 % of frames (biased to eyelid_2/8 confidence), validity gate 99.8 %.
- On the validation animal, 11 strong squint candidates (1–20 s, 67 s total) existed but none passed the 30 s cut → reported as zero severe squinting.

## DLC model (D:\Behavior\Squint_8P_ispi)
- Project `Squint_8P_LargeSep11`, iteration-1, trainset95 shuffle2, resnet_50 (resnet50_gn, GN), 448×448 hybrid crops, AdamW 5e-4 → 1e-4 @90 → 1e-5 @120, 200 epochs, batch 8.
- Training completed 2026-09-30 02:31 (run started 09-28 17:02; two earlier attempts on 09-25/09-26 were interrupted at ~epoch 50 and ~150). Snapshots 050/100/150/200 + `snapshot-best-120.pt` present; `learning_stats.csv` has all 200 epochs.
- Final test metrics (epoch 200): RMSE 7.74 px, RMSE_pcutoff 6.67 px, mAP 97.1, mAR 97.7. Plateau from ~epoch 100 (mAP 97.0–97.3); best-by-mAP snapshot = epoch 120 (97.30).
- `evaluation-results-pytorch/iteration-1/CombinedEvaluation-results.csv` is still from the March model — `evaluate_network` not yet run on the new snapshots.

## Next
- Desktop review tool (trace + synced video, add/remove events → `keep` column), event atlas PNGs.
- Hand-scored clip for precision/recall.
- Run `evaluate_network` on the new iteration-1 shuffle2 model, then analyze Treg videos with it; NTG has phases in filenames already (old snapshot best-60).
