# AQuA2 GUI: setup notes for the afferent pre/post analysis

Working notes from getting AQuA2 (MATLAB) running on the `aqua_input/` tiffs that afferent2p exports. Companion to `claude/afferent-prepost-analysis-pipeline.md`.

## Bugs hit, and their fixes

**1. `aqua_gui` won't start.**

```
Property assignment is not allowed when the object is empty.
Error in ui.com.addCon_dat (line 123)  pMov1.Toolbar.Visible = 'off';
```

AQuA2 hides the per-axes hover toolbar with `X.Toolbar.Visible = 'off'`. MATLAB creates that `AxesToolbar` object lazily, and in some releases/display setups it doesn't exist yet at that point, so `.Toolbar` is an empty placeholder. Fix: replace every such line with

```matlab
if ~isempty(X), X.Toolbar = []; end
```

Assigning `[]` is the supported way to suppress the toolbar and works whether or not the object exists. The same pattern recurs across the `+ui` files, so patch all of them at once (`patch_aqua2_toolbar.m`).

**2. `medfilt2` error during the Pre step.**

```
Error using medfilt2
Expected input number 2, [M N], to be integer-valued.
Error in ui.detect.preProcessRun (line 71)
```

The "Median filter radius" field was set to a non-integer (e.g. 0.5). AQuA2 passes it straight through as `medfilt2(img,[medSize medSize])` with only a `> 0` guard and no rounding. Use `0` (off) or a positive integer. Note the label says "radius" but the value is the **window size** in pixels — `3` means a 3×3 filter.

## Project setup values for this rig

| Field | Value | Why |
|---|---|---|
| Temporal resolution (s/frame) | `0.064` | 15.62 Hz acquisition; afferent2p does no temporal binning |
| Spatial resolution (µm/px) | `0.7862` | Verify Scanbox x vs y are actually square |
| Channel 2 | empty | single green channel |
| Data type preset | Default | no built-in preset targets axonal GCaMP |
| Border exclusion | ≥ max registration shift | read max \|shift\| from afferent2p QC; default 5 px ≈ 4 µm is often too small |

Sanity check: 28000 frames ÷ 15.62 Hz = 1792 s ≈ 30 min, matching the baseline recording length.

## Pre tab

- **Registration: OFF.** afferent2p already registered with Suite2p (rigid + nonrigid); the files are `*_reg.tif`. A second, weaker rigid registration adds interpolation blur and edge artifacts.
- **Photobleach: "Remove bleach globally"**, identical setting for pre and post. Corrects within-run decay only — the pre→post F₀ drop from mCitrine bleaching still needs the opsin-negative control group.
- **Median filter: 0.** Axons are only 1–2 px wide at 0.786 µm/px, so a 3×3 median filter erases them. Raise only if hot pixels are confirmed (bright in the *same* location in every frame — check a minimum-intensity projection).
- **Gaussian filter radius: 0.5** (spatial σ in px). Same caution — few pixels per axon to spend.

## Act tab (active voxel detection)

- **Intensity threshold scaling factor:** start at `3`. Units are noise σ, which makes it relatively robust to the mCitrine bleaching that shifts F₀ between sessions — a useful property for the pre/post design. Drop to 2–2.5 if nothing detects, raise to 4 if everything lights up.
- **Minimum duration:** default `5` frames = 0.32 s, likely too permissive. Axonal GCaMP6 transients span ~8–30 frames at this rate; 8–10 is a better floor. Measure real transient widths with **Check ROI curve** first.
- **Minimum size:** default `20` px ≈ 12 µm². This **excludes individual boutons** (~1 µm ≈ 2 px) and biases toward extended axon-segment events. Lower to 8–10 only if the mean projection shows punctate rather than continuous structure.
- Advanced filters: off initially.

## Workflow rules

- **Validate before tuning.** Use **Check ROI curve** on an ROI drawn over a visible axon (picked from a mean/percentile projection, not a single frame — single frames of this data are shot-noise dominated). If transients aren't visible by eye, the problem is upstream (SNR, registration, z-drift), not in the parameters.
- **Tune on a crop.** The full movie is ~28000 frames ≈ 20 GB as uint16; AQuA2 holds several copies. Cut a ~3000-frame × 256×256 subset for parameter iteration, then batch the full movies.
- **SaveOpts / LoadOpts.** Save the parameter file as soon as detection looks right, and load it for the post movie and every other session. Identical parameters across pre and post is what makes the comparison meaningful.

> Note (Nov 2026): these are the Sep 17 notes for the full-rate suite2p-registered route. The pipeline in use since Sep 19 runs AQuA2 on 2 s/frame binned movies (1.984 s/frame, 0.792 µm/px, border 5) with AQuA2's own registration ON; see `afferent-prepost-analysis-pipeline.md` and the pipeline's GUIDE.md.
