# Mouse Trigeminal Fiber Planner

Source: https://claude.ai/artifact/AjE8wvfrV7xyqx3o3FFvVY (claude.ai artifact, last updated 2026-10-13)

An interactive three.js stereotaxic planner that shows where the mouse trigeminal ganglion (TG) sits relative to the head, cheek and neck, what an optic fiber has to pass through to reach it, and how far the light reaches at a given wavelength, fiber and power.

*Stereotaxic · bregma origin · skull-flat.* Where the trigeminal ganglion sits relative to the head, cheek and neck, and what an optic fiber has to pass through to reach it.

The interactive 3D version (drag to rotate, scroll to zoom, right-drag to pan) is at the source link above. This document records its data, presets, light model, notes and references. The 3D view itself cannot be reproduced here.

## Coordinate convention

- Origin at bregma, skull-flat.
- **+AP rostral · +ML right · −DV ventral** (millimetres).
- "Tip DV" is measured below the bregma plane (DV 0).

## Interactive controls (summary)

**3D viewport.** A schematic mouse head (transparent envelope), brain, both trigeminal ganglia, the V1, V2 and V3 branches and the trigeminal root (mirrored on both sides), whisker-pad follicles, eyes, and reference marks (midline on the skull plane, bregma and lambda dots, a dashed interaural ear-bar line). Floating labels: trigeminal ganglion, V1, "V2 → infraorbital n.", V3, "trigeminal root → pons", whisker pad, bregma, lambda, interaural, neck, cheek, and a "tip" label at the fiber tip. The fiber is drawn as a shaft from 3 mm outside the skin to the tip, with the light volume as a translucent teardrop. It is blue at 470 nm and orange-red at 590/630 nm.

- **View buttons:** 3/4 (default), Sagittal, Coronal, Dorsal.
- **Layer toggles:** Head, Brain, V1–V3 (nerves, whisker pad and their labels), Light (the irradiance volume).
- **Theme button:** switches between light and dark colours.

**Approach presets.** Five buttons (see below). Choosing one loads its fiber and light parameters. Moving any slider clears the preset selection.

**Fiber sliders.**

| Control | Range | Step |
| --- | --- | --- |
| Entry AP | −7 to 12 mm | 0.1 |
| Entry ML | 0 to 5 mm | 0.05 |
| Tip DV (below bregma plane) | −7.5 to −1 mm | 0.05 |
| Tilt in coronal plane | −45° to 45° | 1° |
| Tilt in sagittal plane | −55° to 55° | 1° |

**Light controls.**

| Control | Options / range |
| --- | --- |
| Wavelength | 470 nm, 590 nm, 630 nm |
| Fiber core / NA | 100 µm / NA 0.22; 200 µm / NA 0.39; 400 µm / NA 0.50 |
| Power at tip | 0.5 to 30 mW, step 0.5 |
| Activation threshold | 0.2 to 5 mW/mm², step 0.1 |

**Readout panel.**

- Tip, stereotaxic (AP · ML · DV of the tip)
- Tip to TG centroid (mm)
- Tip position: "in ganglion" or "outside"
- Track length from skin (mm)
- Axial reach to threshold (mm)
- Ganglion above threshold (%), with a bar

**Depth ladder.** A marker moves down the depth-band column (below) to the current tip DV.

**Initial state on load:** AP −1.3, ML 2.0, tip DV −5.6, tilts 0°/0°, 590 nm, 200 µm core / NA 0.39, 10 mW, threshold 1.0 mW/mm². The first preset button is shown as selected, but these values differ from it: that preset's tip DV is −6.05.

## Approach presets

| # | Approach | Note (as shown) | Entry AP | Entry ML | Tip DV | Tilt coronal | Tilt sagittal | Wavelength | Fiber core / NA | Power |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Dorsal vertical to the ganglion | Scaled from the published rat cannula. Shortest path, worst collateral. | −1.3 | 2.0 | −6.05 | 0° | 0° | 590 nm | 200 µm / 0.39 | 10 mW |
| 2 | Stop short, illuminate down | Tip held 0.8–1 mm above the ganglion, off the skull base and its vessels. | −1.3 | 2.0 | −5.1 | 0° | 0° | 630 nm | 400 µm / 0.50 | 20 mW |
| 3 | Lateral-angled, past the midline structures | Enters further out and leans in; shifts the track off VPM but lengthens it. | −1.3 | 3.6 | −5.9 | −16° | 0° | 590 nm | 200 µm / 0.39 | 10 mW |
| 4 | Infraorbital nerve at the whisker pad | Ferrule on the skull, fiber tunnelled under the skin to the cheek. No brain track at all. | 7.8 | 2.5 | −4.6 | 38° | −14° | 590 nm | 400 µm / 0.50 | 15 mW |
| 5 | Open ganglion window (terminal) | Cortex aspirated, glass cuboid on the ganglion — light and objective straight down. | −1.2 | 1.95 | −5.35 | 0° | 0° | 590 nm | 400 µm / 0.50 | 5 mW |

The presets do not change the activation threshold, which stays at its current value (1.0 mW/mm² by default). Entry AP/ML are given on the bregma plane (DV 0). The tip is reached by travelling along the tilted fiber axis down to the tip DV.

## The numbers behind the model

Mouse TG coordinates are not standardised the way a cortical target is: the ganglion sits below the brain, off most atlas plates. These are the published anchors the geometry was fitted to. Treat them as a starting window to verify with a dye or electrode track in your own animals, not as a settled target.

### Landmark coordinates (mm from bregma, skull-flat)

| Landmark | AP | ML | DV | Source |
| --- | --- | --- | --- | --- |
| **TG centroid** (model) | −1.2 | ±1.9 | −6.0 | fitted to rows below |
| TG, V2 region (mouse) | 0.0 | 0.9 | −6.3 to −6.5 | Addgene AAV.PHP.S datasheet, electrode-verified |
| TG body, rat → mouse ×0.52 | −1.8 | 1.9 | −5.9 | KC et al. 2022 (rat −3.5 / 3.6 / −10.4) |
| V3 root, rat → mouse ×0.52 | −1.8 | 1.9 | −6.8 | Tang et al. 2014 (rat −3.5 / 3.6 / −12.0) |
| TG rostral pole (V1/V2 exit) | +0.3 | 1.0 | −5.8 | model |
| TG caudal pole (root to pons) | −2.6 | 2.4 | −6.3 | model |
| Infraorbital foramen (ION exit) | +7.8 | 2.5 | −4.9 | model, skull morphometry |
| Whisker pad centre | +10.5 | 3.0 | −4.8 | model |
| Eye centre | +4.3 | 3.2 | −2.9 | model |
| Lambda | −4.2 | 0 | 0 | Paxinos & Franklin, skull-flat |
| Interaural line (ear bars) | ≈ −3.9 | — | ≈ −5.0 | Paxinos & Franklin mouse atlas |

The "Tip to TG centroid" readout uses the model's centroid point at AP −1.15, ML 1.95, DV −6.05. The table rounds this to −1.2 / ±1.9 / −6.0.

### Depth bands: vertical track at ML ≈ 2.0, AP ≈ −1.3

| DV range (mm) | Structure |
| --- | --- |
| 0.0 … −0.35 | Skull · dura |
| −0.35 … −1.70 | Cortex — S1/S2 barrel field |
| −1.70 … −2.20 | Corpus callosum / external capsule |
| −2.20 … −3.30 | Hippocampus — CA1 to DG |
| **−3.30 … −4.50** | **Thalamus — VPM sits here** (highlighted) |
| −4.50 … −5.40 | Zona incerta / cerebral peduncle |
| −5.40 … −5.90 | Ventral brain surface |
| **−5.90 … −6.70** | **Trigeminal ganglion · Meckel's cave** (highlighted) |
| −6.70 … −7.40 | Skull base · petrous bone, carotid |

These are schematic depth bands for a straight dorsal penetration, not an atlas lookup. The point of the column is the thalamic band: a vertical track to the ganglion passes through VPM, the trigeminal relay.

## Notes

### What a vertical track actually costs (flagged)

From the bregma plane the ganglion is about 6 mm down. A straight 200 µm fiber to that depth passes through barrel cortex, the hippocampus and the VPM thalamus before it arrives, then sits against the petrous bone beside the internal carotid and the cavernous sinus. If the readout is trigeminal, the track runs through the relay being measured.

### Red-shifted opsins buy depth, not much depth

Switching the wavelength from 470 to 590–630 nm roughly halves the scattering coefficient and extends the axial reach by 30–50% at matched power: a fiber that reaches 1.3 mm in blue reaches about 1.7–2.0 mm in red. This is useful for stopping a millimetre short of the skull base, but not enough to illuminate the ganglion from the cortical surface.

### The periphery gives the same spikes

Ganglion somata and their peripheral terminals are the same neurons. Illuminating the infraorbital nerve where it leaves the infraorbital foramen (AP ≈ +7.8, ML ≈ 2.5, 1–2 mm under the skin) or the whisker pad drives the same afferents. It needs no brain track, works in an awake animal, and lets the receptive field define which afferents you recruit.

### Confirm the track before you trust the coordinates

The ganglion window preparation (cortex aspirated, glass cuboid on the ganglion) is terminal, but it is the cheapest way to see where your coordinates actually land. Failing that, run a dye-coated pipette to depth in a pilot animal and recover the ganglion.

## Model and caveats

Head, brain and ganglion are lofted schematic surfaces fitted to the landmarks above for a ~25 g adult C57BL/6. They are not segmented from an atlas volume, so read shapes as topology and spatial relations rather than as morphometry. Light spread uses the conical-emission plus Kubelka–Munk tissue model (Aravanis et al. 2007) with scattering coefficients of 10.3 mm⁻¹ at 470 nm and ~5.5 mm⁻¹ in the red. This is a first-order estimate that ignores absorption by blood and the bone at the skull base. Ganglion coverage is sampled over the modelled ganglion volume.

### Light model as implemented

Implementation details, taken from the page's code:

- Scattering coefficient S (mm⁻¹): **470 nm → 10.3, 590 nm → 6.0, 630 nm → 5.0**. The caveat text above gives "~5.5 mm⁻¹ in the red", the midpoint of the two red values.
- Core radius r = core diameter / 2 (mm). Emission half-angle θ = asin(NA / 1.36), using a tissue refractive index of 1.36 and capping the ratio at 0.99.
- ρ = r / tan θ.
- Irradiance at the fiber face: I₀ = P / (π r²) (mW/mm²).
- Irradiance on axis at distance z (mm) from the tip (conical spread × Kubelka–Munk scattering):

  ```
  I(z) = I0 * rho^2 / ( (S*z + 1) * (z + rho)^2 )
  ```

- **Axial reach** is the largest z in 0–12 mm with I(z) ≥ threshold, found by 48-step bisection. It is 0 if I at the face is already below threshold.
- **Light volume** (drawn only if reach > 0.02 mm) is a lathe whose radius at depth z is (r + z·tan θ)·√(1 − (z/z_max)³).
- **Ganglion coverage** (right ganglion only) uses a sample cloud of 31 stations along the ganglion midline × 10 angles × 4 radial shells (√((b+0.5)/4)), giving 1,240 points. For each point, z is the projection on the fiber axis. Points with z < −0.05 mm are skipped. Inside the cone radius r + z·tan θ the axial z is used; outside it the straight-line distance is used. A point counts if I ≥ threshold. Coverage = % of all samples.
- **"In ganglion"** means the nearest ganglion sample is within 0.14 mm of the tip.
- **Fiber direction** from the two tilts a (coronal) and b (sagittal): dML = sin a, dAP = sin b, dDV = −√max(0.02, 1 − dML² − dAP²), then normalised.
- **Track geometry.** The entry point is on the bregma plane (DV 0) at the given AP/ML. The tip is where the fiber axis reaches the set DV. The skin entry is found by ray-casting along the fiber axis against the head envelope, and "track length from skin" is the skin-to-tip distance.

### Model geometry (schematic, from the page's code)

Head envelope: elliptical cross-sections lofted along AP. w = half-width (ML), cz = section centre (DV), h = half-height (DV), all in mm.

| AP | 14.2 | 12.0 | 10.0 | 8.0 | 6.0 | 4.0 | 2.0 | 0.0 | −2.0 | −4.0 | −6.0 | −8.0 | −10.5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| w | 0.8 | 1.6 | 2.4 | 3.1 | 3.8 | 4.6 | 5.3 | 5.7 | 5.9 | 5.7 | 5.2 | 4.6 | 4.4 |
| cz | −5.0 | −4.9 | −4.8 | −4.8 | −4.6 | −4.3 | −3.9 | −3.6 | −3.4 | −3.3 | −3.4 | −3.7 | −4.2 |
| h | 0.8 | 1.6 | 2.3 | 3.0 | 3.6 | 4.2 | 4.6 | 4.7 | 4.7 | 4.6 | 4.3 | 4.0 | 3.9 |

Brain envelope (same parameters):

| AP | 5.3 | 4.5 | 3.5 | 2.6 | 2.0 | 1.0 | 0.0 | −1.5 | −3.0 | −4.5 | −6.0 | −7.5 | −9.2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| w | 0.6 | 1.4 | 1.9 | 1.9 | 2.6 | 3.4 | 4.1 | 4.7 | 5.0 | 4.9 | 4.3 | 2.6 | 1.6 |
| cz | −2.0 | −1.9 | −1.9 | −2.1 | −2.5 | −2.8 | −3.0 | −3.1 | −3.2 | −3.2 | −3.4 | −4.2 | −4.6 |
| h | 0.6 | 1.2 | 1.5 | 1.6 | 2.1 | 2.4 | 2.6 | 2.7 | 2.8 | 2.8 | 2.6 | 1.6 | 1.3 |

Trigeminal ganglion. The midline runs rostromedial → caudal through the points below (ML, AP, DV) and is swept with a flat elliptical cross-section. Half-width = 0.52 + 0.42·sin(π·clamp(0.95u + 0.05)) mm and half-thickness = 0.20 + 0.16·sin(π·clamp(u)) mm, where u = 0…1 along the midline. The ganglion is mirrored to the left side.

| Point | ML | AP | DV |
| --- | --- | --- | --- |
| 1 | 1.0 | 0.3 | −5.80 |
| 2 | 1.5 | −0.3 | −5.90 |
| 3 | 2.1 | −0.9 | −6.02 |
| 4 | 2.6 | −1.5 | −6.12 |
| 5 | 2.75 | −2.1 | −6.22 |
| 6 | 2.4 | −2.6 | −6.28 |

Nerve paths are tubes through these control points (ML, AP, DV) and are mirrored to the left side. Tube radius is 0.15 mm for V1/V2, 0.17 mm for V3 and 0.22 mm for the root.

| Branch | Control points (ML, AP, DV) |
| --- | --- |
| V1 | (1.05, 0.35, −5.78) → (1.5, 1.6, −5.0) → (2.1, 3.0, −4.0) → (2.9, 4.2, −3.3) |
| V2 (→ infraorbital n.) | (1.25, 0.1, −5.86) → (2.0, 1.8, −5.6) → (2.6, 3.6, −5.3) → (2.75, 5.6, −5.1) → (2.5, 7.8, −4.9) → (2.9, 9.6, −4.8) → (3.1, 11.0, −4.7) |
| V3 | (2.7, −1.6, −6.14) → (3.3, −1.7, −6.8) → (3.9, −1.3, −7.5) → (4.1, −0.2, −8.0) |
| Root (→ pons) | (2.4, −2.6, −6.28) → (1.9, −3.7, −6.0) → (1.35, −4.8, −5.6) |

Other features:

- **Whisker-pad follicles:** 5 rows (r) × 6 columns (c) per side, each a 0.16 mm sphere, at AP = 9.0 + 0.72c, DV = −3.5 − 0.62r − 0.06c, ML = 2.55 + 0.38·sin(0.5c) − 0.03r.
- **Eyes:** 1.55 mm radius spheres centred at ML ±3.2, AP 4.3, DV −2.9.
- **Reference marks:** a midline from AP 13.5 to −8 on the skull plane; dots at bregma (0, 0, 0) and lambda (AP −4.2); an interaural line at AP −3.9, DV −5.0, spanning ML ±7.5.

## References

1. KC E, Islam J, Kim S, Kim HK, Park YS (2022) *Pain relief in a trigeminal neuralgia model via optogenetic inhibition on trigeminal ganglion itself with flexible optic fiber cannula.* Front Cell Neurosci 16:880369. Rat TG cannula at AP −3.5, ML 3.6, DV −10.4; 250 µm core, NA 0.66, 473/589 nm at 10 mW.
2. Tang XF, Bi RY, Meng Z, Gan YH (2014) *Stereotaxis of mandibular nerve initial point of trigeminal ganglion in rats.* Chin J Dent Res 17:99–104. 3.5 mm posterior, 3.6 mm lateral, 12.0 mm below skull surface.
3. Addgene AAV Data Hub, *Labeling trigeminal ganglion afferents with PNS-trophic AAV* (AAV.PHP.S). Mouse TG V2 at AP 0, ML 0.9, DV 6.3–6.5, confirmed by tungsten electrode recording.
4. Bio-protocol (2025) *Hemispherectomy-based optical window for in vivo visualization of trigeminal ganglion neurons in mice.* Glass cuboid window on the TG after cortical aspiration.
5. Aravanis AM et al. (2007) J Neural Eng 4:S143. Light transmission model used for the irradiance volume.

The artifact lists these references without hyperlinks.
