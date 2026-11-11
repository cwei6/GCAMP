# Synchronized optostim + 2p imaging: bring-up plan

**Goal:** stimulate ReaChR (625 nm, fiber-coupled LED) *during* Scanbox 2p imaging of jGCaMP8s afferents, with PMTs protected and every stim pulse timestamped against frames — replacing the lift-the-objective, no-imaging 1-h protocol, and enabling stim-triggered ΔF/F analysis.

**Status (Oct 6 2026, evening):** software side fully verified and exonerated; **no PMT gating observed in any test**. Grab test 2 (opto_test_261006_001, 150 frames, LED 5 / GATE 40 / Count 50 committed with Enter, gate-mode options 0/1/Pulse all tried) was scanned frame-by-frame in both channels: zero gating signature; only a frame-0 recording-start artifact (93 lines ≈ 12 ms, one frame — wrong size/count to be the gate). Gate-check figure saved next to the run. **Next: probe LED_ENABLE electrically, then Neurolabware email.**

## Decoded from the source (Scanbox_Tower_v3-main)

- Each Stim field pushes its value to PSoC2 only via its own callback (Enter/focus change): LED `ttl_on` → `sb_optostim_period` (ms, max 50, floats OK); GATE `ttl_margin` → `sb_optostim_period2`; Count `edit83` → `sb_optostim_frames` (1–255). Nothing pushes at startup or Grab.
- Trigger button → `sb_optostim_trigger` (cmd `0xa2`), **never logged anywhere** (no info field, no event). Network `'e'` presses the same button. Timestamps must come from an EVENT-line wire or gated-band detection.
- `Trigger [2]DIN0` checkbox → `sb_optostim_extint` (arms external TTL trigger). Dropdown (`optostim_gate`, options 0/1/Pulse) → `sb_optostim_gate(0/1/2)`, a deadband-related gate mode; not required for triggering; semantics undocumented.
- Commands verified to go to the Slave port (COM9 = PSoC2, fw 2.14) via `sb_open`. `info.otparam` is optotune-only ([0 0 NaN] = single plane; no opto info).
- DC2200 rows need `led_id` set + USB; with `led_id=''` they throw the harmless `ledstate_Callback` errors.

## Test history

- Test 1 (261005_001, 300 fr): values likely uncommitted; no events in .mat (by design); .sbx too big to transfer (466 MB > 400 MB cap).
- Live tests: LED 5/GATE 40/Count 50 with Enter, all three gate modes — no band on-screen.
- Test 2 (261006_001, 150 fr): full quantitative scan, both channels, sensitive threshold — **no gating anywhere**. A working 40 ms × 50-frame gate would be a ~320-line × 50-frame block; largest real feature is the frame-0 artifact.

## Current fork — settle with the electrical probe

Wire PSoC2 **DIGITAL OUT 0** (SMA→BNC) → NI BNC-2110 **AI0**; NI MAX → device → Test Panels → Analog Input (continuous). Press Trigger in Scanbox **while Focus is scanning**.

- **Pulses visible** (≈3.3 V, 5 ms, ~15.6 Hz, ~3 s): stim engine works; the LED side of the project can proceed immediately (wire the driver!). The missing piece is only PMT gating → ask Neurolabware whether these PMT modules support gating / whether the gate line is wired. While at it, flip the dropdown 0 → 1 → Pulse between triggers and note how the waveform changes (undocumented otherwise).
- **No pulses**: engine/firmware question for Neurolabware.

Email draft: *"We have a Tower (Scanbox_Tower_v3, Master v1.14.1 / Slave v2.14.0) and are setting up the Optogenetics panel. Stim [LED|GATE|Count] values are committed and Trigger pressed during scanning, all three gate-mode options tried, but no PMT gating band appears in the image [and LED_ENABLE on PSoC2 DOUT0 does / does not pulse on trigger]. Does PMT gating require a specific PMT model or wiring on our unit, and what do gate modes 0/1/Pulse mean?"*

## Fallback if PMTs turn out non-gateable

The sync project is not blocked: with 625 nm stim and the green emission filter, running **ungated** at moderate LED power is viable (ch2 gain = 0 during opto runs; expect artifact lines during pulses instead of a clean dark band). Frame-locked timing + EVENT-line timestamps work identically. Gating remains preferable (protection + clean frames) — pursue with Neurolabware, don't wait on it.

## Timestamp logging (settled)

T-split trigger TTL → PSoC2 DIN0 (arm checkbox) and → PSoC1 EVENT 1; `evt_names{2}='opto'`, `imask=bin2dec('00001110')`. Trigger source: NI BNC-2110 (best fit) or network `'e'` from the stim PC — wire EVENT 1 either way.

## DC2200 hookup (unchanged)

Locate DC2200 (serial M00476084) + 625 nm fiber LED → USB → restore `sbconfig.led_id='USB0::0x1313::0x80C8::M00476084::INSTR'` (delete the blanking line; keep `led_stim_com=''`). Panel mode TTL, `TTL [H]` = safe current. PSoC2 DOUT0 → DC2200 mod input. Bench test away from scope first.

## Light delivery during imaging (unchanged)

Oblique fiber at window outside the 16x cone (400 µm, NA 0.39–0.5, 30–45°, printed holder on headbar posts); ~1–10 mW/mm² at window (prior protocol 3.3 mW @ 2 Hz reference). Upgrade: dichroic above objective (reflect 615–635, transmit 920 + green). Ch2 sees 625 nm directly.

## Analysis hooks

Stim frames from EVENT log (+ band cross-check if gating works); exclude gated lines on stim frames; keep the artifact flagger from mislabeling gate bands. Primary readout: stim-triggered responses.

## Open items

- [x] Firmware/panel check → 14.x, panel present
- [x] Decode trigger/dropdown/logging from source
- [x] Big-GATE test + dropdown sweep → no band
- [x] 150-frame Grab scanned → no gating anywhere (figure saved by the run)
- [ ] **NI probe of LED_ENABLE (PSoC2 DOUT0) — decides the fork**
- [ ] Neurolabware email (fill in probe result)
- [ ] Locate DC2200 + 625 LED; USB + `led_id`; SMA→BNC cable
- [ ] EVENT-1 T-split + config edit for timestamps
- [ ] Power at fiber tip/window; printed fiber holder; stim paradigm after pilot
