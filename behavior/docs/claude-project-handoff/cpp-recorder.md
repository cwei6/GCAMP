> Copied on 10 Nov 2026 from the Claude project doc `claude/cpp-recorder.md` (last written 12 Oct 2026).

# CPP Recorder — status (2026-10-12)

## Decision (2026-10-12)
- Webcam check done: Chrome on the Samsung tablet lists only the built-in cameras, not the Logitech HD 720p. The tablet doesn't pass USB webcams to apps or the browser, and Chrome blocks web pages from using webcams over USB, so the web recorder can't use the Logitech on this tablet.
- Logitech has no Android app for its webcams. Its only mobile app (Logi Tune) is for Zone headsets. Its webcam software (Logi Tune, G HUB) is Windows/Mac only.
- Chosen setup: Logitech + the "USB Camera" app by ShenYao (Play Store, com.shenyaocn.android.usbcamera; Pro version removes ads). The app has its own webcam driver. The user starts/stops recordings in the app, and names, times and logs sessions by hand.
- Status: the app works on the tablet and shows the chamber (2026-10-12).
- Two issues seen in the first frame the user sent:
  - It was 320x240, so the app was still on a low default resolution.
  - The camera looked at the box from the front at a low angle.
- Advised for that frame:
  - Switch to 1280x720 MJPG.
  - Mount the camera overhead, centred over the corridor, roughly 55–60 cm above the box floor (estimate for a 720p Logitech with ~55–60° diagonal field of view and a box ~50 cm long). Adjust until both compartments and the corridor just fit.
  - Watch for glare from the lids.
- Settings recommended for CPP (exact menu names not verified):
  - 1280x720 in MJPG/MJPEG format at 30 fps. YUY2 at 720p only reaches about 7–10 fps over USB.
  - Audio off (no audio needed).
  - Timestamp overlay off, unless it sits outside the box.
  - Auto-delete of old videos and auto-segment/loop recording off, so each session stays one complete file.
- The web recorder stays usable on a PC (double-click D:\Behavior\CPP_recorder\index.html) or as a GitHub Pages copy. GitHub Pages was not set up.
- Possible next step offered: a small session timer + log page for the tablet (animal/phase/day, start times, end-of-session alarm, CSV) to run next to the camera app.

## Web recorder (built 2026-10-06 to 10-12)
Files in D:\Behavior\CPP_recorder:
- Source: index.html, store.js, camera.js, recorder.js, panels.js, ui.js, sw.js, manifest.webmanifest, icons, README.md.
- One-file version for hosting: github_pages\index.html and github_pages\sw.js. Rebuild after any edit with `python build_single.py`.

What it does:
- Live preview, camera picker (prefers the USB/external camera), 720p/30 fps, 1 / 2.5 / 5 Mbps, video only (never uses the microphone).
- Session fields: Animal ID, Phase (Hab, Pre, CondDrug, CondVeh, Post, Ext, Other), Day, Compartment (Dots / Stripes / both), Notes.
- Timed auto-stop (remembered per phase) or manual stop (two taps).
- File name: CPP_<YYYYMMDD>_<HHMMSS>_<animal>_<phase>[_D<day>][_Dots|_Stripes].mp4. On stop it downloads to Downloads.
- Crash-safe: 1 s chunks go to IndexedDB with a keyframe every second. After a crash or a closed tab, the recording is kept up to the last second and listed as Interrupted.
- Session log (IndexedDB), exportable as CSV: timings, stop reason, measured fps, worst 5-s fps, interruptions.
- Wake lock, dim-screen mode, framing guides, snapshot PNG, camera settings (exposure / WB / focus, remembered per camera), diagnostics.
- Single-instance lock, offline copy via service worker.

Constraints learned:
- Camera access needs https (or localhost / file://).
- Claude-hosted pages can't use a camera; this account has no camera capability for artifacts (checked 2026-10-12).
- The tablet can't be linked to a Claude session like the PC. The tablet app only gives Claude calendar, alarm, timer and location access.

Testing (headless Chromium + simulated 30 fps CPP scene):
- Timed and manual recordings decode cleanly at 1280x720 / 30 fps.
- Hard kill mid-recording kept all but about the last second.
- Two-tab lock, file:// mode, offline reload, camera settings restore and 400 / 800 / 1280 px layouts all pass.
- An independent code review found 7 issues; all are fixed and re-tested.
- The one-file build works as hosted, offline and installable.
- Never tested with the real webcam.
