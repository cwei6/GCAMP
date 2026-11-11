r"""report — the MAIN self-contained HTML report per animal (GUI tab 8).

    python aqua_report.py <mouse folder>            -> <mouse>\aqua_analysis\<mouse>_report.html
    python aqua_report.py <mouse folder> --qc       -> <mouse>\aqua_analysis\<mouse>_qc.html   (afferent2p.qc_report)

Order of the main report (v32):
  QC paragraph + link to the QC report
  1. general representation (the lab's presentation measurements: sessions on
     the day's clock, activity by running / resting state with event overlays,
     ECDFs, acute / chronic summary of every measure)
  2. heatmap activity comparison (CFU x time, up / down / unchanged, tile maps)
  3. locomotion-evoked activity (trials, responses, responsive CFUs, rest-only)
  4. pupil measures   5. coupling GLMs
  6. behaviour-controlled change in event rate (table, verdicts, one scale)
  7. time course on the day's clock, matched speed
  8. per-experiment details   9. settings and tuning   10. methods
Everything comes from aqua_analysis\ (the analysis step) — nothing is
recomputed here. The per-recording EVIDENCE (detection maps, rasters, clips
and a montage of every event) is the separate QC report (qc_report.py), so
this file stays readable.

Figures are embedded (base64), so the file can be mailed or opened anywhere;
click a figure to zoom.
"""
from __future__ import annotations
import base64
import html as _html
import io
import json
import time
from pathlib import Path

import numpy as np

from .aqua_io import aqua_files, cfu_rules_of, read_cfu, read_curation, read_feature_table
from .timecourse import read_day_labels, day_of, cond_label

# colours (dataviz reference palette; categorical in fixed order)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, MUTED = "#0b0b0b", "#52514e", "#e6e5e1", "#9a9aa2"
KEPT, DROPPED, EXCL, CFU_COL = "#2a78d6", "#e34948", "#9a9aa2", "#eda100"


# ----------------------------------------------------------------------------- small helpers

def _esc(s) -> str:
    return _html.escape("" if s is None else str(s))


def _fig_uri(fig, dpi: int = 100, kind: str = "png", quality: int = 82) -> str:
    """matplotlib figure -> data URI; closes the figure."""
    import matplotlib.pyplot as plt
    buf = io.BytesIO()
    if kind == "jpg":
        fig.savefig(buf, format="jpg", dpi=dpi, pil_kwargs={"quality": quality, "optimize": True})
        mime = "image/jpeg"
    else:
        fig.savefig(buf, format="png", dpi=dpi)
        mime = "image/png"
    plt.close(fig)
    return f"data:{mime};base64," + base64.b64encode(buf.getvalue()).decode()


def _file_uri(path: Path) -> str | None:
    path = Path(path)
    if not path.exists():
        return None
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg", "gif": "image/gif"}.get(path.suffix.lower()[1:], "image/png")
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()


def _figure(uri: str | None, caption: str, cls: str = "", missing: str = "") -> str:
    if not uri:
        return f"<p class='small missing'>{_esc(missing or 'figure not available')}</p>" if (missing or caption) else ""
    return (f"<figure class='{cls}'><a class='zoom' href='#'><img loading='lazy' src='{uri}' alt='{_esc(caption)}'></a>"
            f"<figcaption>{caption} <span class='hint'>(click to zoom)</span></figcaption></figure>")


def _ratio(v, nd: int = 2) -> str:
    try:
        v = float(v)
    except Exception:
        return "–"
    if not np.isfinite(v):
        return "–"
    return f"×{v:.{nd}f}"


def _num(v, nd: int = 2, suffix: str = "") -> str:
    try:
        v = float(v)
    except Exception:
        return "–"
    if not np.isfinite(v):
        return "–"
    return f"{v:.{nd}f}{suffix}"


def _pct(v) -> str:
    try:
        v = float(v)
    except Exception:
        return "–"
    return f"{100 * v:.0f} %" if np.isfinite(v) else "–"


def _mmss(sec: float) -> str:
    sec = max(0.0, float(sec))
    return f"{int(sec // 60):02d}:{int(sec % 60):02d}"


def _short_run(run: str) -> str:
    import re
    m = re.search(r"(\d{6}_\d{3})", run)
    return m.group(1) if m else run


def _style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def _stretch(img: np.ndarray, lo: float = 1.0, hi: float = 99.7):
    a = np.asarray(img, np.float32)
    v0, v1 = np.percentile(a, lo), np.percentile(a, hi)
    if v1 <= v0:
        v1 = v0 + 1
    return v0, v1


# ----------------------------------------------------------------------------- per-run data

def _open_movie(tif: Path):
    import tifffile as tf
    try:
        arr = tf.memmap(str(tif), mode="r")
    except Exception:
        arr = tf.imread(str(tif))
    if arr.ndim == 2:
        arr = arr[None]
    return arr


def _mean_image(tif: Path, movie=None, n: int = 120) -> np.ndarray | None:
    import tifffile as tf
    mp = tif.with_name(tif.stem + "_meanproj.tif")
    if mp.exists():
        try:
            m = np.asarray(tf.imread(str(mp)), np.float32)
            if m.ndim == 2:
                return m
        except Exception:
            pass
    if movie is None:
        return None
    idx = np.unique(np.linspace(0, movie.shape[0] - 1, min(n, movie.shape[0])).astype(int))
    return np.mean(np.stack([np.asarray(movie[int(i)], np.float32) for i in idx]), axis=0)


def _crop_offset(raw_mean, exported_mean):
    """(dy, dx) of AQuA2's cropped mean image in the raw frame; template match, else centred."""
    if raw_mean is None or exported_mean is None:
        return None
    H, W = raw_mean.shape; h, w = exported_mean.shape
    if h > H or w > W:
        return (0, 0)
    try:
        from .tune import crop_offset
        off = crop_offset(raw_mean, exported_mean)
        if off is not None:
            return (int(off[0]), int(off[1]))
    except Exception:
        pass
    return ((H - h) // 2, (W - w) // 2)


def load_run(tif: Path, kept_events: "dict | None", log=print) -> dict:
    """Everything the QC section needs for one recording. kept_events: {event id:
    row of events_<exp>.csv} for this run (None = the analysis has not run)."""
    import pandas as pd
    import tifffile as tf
    tif = Path(tif)
    f = aqua_files(tif)
    d: dict = {"tif": tif, "stem": tif.stem, "short": _short_run(tif.stem), "files": f, "notes": []}
    meta = {}
    if f["json"].exists():
        try:
            meta = json.loads(f["json"].read_text())
        except Exception:
            meta = {}
    d["spf"] = float(meta.get("seconds_per_frame", 2.0) or 2.0)
    d["prereg"] = meta.get("registration") if isinstance(meta.get("registration"), dict) else None   # unpack --register (v33)
    gc = meta.get("global_correction") or {}
    d["frames"] = int(gc.get("frames_in_movie") or meta.get("frames_out") or 0)
    d["excluded_1based"] = [int(v) for v in (gc.get("excluded_frames_1based") or [])]
    d["flagged_1based"] = [int(v) for v in (gc.get("artefact_frames_1based") or gc.get("flagged_frames_1based") or [])]
    d["bad_frames_mode"] = gc.get("bad_frames") or ("filled" if gc.get("replaced_flagged") else "kept")
    # frames csv: recorded frame <-> tiff frame
    d["frames_csv"] = None
    fcsv = Path(str(tif.with_suffix("")) + "_frames.csv")
    if fcsv.exists():
        try:
            from .globalfix import read_frames_csv
            d["frames_csv"] = read_frames_csv(fcsv)
        except Exception:
            d["frames_csv"] = None
    # behaviour
    d["beh"] = None
    if f["behavior"].exists():
        try:
            d["beh"] = pd.read_csv(f["behavior"])
        except Exception:
            d["beh"] = None
    # detected events with location (CFU batch export) — the feature table has no location
    ev_csv = Path(str(tif.with_suffix("")) + "_AQuA2_events.csv")
    d["events"] = None
    if ev_csv.exists():
        try:
            e = pd.read_csv(ev_csv)
            if len(e) and "event" in e.columns:
                e["event"] = e["event"].astype(int)
                d["events"] = e
            else:
                d["events"] = e
        except Exception as ex:
            d["notes"].append(f"could not read {ev_csv.name}: {ex}")
    elif f["features"].exists():
        try:
            ft = read_feature_table(f["features"])
            d["events"] = pd.DataFrame({"event": ft["event"].astype(int), "t0_frame": ft.get("t0_frame", np.nan)})
            d["notes"].append("no per-event locations (_AQuA2_events.csv missing: run run_cfu_batch.m again) — no footprint boxes, montages or clips")
        except Exception:
            d["events"] = None
    # CFUs and curation
    d["cfu"] = None; d["curation"] = None; d["event_cfu"] = {}
    if f["cfu"].exists():
        try:
            c = read_cfu(f["cfu"])
            d["cfu"] = c
            for i, evs in enumerate(c["events"]):
                for ev in evs:
                    d["event_cfu"].setdefault(int(ev), i + 1)
        except Exception as ex:
            d["notes"].append(f"CFU file unreadable: {ex}")
        try:
            d["curation"] = read_curation(tif)
        except Exception:
            d["curation"] = None
    rejected = set(int(i) for i in (d["curation"] or {}).get("rejected", [])) if d["curation"] and not d["curation"].get("stale") else set()
    d["rejected_cfu"] = rejected
    # kept / dropped per detected event
    d["kept"] = kept_events
    if d["events"] is not None and len(d["events"]):
        why = []
        for ev, t0 in zip(d["events"]["event"].to_numpy(int), d["events"].get("t0_frame", pd.Series([np.nan] * len(d["events"]))).to_numpy(float)):
            cfu = d["event_cfu"].get(int(ev))
            if kept_events is not None and int(ev) in kept_events:
                why.append("kept")
            elif cfu is not None and cfu in rejected:
                why.append("noise CFU")
            elif np.isfinite(t0) and d["frames_csv"] is not None and _in_excluded(int(t0), d):
                why.append("excluded frames")
            elif kept_events is None:
                why.append("kept")               # no analysis yet: nothing was dropped
            else:
                why.append("dropped")
        d["events"]["status"] = why
    # images: AQuA2's exported mean image (its crop) and the raw mean projection
    d["meanimg"] = None
    mi = Path(str(tif.with_suffix("")) + "_AQuA2_meanimg.tif")
    if mi.exists():
        try:
            d["meanimg"] = np.asarray(tf.imread(str(mi)), np.float32)
        except Exception:
            d["meanimg"] = None
    # the movie: AQuA2's own registered movie out of the _AQuA2.mat when it can be read (the frames the
    # events were detected on, in the CFU maps' coordinates), else the binned tiff as recorded
    d["movie"] = None; d["movie_src"] = None
    try:
        from .aqua_io import open_movie
        if (f["res"].exists() and f["res"].stat().st_size > 0) or (tif.exists() and tif.stat().st_size > 0):
            d["movie"], d["movie_src"] = open_movie(tif, prefer_registered=True, log=lambda m: d["notes"].append(m.strip()))
    except Exception as ex:
        d["notes"].append(f"movie not readable: {ex}")
    d["raw_mean"] = _mean_image(tif, d["movie"] if d["movie_src"] == "raw" else None)
    if d["meanimg"] is None and d["cfu"] is not None and d["cfu"].get("datPro") is not None and d["cfu"]["datPro"].ndim == 2:
        d["meanimg"] = np.asarray(d["cfu"]["datPro"], np.float32)
    if d["movie_src"] == "registered":
        d["offset"] = (0, 0)                       # same frame as the CFU maps and the event boxes
        d["crop_in_raw"] = (_crop_offset(d["raw_mean"], d["meanimg"]) if d["meanimg"] is not None else None)
    else:
        d["offset"] = (_crop_offset(d["raw_mean"], d["meanimg"]) if d["meanimg"] is not None else None) or (0, 0)
        d["crop_in_raw"] = d["offset"]
    return d


def _in_excluded(t0_tiff: int, d: dict) -> bool:
    fr = d.get("frames_csv")
    if fr is None:
        return False
    tf_ = fr["tiff_frame"]; ex = fr["excluded"]
    hit = np.where(tf_ == t0_tiff)[0]
    if hit.size == 0:
        return False
    i = int(hit[0])
    lo, hi = max(0, i - 1), min(len(ex), i + 2)
    return bool(ex[lo:hi].any())


def _tiff_to_recorded_time(t_tiff: int, d: dict) -> float:
    """Seconds since the recording start for a tiff frame (1-based)."""
    fr = d.get("frames_csv")
    if fr is not None:
        hit = np.where(fr["tiff_frame"] == int(t_tiff))[0]
        if hit.size:
            return float(fr["frame"][int(hit[0])] - 1) * d["spf"]
    return (int(t_tiff) - 1) * d["spf"]


# ----------------------------------------------------------------------------- per-run figures

def detection_map(d: dict) -> str | None:
    """Mean image with every detected event's box and the CFU outlines."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    img = d.get("meanimg")
    if img is None:
        img = d.get("raw_mean")
    if img is None:
        return None
    ev = d.get("events")
    H, W = img.shape
    fig, ax = plt.subplots(figsize=(9.2, 9.2 * H / W + 0.25))
    v0, v1 = _stretch(img, 0.5, 99.8)
    ax.imshow(img, cmap="gray", vmin=v0, vmax=v1, interpolation="nearest")
    ax.set_axis_off()
    # CFU outlines
    c = d.get("cfu")
    if c is not None and c.get("maps") is not None and c["maps"].ndim == 3 and c["maps"].shape[1:] == img.shape:
        for i in range(c["n"]):
            m = c["maps"][i]
            if not np.any(m > 0):
                continue
            thr = 0.3 * float(m.max())
            try:
                ax.contour(m, levels=[thr], colors=[DROPPED if (i + 1) in d["rejected_cfu"] else CFU_COL], linewidths=0.6,
                           linestyles=["dashed" if (i + 1) in d["rejected_cfu"] else "solid"])
            except Exception:
                pass
    if ev is not None and len(ev) and {"y0", "y1", "x0", "x1"} <= set(ev.columns):
        for r in ev.itertuples(index=False):
            st = getattr(r, "status", "kept")
            col = KEPT if st == "kept" else (DROPPED if st == "noise CFU" else EXCL)
            x0, x1, y0, y1 = float(r.x0) - 1, float(r.x1) - 1, float(r.y0) - 1, float(r.y1) - 1
            ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0 + 1, y1 - y0 + 1, fill=False, ec=col, lw=0.9,
                                   ls="-" if st == "kept" else "--"))
            ax.text(x1 + 1.5, y0 - 1.5, str(int(r.event)), color=col, fontsize=8, fontweight="bold", va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
    ax.set_title(f"{d['short']} — mean image with the detected events (boxes, numbered) and the CFU outlines", fontsize=9, loc="left", color=INK)
    fig.tight_layout(pad=0.3)
    return _fig_uri(fig, dpi=100)


def raster_figure(d: dict, trials: list | None = None) -> str | None:
    """Running speed, pupil and event onsets on one time axis, excluded frames shaded;
    trials = [(onset recorded frame 0-based, 'strict' | 'relaxed')] marked on the speed axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ev = d.get("events")
    beh = d.get("beh")
    spf = d["spf"]
    n = d["frames"] or (int(beh["frame"].max()) + 1 if beh is not None else 0)
    if not n:
        return None
    minutes = n * spf / 60.0
    lane_guess = 1
    if ev is not None and len(ev) and {"t0_frame", "t1_frame"} <= set(ev.columns):
        t0s = ev["t0_frame"].to_numpy(float); t1s = ev["t1_frame"].to_numpy(float)
        lane_guess = int(min(8, max(1, np.max([np.sum((t0s <= t) & (t1s >= t)) for t in t0s[np.isfinite(t0s)]] or [1]))))
    fig, axes = plt.subplots(3, 1, figsize=(11, 3.2 + 0.28 * lane_guess), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.2, 1.2, 0.6 + 0.3 * lane_guess], hspace=0.12))
    ax_v, ax_p, ax_e = axes
    for ax in axes:
        _style(ax)
    # excluded frames (recorded) as spans
    fr = d.get("frames_csv")
    spans = []
    if fr is not None:
        exc = np.where(fr["excluded"])[0]
        if exc.size:
            start = prev = int(exc[0])
            for i in exc[1:]:
                if int(i) != prev + 1:
                    spans.append((start, prev)); start = int(i)
                prev = int(i)
            spans.append((start, prev))
    for a, b in spans:
        for ax in axes:
            ax.axvspan(a * spf / 60, (b + 1) * spf / 60, color="#fde2e0", lw=0, zorder=0)
    if beh is not None and "velocity_mean" in beh.columns:
        t = beh["time_s"].to_numpy(float) / 60.0 if "time_s" in beh.columns else (beh["frame"].to_numpy(float) + 0.5) * spf / 60
        ax_v.plot(t, beh["velocity_mean"].to_numpy(float), color=SERIES[0], lw=1.0)
        ax_v.set_ylabel("speed", fontsize=8, color=INK2)
        if trials:
            vmax = float(np.nanmax(beh["velocity_mean"].to_numpy(float))) if np.isfinite(beh["velocity_mean"].to_numpy(float)).any() else 1.0
            for onset, kind in trials:
                ax_v.plot((onset + 0.5) * spf / 60, vmax * 1.08, marker="v", ms=6 if kind == "strict" else 4.5,
                          color=INK if kind == "strict" else MUTED, lw=0, clip_on=False)
            ax_v.set_ylim(top=vmax * 1.2)
        if "pupil_diam_px_mean" in beh.columns and np.isfinite(beh["pupil_diam_px_mean"].to_numpy(float)).any():
            ax_p.plot(t, beh["pupil_diam_px_mean"].to_numpy(float), color=SERIES[2], lw=1.0)
            ax_p.set_ylabel("pupil px", fontsize=8, color=INK2)
        else:
            ax_p.text(0.01, 0.5, "no pupil data", transform=ax_p.transAxes, fontsize=8, color=MUTED, va="center")
    else:
        ax_v.text(0.01, 0.5, "no behaviour csv", transform=ax_v.transAxes, fontsize=8, color=MUTED, va="center")
    n_kept = n_drop = 0
    lanes_end: list[float] = []
    n_lanes = 1
    if ev is not None and len(ev) and "t0_frame" in ev.columns:
        items = []
        for r in ev.itertuples(index=False):
            t0 = float(r.t0_frame)
            if not np.isfinite(t0):
                continue
            t1v = float(getattr(r, "t1_frame", np.nan))
            t1 = t1v if np.isfinite(t1v) else t0
            tsec = _tiff_to_recorded_time(int(t0), d)
            items.append((tsec / 60, tsec / 60 + max(t1 - t0 + 1, 1) * spf / 60, getattr(r, "status", "kept"), int(r.event)))
        items.sort()
        for a, b, st, eid in items:
            lane = next((i for i, e in enumerate(lanes_end) if e + 0.05 < a), None)
            if lane is None:
                lanes_end.append(b); lane = len(lanes_end) - 1
            else:
                lanes_end[lane] = b
            col = KEPT if st == "kept" else (DROPPED if st == "noise CFU" else EXCL)
            if st == "kept":
                n_kept += 1
            else:
                n_drop += 1
            y = lane + 0.5
            ax_e.plot([a, b], [y, y], color=col, lw=3.5, solid_capstyle="butt", alpha=0.85)
            ax_e.plot([a], [y], marker="|", color=col, ms=9, mew=1.4)
            if len(items) <= 60:
                ax_e.text(a, y + 0.22, str(eid), fontsize=6.5, color=col, ha="center", va="bottom")
        n_lanes = max(1, len(lanes_end))
    ax_e.set_ylim(0, n_lanes + 0.3); ax_e.set_yticks([])
    ax_e.set_ylabel("events", fontsize=8, color=INK2)
    ax_e.set_xlabel("minutes", fontsize=8, color=INK2)
    ax_e.set_xlim(0, minutes)
    ax_e.grid(visible=False)
    ax_v.set_title(f"{d['short']} — running speed, pupil size and the detected events in time ({n_kept} kept, {n_drop} dropped; "
                   f"excluded frames shaded)", fontsize=9, loc="left", color=INK)
    return _fig_uri(fig, dpi=100)


def _box_raw(r, d: dict, movie_shape):
    """Event box in raw-frame 0-based pixel indices (inclusive), clipped."""
    dy, dx = d.get("offset") or (0, 0)
    H, W = movie_shape[-2], movie_shape[-1]
    y0 = int(round(float(r.y0))) - 1 + dy; y1 = int(round(float(r.y1))) - 1 + dy
    x0 = int(round(float(r.x0))) - 1 + dx; x1 = int(round(float(r.x1))) - 1 + dx
    return (int(np.clip(y0, 0, H - 1)), int(np.clip(y1, 0, H - 1)), int(np.clip(x0, 0, W - 1)), int(np.clip(x1, 0, W - 1)))


def _crop_window(box, movie_shape, min_half: int = 26, max_half: int = 90):
    y0, y1, x0, x1 = box
    H, W = movie_shape[-2], movie_shape[-1]
    cy, cx = (y0 + y1) / 2, (x0 + x1) / 2
    half = int(np.clip(max(y1 - y0, x1 - x0) * 1.1, min_half, max_half))
    cy0 = int(np.clip(round(cy) - half, 0, max(H - 2 * half, 0))); cx0 = int(np.clip(round(cx) - half, 0, max(W - 2 * half, 0)))
    return cy0, min(cy0 + 2 * half, H), cx0, min(cx0 + 2 * half, W)


def _box_trace(movie, box) -> np.ndarray:
    y0, y1, x0, x1 = box
    sub = movie[:, y0:y1 + 1, x0:x1 + 1]
    return np.asarray(sub, np.float32).reshape(sub.shape[0], -1).mean(axis=1)


def event_montage(d: dict, r, n_frames: int = 8) -> tuple[str | None, dict]:
    """Movie frames around one event with its box, and the box's dF/F trace with
    the running speed underneath. Returns (jpg data URI, metrics)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    movie = d["movie"]
    T = movie.shape[0]
    t0 = int(round(float(r.t0_frame))); t1 = int(round(float(getattr(r, "t1_frame", t0))))
    t1 = max(t1, t0)
    box = _box_raw(r, d, movie.shape)
    cy0, cy1, cx0, cx1 = _crop_window(box, movie.shape)
    # frames to show: before, onset, inside, end, after
    inside = np.unique(np.round(np.linspace(t0, t1, min(n_frames - 3, max(2, t1 - t0 + 1)))).astype(int))
    want = [t0 - 3, t0 - 1] + list(inside) + [t1 + 2, t1 + 6]
    frames = sorted({int(np.clip(t, 1, T)) for t in want})
    if len(frames) > n_frames + 1:
        frames = frames[:n_frames + 1]
    crops = [np.asarray(movie[t - 1, cy0:cy1, cx0:cx1], np.float32) for t in frames]
    stack = np.stack(crops)
    v0, v1 = _stretch(stack, 1.0, 99.8)
    # trace
    trace = _box_trace(movie, box)
    f0 = float(np.percentile(trace, 20))
    dff = (trace - f0) / max(f0, 1e-6)
    w0, w1 = max(1, t0 - 30), min(T, t1 + 30)
    tt = np.arange(w0, w1 + 1)
    spf = d["spf"]
    fig = plt.figure(figsize=(11, 4.3))
    gs = fig.add_gridspec(2, len(frames), height_ratios=[1.0, 1.15], hspace=0.3, wspace=0.08, left=0.045, right=0.99, top=0.91, bottom=0.12)
    for k, (t, cr) in enumerate(zip(frames, crops)):
        ax = fig.add_subplot(gs[0, k])
        ax.imshow(cr, cmap="gray", vmin=v0, vmax=v1, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(GRID)
        if t0 <= t <= t1:
            y0, y1, x0, x1 = box
            ax.add_patch(Rectangle((x0 - cx0 - 0.5, y0 - cy0 - 0.5), x1 - x0 + 1, y1 - y0 + 1, fill=False, ec=SERIES[1], lw=1.2))
        rel = (t - t0) * spf
        ax.set_title(("onset" if t == t0 else f"{rel:+.0f} s"), fontsize=8, color=INK2 if t != t0 else SERIES[1], pad=3)
    gs_b = gs[1, :].subgridspec(2, 1, height_ratios=[1.0, 0.28], hspace=0.08)
    ax_t = fig.add_subplot(gs_b[0]); ax_s = fig.add_subplot(gs_b[1], sharex=ax_t)
    _style(ax_t); _style(ax_s)
    ax_t.axvspan((t0 - 1) * spf, t1 * spf, color="#fdebe3", lw=0, zorder=0)
    ax_s.axvspan((t0 - 1) * spf, t1 * spf, color="#fdebe3", lw=0, zorder=0)
    ax_t.plot((tt - 1) * spf, dff[w0 - 1:w1], color=INK, lw=1.2)
    for t in frames:
        ax_t.axvline((t - 1) * spf, color=GRID, lw=0.7, zorder=0)
    ax_t.set_ylabel("dF/F (box)", fontsize=8, color=INK2)
    ax_t.tick_params(labelbottom=False)
    ax_s.set_xlabel("seconds since the recording start (movie frames)", fontsize=8, color=INK2)
    ax_t.set_xlim((w0 - 1) * spf, (w1 - 1) * spf)
    beh = d.get("beh")
    drawn = False
    if beh is not None and "velocity_mean" in beh.columns:
        tb = beh["time_s"].to_numpy(float) if "time_s" in beh.columns else (beh["frame"].to_numpy(float) + 0.5) * spf
        vb = beh["velocity_mean"].to_numpy(float)
        sel = (tb >= (w0 - 1) * spf) & (tb <= (w1 - 1) * spf)
        if sel.any():
            ax_s.fill_between(tb[sel], 0, vb[sel], color=SERIES[0], alpha=0.5, lw=0, step="mid")
            drawn = True
    ax_s.set_ylabel("speed", fontsize=7, color=INK2)
    ax_s.set_yticks([])
    if not drawn:
        ax_s.text(0.005, 0.5, "no behaviour csv", transform=ax_s.transAxes, fontsize=7, color=MUTED, va="center")
    metrics = {"t0": t0, "t1": t1, "onset_s": _tiff_to_recorded_time(t0, d), "dur_s": (t1 - t0 + 1) * spf,
               "peak_dff_box": float(np.nanmax(dff[t0 - 1:t1])) if t1 >= t0 else float("nan")}
    return _fig_uri(fig, dpi=80, kind="jpg", quality=80), metrics


def event_clip(d: dict, r, fps: float = 6.0, max_frames: int = 40, pad_frames: int = 5, scale_to: int = 220) -> str | None:
    """Animated clip (GIF data URI): the movie around one event, cropped, with the detection box."""
    from PIL import Image, ImageDraw
    movie = d["movie"]
    T = movie.shape[0]
    t0 = int(round(float(r.t0_frame))); t1 = max(int(round(float(getattr(r, "t1_frame", t0)))), t0)
    box = _box_raw(r, d, movie.shape)
    cy0, cy1, cx0, cx1 = _crop_window(box, movie.shape, min_half=30, max_half=100)
    a, b = max(1, t0 - pad_frames), min(T, t1 + pad_frames)
    frames = list(range(a, b + 1))
    if len(frames) > max_frames:
        keep = np.unique(np.round(np.linspace(0, len(frames) - 1, max_frames)).astype(int))
        frames = [frames[i] for i in keep]
    crops = np.stack([np.asarray(movie[t - 1, cy0:cy1, cx0:cx1], np.float32) for t in frames])
    v0, v1 = _stretch(crops, 1.0, 99.8)
    k = max(1, int(round(scale_to / max(crops.shape[1], crops.shape[2]))))
    spf = d["spf"]
    ims = []
    y0, y1, x0, x1 = box
    for t, cr in zip(frames, crops):
        g = np.clip((cr - v0) / (v1 - v0), 0, 1)
        im = Image.fromarray((g * 255).astype(np.uint8)).convert("RGB")
        if k > 1:
            im = im.resize((im.width * k, im.height * k), Image.NEAREST)
        dr = ImageDraw.Draw(im)
        if t0 <= t <= t1:
            dr.rectangle([(x0 - cx0) * k, (y0 - cy0) * k, (x1 - cx0 + 1) * k - 1, (y1 - cy0 + 1) * k - 1], outline=(235, 104, 52), width=2)
        rel = (t - t0) * spf
        label = f"#{int(r.event)}  {rel:+.0f} s" + ("  onset" if t == t0 else "")
        dr.rectangle([0, 0, 8 + 6 * len(label), 13], fill=(0, 0, 0))
        dr.text((3, 1), label, fill=(255, 255, 255))
        ims.append(im.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
    buf = io.BytesIO()
    ims[0].save(buf, format="GIF", save_all=True, append_images=ims[1:], duration=int(1000 / fps), loop=0, optimize=True, disposal=1)
    return "data:image/gif;base64," + base64.b64encode(buf.getvalue()).decode()


def choose_clip_events(ev, n: int) -> list[int]:
    """Which events get a clip: kept events spread over the amplitude range (largest first)."""
    if ev is None or len(ev) == 0 or n <= 0:
        return []
    kept = ev[ev.get("status", "kept") == "kept"] if "status" in ev.columns else ev
    if len(kept) == 0:
        kept = ev
    if "dff_max" in kept.columns and np.isfinite(kept["dff_max"].to_numpy(float)).any():
        order = kept.sort_values("dff_max", ascending=False)
    else:
        order = kept
    if len(order) <= n:
        return [int(v) for v in order["event"]]
    idx = np.unique(np.round(np.linspace(0, len(order) - 1, n)).astype(int))
    return [int(order["event"].iloc[int(i)]) for i in idx]


# ----------------------------------------------------------------------------- sections

def _read_csv(p: Path):
    import pandas as pd
    try:
        return pd.read_csv(p) if p.exists() else None
    except Exception:
        return None


def _cond_name(c: str, day_labels: dict) -> str:
    c = str(c)
    lab = cond_label(c, day_labels) if c.startswith("d") else ""
    return f"{c} ({lab})" if lab else c


def _f(r, k, default=np.nan):
    try:
        v = float(r.get(k, default))
    except (TypeError, ValueError):
        return default
    return v


def _verdict(r) -> str:
    """One plain sentence per experiment row of experiments_summary.csv, scale-aware:
    counts and the exact interval decide what may be said, not the bare ratio."""
    parts = []
    rr = _f(r, "ratio_raw"); rb, rc = _f(r, "rate_baseline"), _f(r, "rate_condition")
    nb, nc = _f(r, "n_events_baseline"), _f(r, "n_events_condition")
    lo, hi = _f(r, "ratio_ci_lo"), _f(r, "ratio_ci_hi")
    d = _f(r, "delta_rate")
    low = str(r.get("low_counts", "")).lower() in ("true", "1", "1.0")
    nb_s, nc_s = (f"{int(nb)}" if np.isfinite(nb) else "?"), (f"{int(nc)}" if np.isfinite(nc) else "?")
    ci = ""
    if np.isfinite(lo):
        ci = f"; 95 % interval of the ratio {lo:.2f}–{'∞' if not np.isfinite(hi) else f'{hi:.2f}'}"
    up = np.isfinite(lo) and lo > 1
    down = np.isfinite(hi) and hi < 1
    if np.isfinite(nb) and np.isfinite(nc) and nb + nc == 0:
        parts.append("no events in either condition")
    elif low:
        head = (f"{nb_s} → {nc_s} events ({rb:.2f} → {rc:.2f} events/min, difference {d:+.2f}/min{ci})")
        if up or down:
            parts.append(head + f": the {'increase' if up else 'decrease'} itself is unlikely to be chance, but on so few events its SIZE "
                                f"cannot be compared with a well-sampled experiment — read the rate difference, not the ratio")
        else:
            parts.append(head + ": too few events for a conclusion")
    else:
        if np.isfinite(rr) and rb > 0:
            word = "up" if rr >= 1.25 else ("down" if rr <= 0.8 else "about the same")
            parts.append(f"rate {word} ×{rr:.2f} ({rb:.2f} → {rc:.2f} events/min, {nb_s} → {nc_s} events, difference {d:+.2f}/min{ci})")
        else:
            parts.append(f"{nb_s} → {nc_s} events ({rb:.2f} → {rc:.2f} events/min{ci})")
        if up:
            parts.append("an increase the interval supports")
        elif down:
            parts.append("a decrease the interval supports")
        else:
            parts.append("not distinguishable from no change (the interval includes 1)")
    sb, sc = _f(r, "still_frac_baseline"), _f(r, "still_frac_condition")
    comp = r.get("behaviour_comparable", None)
    if np.isfinite(sb) and np.isfinite(sc) and abs(sb - sc) >= 0.3:
        parts.append(f"the mouse was still {100 * sb:.0f} % of the baseline but {100 * sc:.0f} % of the condition — a behaviour difference, "
                     f"so trust the matched columns over the raw numbers")
    if comp is not None and str(comp).lower() in ("false", "0", "0.0"):
        parts.append("at matched running speed the conditions are NOT COMPARABLE (too little overlap in behaviour)")
    ra = _f(r, "ratio_adjusted"); alo = _f(r, "ratio_adjusted_ci_lo"); ahi = _f(r, "ratio_adjusted_ci_hi")
    if np.isfinite(ra):
        aci = f" (95 % CI {alo:.2f}–{ahi:.2f})" if np.isfinite(alo) and np.isfinite(ahi) else ""
        if np.isfinite(alo) and np.isfinite(ahi) and (alo > 1 or ahi < 1):
            parts.append(f"behaviour-adjusted ×{ra:.2f}{aci}: the change holds with the behaviour taken out")
        else:
            parts.append(f"behaviour-adjusted ×{ra:.2f}{aci}: not distinguishable from no change")
    elif not low:
        parts.append("behaviour-adjusted ratio not computed (too few events or bins)")
    # how the two conditions were detected (aqua_tune v36: matched sensitivity / one common setting), and the ratio with
    # the detector's own expected false alarms subtracted
    rule = str(r.get("detection_rule", "") or "")
    if rule and rule != "nan":
        sbase, scond = str(r.get("setting_baseline", "") or ""), str(r.get("setting_condition", "") or "")
        rcb, rcc = _f(r, "recall_baseline"), _f(r, "recall_condition")
        fb, fc = _f(r, "fp_rate_baseline"), _f(r, "fp_rate_condition")
        comparable = str(r.get("detection_comparable", "True")).lower() not in ("false", "0", "0.0")
        if sbase and sbase != "nan":
            how = "matched sensitivity" if rule == "matched" else "one common setting"
            txt = (f"detection ({how}): {sbase} → {scond}"
                   + (f", recall {rcb:.2f} / {rcc:.2f}" if np.isfinite(rcb) and np.isfinite(rcc) else "")
                   + (f", false alarms {fb:.2f} / {fc:.2f} per min" if np.isfinite(fb) and np.isfinite(fc) else ""))
            rf = _f(r, "ratio_fpcorr")
            if np.isfinite(rf):
                txt += f"; false-alarm-subtracted ratio ×{rf:.2f}"
            parts.append(txt)
        if not comparable:
            parts.append("DETECTION NOT COMPARABLE: the two conditions could not be detected with the same sensitivity at any clean "
                         "setting — this change is not interpretable as it stands")
    return "; ".join(parts) + "."


def _exp_tag(exp: str, exp_meta: dict, day_labels: dict) -> str:
    """'260917 (sham) · FOV1' for an optostim experiment, 'FOV1' for a Treg one."""
    m = exp_meta.get(exp, {})
    days = m.get("days", [])
    fov = f"FOV{m.get('fov')}" if m.get("fov") is not None else ""
    if len(days) == 1:
        d = days[0]
        return f"{d}{' (' + day_labels[d] + ')' if day_labels.get(d) else ''} · {fov}"
    return fov


def section_comparison(mouse: str, paradigm_of: dict, A: Path, day_labels: dict, out_uris: dict, exp_meta: dict | None = None) -> str:
    """The result: one table per paradigm layout + the pipeline's figures."""
    exp_meta = exp_meta or {}
    es = _read_csv(A / "experiments_summary.csv")
    cs = _read_csv(A / "condition_summary.csv")
    h = ["<h2 id='result'>6. Behaviour-controlled change in event rate</h2>"]
    if es is None or len(es) == 0:
        h.append("<p class='missing'>No aqua_analysis\\experiments_summary.csv yet — run the analysis (tab 7) first; this section is empty.</p>")
        return "\n".join(h)
    spont = any(p == "spontaneous" for p in paradigm_of.values())
    h.append("<p>Each experiment compares a condition with its baseline — " + ("the days of a FOV with its first day (Treg)" if spont else
             "post with pre, per day and FOV (optostim)") + ". The raw ratio of event rates is shown first, then the same change with the "
             "behaviour controlled: on still frames only, on running frames only, at matched pupil size, at matched running speed "
             "(only where both conditions visited the same speeds) and by a regression that takes speed, running and pupil out "
             "(behaviour-adjusted, with a 95 % confidence interval). A change that survives these controls is not a change of "
             "behavioural state.</p>")
    head = ("<tr><th>Experiment</th><th>Day · FOV</th><th>Baseline</th><th>Condition</th><th>Baseline events / min</th><th>Condition events / min</th>"
            "<th>Difference /min</th><th>Raw ratio (95 % interval)</th><th>Modulation index</th><th>Still only</th><th>Running only</th><th>Pupil-matched</th>"
            "<th>Speed-matched</th><th>Adjusted (95 % CI)</th><th>Common region</th></tr>")
    rows = []
    verdicts = []
    cond_lookup = {}
    if cs is not None:
        for r in cs.to_dict("records"):
            cond_lookup[(r["experiment"], r["condition"])] = r
    for r in es.to_dict("records"):
        b = cond_lookup.get((r["experiment"], r["baseline"]), {}); c = cond_lookup.get((r["experiment"], r["condition"]), {})
        adj = _ratio(r.get("ratio_adjusted"))
        lo, hi = r.get("ratio_adjusted_ci_lo", np.nan), r.get("ratio_adjusted_ci_hi", np.nan)
        if adj != "–" and np.isfinite(float(lo)) and np.isfinite(float(hi)):
            adj += f" <span class='small'>({float(lo):.2f}–{float(hi):.2f})</span>"
        comp = r.get("behaviour_comparable", None)
        sm = _ratio(r.get("ratio_speed_matched"))
        if comp is not None and str(comp).lower() in ("false", "0", "0.0"):
            sm = "<span class='warn'>not comparable</span>"
        lo, hi = _f(r, "ratio_ci_lo"), _f(r, "ratio_ci_hi")
        raw = _ratio(r.get("ratio_raw"))
        if raw == "–" and np.isfinite(_f(r, "log2_ratio_pc")) and _f(r, "n_events_condition") > 0:
            raw = "∞ (no baseline events)"
        if np.isfinite(lo):
            raw += f" <span class='small'>({lo:.2f}–{'∞' if not np.isfinite(hi) else f'{hi:.2f}'})</span>"
        low = str(r.get("low_counts", "")).lower() in ("true", "1", "1.0")
        if low:
            raw += " <span class='st few'>few events</span>"
        dlt = _num(r.get("delta_rate"), 3)
        if dlt != "–" and _f(r, "delta_rate") > 0:
            dlt = "+" + dlt
        rows.append(f"<tr><td>{_esc(r['experiment'])}</td><td>{_esc(_exp_tag(r['experiment'], exp_meta, day_labels))}</td>"
                    f"<td>{_esc(_cond_name(r['baseline'], day_labels))}</td>"
                    f"<td>{_esc(_cond_name(r['condition'], day_labels))}</td>"
                    f"<td>{int(b.get('n_events', 0)) if b else '–'} in {_num(r.get('minutes_baseline'), 1)} min = <b>{_num(r.get('rate_baseline'), 3)}</b></td>"
                    f"<td>{int(c.get('n_events', 0)) if c else '–'} in {_num(r.get('minutes_condition'), 1)} min = <b>{_num(r.get('rate_condition'), 3)}</b></td>"
                    f"<td><b>{dlt}</b></td><td>{raw}</td><td>{_num(r.get('modulation_index'), 2)}</td>"
                    f"<td>{_ratio(r.get('ratio_still'))}</td><td>{_ratio(r.get('ratio_running'))}</td>"
                    f"<td>{_ratio(r.get('ratio_pupil_matched'))}</td><td>{sm}</td><td>{adj}</td><td>{_ratio(r.get('ratio_common_region'))}</td></tr>")
        verdicts.append(f"<li><b>{_esc(r['experiment'])}</b> ({_esc(_exp_tag(r['experiment'], exp_meta, day_labels))}) — "
                        f"{_esc(_cond_name(r['condition'], day_labels))} vs "
                        f"{_esc(_cond_name(r['baseline'], day_labels))}: {_esc(_verdict(r))}</li>")
    h.append("<table class='cmp'><thead>" + head + "</thead><tbody>" + "".join(rows) + "</tbody></table>")
    h.append("<p class='small'>Rates are events per minute of usable recording. <b>Difference</b> = condition − baseline in events/min — the one "
             "measure on the same scale for every experiment. <b>Raw ratio</b> = condition / baseline with its exact 95 % interval (two Poisson "
             "rates): a ratio built on a handful of events has an interval spanning orders of magnitude and is flagged <i>few events</i> "
             "(fewer than 10 in total or fewer than 3 in a condition); it says something about the direction, nothing about the size. "
             "<b>Modulation index</b> = (condition − baseline) / (condition + baseline), −1 … +1, bounded so that a tiny baseline cannot "
             "inflate it. The other ratios are the same change with the behaviour controlled; '–' = not computable. Common region = the ratio "
             "inside the part of the field present in every run after alignment.</p>")
    h.append("<h3>In words</h3><ul class='verdicts'>" + "".join(verdicts) + "</ul>")
    h.append("<h3>On one scale</h3>")
    h.append("<p>The ratios above live on different scales: a day with 1 baseline event and 10 after reads ×10, a day with 51 and 144 reads "
             "×2.8, and a day with 0 and 8 has no ratio at all — yet in events per minute the first two changes are +0.3 and +3.0 and the third "
             "+0.26. The figure below puts every condition of every experiment on one rate axis with its exact Poisson interval, and the change "
             "underneath as a log ratio with a pseudo-count so that experiments without baseline events keep a place; hollow markers are the "
             "ones with too few events to size the change. This is the panel to compare days or FOVs against each other.</p>")
    h.append(_figure(out_uris.get("rates"), "Every condition on one rate axis (top; dots = conditions with exact 95 % Poisson intervals, small dots "
                     "= single runs, numbers = events) and the change vs baseline on one axis (bottom; log2 ratio with pseudo-count 0.5, bars = "
                     "exact 95 % interval of the rate ratio, dashed = open-ended because a condition has no events, hollow = too few events).",
                     missing="no rates figure (run the analysis with v31 or later)"))
    h.append(_figure(out_uris.get("effects"), "Change in event rate per experiment (log scale, 1 = no change): raw with its exact 95 % interval "
                     "(arrows = beyond the axis; hollow = pseudo-count ratio where the baseline had no events) and with each behaviour control; "
                     "the adjusted estimate carries its 95 % CI; rows with too few events are labelled.", missing="no effects figure (run the analysis)"))
    return "\n".join(h)


def section_timecourse(A: Path, out_uris: dict) -> str:
    """The later panel: time course on the day's clock and the matched-speed bands."""
    h = ["<h2 id='timecourse'>7. Time course on the day's clock and matched running speed</h2>"]
    tc = A / "timecourse_summary.csv"
    if not out_uris.get("timecourse") and not out_uris.get("speed"):
        h.append("<p class='small'>No time-course or speed-band figure (spontaneous paradigm, or the analysis has not run).</p>")
        return "\n".join(h)
    if out_uris.get("timecourse"):
        h.append("<h3>Time course on the day's clock (optostim days)</h3>")
        h.append(_figure(out_uris["timecourse"], "Every run of a day on the session clock: event rate per 2-min bin and per run, pre (blue) and post "
                         "(orange) with the stimulation gap; below, the post rate relative to the same field's pre rate against minutes since the first "
                         "post recording, with the fitted slope (does the change decay, build up or stay flat)."))
        ts = _read_csv(tc)
        if ts is not None and len(ts):
            keep = [c for c in ts.columns if c in ("day", "label", "gap_min", "n_post_events", "slope_per_10min", "slope_ci_lo", "slope_ci_hi",
                                                     "slope_adjusted_per_10min", "step_ratio", "within_drift_ratio", "second_vs_first_post")]
            if keep:
                h.append("<table><thead><tr>" + "".join(f"<th>{_esc(c)}</th>" for c in keep) + "</tr></thead><tbody>"
                         + "".join("<tr>" + "".join(f"<td>{_esc(_num(r[c], 2) if isinstance(r[c], float) else r[c])}</td>" for c in keep) + "</tr>"
                                   for r in ts.to_dict('records')) + "</tbody></table>")
    if out_uris.get("speed"):
        h.append("<h3>At matched running speed</h3>")
        h.append(_figure(out_uris["speed"], "Event rate per speed band (still, slow, medium, fast running) and condition, with the minutes spent "
                         "in each band; the matched-speed ratio averages the band rates with fixed weights over the bands both conditions visited."))
    return "\n".join(h)


def section_experiment(exp: str, info: dict, A: Path, rs, day_labels: dict, summary_txt: str) -> str:
    h = [f"<h3 id='exp-{_esc(exp)}'>{_esc(exp)} <span class='cond'>{_esc(info.get('paradigm', ''))}, FOV {info.get('fov', '')}</span></h3>"]
    h.append(f"<p class='small'>Figures of this experiment above: <a href='#gen-{_esc(exp)}'>general</a> · <a href='#heat-{_esc(exp)}'>heatmap</a> · "
             f"<a href='#loco-{_esc(exp)}'>locomotion trials</a>.</p>")
    runs = info.get("runs", [])
    days = list(dict.fromkeys(day_of(r["run"]) for r in runs if day_of(r["run"])))
    labs = [f"{d}: {day_labels[d]}" for d in days if day_labels.get(d)]
    if labs:
        h.append("<p class='small'>" + "; ".join(_esc(x) for x in labs) + "</p>")
    if rs is not None and len(rs):
        sub = rs[rs["experiment"] == exp]
        if len(sub):
            rows = []
            for r in sub.to_dict("records"):
                still = r.get("minutes_still", np.nan); run_ = r.get("minutes_running", np.nan)
                sfrac = still / (still + run_) if np.isfinite(still) and np.isfinite(run_) and (still + run_) > 0 else np.nan
                cur = "yes" if str(r.get("curated", "")).lower() in ("true", "1") else "no"
                rows.append(f"<tr><td>{_esc(_short_run(r['run']))}</td><td>{_esc(_cond_name(r['condition'], day_labels))}</td>"
                            f"<td>{int(r['n_events'])} / {int(r['n_events_detected'])}</td><td>{_num(r['minutes_used'], 1)}</td>"
                            f"<td><b>{_num(r['events_per_min'], 3)}</b></td><td>{int(r.get('excluded_frames', 0))} ({int(r.get('flagged_frames', 0))})</td>"
                            f"<td>{_pct(sfrac)}</td><td>{_num(r.get('velocity_median'), 0)}</td><td>{_num(r.get('pupil_median_px'), 1)}</td>"
                            f"<td>{_num(r.get('dff_max_median'), 2)}</td><td>{_num(r.get('area_um2_median'), 0)}</td><td>{_num(r.get('duration_s_median'), 1)}</td><td>{cur}</td></tr>")
            h.append("<table class='runs'><thead><tr><th>Run</th><th>Condition</th><th>Events kept / detected</th><th>Minutes used</th><th>Events / min</th>"
                     "<th>Excluded frames (flagged)</th><th>Time still</th><th>Median speed</th><th>Pupil px</th><th>Median peak dF/F</th>"
                     "<th>Median area µm²</th><th>Median duration s</th><th>Curated</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")
    for kind, cap in (("overview", "event rate per run (y axis shared by every experiment of this animal, bars = 95 % Poisson intervals) "
                                   "and the event-property distributions per condition"),
                      ("cfu_map", "CFU maps per run, colour = matched identity across runs"),
                      ("maps", "where the activity is: aligned runs, fibre atlas, tile maps, fibre segments"),
                      ("pupil", "pupil size and afferent activity")):
        p = A / f"{exp}_{kind}.png"
        if p.exists():
            h.append(_figure(_file_uri(p), f"{_esc(exp)} — {cap}.", cls="half" if kind in ("cfu_map", "pupil") else ""))
    if summary_txt:
        h.append(f"<details><summary>Full text of summary.txt for this experiment</summary><pre>{_esc(summary_txt)}</pre></details>")
    return "\n".join(h)


def _summary_chunks(A: Path) -> dict:
    """summary.txt split by experiment header lines."""
    p = A / "summary.txt"
    if not p.exists():
        return {}
    chunks: dict = {}
    cur = None
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if line and not line.startswith(" ") and "[" in line and "FOV" in line and "baseline" in line:
            cur = line.split("  ")[0].strip()
            chunks[cur] = [line]
        elif cur:
            chunks[cur].append(line)
    return {k: "\n".join(v).rstrip() for k, v in chunks.items()}


# ----------------------------------------------------------------------------- v32 sections: general representation, heatmaps, locomotion, pupil, GLM

def _mean_sem(v) -> str:
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if v.size == 0:
        return "–"
    if v.size == 1:
        return f"{v[0]:.2f}"
    return f"{v.mean():.2f} ± {v.std(ddof=1) / np.sqrt(v.size):.2f} <span class='small'>(n={v.size})</span>"


def _table(head: list, rows: list, cls: str = "") -> str:
    return (f"<table class='{cls}'><thead><tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead><tbody>"
            + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</tbody></table>")


def _exp_figs(exp: str, A: Path, kinds: list, anchor: str, day_labels: dict, exp_meta: dict) -> list:
    """Figures <exp>_<kind>.png of one experiment, with a small heading."""
    h = [f"<h3 id='{anchor}-{_esc(exp)}'>{_esc(exp)} <span class='cond'>{_esc(_exp_tag(exp, exp_meta, day_labels))}</span></h3>"]
    found = False
    for kind, cap, cls in kinds:
        p = A / f"{exp}_{kind}.png"
        if p.exists():
            found = True
            h.append(_figure(_file_uri(p), f"{_esc(exp)} — {cap}", cls=cls))
    if not found:
        h.append("<p class='small missing'>no figure for this experiment (run the analysis with v32 or later)</p>")
    return h


def section_qc_summary(mouse_dir: Path, A: Path, rs, qc_path: Path | None) -> str:
    """One paragraph on the state of the detection + a link to the QC report."""
    h = ["<h2 id='qc'>Quality control</h2>"]
    if rs is None or not len(rs):
        h.append("<p class='small'>No run_summary.csv yet.</p>")
    else:
        n = len(rs)
        n_cur = int(sum(1 for v in rs["curated"] if str(v).lower() in ("true", "1"))) if "curated" in rs.columns else 0
        exc = rs["excluded_frames"].to_numpy(float) if "excluded_frames" in rs.columns else np.zeros(n)
        tot = rs["n_events_detected"].sum() if "n_events_detected" in rs.columns else 0
        kept = rs["n_events"].sum() if "n_events" in rs.columns else 0
        zero = int((rs["n_events"] == 0).sum()) if "n_events" in rs.columns else 0
        drop_cur = int(rs["n_events_dropped_by_curation"].fillna(0).sum()) if "n_events_dropped_by_curation" in rs.columns else 0
        h.append(f"<p>{n} recordings; {int(kept)} of {int(tot)} detected events kept ({drop_cur} dropped by the curation, the rest in excluded frames); "
                 f"{n_cur} of {n} recordings curated; {int(np.nansum(exc))} excluded frames in total"
                 + (f" (max {int(np.nanmax(exc))} in one recording)" if np.isfinite(exc).any() and np.nanmax(exc) > 0 else "")
                 + (f"; {zero} recording(s) without any kept event" if zero else "") + ".</p>")
    mo = _read_csv(A / "motion_summary.csv")
    if mo is not None and len(mo):
        v = mo["verdict"].astype(str)
        n_w, n_c = int((v == "warn").sum()), int((v == "check").sum())
        med = mo["shift_median_px"].to_numpy(float); jp = mo["jumps_ge2px_per_min"].to_numpy(float)
        near = mo["onsets_near_motion"].fillna(0).sum(); nev = mo["n_event_onsets"].fillna(0).sum(); expd = mo["onsets_near_motion_expected"].fillna(0).sum()
        defm = mo["deformation_median_px"].to_numpy(float) if "deformation_median_px" in mo.columns else np.array([np.nan])
        h.append(f"<p><b>Motion and registration</b> ({len(mo)} recordings): AQuA2 detects on its own registered movie (one whole-pixel translation per "
                 f"frame against the mean of the first 10 frames, border cropped); the raw frames sat {np.nanmedian(med):.1f} px (median of the "
                 f"recordings' medians, up to {np.nanmax(mo['shift_max_px'].to_numpy(float)):.0f} px) from that reference, with "
                 f"{np.nanmedian(jp):.1f} jumps of ≥ 2 px per minute"
                 + (f"; what the registration left: a deformation of {np.nanmedian(defm):.1f} px (median) that one translation cannot remove" if np.isfinite(defm).any() else "")
                 + f". Event onsets on moving frames: {int(near)} of {int(nev)} ({expd:.1f} expected by chance)"
                 + (f" — <span class='warn'>{n_w} recording(s) where the events follow the movement or the field moved a lot</span>" if n_w else "")
                 + (f", {n_c} to check" if n_c else "")
                 + (", nothing to flag" if not (n_w or n_c) else "") + ". Per recording: the QC report.</p>")
    if qc_path is not None and qc_path.exists():
        h.append(f"<p>The evidence — every recording's motion check, detection map, raster, clips and a montage of every event, with a QC table and warnings — is the "
                 f"separate <a href='{_esc(qc_path.name)}'>QC report</a> ({qc_path.stat().st_size / 1e6:.0f} MB, same folder).</p>")
    else:
        h.append("<p class='small'>The per-recording evidence (detection maps, rasters, clips, montages) is a separate QC report: build it with "
                 "<i>Build QC report</i> on the Report tab (or <code>aqua_report.py --qc</code>).</p>")
    return "\n".join(h)


def section_general(mouse: str, exps: list, A: Path, day_labels: dict, exp_meta: dict, spont: bool) -> str:
    """1. The general representation: session overviews, rates by state, ECDFs, acute / chronic summary."""
    h = ["<h2 id='general'>1. General representation of the data</h2>",
         "<p>The measurements of the lab's standard presentation, per experiment: the whole session on the day's clock (events per minute, "
         "running with its bouts and locomotion trials, pupil), the activity split by behavioural state (running / resting) with the event "
         "footprints on the field, and the distributions of event amplitude, duration and area. Within one animal the dots are recordings, "
         "trials or CFUs; the unit for a study-level statistic is the animal.</p>"]
    sf = _read_csv(A / "session_figures.csv")
    if sf is not None and len(sf):
        h.append("<h3>The sessions on the clock</h3>")
        for r in sf.to_dict("records"):
            p = A / str(r["png"])
            if not p.exists():
                continue
            est = str(r.get("estimated_clock", "")).lower() in ("true", "1")
            cap = (f"Day {r['key']}" if r.get("kind") == "day" else f"{r['key']}") + ": events per 1-min bin (baseline grey, condition coloured; "
            cap += "hatched = the stimulation between the recordings; " if r.get("kind") == "day" else "one block per day; "
            cap += "orange = excluded frames), running speed with the running state shaded, bouts (bars) and locomotion trials (▼), and the pupil."
            cap += " Clock estimated from the run order (no recording times found)." if est else ""
            h.append(_figure(_file_uri(p), cap))
    else:
        h.append("<p class='small missing'>no session figures (run the analysis with v32 or later)</p>")
    # rates by state: a compact table from state_rates.csv
    st = _read_csv(A / "state_rates.csv")
    h.append("<h3>Activity by behavioural state</h3>")
    if st is not None and len(st):
        rows = []
        for exp in [e["name"] for e in exps]:
            sub = st[st["experiment"] == exp]
            for c in list(dict.fromkeys(sub["condition"])):
                s = sub[sub["condition"] == c]
                rows.append([_esc(exp), _esc(_cond_name(c, day_labels)), f"{len(s)}", _mean_sem(s["pct_running"]),
                             f"{int(s['n_events_running'].sum())} / {int(s['n_events_resting'].sum())}",
                             _mean_sem(s["events_per_min_running"]), _mean_sem(s["events_per_min_resting"]),
                             _mean_sem(s["dff_max_median_running"]), _mean_sem(s["dff_max_median_resting"])])
        h.append(_table(["Experiment", "Condition", "Runs", "% running", "Events running / resting", "Events/min running", "Events/min resting",
                         "Max dF/F running", "Max dF/F resting"], rows, cls="cmp"))
        h.append("<p class='small'>Mean ± sem over the runs of the condition (running = still fraction &lt; 0.5 in the frame; a rate needs ≥ 15 s in the "
                 "state). The afferents of these animals fire mostly while the mouse runs, and the animals run most of the time — compare the "
                 "resting rates with that in mind.</p>")
    for e in exps:
        h += _exp_figs(e["name"], A, [("state", "event rate and amplitude by state (dot = run), running %, and the event footprints (boxes) that started "
                                                "while running vs resting, per condition.", ""),
                                      ("ecdf", "cumulative distributions of event amplitude (÷ the FOV's 99th percentile), duration and area, "
                                               "by condition (colour) and state at onset (solid running, dashed resting).", "")],
                       "gen", day_labels, exp_meta)
    # acute / chronic
    ac = _read_csv(A / "acute_chronic.csv")
    h.append("<h3 id='acute'>Acute and chronic change of every measure</h3>")
    if ac is not None and len(ac):
        h.append(_figure(_file_uri(A / f"{mouse}_acute_chronic.png"),
                         "Per FOV and day: acute change = post / pre of the day (solid; RearChr) and chronic change = the pre value of the day / the pre "
                         "value of the first day (dashed) for event rates, amplitude, locomotion-trial responses, behaviour, pupil and coupling "
                         "(1 = no change, log axis). Ratios built on a handful of events are as fragile here as in section 6 — read them with the "
                         "counts of that section."))
        keys = [k for k in ac["measure"].unique()]
        cols = [c for c in ac[["experiment", "condition"]].drop_duplicates().itertuples(index=False)]
        rows = []
        for k in keys:
            sub = ac[ac["measure"] == k]
            lab = sub["label"].iloc[0]
            cells = [_esc(lab)]
            for exp, c in cols:
                r = sub[(sub["experiment"] == exp) & (sub["condition"] == c)]
                if not len(r):
                    cells.append("–"); continue
                r = r.iloc[0]
                v = _num(r["value"], 2)
                a = _ratio(r["acute_ratio"]); ch = _ratio(r["chronic_ratio"])
                tag = (f" <span class='small'>acute {a}</span>" if a != "–" else "") + (f" <span class='small'>chronic {ch}</span>" if ch != "–" else "")
                cells.append(v + tag)
            rows.append(cells)
        h.append(_table(["Measure"] + [f"{_esc(exp)}<br>{_esc(_cond_name(c, day_labels))}" for exp, c in cols], rows, cls="cmp"))
        h.append("<p class='small'>Values per experiment × condition (acute_chronic.csv); acute = condition / baseline within the experiment "
                 "(RearChr), chronic = this day's baseline / the first day's baseline of the same FOV (both paradigms; for Treg the days are the "
                 "conditions, so chronic = day / day 1).</p>")
    else:
        h.append("<p class='small missing'>no acute_chronic.csv (run the analysis with v32 or later)</p>")
    return "\n".join(h)


def section_heatmaps(exps: list, A: Path, day_labels: dict, exp_meta: dict) -> str:
    h = ["<h2 id='heatmaps'>2. Heatmap activity comparison</h2>",
         "<p>CFU × time per experiment, the runs concatenated in clock order (white dashed lines = run boundaries, the thick one = the "
         "stimulation or the day change). Rows are the experiment's CFU identities — the same unit followed across runs by map overlap — and "
         "each row is read from every run's movie under that identity's footprint (dF/F, F0 = 20th percentile), then scaled 0–1 per CFU "
         "over the whole experiment; rows are grouped by k-means (k = 3) on their per-condition means and ordered by their change. Underneath: "
         "the event onsets on the same rows and the running speed. The pies count the CFUs whose 20-s bin means are higher (up) or lower (down) in "
         "the condition than in the baseline (Mann-Whitney p &lt; 0.05); the scatter is the per-CFU event rate.</p>"]
    hs = _read_csv(A / "heatmap_summary.csv")
    if hs is not None and len(hs) and str(hs["source"].iloc[0]) != "movie":
        h.append("<p class='small warn'>The binned movies were not readable when the analysis ran, so the rows show AQuA2's own CFU dF/F where the "
                 "CFU was detected and grey elsewhere; the up / down test needs the movies (run the analysis on the computer that holds them).</p>")
    for e in exps:
        h += _exp_figs(e["name"], A, [("heatmap", "CFU identities × time (0–1 per CFU), event onsets and running speed.", ""),
                                      ("updown", "CFUs up / down / unchanged per condition vs baseline (Total = N identities) and their event rates.", "half"),
                                      ("maps", "where the activity is (25-µm tiles): aligned runs, fibre atlas, event density per condition and "
                                               "the log2 ratio with the behaviour factor taken out, fibre segments.", "")],
                       "heat", day_labels, exp_meta)
    return "\n".join(h)


def _trigger_table(A: Path, day_labels: dict, exp_meta: dict) -> list:
    """The three-trigger comparison (trigger_comparison.csv) as one table: per experiment and condition, running onsets /
    still pupil dilations / rest with the evoked rate, its ratio interval, the responsive fraction, the amplitude and the verdict."""
    tc = _read_csv(A / "trigger_comparison.csv")
    if tc is None or not len(tc):
        return ["<p class='small missing'>no trigger_comparison.csv (run the analysis with v34 or later)</p>"]
    rows = []
    for r in tc.to_dict("records"):
        trig = str(r.get("trigger", ""))
        testable = str(r.get("testable", "")).lower() in ("true", "1")
        if trig == "rest":
            n = f"{_num(r.get('minutes_base'), 1)} / {_num(r.get('minutes_cond'), 1)} min"
        else:
            n = f"{int(_f(r, 'n_base', 0))} / {int(_f(r, 'n_cond', 0))} trials"
        rate = f"{_num(r.get('rate_base'), 3)} → {_num(r.get('rate_cond'), 3)}"
        lo, hi = _f(r, "rate_ratio_ci_lo"), _f(r, "rate_ratio_ci_hi")
        if testable and np.isfinite(lo):
            rate += f" <span class='small'>(×{_num(r.get('rate_ratio'), 2)}, CI {_num(lo, 2)}–{'∞' if hi == np.inf else _num(hi, 2)})</span>"
        frac = "–"
        if np.isfinite(_f(r, "frac_baseline")):
            frac = (f"{_num(100 * _f(r, 'frac_baseline'), 0)} → {_num(100 * _f(r, 'frac_condition'), 0)} %"
                    + (f" <span class='small'>(p = {_num(r.get('frac_p_fisher'), 2)})</span>" if np.isfinite(_f(r, "frac_p_fisher")) else ""))
        amp = "–"
        if np.isfinite(_f(r, "amp_delta")):
            amp = (f"{_num(r.get('amp_baseline'), 3)} → {_num(r.get('amp_condition'), 3)}"
                   + (f" <span class='small'>(p = {_num(r.get('amp_p_mwu'), 2)})</span>" if np.isfinite(_f(r, "amp_p_mwu")) else ""))
        verdict = _esc(str(r.get("verdict", "")))
        cls = "ok" if testable else "nt"
        rows.append([f"<span class='cond'>{_esc(_exp_tag(r['experiment'], exp_meta, day_labels))}</span> {_esc(r['experiment'])}",
                     f"{_esc(_cond_name(str(r.get('condition', '')), day_labels))} vs {_esc(str(r.get('baseline', '')))}",
                     _esc(str(r.get("label", trig))), n, rate, frac, amp, f"<span class='verdict {cls}'>{verdict}</span>"])
    n_t = int(sum(1 for r in tc.to_dict("records") if str(r.get("testable", "")).lower() in ("true", "1")))
    return [_table(["Experiment", "Comparison", "Trigger", "N", "Events/min (base → cond)", "Responsive trials", "Response dF/F (peak − base)", "Verdict"],
                   rows, cls="cmp"),
            f"<p class='small'>{n_t} of {len(tc)} comparisons testable (≥ 3 trials in both conditions; rest: ≥ 1 min in both). Events/min = kept "
            "events starting in the 30-s response windows per minute of window (rest: per minute of rest), with the exact Poisson interval of "
            "the ratio; responsive = z rule (peak − baseline ≥ 2 × the still-period noise), Fisher's exact p; response = population dF/F peak "
            "− baseline per trial, Mann-Whitney p. A comparison that is not testable says why: when the animal never ran, the running trigger "
            "has no trials and the pupil and rest rows carry the comparison; when it ran all the time, the rest row is empty and the running "
            "row carries it.</p>"]


def section_loco(exps: list, A: Path, day_labels: dict, exp_meta: dict, spont: bool) -> str:
    h = ["<h2 id='loco'>3. Triggered activity: running onsets, pupil dilations, rest</h2>",
         "<p>The afferent activity in these recordings is spontaneous and running is its main trigger, so the effect of the stimulation "
         + ("across days" if spont else "(acute, post vs pre)") + " is read on the triggered responses, as in the lab's astrocyte "
         "pipeline — on three triggers side by side, so that a session in which the animal does not run is not left without a comparison: "
         "<b>running onsets</b> (the locomotion trials below), <b>pupil dilations during stillness</b> (the same windows around a rise of the "
         "pupil ≥ 20 % of the run's range within 10 s after a stable, still 10-s baseline, with no running in the response window either — "
         "arousal without movement; whether a dilation triggers afferent activity is an open question this answers per session, and the pupil "
         "remains a control elsewhere: section 4's matched-pupil rates and the GLM), and <b>resting activity</b> (events in still frames ≥ 20 s "
         "from any running — the whole run when the animal never runs). Each trigger carries the same scale-aware statistics as section 1 "
         "(rate difference, exact interval of the ratio, modulation index) plus the responsive-trial fraction and the response amplitude.</p>"]
    h += _trigger_table(A, day_labels, exp_meta)
    h.append("<h3>The locomotion trials in detail</h3>"
             "<p>A <b>trial</b> is one minute around a running onset: a pre-bout window of 10 s (5 frames) with no running, a running window "
             "of 30 s (15 frames) and a post-bout window of 20 s (10 frames) with no running — the <b>strict</b> set; because the animals run "
             "most of the time, a <b>relaxed</b> set (≥ 10 s still before the onset, bout ≥ 4 s, nothing required after) is scored too, and every "
             "number carries its N. Responses per trial: the population dF/F (mean of the kept CFUs), the event-count trace and each CFU; baseline "
             "= pre-bout mean, peak = maximum in the running window.</p>"
             "<p><b>Responsive trial</b>, three rules side by side: <i>yu</i> = peak ≥ 1.5 × baseline on F/F0 (the slide's definition — with afferent "
             "dF/F baselines near 0 it needs a peak of +0.5 dF/F and is rarely met), <i>z</i> = peak − baseline ≥ 2 × the run's still-period noise "
             "(the default in the figures and in the trigger table), <i>ev</i> = at least one kept event starts in the running window. Please tell "
             "the pipeline which one the lab uses.</p>")
    lc = _read_csv(A / "loco_conditions.csv")
    if lc is not None and len(lc):
        rows = []
        for r in lc.to_dict("records"):
            rows.append([_esc(r["experiment"]), _esc(_cond_name(r["condition"], day_labels)), f"{int(r.get('n_runs', 0))}",
                         f"{int(r.get('n_trials_strict', 0))} / {int(r.get('n_trials_relaxed', 0))}",
                         f"{_num(r.get('pct_responsive_z_relaxed'), 0)} / {_num(r.get('pct_responsive_yu_relaxed'), 0)} / {_num(r.get('pct_responsive_ev_relaxed'), 0)}",
                         _num(r.get("response_mean_relaxed"), 3) + (f" ± {_num(r.get('response_sem_relaxed'), 3)}" if np.isfinite(_f(r, 'response_sem_relaxed')) else ""),
                         _num(r.get("peak_responsive_mean_relaxed"), 3),
                         _num(r.get("events_per_trial_relaxed"), 2),
                         (f"{_num(r.get('pct_cfu_loco_responsive'), 0)} % of {int(_f(r, 'n_cfu_tested'))}" if np.isfinite(_f(r, "n_cfu_tested")) and _f(r, "n_cfu_tested") > 0 else "–"),
                         (f"{_num(r.get('rest_events_per_min'), 3)} <span class='small'>({int(_f(r, 'rest_n_events', 0))} in {_num(r.get('rest_minutes'), 1)} min)</span>"
                          if np.isfinite(_f(r, "rest_minutes")) else "–"),
                         _num(r.get("pct_running_mean"), 0), _num(r.get("n_bouts_mean"), 1),
                         _num(r.get("dev_l2f_population_mean"), 2)])
        h.append(_table(["Experiment", "Condition", "Runs", "Trials strict / relaxed", "Responsive % z / yu / ev (relaxed)", "Response dF/F (peak − base)",
                         "Peak of responsive trials", "Events per trial", "Loco-responsive CFUs", "Rest-only events/min", "% running", "Bouts / run",
                         "GLM dev. expl. (speed)"], rows, cls="cmp"))
        h.append("<p class='small'>Trials of all runs of the condition pooled. Loco-responsive CFUs: dF/F in the running windows above the pre-bout "
                 "windows beyond 95 % of 500 circular shifts of the onsets. Rest-only = events in still frames ≥ 20 s away from any running — the "
                 "activity that running did not trigger. GLM = held-out deviance explained of a ridge regression of the population dF/F on the "
                 "running speed (lags −5…10 s) over the trial frames.</p>")
    else:
        h.append("<p class='small missing'>no loco_conditions.csv (run the analysis with v32 or later)</p>")
    for e in exps:
        h += _exp_figs(e["name"], A, [("triggers", "the three triggers: onset-aligned population dF/F for running onsets and still pupil dilations, "
                                                   "resting and evoked event rates with exact intervals, and the verdicts.", ""),
                                      ("loco", "the locomotion trials: onset-aligned population dF/F and event-count traces, responsive-trial pies, trials per run, "
                                               "per-trial peaks and event counts, locomotion-responsive CFUs, rest-only activity, the behaviour "
                                               "check (running, bouts, pupil dilation) and the GLM fit.", ""),
                                      ("loco_cfu", "peri-onset activity of every CFU (trial-averaged, 0–1 per CFU), per condition.", "half")],
                       "loco", day_labels, exp_meta)
    h.append("<p class='small'><b>Caveats.</b> Afferent events are sparse: most runs give a few trials and few of them responsive — read the N; the "
             "event-count trace is often more robust than the dF/F trace. Dots are trials, runs or CFUs within one animal, as in the slides; the "
             "study-level unit is the animal. Trials whose windows contain excluded frames are dropped.</p>")
    return "\n".join(h)


def section_pupil(exps: list, A: Path, day_labels: dict, exp_meta: dict) -> str:
    h = ["<h2 id='pupil'>4. Arousal — pupil measures</h2>",
         "<p>Pupil diameter normalised 0–1 per FOV and day (1st–99th percentile over that day's runs of the FOV, not per run, so that a change "
         "between runs is kept). Pupil variation = SD of the normalised trace; resting pupil = its mean over resting frames (valid when ≥ 10 % "
         "of the frames are resting); locomotion-related dilation = peak in the running window − pre-bout mean per locomotion trial; spontaneous "
         "dilations = rises of ≥ 10 % of the range during rest ≥ 20 s away from running.</p>"]
    pm = _read_csv(A / "pupil_measures.csv")
    if pm is not None and len(pm) and "pupil_variation" in pm.columns:
        rows = []
        for exp in [e["name"] for e in exps]:
            sub = pm[pm["experiment"] == exp]
            for c in list(dict.fromkeys(sub["condition"])):
                s = sub[sub["condition"] == c]
                rows.append([_esc(exp), _esc(_cond_name(c, day_labels)), f"{int(s['has_pupil'].astype(bool).sum())} / {len(s)}",
                             _mean_sem(s["pupil_variation"]), _mean_sem(s["resting_pupil"]) + (f" <span class='small'>({int(s['resting_pupil_valid'].fillna(0).sum())} valid)</span>" if "resting_pupil_valid" in s.columns else ""),
                             _mean_sem(s["running_pupil"]) if "running_pupil" in s.columns else "–",
                             _mean_sem(s["loco_dilation_norm_mean"]) if "loco_dilation_norm_mean" in s.columns else "–",
                             _mean_sem(s["spont_dilations_per_min"]) if "spont_dilations_per_min" in s.columns else "–",
                             _mean_sem(s["spont_dilation_amp_norm"]) if "spont_dilation_amp_norm" in s.columns else "–"])
        h.append(_table(["Experiment", "Condition", "Runs with pupil", "Pupil variation (SD)", "Resting pupil (0–1)", "Running pupil (0–1)",
                         "Loco-related dilation (0–1)", "Spontaneous dilations / min", "Their amplitude (0–1)"], rows, cls="cmp"))
    else:
        h.append("<p class='small'>No pupil data for this animal (no DLC eye csv next to the movies), or the analysis predates v32.</p>")
    for e in exps:
        p = A / f"{e['name']}_pupil.png"
        if p.exists():
            h.append(_figure(_file_uri(p), f"{_esc(e['name'])} — pupil size and afferent activity: coupling, lags and the pre/post comparison at matched pupil size.", cls="half"))
    return "\n".join(h)


def section_glm(exps: list, A: Path, day_labels: dict, exp_meta: dict) -> str:
    h = ["<h2 id='glm'>5. Coupling between behaviour and activity (GLMs)</h2>",
         "<p>Ridge regressions of the dF/F on lagged behaviour, as in the slides: <b>Loco2Fluor</b> — running speed with lags −5…10 s on the "
         "trial frames (train 75 %, 50 random splits, penalty by 5-fold CV); <b>LocoPupil2Fluor</b> — speed and pupil (0–1) with lags "
         "−20…20 s on the whole recording (20 splits, 10-fold CV). Held-out deviance explained ≥ 0.2 = well fitted; leave-one-predictor-out "
         "gives the relative value of speed and pupil. Fitted on the population dF/F and on every CFU with ≥ 3 events (afferent CFUs mostly "
         "have one event, so the population trace is the main readout).</p>"]
    import pandas as pd
    parts, kparts = [], []
    for e in exps:
        g = _read_csv(A / f"loco_glm_{e['name']}.csv")
        if g is not None and len(g):
            parts.append(g)
        k = _read_csv(A / f"loco_glm_kernels_{e['name']}.csv")
        if k is not None and len(k):
            k = k.assign(experiment=e["name"]); kparts.append(k)
    if not parts:
        h.append("<p class='small missing'>no GLM output (run the analysis with v32 or later, without --no-glm)</p>")
        return "\n".join(h)
    g = pd.concat(parts, ignore_index=True)
    rows = []
    for exp in [e["name"] for e in exps]:
        sub = g[(g["experiment"] == exp) & (g["unit"] == "population")]
        for c in list(dict.fromkeys(sub["condition"])):
            s = sub[sub["condition"] == c]
            l2f = s[s["model"] == "loco2fluor"]; lp = s[s["model"] == "locopupil2fluor"]
            cf = g[(g["experiment"] == exp) & (g["condition"] == c) & (g["unit"] != "population") & (g["model"] == "loco2fluor")]
            rows.append([_esc(exp), _esc(_cond_name(c, day_labels)), _mean_sem(l2f["dev_explained"]), _mean_sem(l2f["peak_coef_speed"]),
                         _mean_sem(lp["dev_explained"]), _mean_sem(lp["rel_value_speed"]) if "rel_value_speed" in lp.columns else "–",
                         _mean_sem(lp["rel_value_pupil"]) if "rel_value_pupil" in lp.columns else "–",
                         f"{int(cf['well_fitted'].sum())} / {len(cf)}" if len(cf) else "0 / 0"])
    h.append(_table(["Experiment", "Condition", "Loco2Fluor dev. explained", "Peak speed coefficient", "LocoPupil2Fluor dev. explained",
                     "Relative value: speed", "Relative value: pupil", "CFUs well fitted / fitted"], rows, cls="cmp"))
    h.append("<p class='small'>Mean ± sem over runs (population trace); coefficients are per SD of the predictor; the relative value is "
             "(dev<sub>full</sub> − dev<sub>without</sub>) / dev<sub>full</sub>, clipped to 0–1. A drop of the deviance explained or of the "
             "speed coefficient after the stimulation = weaker coupling between running and afferent activity.</p>")
    if kparts:
        try:
            h.append(_figure(_kernel_figure(pd.concat(kparts, ignore_index=True), [e["name"] for e in exps], day_labels),
                             "Kernels of the population GLMs (coefficient vs lag of the predictor; mean over runs per condition): speed in the "
                             "Loco2Fluor model (left) and speed and pupil in the LocoPupil2Fluor model (right)."))
        except Exception as ex:
            h.append(f"<p class='small warn'>kernel figure failed: {_esc(ex)}</p>")
    return "\n".join(h)


def _kernel_figure(k, exp_names: list, day_labels: dict) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    k = k[k["unit"] == "population"]
    n = len(exp_names)
    fig, axes = plt.subplots(n, 3, figsize=(13, 2.6 * n + 0.6), squeeze=False)
    for i, exp in enumerate(exp_names):
        sub = k[k["experiment"] == exp]
        conds = list(dict.fromkeys(sub["condition"]))
        for j, (model, pred, ttl) in enumerate((("loco2fluor", "speed", "Loco2Fluor: speed"), ("locopupil2fluor", "speed", "LocoPupil2Fluor: speed"),
                                                ("locopupil2fluor", "pupil", "LocoPupil2Fluor: pupil"))):
            ax = axes[i, j]
            for ci, c in enumerate(conds):
                s = sub[(sub["condition"] == c) & (sub["model"] == model) & (sub["predictor"] == pred)]
                if not len(s):
                    continue
                m = s.groupby("lag_s")["coef"].mean()
                ax.plot(m.index.to_numpy(), m.to_numpy(), color=SERIES[ci % len(SERIES)], lw=1.4, label=_cond_name(c, day_labels))
            ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
            ax.set_title(f"{_short_exp(exp)} — {ttl}", fontsize=8); ax.set_xlabel("lag (s; > 0 = predictor before dF/F)", fontsize=7)
            ax.set_ylabel("coefficient", fontsize=7); _style(ax)
            if j == 0:
                ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    return _fig_uri(fig, dpi=100)


def _short_exp(exp: str) -> str:
    import re
    m = re.search(r"(\d{6}_FOV\d+|FOV\d+_days)", exp)
    return m.group(1) if m else exp


def section_settings(man: dict, mouse_dir: Path, exps: list | None = None) -> str:
    """Which AQuA2 settings each experiment was analysed with, the CFU rule the results were made with, and the tuning sweep if there was one."""
    h = ["<h2 id='settings'>9. Detection settings and the tuning sweep</h2>"]
    presets = {p["preset"]: p for p in man.get("presets", [])}
    rows = []
    keys = ["thrARScl", "minSize", "minDur", "maxSize", "smoXY", "medSmo", "sigThr", "sourceSensitivity", "registrateCorrect", "bleachCorrect", "regMaskGap"]
    for e in man.get("experiments", []):
        if not e.get("runs"):
            continue
        pidx = e["runs"][0].get("preset")
        p = presets.get(pidx, {})
        vals = p.get("values", {})
        src = Path(str(e.get("opt_file") or "")).name or "–"
        rows.append(f"<tr><td>{_esc(e['name'])}</td><td>{_esc(src)}</td>" + "".join(f"<td>{_esc(vals.get(k, '–'))}</td>" for k in keys) + "</tr>")
    if rows:
        h.append("<table><thead><tr><th>Experiment</th><th>Settings file</th>" + "".join(f"<th>{k}</th>" for k in keys) + "</tr></thead><tbody>"
                 + "".join(rows) + "</tbody></table>")
        h.append("<p class='small'>One settings file per FOV for all its runs (the batch enforces this): thrARScl = active-voxel threshold, minSize / minDur = "
                 "smallest footprint and duration, smoXY / medSmo = spatial smoothing and median filter, sigThr = significance of a temporal segment, "
                 "sourceSensitivity = spatial splitting, registrateCorrect / bleachCorrect = motion and bleach correction modes, regMaskGap = border excluded. "
                 "A file in aqua_opts\\ was written by the tuning sweep; FOV&lt;k&gt;_Opt.csv by hand in the AQuA2 GUI; default_Opt.csv is the pipeline's default.</p>")
    # the CFU rule, read from the result files themselves (their saved cfuOpts), not from the batch settings
    found: dict[str, int] = {}
    for e in exps or []:
        for r in e.get("runs", []):
            f = aqua_files(Path(r["tif"]))["cfu"] if r.get("tif") else None
            if f is not None and f.exists():
                rule = cfu_rules_of(f) or "unreadable cfuOpts"
                found[rule] = found.get(rule, 0) + 1
    if found:
        listing = "; ".join(f"<b>{_esc(k)}</b> ({v} run{'s' if v != 1 else ''})" for k, v in sorted(found.items(), key=lambda kv: -kv[1]))
        newer = any("post-merge" in k for k in found)
        h.append(f"<p class='small'><b>CFU rule.</b> The units of this report were made with: {listing}. "
                 "AQuA2's published consensus-functional-unit rule: events whose duration-weighted footprints overlap by at least the overlap "
                 "threshold (weighted intersection over union, single-linkage clustering) form one unit, and a unit needs 'min events' events. "
                 "AQuA2's own defaults are 0.5 and 3; the pipeline's 0.3 / 1 were optimised on all recordings of the three RearChr animals "
                 "(bouton events of different extent at one spot overlap only 0.25–0.45, so at 0.5 a repeatedly firing spot split into "
                 "one-event units; 0.3 is where the clustering best reproduces same-location repeats without chaining distinct places)."
                 + (" <b>Some results carry a post-merge entry: they were made by an AQuA2 copy of 28 Aug 2026 or later, whose newer rules "
                    "(basin splitting + post-merge of correlated neighbours) differ from the pipeline's.</b>" if newer else "")
                 + (" <b>The runs do not all share one rule — rebuild the batch and run the CFU step again before comparing them.</b>"
                    if len(found) > 1 else "") + "</p>")
    else:
        cfu = man.get("cfu") or {}
        h.append("<p class='small'><b>CFU rule.</b> No CFU result file was found; the batch asks for AQuA2's published overlap clustering with "
                 f"overlap ≥ {cfu.get('overlap', '0.5 (batch built before v41)')} and min events {cfu.get('min_events', '?')}.</p>")
    chk = mouse_dir / "aqua_batch" / "cfu_overlap_check.txt"
    if chk.exists():
        txt = chk.read_text(encoding="utf-8", errors="replace")
        verdict = next((ln for ln in txt.splitlines() if ln.startswith("Verdict:") or ln.startswith("CHECK means")), "")
        pooled = next((ln for ln in txt.splitlines() if ln.startswith("POOLED:")), "")
        h.append("<p class='small'><b>Overlap-threshold check</b> (tab 4, v42: AQuA2's pairwise event overlaps re-clustered at every threshold and "
                 f"compared with overlap-free same-place / two-place labels, per day × FOV and pooled): {_esc(pooled)} {_esc(verdict)}</p>")
        h.append(_figure(_file_uri(mouse_dir / "aqua_batch" / "cfu_overlap_check.png"),
                         "Left: how many same-place / two-place pairs end up in one unit vs the threshold; middle: their agreement (F1) pooled; "
                         "right: per experiment, dashed = CHECK.", missing=""))
        h.append(f"<details><summary>Full check</summary><pre>{_esc(txt)}</pre></details>")
    tune_dir = mouse_dir / "aqua_tune"
    rec = tune_dir / "tune_recommendation.txt"
    if rec.exists():
        txt = rec.read_text(encoding="utf-8", errors="replace")
        chosen = [ln.strip() for ln in txt.splitlines() if ln.strip().startswith("->")]
        h.append("<h3>Tuning sweep (tab 3, Detection settings)</h3><p>Settings were chosen by measurement: every grid setting was run on the real movie, on a null copy "
                 "(each pixel's time series shifted by its own random offset — every event found there is false) and on a copy with synthetic transients "
                 "of known amplitude, size and duration injected on the fibres. The chosen setting has the highest recall of the injected transients among "
                 "those whose false-alarm share on the fibres stays low.</p>")
        if chosen:
            h.append("<ul class='small'>" + "".join(f"<li>{_esc(c)}</li>" for c in chosen) + "</ul>")
        fig = tune_dir / f"{man.get('mouse', mouse_dir.name)}_tune.png"
        h.append(_figure(_file_uri(fig), "Per tuned movie: recall of the injected transients vs threshold, recall by amplitude and duration, false-alarm share "
                         "on the fibres, and the real movie's event rate on and off the fibres; star = the chosen setting.", missing=""))
        h.append(f"<details><summary>Full tuning report</summary><pre>{_esc(txt)}</pre></details>")
    else:
        h.append("<p class='small'>No tuning sweep for this animal: the settings above were set by hand or are the pipeline's default.</p>")
    return "\n".join(h)


def section_methods(man: dict, spont: bool) -> str:
    return """<h2 id='methods'>10. How the numbers were made</h2>
<p><b>Recording and binning.</b> Scanbox movies are binned 31 frames → 1 frame (1.984 s per frame, 0.792 µm per pixel). The brightness step
flattens slow global level changes and flags frames with artefacts (brightness jumps, state changes); frames decided as excluded in Clean &amp; check
never count — events starting in or next to them are dropped and their time leaves the denominator.</p>
<p><b>Detection (AQuA2).</b> Every recording of a FOV is analysed with the same settings (one file per FOV, section 9), so event rates are comparable
within a FOV. AQuA2 registers the movie, finds active voxels above <i>thrARScl</i> × the noise, grows them into events with a minimum footprint and
duration, and reports per event the onset, duration, footprint and peak dF/F. The CFU step groups events that recur at the same place into units
(CFUs) by AQuA2's published rule — single-linkage clustering of events whose weighted footprints overlap at least the overlap threshold (section 9
gives the values used); the curation window lets you mark units that are not afferents as noise — their events are dropped here — and add missed
regions.</p>
<p><b>Behavioural state.</b> Running = the frame's still fraction below 0.5 (the wheel encoder, still = |speed| &lt; 5 counts/s for ≥ 5 s). Bouts bridge
single-frame gaps and last ≥ 2 frames. Locomotion trials, their responses and the three responsive rules are defined in section 3; rest-only
activity counts the events in still frames ≥ 20 s away from any running. Pupil measures (section 4) use the diameter from the DLC eye tracking,
normalised 0–1 per FOV and day.</p>
<p><b>Heatmaps.</b> CFU identities are matched across the runs of an experiment by overlap of their maps after aligning the mean images
(rotation + shift); each identity's trace is read from every run's binned movie under its footprint (mapped with the alignment and AQuA2's crop
offset), dF/F with the 20th percentile as F0, excluded frames blanked, then scaled 0–1 per CFU; k-means (k = 3) on the per-condition means orders the
rows. Up / down: 20-s bin means of the 0–1 trace, condition vs baseline, Mann-Whitney p &lt; 0.05.</p>
<p><b>GLMs.</b> Gaussian ridge regression of the dF/F on lagged predictors (standardised), penalty chosen by K-fold cross-validation inside each
random train/test split; the reported deviance explained is 1 − SS<sub>res</sub>/SS<sub>tot</sub> on the held-out quarter, averaged over the splits;
leave-one-predictor-out relative value = (dev<sub>full</sub> − dev<sub>without</sub>) / dev<sub>full</sub>.</p>
<p><b>Behaviour-controlled comparison (section 6).</b> For each experiment (""" + ("one FOV across days; the earliest day is the baseline" if spont else
                                                      "one day and FOV; the pre runs are the baseline") + """) the event rate per minute of usable
recording is compared between conditions: raw; on still frames only and on running frames only (locomotion matched); within pupil-size terciles
(arousal matched); at matched running speed (direct standardisation over the speed bands both conditions visited, with the overlap reported);
and by a Poisson regression on 20-s bins with speed, running and pupil as covariates (behaviour-adjusted ratio with a block-bootstrap 95 % interval).
Where the activity is: the runs are aligned, and the ratio is also taken inside the region present in every run, per 25-µm tile and per fibre segment.
Event-level p values treat events as independent and are descriptive; the statistic for a study is one value per animal (or per FOV, stated as such).</p>
<p><b>This report</b> restates aqua_analysis\\ (the csvs and figures written by the analysis); it recomputes nothing. The per-recording evidence
(detection maps, rasters, clips, montages) is the separate QC report.</p>"""


# ----------------------------------------------------------------------------- CSS / JS

CSS = """
:root { --ink:#0b0b0b; --muted:#52514e; --line:#e6e5e1; --accent:#2a78d6; --bg:#fcfcfb; --soft:#f3f3f1; --warn:#b3261e; }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--ink); font: 15px/1.5 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 26px 24px 80px; }
h1 { font-size: 25px; line-height:1.25; margin: 0 0 4px; }
h2 { font-size: 20px; margin: 42px 0 10px; padding-top: 14px; border-top: 1px solid var(--line); }
h3 { font-size: 16px; margin: 26px 0 6px; }
p { margin: 0 0 10px; }
.byline { color: var(--muted); margin-bottom: 16px; }
.lead { background: var(--soft); border-left: 4px solid var(--accent); padding: 10px 14px; margin: 14px 0 20px; }
figure { margin: 14px 0 18px; }
figure img { width: 100%; height: auto; border: 1px solid var(--line); border-radius: 4px; cursor: zoom-in; display:block; background:#fff; }
figure.half img { max-width: 860px; }
figcaption { font-size: 13px; color: var(--muted); margin-top: 5px; }
figure.montage { margin: 8px 0 14px; } figure.montage img { max-width: 1100px; }
figure.montage figcaption { margin-top: 2px; }
.clips { display:flex; flex-wrap: wrap; gap: 12px; margin: 6px 0 4px; }
figure.clip { margin: 0; width: 230px; } figure.clip img { width: 100%; cursor: default; }
figure.clip figcaption { font-size: 12px; }
.hint { color:#9a9aa2; }
table { border-collapse: collapse; width: 100%; font-size: 13px; margin: 8px 0 14px; }
th, td { border-bottom: 1px solid var(--line); padding: 4px 7px; text-align: left; vertical-align: top; }
th { background: var(--soft); font-weight: 600; }
table.cmp td:nth-child(n+7), table.runs td:nth-child(n+3) { white-space: nowrap; } table.cmp { font-size: 12.5px; }
.small { font-size: 13px; color: var(--muted); }
.warn { color: var(--warn); }
.missing { color: var(--warn); font-style: italic; }
nav.toc { font-size: 14px; margin: 8px 0 6px; line-height: 1.9; } nav.toc a { margin-right: 14px; color: var(--accent); text-decoration: none; }
details { margin: 10px 0; } summary { cursor: pointer; font-weight: 600; }
pre { font: 12px/1.4 Consolas, Menlo, monospace; background: var(--soft); padding: 10px; overflow-x: auto; white-space: pre; }
.cond { display:inline-block; padding: 1px 7px; border-radius: 10px; background: var(--soft); border:1px solid var(--line); font-size: 12.5px; font-weight: normal; color: var(--muted); }
.st { display:inline-block; padding: 0 6px; border-radius: 8px; font-size: 12px; background:#e3eefc; color:#1c5cab; }
.st.noise-CFU { background:#fde2e0; color:#b3261e; } .st.excluded-frames, .st.dropped { background:#ececea; color:#52514e; }
.st.few { background:#fff3d6; color:#8a5a00; }
ul.verdicts li { margin-bottom: 6px; }
span.verdict.nt { color: var(--muted); font-style: italic; } span.verdict.ok { color: inherit; }
table.cmp td:nth-child(n+7) { white-space: normal; }
.glance { display:flex; flex-wrap:wrap; gap: 10px; margin: 10px 0 6px; }
.tile { background: var(--soft); border: 1px solid var(--line); border-radius: 6px; padding: 8px 12px; min-width: 150px; }
.tile .n { font-size: 22px; font-weight: 600; } .tile .l { font-size: 12.5px; color: var(--muted); }
#lb { position: fixed; inset: 0; background: rgba(10,10,12,0.92); display: none; z-index: 1000; overflow: hidden; }
#lb.open { display: block; }
#lb img { position: absolute; left: 0; top: 0; transform-origin: 0 0; user-select: none; -webkit-user-drag: none; cursor: grab; max-width:none; }
#lb .bar { position: absolute; top: 10px; right: 14px; color: #eee; font-size: 14px; background: rgba(0,0,0,0.5); padding: 6px 10px; border-radius: 6px; }
#lb .bar button { margin-left: 8px; font-size: 13px; }
@media print { figure img { cursor: default; } .hint, nav.toc { display: none; } details { display:block; } }
"""

JS = """
(function(){
  const lb = document.getElementById('lb'), im = document.getElementById('lbimg'), zl = document.getElementById('zl');
  let s = 1, tx = 0, ty = 0, drag = null;
  function apply(){ im.style.transform = `translate(${tx}px, ${ty}px) scale(${s})`; zl.textContent = Math.round(s*100) + ' %'; }
  function fit(){ const W = lb.clientWidth, H = lb.clientHeight; s = Math.min(W / im.naturalWidth, H / im.naturalHeight, 1);
    tx = (W - im.naturalWidth * s) / 2; ty = (H - im.naturalHeight * s) / 2; apply(); }
  function one(){ const W = lb.clientWidth, H = lb.clientHeight; s = 1; tx = (W - im.naturalWidth) / 2; ty = (H - im.naturalHeight) / 2; apply(); }
  function open(src){ im.src = src; lb.classList.add('open'); if (im.complete) fit(); else im.onload = fit; }
  window.openLightbox = open;
  document.querySelectorAll('a.zoom').forEach(a => a.addEventListener('click', e => { e.preventDefault(); open(a.querySelector('img').src); }));
  function close(){ lb.classList.remove('open'); }
  document.getElementById('lbclose').addEventListener('click', close);
  document.getElementById('lbfit').addEventListener('click', fit);
  document.getElementById('lbone').addEventListener('click', one);
  document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); if (lb.classList.contains('open') && (e.key === '0')) fit(); if (lb.classList.contains('open') && (e.key === '1')) one(); });
  lb.addEventListener('wheel', e => { e.preventDefault(); const k = e.deltaY < 0 ? 1.15 : 1/1.15; const ns = Math.min(Math.max(s * k, 0.1), 12);
    const r = lb.getBoundingClientRect(); const mx = e.clientX - r.left, my = e.clientY - r.top;
    tx = mx - (mx - tx) * (ns / s); ty = my - (my - ty) * (ns / s); s = ns; apply(); }, {passive:false});
  lb.addEventListener('mousedown', e => { if (e.target === im) { drag = {x: e.clientX - tx, y: e.clientY - ty}; im.style.cursor = 'grabbing'; e.preventDefault(); } });
  window.addEventListener('mousemove', e => { if (drag) { tx = e.clientX - drag.x; ty = e.clientY - drag.y; apply(); } });
  window.addEventListener('mouseup', () => { drag = null; im.style.cursor = 'grab'; });
  lb.addEventListener('dblclick', e => { if (e.target === im) { if (s === 1) fit(); else one(); } });
  lb.addEventListener('click', e => { if (e.target === lb) close(); });
})();
"""


# ----------------------------------------------------------------------------- main

def report_context(mouse_dir: Path) -> dict:
    """What both reports need: manifest, experiments, labels, tables."""
    mouse_dir = Path(mouse_dir)
    A = mouse_dir / "aqua_analysis"
    A.mkdir(exist_ok=True)
    mp = mouse_dir / "aqua_batch" / "manifest.json"
    if not mp.exists():
        raise SystemExit(f"{mp} not found — build the batch first (tab 4, aqua_batch.py)")
    man = json.loads(mp.read_text())
    day_labels = read_day_labels(mouse_dir)
    exps = [e for e in man.get("experiments", []) if e.get("runs")]
    paradigm_of = {e["name"]: e.get("paradigm", "") for e in exps}
    exp_meta = {e["name"]: {"fov": e.get("fov"), "paradigm": e.get("paradigm", ""),
                            "days": list(dict.fromkeys(day_of(r["run"]) for r in e["runs"] if day_of(r["run"])))} for e in exps}
    return {"mouse_dir": mouse_dir, "A": A, "man": man, "day_labels": day_labels, "exps": exps, "paradigm_of": paradigm_of,
            "exp_meta": exp_meta, "spont": any(p == "spontaneous" for p in paradigm_of.values()),
            "rs": _read_csv(A / "run_summary.csv"), "es": _read_csv(A / "experiments_summary.csv")}


def page_html(title: str, body: str, extra_css: str = "", extra_js: str = "") -> str:
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
            f"<title>{_esc(title)}</title><style>{CSS}{extra_css}</style></head><body><main>{body}</main>"
            f"<div id='lb'><img id='lbimg' alt=''><div class='bar'><span id='zl'>100 %</span><button id='lbfit'>fit</button><button id='lbone'>100 %</button>"
            f"<button id='lbclose'>close (Esc)</button></div></div><script>{JS}{extra_js}</script></body></html>")


def build_report(mouse_dir: Path, out: Path | None = None, max_events: int = 0, clips_per_run: int = 6, no_clips: bool = False,
                 no_qc: bool = False, log=print) -> Path:
    """The MAIN report: general representation, heatmaps, locomotion-evoked activity, pupil, GLMs, the behaviour-controlled change, the time
    course, per-experiment details, settings, methods. The per-recording evidence is the separate QC report (qc_report.build_qc_report);
    max_events / clips_per_run / no_clips / no_qc are accepted for compatibility and only matter there."""
    import matplotlib
    matplotlib.use("Agg")
    ctx = report_context(mouse_dir)
    mouse_dir, A, man, day_labels, exps = ctx["mouse_dir"], ctx["A"], ctx["man"], ctx["day_labels"], ctx["exps"]
    paradigm_of, exp_meta, spont, rs, es = ctx["paradigm_of"], ctx["exp_meta"], ctx["spont"], ctx["rs"], ctx["es"]
    out = Path(out) if out else A / f"{mouse_dir.name}_report.html"
    t_start = time.time()
    log(f"report for {mouse_dir.name}: {len(exps)} experiment(s), {sum(len(e['runs']) for e in exps)} recording(s)"
        + ("" if es is not None else " — no analysis output yet, most sections will be empty"))
    uris = {"effects": _file_uri(A / f"{mouse_dir.name}_effects.png"), "timecourse": _file_uri(A / f"{mouse_dir.name}_timecourse.png"),
            "speed": _file_uri(A / f"{mouse_dir.name}_speed_tuning.png"), "rates": _file_uri(A / f"{mouse_dir.name}_rates.png")}
    # --- at a glance
    n_runs = sum(len(e["runs"]) for e in exps)
    days = sorted({day_of(r["run"]) for e in exps for r in e["runs"] if day_of(r["run"])})
    n_events = int(rs["n_events"].sum()) if rs is not None and "n_events" in rs.columns else None
    n_det = int(rs["n_events_detected"].sum()) if rs is not None and "n_events_detected" in rs.columns else None
    n_cur = int(sum(1 for v in rs["curated"] if str(v).lower() in ("true", "1"))) if rs is not None and "curated" in rs.columns else 0
    lc = _read_csv(A / "loco_summary.csv")
    n_trials = int(lc["n_trials_relaxed"].fillna(0).sum()) if lc is not None and "n_trials_relaxed" in lc.columns else None
    n_strict = int(lc["n_trials_strict"].fillna(0).sum()) if lc is not None and "n_trials_strict" in lc.columns else None
    tiles = [("recordings", n_runs), ("days", len(days)), ("FOVs", len({e.get('fov') for e in exps})),
             ("events kept / detected", f"{n_events} / {n_det}" if n_events is not None else "–"),
             ("locomotion trials (strict / relaxed)", f"{n_strict} / {n_trials}" if n_trials is not None else "–"),
             ("recordings curated", f"{n_cur} / {n_runs}"),
             ("tuned by measurement", "yes" if (mouse_dir / "aqua_tune" / "tune_recommendation.txt").exists() else "no")]
    glance = "<div class='glance'>" + "".join(f"<div class='tile'><div class='n'>{_esc(v)}</div><div class='l'>{_esc(k)}</div></div>" for k, v in tiles) + "</div>"
    para_txt = ("Treg / spontaneous paradigm: one experiment per FOV, its days are the conditions and the earliest day is the baseline" if spont
                else "Optostim (RearChr) paradigm: one experiment per day and FOV, pre vs post the stimulation")
    labs = "; ".join(f"{d} = {day_labels[d]}" for d in days if day_labels.get(d))
    qc_path = A / f"{mouse_dir.name}_qc.html"
    parts = [f"<h1>{_esc(mouse_dir.name)} — afferent activity report</h1>",
             f"<p class='byline'>{_esc(para_txt)}. Generated {time.strftime('%Y-%m-%d %H:%M')} by afferent2p (aqua_report.py) from {_esc(str(mouse_dir))}.</p>",
             (f"<p class='small'>Days: {_esc(labs)}</p>" if labs else ""), glance,
             "<nav class='toc'><a href='#qc'>QC</a><a href='#general'>1 General representation</a><a href='#heatmaps'>2 Heatmaps</a>"
             "<a href='#loco'>3 Triggers</a><a href='#pupil'>4 Pupil</a><a href='#glm'>5 Coupling GLMs</a>"
             "<a href='#result'>6 Behaviour-controlled change</a><a href='#timecourse'>7 Time course</a><a href='#experiments'>8 Per experiment</a>"
             "<a href='#settings'>9 Settings</a><a href='#methods'>10 Methods</a></nav>"]
    parts.append(section_qc_summary(mouse_dir, A, rs, qc_path))
    parts.append(section_general(mouse_dir.name, exps, A, day_labels, exp_meta, spont))
    parts.append(section_heatmaps(exps, A, day_labels, exp_meta))
    parts.append(section_loco(exps, A, day_labels, exp_meta, spont))
    parts.append(section_pupil(exps, A, day_labels, exp_meta))
    parts.append(section_glm(exps, A, day_labels, exp_meta))
    parts.append(section_comparison(mouse_dir.name, paradigm_of, A, day_labels, uris, exp_meta))
    parts.append(section_timecourse(A, uris))
    parts.append("<h2 id='experiments'>8. Per experiment — runs, properties, where the activity is</h2>")
    chunks = _summary_chunks(A)
    for e in exps:
        parts.append(section_experiment(e["name"], e, A, rs, day_labels, chunks.get(e["name"], "")))
    parts.append(section_settings(man, mouse_dir, exps))
    parts.append(section_methods(man, spont))
    out.write_text(page_html(f"{mouse_dir.name} — afferent activity report", "\n".join(parts)), encoding="utf-8")
    log(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB) in {(time.time() - t_start) / 60:.1f} min")
    return out
