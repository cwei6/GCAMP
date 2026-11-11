#!/usr/bin/env python3
"""reg_demo — a movie that shows what the registration does to one run.

Why: shift traces and sharpness numbers say that the field drifts by a few
pixels over a run, that AQuA2's own registration takes that drift out and that
the legacy MATLAB chain leaves it in — a movie shows it. This tool puts the
movies side by side, zoomed in on the region with the most structure, and
draws a FIXED outline of each movie's reference (the mean of its first 10
frames — AQuA2's reference) on every frame: whatever slides away from the
outline is motion that the route left in; whatever stays on it was corrected
(or never moved).

Panels (those that exist for the run):
  as recorded     the binned tiff AQuA2 receives (<run>_bin31f_1p98s.tif)
  after AQuA2     res.datOrg1 of <run>..._AQuA2.mat — the movie AQuA2 detected
                  the events on (its integer shift per frame against the mean
                  of the first 10 frames, border crop). When the .mat is not
                  there or not readable, the same step is emulated in numpy
                  (afferent2p.motion.aqua2_style_register) and the panel says so.
  legacy chain    reg_compare\\<run>\\D_legacy\\<stem>.tif — the lab's MATLAB
                  registration (legacy_register_run.m), when the run went
                  through reg_compare.py prepare --add-legacy
  Python pre-reg  reg_compare\\<run>\\B_prereg\\<stem>.tif (unpack --register)

Below the panels: the whole field (as recorded) with the zoom region marked,
the shift of every frame of every movie against its own reference (phase
correlation, afferent2p.motion.frame_shifts) and the running speed, with a
cursor at the current frame. Each panel prints its frame's shift.

Output: <mouse>\\reg_compare\\<run>_reg_demo.mp4 (when reg_compare exists,
else next to the tiff), a 6-s GIF excerpt of the most mobile stretch
(_reg_demo_excerpt.gif) and one PNG frame. Needs imageio + imageio-ffmpeg
(pip install imageio imageio-ffmpeg) for the mp4; the GIF and PNG need only
pillow.

Usage:
  python reg_demo.py "D:\\2photon\\Chao\\GCAMP_affrerents\\NavRearChrTIGRE-2608-1" --run 260921_006
  python reg_demo.py <mouse_dir> --run 260921_006 --variants raw,aqua,legacy --fps 20 --scale 2.5
  python reg_demo.py <mouse_dir> --run 260921_006 --region 200,300,160,240   (y0,x0,height,width in the
                                                                             as-recorded frame; default: auto)
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from afferent2p import motion                                   # noqa: E402
from afferent2p.aqua_io import aqua_files, open_registered, pick_binned_tif   # noqa: E402

VARIANTS = ("raw", "aqua", "legacy", "prereg")
LABEL = {"raw": "as recorded (what AQuA2 receives)",
         "aqua": "after AQuA2's own registration (route A)",
         "legacy": "legacy MATLAB chain (route D)",
         "prereg": "Python pre-registration (route B)"}
COLOR = {"raw": (235, 235, 235), "aqua": (86, 160, 236), "legacy": (96, 200, 96), "prereg": (255, 170, 80)}
OUTLINE = (0, 235, 235)
BG = (18, 18, 22)
REF_FRAMES = motion.REF_FRAMES


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------------- inputs

def find_run_tif(mouse_dir: Path, run: str) -> Path:
    """The as-recorded binned tiff of run '260921_006' (or 'NavRear..._260921_006')."""
    mouse_dir = Path(mouse_dir)
    name = run if run.startswith(mouse_dir.name) else f"{mouse_dir.name}_{run}"
    folder = mouse_dir / name
    if not folder.is_dir():
        raise SystemExit(f"run folder not found: {folder}")
    tif = pick_binned_tif(folder, 2.0, log=log)
    if tif is None:
        raise SystemExit(f"no binned tiff in {folder}")
    return tif


def open_tif(path: Path):
    import tifffile as tf
    try:
        arr = tf.memmap(str(path), mode="r")
    except Exception:
        arr = tf.imread(str(path))
    return arr[None] if arr.ndim == 2 else arr


def load_variants(mouse_dir: Path, tif: Path, wanted: list[str]) -> list[dict]:
    """[{key, label, movie, note}] in the order asked, for the variants that exist."""
    out = []
    run_name = tif.parent.name
    rc = Path(mouse_dir) / "reg_compare" / run_name
    for key in wanted:
        if key == "raw":
            out.append({"key": key, "label": LABEL[key], "movie": open_tif(tif), "note": ""})
        elif key == "aqua":
            reg = open_registered(tif, log=log)
            if reg is not None:
                out.append({"key": key, "label": LABEL[key], "movie": reg,
                            "note": f"AQuA2's registered movie: {aqua_files(tif)['res'].name}"})
            else:
                log("  the AQuA2 result with its registered movie is not here — emulating AQuA2's "
                    "registration (integer shift per frame against the mean of the first 10 frames) ...")
                raw = open_tif(tif)
                t0 = time.time()
                mv, info = motion.aqua2_style_register(raw, gap=5, max_shift=40, log=log)
                log(f"  emulated in {time.time() - t0:.0f} s -> {mv.shape[1]} x {mv.shape[2]}")
                out.append({"key": key, "label": LABEL[key] + " — emulated", "movie": mv,
                            "note": "emulated (afferent2p.motion.aqua2_style_register): AQuA2's own result was not readable here"})
        else:
            sub = {"legacy": "D_legacy", "prereg": "B_prereg"}[key]
            p = rc / sub / tif.name
            if p.exists() and p.stat().st_size > 0:
                out.append({"key": key, "label": LABEL[key], "movie": open_tif(p), "note": str(p)})
            else:
                log(f"  {key}: {p} not found — panel left out")
    if not out:
        raise SystemExit("nothing to show")
    return out


def read_behavior(tif: Path, T: int):
    import pandas as pd
    p = aqua_files(tif)["behavior"]
    if not p.exists():
        return None
    b = pd.read_csv(p)
    if "velocity_mean" not in b.columns:
        return None
    b = b.iloc[:T]
    return b


# ----------------------------------------------------------------------------- geometry

def reference(movie) -> np.ndarray:
    n = min(REF_FRAMES, int(movie.shape[0]))
    return np.mean(np.stack([np.asarray(movie[i], np.float32) for i in range(n)]), axis=0)


def auto_region(ref: np.ndarray, h: int, w: int, margin: int = 16) -> tuple[int, int]:
    """Top-left corner of the h x w window of `ref` with the most edge energy
    (the fibres' outlines are what shows motion), kept `margin` px off the borders."""
    from scipy.ndimage import gaussian_filter, uniform_filter
    f = gaussian_filter(ref, 1.5)
    gy, gx = np.gradient(f)
    g = np.hypot(gy, gx)
    g = g - np.percentile(g, 50)                    # structure above the noise floor only
    g[g < 0] = 0
    e = uniform_filter(g, size=(h, w), mode="constant")       # mean over the window centred on each pixel
    H, W = ref.shape
    cy0, cy1 = margin + h // 2, H - margin - (h - h // 2)
    cx0, cx1 = margin + w // 2, W - margin - (w - w // 2)
    sub = e[cy0:cy1, cx0:cx1]
    iy, ix = np.unravel_index(int(np.argmax(sub)), sub.shape)
    return int(iy + cy0 - h // 2), int(ix + cx0 - w // 2)


def offset_to(ref_raw: np.ndarray, ref_var: np.ndarray) -> tuple[int, int]:
    """Where the as-recorded reference's pixel (y, x) sits in the variant's
    reference: (y + dy, x + dx). Whole pixels (the variants are whole-pixel
    crops / shifts of the same frames)."""
    h = min(ref_raw.shape[0], ref_var.shape[0]); w = min(ref_raw.shape[1], ref_var.shape[1])
    M = motion._Matcher(ref_raw[:h, :w], max_shift=48, down=1, sigma=1.0, hp=30.0)   # same size: the top-left corners
    dy, dx, _ = M.shift(ref_var[:h, :w])
    return int(round(dy)), int(round(dx))


def contours_of(img: np.ndarray, level_pct: float = 78.0):
    from scipy.ndimage import gaussian_filter
    from skimage import measure
    f = gaussian_filter(np.asarray(img, np.float32), 1.5)
    lvl = float(np.percentile(f, level_pct))
    return [c for c in measure.find_contours(f, lvl) if len(c) >= 12]


# ----------------------------------------------------------------------------- drawing

def _font(size: int):
    from PIL import ImageFont
    import matplotlib
    try:
        return ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf"), size)
    except Exception:
        return ImageFont.load_default()


def to8(crop: np.ndarray, lo: float, hi: float, gamma: float = 0.75) -> np.ndarray:
    v = (np.asarray(crop, np.float32) - lo) / max(hi - lo, 1e-6)
    v = np.clip(v, 0, 1) ** gamma
    return (v * 255 + 0.5).astype(np.uint8)


def trace_strip(width: int, height: int, shifts: dict, beh, T: int, spf: float, colors: dict, labels: dict):
    """The static trace image (RGB array) and the pixel box of the time axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dpi = 100
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi, facecolor=np.array(BG) / 255)
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1], left=0.075, right=0.985, top=0.9, bottom=0.2, hspace=0.12)
    t = np.arange(T) * spf / 60.0
    ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1], sharex=ax)
    for a in (ax, ax2):
        a.set_facecolor(np.array(BG) / 255)
        for s in a.spines.values():
            s.set_color("0.5")
        a.tick_params(colors="0.85", labelsize=8)
    if beh is not None:
        run = (beh["still_fraction"].values < 0.5)
        for a in (ax, ax2):
            a.fill_between(t, 0, 1, where=run, transform=a.get_xaxis_transform(), color="0.35", alpha=0.45, lw=0)
        ax2.fill_between(t, 0, beh["velocity_mean"].values, color="0.75", lw=0)
        ax2.set_ylabel("speed", color="0.85", fontsize=8)
    ymax = 1.0
    for k, s in shifts.items():
        a = np.hypot(s["dy"], s["dx"])
        ax.plot(t, a, color=np.array(colors[k]) / 255, lw=1.0, label=labels[k].split(" (")[0].split(" —")[0])
        ymax = max(ymax, float(np.percentile(a, 99.5)) * 1.15)
    ax.set_ylim(0, ymax); ax.set_xlim(0, t[-1])
    ax.set_ylabel("|shift| vs the movie's\nreference (px)", color="0.85", fontsize=8)
    ax.legend(loc="upper right", fontsize=7.5, ncol=len(shifts), frameon=False, labelcolor="0.9")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax2.set_xlabel("min  (grey bands: running)", color="0.85", fontsize=8)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
    bb = ax.get_window_extent(); bb2 = ax2.get_window_extent()
    box = (float(bb.x0), float(height - bb.y1), float(bb.x1), float(height - bb2.y0))   # x0, y_top, x1, y_bottom in image px
    plt.close(fig)
    return img, box, float(t[-1])


def render(run_label: str, variants: list[dict], shifts: dict, beh, region: tuple, scale: float, fps: float,
           spf: float, out_mp4: Path, gif_seconds: float = 5.0, bitrate: str = "8000k", max_frames: int = 0, log=log):
    from PIL import Image, ImageDraw
    y0, x0, h, w = region
    T = min(int(v["movie"].shape[0]) for v in variants)
    if max_frames:
        T = min(T, int(max_frames))
    shifts = {k: {kk: np.asarray(vv)[:T] for kk, vv in sv.items()} for k, sv in shifts.items()}   # the movies may differ in length
    if beh is not None:
        beh = beh.iloc[:T]
    ph, pw = int(round(h * scale)), int(round(w * scale))
    gap = 18
    n = len(variants)
    title_h, status_h, strip_h, foot_h = 40, 36, 250, 10
    W = max(n * pw + (n + 1) * gap, 1200)
    W += (-W) % 16
    H = title_h + status_h + 24 + ph + gap + strip_h + foot_h
    H += (-H) % 16
    thumb_h = strip_h - 20
    raw_ref = variants[0]["ref"]
    tw = int(round(raw_ref.shape[1] * thumb_h / raw_ref.shape[0]))
    strip_w = W - tw - 3 * gap
    strip_img, box, tmax = trace_strip(strip_w, strip_h, shifts, beh, T, spf, {v["key"]: COLOR[v["key"]] for v in variants},
                                       {v["key"]: v["label"] for v in variants})
    f_title, f_lab, f_txt, f_small = _font(19), _font(17), _font(16), _font(13)
    # per-variant prepared things
    for v in variants:
        oy, ox = v["offset"]
        v["box"] = (y0 + oy, x0 + ox)
        ref_crop = v["ref"][y0 + oy:y0 + oy + h, x0 + ox:x0 + ox + w]
        v["lo"], v["hi"] = float(np.percentile(ref_crop, 0.5)), float(np.percentile(ref_crop, 99.8))
        v["contours"] = [[(float(c[1] * scale), float(c[0] * scale)) for c in cc] for cc in contours_of(ref_crop)]
    writer = None
    try:
        import imageio
        writer = imageio.get_writer(str(out_mp4), fps=fps, codec="libx264", quality=None, bitrate=bitrate,
                                    pixelformat="yuv420p", macro_block_size=16)
    except Exception as e:                                             # noqa: BLE001
        log(f"  mp4 writer not available ({e}) — writing the GIF excerpt and the PNG only")
    strip_x = tw + 2 * gap
    strip_y = title_h + status_h + 24 + ph + gap
    gif_frames = []
    # the excerpt: the 6 s (at `fps`) in which the as-recorded frames move the most
    nex = int(round(gif_seconds * fps))
    s0 = shifts[variants[0]["key"]]
    mob = np.hypot(np.diff(s0["dy"], prepend=s0["dy"][0]), np.diff(s0["dx"], prepend=s0["dx"][0]))
    if T > nex:
        csum = np.convolve(mob, np.ones(nex), mode="valid")
        ex0 = int(np.argmax(csum))
    else:
        ex0 = 0
    ex1 = min(T, ex0 + nex)
    png_frame = ex0 + nex // 2
    gif_fps = min(fps, 10.0)
    gif_step = max(1, int(round(fps / gif_fps)))
    for t in range(T):
        canvas = Image.new("RGB", (W, H), BG)
        d = ImageDraw.Draw(canvas)
        d.text((gap, 10), f"{run_label} — 2-s binned frames; zoom x{scale:g} on the {w} x {h} px region with the most structure; "
                          f"cyan: outline of each movie's reference (mean of its first {REF_FRAMES} frames)", fill=(230, 230, 230), font=f_title)
        sp = ""
        if beh is not None and t < len(beh):
            vel = float(beh["velocity_mean"].values[t]); st = float(beh["still_fraction"].values[t])
            sp = f"   speed {vel:.0f}   {'RUNNING' if st < 0.5 else 'still'}"
        secs = t * spf
        d.text((gap, title_h + 6), f"frame {t + 1} / {T}   {int(secs // 60):d}:{secs % 60:04.1f} min" + sp, fill=(200, 200, 200), font=f_lab)
        for i, v in enumerate(variants):
            by, bx = v["box"]
            fr = np.asarray(v["movie"][t], np.float32)
            crop = fr[by:by + h, bx:bx + w]
            img = Image.fromarray(to8(crop, v["lo"], v["hi"])).resize((pw, ph), Image.BICUBIC).convert("RGB")
            dd = ImageDraw.Draw(img)
            for poly in v["contours"]:
                dd.line(poly, fill=OUTLINE, width=1)
            s = shifts[v["key"]]
            dy, dx = float(s["dy"][t]), float(s["dx"][t])
            amp = float(np.hypot(dy, dx))
            word = "residual" if v["key"] in ("aqua",) else "shift vs reference"
            dd.rectangle([0, ph - 26, pw, ph], fill=(0, 0, 0))
            dd.text((6, ph - 24), f"{word}: dy {dy:+.1f}  dx {dx:+.1f}  ({amp:.1f} px)", fill=COLOR[v["key"]], font=f_txt)
            px = gap + i * (pw + gap)
            py = title_h + status_h + 24
            canvas.paste(img, (px, py))
            d.rectangle([px - 1, py - 1, px + pw, py + ph], outline=COLOR[v["key"]], width=2)
            d.text((px + 4, py - 23), v["label"], fill=COLOR[v["key"]], font=f_lab)
        # whole field (as recorded) with the zoom region
        raw = variants[0]
        fr = np.asarray(raw["movie"][t], np.float32)
        th = Image.fromarray(to8(fr, raw["lo"], raw["hi"] * 1.3)).resize((tw, thumb_h), Image.BILINEAR).convert("RGB")
        ImageDraw.Draw(th).rectangle([x0 * tw / fr.shape[1], y0 * thumb_h / fr.shape[0],
                                      (x0 + w) * tw / fr.shape[1], (y0 + h) * thumb_h / fr.shape[0]], outline=(255, 230, 0), width=2)
        canvas.paste(th, (gap, strip_y + 10))
        d.text((gap, strip_y + thumb_h + 10), "as recorded, whole field; box = the zoom", fill=(170, 170, 170), font=f_small)
        canvas.paste(Image.fromarray(strip_img), (strip_x, strip_y))
        # cursor
        bx0, byt, bx1, byb = box
        cx = strip_x + bx0 + (bx1 - bx0) * (t * spf / 60.0) / max(tmax, 1e-9)
        d.line([(cx, strip_y + byt), (cx, strip_y + byb)], fill=(255, 80, 80), width=2)
        arr = np.asarray(canvas)
        if writer is not None:
            writer.append_data(arr)
        if ex0 <= t < ex1 and (t - ex0) % gif_step == 0:
            gif_frames.append(canvas.resize((int(W * 0.45) // 2 * 2, int(H * 0.45) // 2 * 2), Image.BILINEAR)
                              .quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE))
        if t == png_frame:
            canvas.save(str(out_mp4.with_suffix(".png")))
        if t % 100 == 99:
            log(f"  frame {t + 1} / {T}")
    if writer is not None:
        writer.close()
        log(f"wrote {out_mp4}")
    gif = Path(str(out_mp4.with_suffix("")) + "_excerpt.gif")
    if gif_frames:
        gif_frames[0].save(str(gif), save_all=True, append_images=gif_frames[1:], duration=int(1000 * gif_step / fps), loop=0, optimize=True)
        log(f"wrote {gif} (frames {ex0 + 1}-{ex1}, the most mobile {gif_seconds:g} s of the movie)")
    log(f"wrote {out_mp4.with_suffix('.png')} (frame {png_frame + 1})")
    return {"mp4": out_mp4 if writer is not None else None, "gif": gif, "png": out_mp4.with_suffix(".png"),
            "excerpt_frames": (ex0, ex1)}


# ----------------------------------------------------------------------------- main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mouse_dir")
    ap.add_argument("--run", required=True, help="run id, e.g. 260921_006")
    ap.add_argument("--variants", default=",".join(VARIANTS), help="comma list of raw,aqua,legacy,prereg (those found are shown)")
    ap.add_argument("--region", default="auto", help="y0,x0,height,width in the as-recorded frame, or auto")
    ap.add_argument("--region-size", default="160x240", help="height x width of the auto region (px)")
    ap.add_argument("--scale", type=float, default=2.5)
    ap.add_argument("--fps", type=float, default=20.0)
    ap.add_argument("--gif-seconds", type=float, default=5.0)
    ap.add_argument("--bitrate", default="8000k", help="mp4 bitrate (ffmpeg -b:v)")
    ap.add_argument("--max-frames", type=int, default=0, help="render only the first N frames (a quick test)")
    ap.add_argument("--out", default=None, help="mp4 path (default: <mouse>\\reg_compare\\<run>_reg_demo.mp4)")
    a = ap.parse_args(argv)
    mouse_dir = Path(a.mouse_dir)
    tif = find_run_tif(mouse_dir, a.run)
    log(f"run: {tif}")
    wanted = [k.strip() for k in a.variants.split(",") if k.strip()]
    bad = [k for k in wanted if k not in VARIANTS]
    if bad:
        raise SystemExit(f"unknown variant(s) {bad}; choose from {VARIANTS}")
    if "raw" not in wanted:
        wanted = ["raw"] + wanted
    variants = load_variants(mouse_dir, tif, wanted)
    T = min(int(v["movie"].shape[0]) for v in variants)
    tag = None
    try:
        from afferent2p.aqua_io import parse_bin_tag
        tag = parse_bin_tag(tif)
    except Exception:
        pass
    spf = float(tag[1]) if tag else 2.0
    beh = read_behavior(tif, T)
    # references, offsets, shifts
    for v in variants:
        v["ref"] = reference(v["movie"])
    raw_ref = variants[0]["ref"]
    variants[0]["offset"] = (0, 0)
    for v in variants[1:]:
        v["offset"] = offset_to(raw_ref, v["ref"])
        log(f"  {v['key']}: {v['movie'].shape[1]} x {v['movie'].shape[2]} px, sits at {v['offset']} px from the as-recorded frame")
    shifts = {}
    for v in variants:
        t0 = time.time()
        s = motion.frame_shifts(v["movie"], ref=v["ref"], with_sharpness=False)
        amp = np.hypot(s["dy"], s["dx"])
        shifts[v["key"]] = s
        log(f"  {v['key']}: |shift| vs its reference median {np.median(amp):.2f} px, 95th pct {np.percentile(amp, 95):.2f}, "
            f"max {amp.max():.2f} ({time.time() - t0:.0f} s)")
    hh, ww = (int(x) for x in a.region_size.lower().split("x"))
    if a.region == "auto":
        y0, x0 = auto_region(raw_ref, hh, ww)
        region = (y0, x0, hh, ww)
    else:
        region = tuple(int(x) for x in a.region.split(","))
    log(f"  zoom region (as recorded): y {region[0]}-{region[0] + region[2]}, x {region[1]}-{region[1] + region[3]}")
    # every variant must contain the region after its offset
    for v in variants:
        oy, ox = v["offset"]
        Hv, Wv = int(v["movie"].shape[1]), int(v["movie"].shape[2])
        if not (0 <= region[0] + oy and region[0] + oy + region[2] <= Hv and 0 <= region[1] + ox and region[1] + ox + region[3] <= Wv):
            raise SystemExit(f"the zoom region falls outside the {v['key']} movie ({Hv} x {Wv}, offset {v['offset']}); choose --region inside")
    if a.out:
        out_mp4 = Path(a.out)
    else:
        rc = mouse_dir / "reg_compare"
        out_mp4 = (rc if rc.is_dir() else tif.parent) / f"{tif.parent.name}_reg_demo.mp4"
    out_mp4.parent.mkdir(parents=True, exist_ok=True)
    res = render(tif.parent.name, variants, shifts, beh, region, a.scale, a.fps, spf, out_mp4, a.gif_seconds, a.bitrate,
                 max_frames=a.max_frames, log=log)
    notes = [f"{v['label']}: {v['note']}" for v in variants if v["note"]]
    for n_ in notes:
        log("  " + n_)
    for v in variants:
        m = v["movie"]
        if hasattr(m, "close"):
            m.close()
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
