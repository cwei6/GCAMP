r"""qc_report — the EVIDENCE as a report of its own: "is the pipeline detecting
and analysing the right thing?" (v32; before that it was section 3 of the
main report and made it too long to read).

    python aqua_report.py <mouse folder> --qc      -> <mouse>\aqua_analysis\<mouse>_qc.html

Structure
  * top: a QC summary table, one row per recording — minutes, excluded frames,
    events detected / kept / dropped, CFUs kept / noise, "looks" counts (fibre
    / off-fibre / blob, afferent2p.cfu_shape), curation status, crop match,
    movie found (registered = AQuA2's own movie out of the .mat), exports
    present, locomotion trials, the motion check (afferent2p.motion: shift of
    the raw frames, jumps per minute, deformation left after AQuA2, event
    onsets on moving frames vs chance, verdict) — with warning chips (stale
    curation, > 30 % non-fibre CFUs, missing event export, > 10 % frames
    excluded, no events, < 3 locomotion trials, events follow the motion);
    columns sort on click; every row links to its recording;
  * a sticky left navigation: experiments -> recordings with counts;
  * per recording, collapsed by default: the motion figure (shift and jump
    traces, residual after AQuA2, running speed, event onsets), the detection
    map (CFU outlines coloured by "looks") and the raster with the locomotion
    trials marked, up
    to 6 clips, the event gallery (12 montages per page with "more", filters
    all / kept / dropped / in non-fibre CFU / in a locomotion trial / top
    amplitude / long, sort by onset / amplitude / duration) and the CFU table.
  Each recording's images sit in a <template> until its section is opened, so
  the page opens fast however many events there are. When the whole file would
  exceed ~80 MB, one page per experiment is written plus an index.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from .aqua_io import relocate
from .report import (_esc, _figure, _num, _mmss, _short_run, _cond_name, _read_csv, _fig_uri, _stretch,
                     load_run, raster_figure, event_montage, event_clip, choose_clip_events, _tiff_to_recorded_time,
                     report_context, page_html, KEPT, DROPPED, EXCL, CFU_COL, INK)

LOOK_COL = {"fibre": CFU_COL, "off-fibre": "#e34948", "blob": "#e87ba4"}


# ----------------------------------------------------------------------------- per-run pieces

def _looks_for(d: dict, um_per_px: float) -> dict:
    """{cfu number (1-based): feature dict} from cfu_shape on the run's mean image."""
    c = d.get("cfu")
    if c is None or c.get("maps") is None or c["maps"].ndim != 3 or c["n"] == 0:
        return {}
    img = d.get("meanimg")
    if img is None or img.shape != c["maps"].shape[1:]:
        img = c.get("datPro")
    try:
        from .cfu_shape import cfu_features
        masks = [c["maps"][k] > 0.1 for k in range(c["n"])]
        feats = cfu_features(masks, img, um_per_px)
        return {k + 1: f for k, f in enumerate(feats)}
    except Exception:
        return {}


def detection_map_looks(d: dict, looks: dict) -> str | None:
    """Mean image with the event boxes and the CFU outlines coloured by 'looks'."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    img = d.get("meanimg") if d.get("meanimg") is not None else d.get("raw_mean")
    if img is None:
        return None
    ev = d.get("events")
    H, W = img.shape
    fig, ax = plt.subplots(figsize=(9.2, 9.2 * H / W + 0.25))
    v0, v1 = _stretch(img, 0.5, 99.8)
    ax.imshow(img, cmap="gray", vmin=v0, vmax=v1, interpolation="nearest")
    ax.set_axis_off()
    c = d.get("cfu")
    if c is not None and c.get("maps") is not None and c["maps"].ndim == 3 and c["maps"].shape[1:] == img.shape:
        for i in range(c["n"]):
            m = c["maps"][i]
            if not np.any(m > 0):
                continue
            rejected = (i + 1) in d["rejected_cfu"]
            col = LOOK_COL.get(looks.get(i + 1, {}).get("suggest", "fibre"), CFU_COL)
            try:
                ax.contour(m, levels=[0.3 * float(m.max())], colors=[col], linewidths=0.7 if not rejected else 0.6,
                           linestyles=["dashed" if rejected else "solid"])
                ys, xs = np.nonzero(m > 0.3 * float(m.max()))
                if ys.size:
                    ax.text(xs.min(), ys.min() - 2, f"C{i + 1}", color=col, fontsize=6.5, va="bottom", ha="left")
            except Exception:
                pass
    if ev is not None and len(ev) and {"y0", "y1", "x0", "x1"} <= set(ev.columns):
        for r in ev.itertuples(index=False):
            st = getattr(r, "status", "kept")
            col = KEPT if st == "kept" else (DROPPED if st == "noise CFU" else EXCL)
            x0, x1, y0, y1 = float(r.x0) - 1, float(r.x1) - 1, float(r.y0) - 1, float(r.y1) - 1
            ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0 + 1, y1 - y0 + 1, fill=False, ec=col, lw=0.9, ls="-" if st == "kept" else "--"))
            ax.text(x1 + 1.5, y0 - 1.5, str(int(r.event)), color=col, fontsize=8, fontweight="bold", va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
    ax.set_title(f"{d['short']} — detected events (boxes, numbered) and CFU outlines (C#) coloured by what they look like", fontsize=9, loc="left", color=INK)
    fig.tight_layout(pad=0.3)
    return _fig_uri(fig, dpi=100)


def _motion_for(d: dict, log=print) -> tuple:
    """(summary, per-frame table) of the run's motion check — the cached files, else measured now."""
    from .motion import read_motion, run_motion
    tif = d["tif"]
    try:
        summ = run_motion(tif, beh=d.get("beh"), events=d.get("events"), log=log)   # cached unless the tiff / result is newer
    except Exception as e:
        log(f"    motion check failed ({e})")
        summ = None
    s2, df = read_motion(tif)
    return (summ if summ is not None else s2), df


def _run_numbers(d: dict, trials: list, looks: dict, motion: dict | None = None) -> dict:
    ev = d.get("events")
    n_det = 0 if ev is None else len(ev)
    st = ev["status"] if ev is not None and "status" in ev.columns else None
    n_kept = int((st == "kept").sum()) if st is not None else n_det
    n_noise = int((st == "noise CFU").sum()) if st is not None else 0
    n_exc = int((st == "excluded frames").sum()) if st is not None else 0
    cur = d.get("curation")
    c = d.get("cfu")
    n_cfu = int(c["n"]) if c else 0
    n_rej = len(d.get("rejected_cfu", set()))
    kinds = {"fibre": 0, "off-fibre": 0, "blob": 0}
    for k, f in looks.items():
        if k not in d.get("rejected_cfu", set()):
            kinds[f.get("suggest", "fibre")] = kinds.get(f.get("suggest", "fibre"), 0) + 1
    frames = int(d.get("frames") or 0)
    n_ex_frames = len(d.get("excluded_1based", []))
    return {"minutes": frames * d["spf"] / 60.0 if frames else np.nan, "frames": frames, "excluded_frames": n_ex_frames,
            "excluded_pct": 100.0 * n_ex_frames / frames if frames else np.nan,
            "n_detected": n_det, "n_kept": n_kept, "n_noise": n_noise, "n_excluded": n_exc,
            "n_cfu": n_cfu, "n_cfu_noise": n_rej, "n_cfu_kept": n_cfu - n_rej, "looks_fibre": kinds["fibre"], "looks_off": kinds["off-fibre"],
            "looks_blob": kinds["blob"],
            "curation": ("stale" if cur and cur.get("stale") else "curated" if cur else "not curated"),
            "movie": d.get("movie") is not None, "movie_src": d.get("movie_src"),
            "events_csv": ev is not None and "y0" in (ev.columns if ev is not None else []),
            "offset": d.get("crop_in_raw") or d.get("offset"), "n_trials": len(trials), "n_trials_strict": sum(1 for t in trials if t[1] == "strict"),
            "motion": motion}


def _warnings(n: dict) -> list:
    w = []
    if n["curation"] == "stale":
        w.append(("stale curation", "warn"))
    non = n["looks_off"] + n["looks_blob"]
    if n["n_cfu_kept"] and non / n["n_cfu_kept"] > 0.3:
        w.append((f"{100 * non / n['n_cfu_kept']:.0f} % non-fibre CFUs", "warn"))
    if not n["events_csv"]:
        w.append(("no event export", "warn"))
    if np.isfinite(n["excluded_pct"]) and n["excluded_pct"] > 10:
        w.append((f"{n['excluded_pct']:.0f} % frames excluded", "warn"))
    if n["n_kept"] == 0:
        w.append(("0 events", "note"))
    if n["n_trials"] < 3:
        w.append((f"{n['n_trials']} loco trials", "note"))
    if not n["movie"]:
        w.append(("movie not readable", "note"))
    m = n.get("motion")
    if m:
        v = m.get("verdict", ["ok", ""])
        if v[0] == "warn":
            w.append(("events follow the motion" if (m.get("events", {}).get("ratio_motion") or 0) >= 2 else "strong motion", "warn"))
        elif v[0] == "check":
            w.append(("motion: check", "note"))
    return w


def _cfu_table(d: dict, looks: dict, loco_cfu: dict) -> str:
    c = d.get("cfu")
    if c is None or c["n"] == 0:
        return "<p class='small'>No CFUs in this recording.</p>"
    rows = []
    for k in range(c["n"]):
        f = looks.get(k + 1, {})
        rejected = (k + 1) in d.get("rejected_cfu", set())
        lr = loco_cfu.get(k + 1)
        rows.append(f"<tr><td>C{k + 1}</td><td>{len(c['events'][k])}</td>"
                    f"<td><span class='lk {f.get('suggest', 'fibre').replace(' ', '-')}'>{_esc(f.get('suggest', '?'))}</span></td>"
                    f"<td>{_num(100 * f.get('on_fibre', np.nan), 0)}</td><td>{_num(f.get('bright_ratio', np.nan), 1)}</td><td>{_num(f.get('elongation', np.nan), 1)}</td>"
                    f"<td>{_num(f.get('area_um2', np.nan), 0)}</td>"
                    f"<td>{('yes (p=' + _num(lr['p_dff'], 2) + ')') if lr and lr.get('loco_responsive') else ('no (p=' + _num(lr['p_dff'], 2) + ')' if lr else '–')}</td>"
                    f"<td>{'<span class=warn>noise</span>' if rejected else 'kept'}</td><td class='small'>{_esc(f.get('reason', ''))}</td></tr>")
    return ("<table class='cfus'><thead><tr><th>CFU</th><th>Events</th><th>Looks</th><th>On fibre %</th><th>Bright ratio</th><th>Elongation</th>"
            "<th>Area µm²</th><th>Loco-responsive</th><th>Status</th><th>Why</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _gallery(d: dict, looks: dict, trials: list, max_events: int, log=print) -> tuple[str, int]:
    """Every event's montage as a card with data attributes for the filters; the cards go in a <template>."""
    ev = d.get("events")
    movie = d.get("movie")
    if ev is None or len(ev) == 0:
        return "<p class='small'>No detected events in this recording.</p>", 0
    if movie is None or not {"y0", "y1", "x0", "x1", "t0_frame"} <= set(ev.columns):
        return ("<p class='small warn'>No montages: " + ("the movie is not readable here" if movie is None else
                                                          "no per-event locations (run run_cfu_batch.m again)") + ".</p>", 0)
    order = ev.sort_values("t0_frame")
    note = ""
    if max_events and len(order) > max_events:
        order = order.iloc[:max_events]
        note = f"<p class='small'>Montages for the first {max_events} events (by onset) of {len(ev)}; raise the limit to see all.</p>"
    amps = ev["dff_max"].to_numpy(float) if "dff_max" in ev.columns else np.full(len(ev), np.nan)
    top_amp = np.nanpercentile(amps, 80) if np.isfinite(amps).any() else np.inf
    spf = d["spf"]
    trial_windows = [(t0, t0 + 25) for t0, _ in trials]            # recorded frames (0-based) of the trial windows
    fr = d.get("frames_csv")
    cards = []
    for r in order.itertuples(index=False):
        try:
            uri, m = event_montage(d, r)
        except Exception as e:
            cards.append(f"<p class='small warn'>event #{int(r.event)}: montage failed ({_esc(e)})</p>")
            continue
        st = getattr(r, "status", "kept")
        cfu = d["event_cfu"].get(int(r.event))
        lk = looks.get(cfu, {}).get("suggest", "none") if cfu else "none"
        amp = float(getattr(r, "dff_max", np.nan))
        t0_rec = int(round(m["onset_s"] / spf)) if fr is not None else int(m["t0"]) - 1
        in_trial = int(any(a <= t0_rec < b for a, b in trial_windows))
        cap = (f"<b>#{int(r.event)}</b> <span class='st {st.replace(' ', '-')}'>{_esc(st)}</span> — onset {_mmss(m['onset_s'])} (frame {m['t0']}), "
               f"{_num(m['dur_s'], 0)} s, AQuA2 peak dF/F {_num(amp)}, area {_num(getattr(r, 'area_um2', np.nan), 0)} µm²"
               + (f", CFU C{cfu} <span class='lk {lk.replace(' ', '-')}'>{_esc(lk)}</span>" if cfu else ", no CFU")
               + (" <span class='st trial'>in a locomotion trial</span>" if in_trial else ""))
        cards.append(f"<figure class='montage' data-status='{_esc(st)}' data-looks='{_esc(lk)}' data-trial='{in_trial}' "
                     f"data-amp='{amp if np.isfinite(amp) else -1}' data-top='{int(np.isfinite(amp) and amp >= top_amp)}' "
                     f"data-dur='{m['dur_s']:.1f}' data-onset='{m['onset_s']:.1f}'>"
                     f"<a class='zoom' href='#'><img loading='lazy' src='{uri}' alt='event {int(r.event)}'></a><figcaption>{cap}</figcaption></figure>")
    controls = ("<div class='galctl'><span class='small'>show:</span> "
                + "".join(f"<button data-f='{f}'>{lab}</button>" for f, lab in (("all", "all"), ("kept", "kept"), ("dropped", "dropped"),
                                                                                   ("nonfibre", "in non-fibre CFU"), ("trial", "in a locomotion trial"),
                                                                                   ("top", "top amplitude"), ("long", "long > 60 s")))
                + " <span class='small'>sort:</span> " + "".join(f"<button data-s='{s}'>{lab}</button>" for s, lab in (("onset", "onset"), ("amp", "amplitude"), ("dur", "duration")))
                + " <span class='galn small'></span></div>")
    html = (note + controls + "<div class='galbody'></div><button class='more'>show 12 more</button>"
            + "<template class='gal'>" + "".join(cards) + "</template>")
    return html, len(cards)


def run_section(d: dict, exp: str, cond: str, day_labels: dict, looks: dict, trials: list, loco_cfu: dict, numbers: dict,
                clips_per_run: int, max_events: int, log=print) -> str:
    rid = f"run-{_esc(d['short'])}"
    warns = _warnings(numbers)
    chips = "".join(f"<span class='chip {k}'>{_esc(w)}</span>" for w, k in warns)
    head = (f"<details class='rec' id='{rid}'><summary><span class='rn'>{_esc(d['short'])}</span> <span class='cond'>{_esc(_cond_name(cond, day_labels))}</span> "
            f"<span class='small'>{_num(numbers['minutes'], 1)} min · {numbers['n_kept']} / {numbers['n_detected']} events kept · "
            f"{numbers['n_cfu_kept']} CFUs</span> {chips}</summary><div class='recbody'>")
    cur_txt = {"stale": "curation STALE (older than the result — curate again)", "curated": "curated", "not curated": "not curated yet"}[numbers["curation"]]
    cur = d.get("curation")
    if cur and numbers["curation"] == "curated":
        cur_txt += f" ({len(cur.get('rejected', []))} CFU marked noise, {len(cur.get('added', []))} added)"
    body = [f"<p class='small'>{_num(numbers['minutes'], 1)} min, {numbers['frames']} frames at {d['spf']:.3f} s; {numbers['excluded_frames']} excluded frame(s) "
            f"({len(d['flagged_1based'])} flagged by the brightness step, mode: {_esc(d['bad_frames_mode'])}); detected {numbers['n_detected']} events, kept "
            f"{numbers['n_kept']}" + (f", dropped {numbers['n_noise']} in noise CFUs" if numbers["n_noise"] else "")
            + (f", dropped {numbers['n_excluded']} in excluded frames" if numbers["n_excluded"] else "")
            + f"; {cur_txt}" + (f"; AQuA2's crop at row {d['crop_in_raw'][0]}, column {d['crop_in_raw'][1]} of the raw frame" if d.get("crop_in_raw") else "")
            + (f"; movie shown: {'AQuA2 registered' if d.get('movie_src') == 'registered' else 'binned tiff as recorded'}" if d.get("movie") is not None else "")
            + (f"; raw frames registered before binning (unpack --register: shift median {_num((d['prereg'].get('stats') or {}).get('shift_median_px'), 1)} px, "
               f"border {d['prereg'].get('crop_in_raw_px', 0)} px cropped)" if d.get("prereg") else "")
            + f"; {numbers['n_trials']} locomotion trials ({numbers['n_trials_strict']} strict).</p>"]
    for n in d.get("notes", []):
        body.append(f"<p class='small warn'>{_esc(n)}</p>")
    m = numbers.get("motion")
    if m:
        e_ = m.get("events", {})
        reg = m.get("registered", {}).get("available")
        body.append(f"<p class='small'><b>Motion:</b> the raw frames sit {_num(m['shift_median'], 1)} px (median, up to {_num(m['shift_max'], 1)}) from AQuA2's "
                    f"reference; {_num(m['jumps_ge2_per_min'], 1)} jumps of ≥ {_num(m.get('jump_threshold', 2), 1)} px per minute between consecutive frames"
                    + (f" ({_num(m['jumps_ge2_per_min_still'], 1)} still / {_num(m['jumps_ge2_per_min_run'], 1)} running)" if m.get("jumps_ge2_per_min_run") is not None else "")
                    + (f"; after AQuA2's registration the field is still deformed by {_num(m.get('nonrigid_median'), 1)} px (median) with "
                       f"{_num(m.get('nonrigid_jumps_per_min'), 1)} deformation jumps per minute" if reg else "; the registered movie was not readable, so only the raw motion is measured")
                    + f". Event onsets on moving frames: {e_.get('near_motion', 0)} of {e_.get('n_onsets', 0)}, {_num(e_.get('expected_near_motion'), 1)} expected by chance"
                    + (f" (ratio {_num(e_.get('ratio_motion'), 1)})" if e_.get("n_onsets") else "")
                    + f" — <b>{_esc(m['verdict'][0])}</b>: {_esc(m['verdict'][1])}.</p>")
    imgs = []
    if m and d.get("motion_df") is not None:
        try:
            from .motion import motion_figure
            imgs.append(_figure(motion_figure(d["motion_df"], m, d.get("events"), title=d["short"]),
                                "Motion of this recording: the rigid shift of every raw frame against AQuA2's reference (mean of its first 10 frames), the jump "
                                "between consecutive frames (dashed = the jump that counts as motion here), what is left after AQuA2's registration (purple = how "
                                "far the 3 × 4 blocks of the registered movie moved since the previous frame, orange = their deformation relative to the "
                                "reference) and the running speed; teal lines = event onsets, pink = frames that count as moving."))
        except Exception as e:
            imgs.append(f"<p class='small warn'>motion figure failed: {_esc(e)}</p>")
    try:
        imgs.append(_figure(detection_map_looks(d, looks), "Mean image of the registered movie with every detected event's box and the CFU outlines coloured by what "
                            "they look like (yellow fibre, red off the fibres, pink compact bright blob; dashed = marked noise).", cls="half"))
    except Exception as e:
        imgs.append(f"<p class='small warn'>detection map failed: {_esc(e)}</p>")
    try:
        imgs.append(_figure(raster_figure(d, trials=trials), "Running speed (▼ = locomotion trial onsets: black strict, grey relaxed), pupil size and the "
                            "event onsets (bars = duration) on one clock; shaded = excluded frames."))
    except Exception as e:
        imgs.append(f"<p class='small warn'>raster failed: {_esc(e)}</p>")
    ev = d.get("events")
    if clips_per_run and d.get("movie") is not None and ev is not None and len(ev) and {"y0", "y1", "x0", "x1", "t0_frame"} <= set(ev.columns):
        cards = []
        for eid in choose_clip_events(ev, clips_per_run):
            r = ev[ev["event"] == eid].iloc[0]
            try:
                uri = event_clip(d, r)
                cards.append(f"<figure class='clip'><img src='{uri}' alt='event {eid} clip'><figcaption>#{eid} {_esc(getattr(r, 'status', 'kept'))}, onset "
                             f"{_mmss(_tiff_to_recorded_time(int(r.t0_frame), d))}, {_num((float(r.t1_frame) - float(r.t0_frame) + 1) * d['spf'], 0)} s, "
                             f"peak dF/F {_num(getattr(r, 'dff_max', np.nan))}</figcaption></figure>")
            except Exception as e:
                cards.append(f"<p class='small warn'>clip #{eid} failed: {_esc(e)}</p>")
        src = "AQuA2's registered movie" if d.get("movie_src") == "registered" else "the binned movie as recorded (AQuA2's crop offset applied)"
        imgs.append("<h4>Clips</h4><div class='clips'>" + "".join(cards) + f"</div><p class='small'>{src} around the event (orange box = the "
                    "detected footprint while the event lasts), chosen across the amplitude range of the kept events.</p>")
    gal, n_cards = _gallery(d, looks, trials, max_events, log=log)
    imgs.append(f"<h4>Every detected event ({n_cards})</h4>" + gal)
    imgs.append("<h4>CFUs</h4>" + _cfu_table(d, looks, loco_cfu))
    # the heavy part waits in a template until the section is opened
    body.append("<div class='lazy'></div><template class='sec'>" + "".join(imgs) + "</template>")
    return head + "".join(body) + "</div></details>"


# ----------------------------------------------------------------------------- the page

def _mv(n: dict, key: str, sub: str | None = None):
    m = n.get("motion") or {}
    if sub:
        m = m.get(sub) or {}
    v = m.get(key)
    return v if v is not None else np.nan


def summary_table(rows: list) -> str:
    head = ["Recording", "Experiment", "Condition", "Minutes", "Excluded frames", "Events detected", "Kept", "Dropped (noise / excluded)", "CFUs kept / noise",
            "Looks fibre / off / blob", "Curation", "Crop offset", "Movie", "Event export", "Loco trials (strict)",
            "Shift median / max px", "Jumps per min", "Deformation px", "Onsets on moving frames (expected)", "Motion", "Warnings"]
    body = []
    for r in rows:
        n = r["numbers"]
        chips = "".join(f"<span class='chip {k}'>{_esc(w)}</span>" for w, k in r["warnings"])
        body.append("<tr>" + "".join(f"<td data-v='{v}'>{c}</td>" for v, c in (
            (r["short"], f"<a href='{r['href']}'>{_esc(r['short'])}</a>"), (r["exp"], _esc(r["exp"])), (r["cond"], _esc(r["cond"])),
            (n["minutes"], _num(n["minutes"], 1)), (n["excluded_frames"], f"{n['excluded_frames']} ({_num(n['excluded_pct'], 0)} %)"),
            (n["n_detected"], n["n_detected"]), (n["n_kept"], n["n_kept"]), (n["n_noise"] + n["n_excluded"], f"{n['n_noise']} / {n['n_excluded']}"),
            (n["n_cfu_kept"], f"{n['n_cfu_kept']} / {n['n_cfu_noise']}"), (n["looks_off"] + n["looks_blob"], f"{n['looks_fibre']} / {n['looks_off']} / {n['looks_blob']}"),
            (n["curation"], _esc(n["curation"])), (0, f"({n['offset'][0]}, {n['offset'][1]})" if n.get("offset") else "–"),
            (int(n["movie"]), ("registered" if n.get("movie_src") == "registered" else "raw") if n["movie"] else "no"),
            (int(n["events_csv"]), "yes" if n["events_csv"] else "no"),
            (n["n_trials"], f"{n['n_trials']} ({n['n_trials_strict']})"),
            (_mv(n, "shift_median"), f"{_num(_mv(n, 'shift_median'), 1)} / {_num(_mv(n, 'shift_max'), 1)}" if n.get("motion") else "–"),
            (_mv(n, "jumps_ge2_per_min"), _num(_mv(n, "jumps_ge2_per_min"), 1) if n.get("motion") else "–"),
            (_mv(n, "nonrigid_median"), _num(_mv(n, "nonrigid_median"), 1) if n.get("motion") else "–"),
            (_mv(n, "ratio_motion", "events"), (f"{n['motion']['events'].get('near_motion', 0)} of {n['motion']['events'].get('n_onsets', 0)} "
                                                 f"({_num(n['motion']['events'].get('expected_near_motion'), 1)})") if n.get("motion") else "–"),
            ({"ok": 0, "check": 1, "warn": 2}.get((n.get("motion") or {}).get("verdict", ["", ""])[0], -1),
             f"<span class='chip {'warn' if (n.get('motion') or {}).get('verdict', [''])[0] == 'warn' else 'note' if (n.get('motion') or {}).get('verdict', [''])[0] == 'check' else ''}'>"
             f"{_esc((n.get('motion') or {}).get('verdict', ['–'])[0])}</span>" if n.get("motion") else "–"),
            (len(r["warnings"]), chips))) + "</tr>")
    return ("<div class='scroll'><table class='qc sortable'><thead><tr>" + "".join(f"<th>{h} <span class='sorthint'>⇅</span></th>" for h in head)
            + "</tr></thead><tbody>" + "".join(body) + "</tbody></table></div>")


QC_CSS = """
body { padding-left: 230px; } @media (max-width: 900px) { body { padding-left: 0; } nav.side { display:none; } }
nav.side { position: fixed; left: 0; top: 0; bottom: 0; width: 220px; overflow-y: auto; background: #f3f3f1; border-right: 1px solid var(--line); padding: 14px 10px; font-size: 12.5px; }
nav.side b { display:block; margin: 10px 0 2px; } nav.side a { display:block; color: var(--accent); text-decoration: none; padding: 1px 4px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
nav.side a:hover { background: #e9e9e6; }
details.rec { border: 1px solid var(--line); border-radius: 6px; margin: 10px 0; background: #fff; }
details.rec > summary { padding: 8px 12px; } details.rec .recbody { padding: 4px 14px 14px; }
.rn { font-weight: 600; font-size: 15px; }
.chip { display:inline-block; padding: 0 7px; border-radius: 9px; font-size: 11.5px; margin-left: 4px; background: #ececea; color: #52514e; font-weight: normal; }
.chip.warn { background: #fde2e0; color: #b3261e; } .chip.note { background: #fff3d6; color: #8a5a00; }
.lk { display:inline-block; padding: 0 6px; border-radius: 8px; font-size: 11.5px; background:#fff3d6; color:#8a5a00; }
.lk.off-fibre { background:#fde2e0; color:#b3261e; } .lk.blob { background:#fbe3ee; color:#a3306b; } .lk.none { background:#ececea; color:#52514e; }
.st.trial { background:#e3f6ee; color:#0d6b43; }
table.qc { font-size: 12px; } table.qc th { cursor: pointer; white-space: nowrap; } .sorthint { color:#bbb; font-size: 10px; }
.scroll { overflow-x: auto; margin: 8px 0 14px; } .scroll table { margin: 0; min-width: 1700px; }
table.cfus { font-size: 12px; } table.cfus td.small { max-width: 320px; }
.galctl { margin: 6px 0; } .galctl button { font-size: 12px; margin: 2px 2px; } .galctl button.on { background: #2a78d6; color: #fff; border-color: #2a78d6; }
button.more { font-size: 13px; margin: 6px 0 12px; } button.more[hidden] { display:none; }
h4 { margin: 18px 0 4px; font-size: 14px; }
"""

QC_JS = """
(function(){
  // sortable summary table
  document.querySelectorAll('table.sortable th').forEach((th, i) => th.addEventListener('click', () => {
    const tb = th.closest('table').querySelector('tbody'); const rows = Array.from(tb.rows);
    const asc = th.dataset.asc !== '1'; th.dataset.asc = asc ? '1' : '0';
    rows.sort((a, b) => { const va = a.cells[i].dataset.v, vb = b.cells[i].dataset.v; const na = parseFloat(va), nb = parseFloat(vb);
      const c = (!isNaN(na) && !isNaN(nb)) ? na - nb : String(va).localeCompare(String(vb)); return asc ? c : -c; });
    rows.forEach(r => tb.appendChild(r)); }));
  // lazy sections + galleries
  function wire(det){
    const tpl = det.querySelector('template.sec'); const host = det.querySelector('.lazy');
    if (!tpl || !host || host.dataset.done) return;
    host.appendChild(tpl.content.cloneNode(true)); host.dataset.done = '1'; tpl.remove();
    host.querySelectorAll('a.zoom').forEach(a => a.addEventListener('click', e => { e.preventDefault(); window.openLightbox(a.querySelector('img').src); }));
    const gal = host.querySelector('template.gal'); if (!gal) return;
    const body = host.querySelector('.galbody'), more = host.querySelector('button.more'), ctl = host.querySelector('.galctl'), nlab = host.querySelector('.galn');
    const cards = Array.from(gal.content.children).filter(c => c.tagName === 'FIGURE'); gal.remove();
    let filt = 'all', sort = 'onset', shown = 0, cur = [];
    const pass = c => { const d = c.dataset; switch (filt) {
      case 'kept': return d.status === 'kept'; case 'dropped': return d.status !== 'kept';
      case 'nonfibre': return d.looks === 'off-fibre' || d.looks === 'blob'; case 'trial': return d.trial === '1';
      case 'top': return d.top === '1'; case 'long': return parseFloat(d.dur) > 60; default: return true; } };
    const key = c => sort === 'amp' ? -parseFloat(c.dataset.amp) : sort === 'dur' ? -parseFloat(c.dataset.dur) : parseFloat(c.dataset.onset);
    function render(reset){ if (reset) { body.innerHTML = ''; shown = 0; cur = cards.filter(pass).sort((a, b) => key(a) - key(b)); }
      const next = cur.slice(shown, shown + 12); next.forEach(c => { body.appendChild(c); c.querySelector('a.zoom').onclick = e => { e.preventDefault(); window.openLightbox(c.querySelector('img').src); }; });
      shown += next.length; more.hidden = shown >= cur.length; nlab.textContent = shown + ' of ' + cur.length + ' shown';
      ctl.querySelectorAll('button').forEach(b => b.classList.toggle('on', b.dataset.f === filt || b.dataset.s === sort)); }
    ctl.querySelectorAll('button[data-f]').forEach(b => b.addEventListener('click', () => { filt = b.dataset.f; render(true); }));
    ctl.querySelectorAll('button[data-s]').forEach(b => b.addEventListener('click', () => { sort = b.dataset.s; render(true); }));
    more.addEventListener('click', () => render(false));
    render(true);
  }
  document.querySelectorAll('details.rec').forEach(det => det.addEventListener('toggle', () => { if (det.open) wire(det); }));
  if (location.hash) { const d = document.querySelector(location.hash); if (d && d.tagName === 'DETAILS') { d.open = true; wire(d); } }
  window.addEventListener('hashchange', () => { const d = document.querySelector(location.hash); if (d && d.tagName === 'DETAILS') { d.open = true; wire(d); d.scrollIntoView(); } });
  document.querySelectorAll('a.open-all').forEach(a => a.addEventListener('click', e => { e.preventDefault();
    document.querySelectorAll('details.rec').forEach(d => { d.open = true; wire(d); }); }));
})();
"""


def _nav(exp_rows: dict, pages: dict | None = None) -> str:
    h = ["<nav class='side'><b>QC report</b><a href='#top'>summary table</a>"]
    for exp, rows in exp_rows.items():
        h.append(f"<b>{_esc(exp)}</b>")
        for r in rows:
            href = r["href"] if pages is None else pages.get(exp, "") + "#" + r["href"].lstrip("#")
            h.append(f"<a href='{href}'>{_esc(r['short'])} · {r['numbers']['n_kept']} ev</a>")
    h.append("</nav>")
    return "".join(h)


def build_qc_report(mouse_dir: Path, out: Path | None = None, max_events: int = 0, clips_per_run: int = 6, no_clips: bool = False,
                    um_per_px: float = 0.792, split_mb: float = 80.0, log=print) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    ctx = report_context(mouse_dir)
    mouse_dir, A, day_labels, exps, rs = ctx["mouse_dir"], ctx["A"], ctx["day_labels"], ctx["exps"], ctx["rs"]
    out = Path(out) if out else A / f"{mouse_dir.name}_qc.html"
    t_start = time.time()
    # kept events, locomotion trials and per-CFU responsiveness per run from the analysis tables
    kept_by_run, trials_by_run, loco_cfu_by_run = {}, {}, {}
    for e in exps:
        t = _read_csv(A / f"events_{e['name']}.csv")
        if t is not None:
            for r in t.to_dict("records"):
                kept_by_run.setdefault(r["run"], {})[int(r["event"])] = r
        t = _read_csv(A / f"loco_trials_{e['name']}.csv")
        if t is not None:
            for r in t.to_dict("records"):
                trials_by_run.setdefault(r["run"], []).append((int(r["onset"]), str(r["kind"])))
        t = _read_csv(A / f"loco_cfu_{e['name']}.csv")
        if t is not None:
            for r in t.to_dict("records"):
                loco_cfu_by_run.setdefault(r["run"], {})[int(r["cfu"])] = r
    have_analysis = rs is not None
    log(f"QC report for {mouse_dir.name}: {len(exps)} experiment(s), {sum(len(e['runs']) for e in exps)} recording(s)")
    exp_rows: dict = {}
    sections: dict = {}
    n_rec = sum(len(e["runs"]) for e in exps)
    k_rec = 0
    for e in exps:
        exp_rows[e["name"]] = []
        secs = []
        for r in e["runs"]:
            tif = relocate(r["tif"], mouse_dir)
            k_rec += 1
            log(f"  QC: {_short_run(r['run'])} …  [{k_rec}/{n_rec}]")
            try:
                d = load_run(Path(tif), kept_by_run.get(r["run"]) if have_analysis else None, log=log)
                looks = _looks_for(d, um_per_px)
                trials = trials_by_run.get(r["run"], [])
                msumm, mdf = _motion_for(d, log=log)
                d["motion_df"] = mdf
                numbers = _run_numbers(d, trials, looks, motion=msumm)
                row = {"short": d["short"], "exp": e["name"], "cond": _cond_name(r.get("condition", ""), day_labels), "numbers": numbers,
                       "warnings": _warnings(numbers), "href": f"#run-{d['short']}"}
                exp_rows[e["name"]].append(row)
                secs.append(run_section(d, e["name"], r.get("condition", ""), day_labels, looks, trials, loco_cfu_by_run.get(r["run"], {}),
                                        numbers, 0 if no_clips else clips_per_run, max_events, log=log))
                try:
                    if hasattr(d.get("movie"), "close"):
                        d["movie"].close()
                except Exception:
                    pass
                d["movie"] = None; d["motion_df"] = None
            except Exception as ex:
                secs.append(f"<h3>{_esc(_short_run(r['run']))}</h3><p class='warn'>could not build the QC section: {_esc(ex)}</p>")
                log(f"    failed: {ex}")
        sections[e["name"]] = secs
    all_rows = [r for rows in exp_rows.values() for r in rows]
    intro = (f"<h1 id='top'>{_esc(mouse_dir.name)} — QC report: the detection on the real frames</h1>"
             f"<p class='byline'>Generated {time.strftime('%Y-%m-%d %H:%M')} by afferent2p (aqua_report.py --qc). The results are in the main report "
             f"<a href='{_esc(mouse_dir.name)}_report.html'>{_esc(mouse_dir.name)}_report.html</a>.</p>"
             "<p>One row per recording; click a column to sort, a recording to open its section (sections are collapsed and load when opened; "
             "<a class='open-all' href='#'>open all</a>). Per recording: the motion check (how far the raw frames moved, what AQuA2's own "
             "registration left behind, and whether the event onsets fall on moving frames more often than chance — <b>warn</b> when they do), "
             "the mean image with every detected event's box and the CFU outlines "
             "coloured by what they look like, the raster against running and pupil with the locomotion trials marked, a few clips, a montage of "
             "<b>every</b> detected event (frames before, at and after the onset with the detection box and the box's dF/F — a real afferent transient "
             "rises over one or two frames on a fibre-shaped footprint and decays; a noise detection is a flicker, a motion smear or a bright "
             "non-fibre structure), and the CFU table.</p>")
    per_exp_html = {exp: "".join(secs) for exp, secs in sections.items()}
    total_mb = sum(len(h) for h in per_exp_html.values()) / 1e6
    if total_mb > split_mb and len(exps) > 1:
        pages = {exp: f"{mouse_dir.name}_qc_{exp}.html" for exp in per_exp_html}
        for exp, html in per_exp_html.items():
            body = (_nav(exp_rows, pages) + f"<h1 id='top'>{_esc(mouse_dir.name)} — QC: {_esc(exp)}</h1>"
                    f"<p class='byline'><a href='{_esc(out.name)}'>back to the QC summary</a></p>" + summary_table(exp_rows[exp]) + html)
            (out.parent / pages[exp]).write_text(page_html(f"{mouse_dir.name} QC — {exp}", body, QC_CSS, QC_JS), encoding="utf-8")
        for rows in exp_rows.values():
            for r in rows:
                r["href"] = pages[r["exp"]] + r["href"]
        body = (_nav(exp_rows, {e: "" for e in pages}) + intro
                + f"<p class='small'>The per-recording sections are on one page per experiment ({total_mb:.0f} MB in all): "
                + ", ".join(f"<a href='{p}'>{_esc(exp)}</a>" for exp, p in pages.items()) + ".</p>" + summary_table(all_rows))
        out.write_text(page_html(f"{mouse_dir.name} — QC report", body, QC_CSS, QC_JS), encoding="utf-8")
        log(f"wrote {out} + {len(pages)} experiment page(s) ({total_mb:.0f} MB) in {(time.time() - t_start) / 60:.1f} min")
        return out
    body = _nav(exp_rows) + intro + summary_table(all_rows)
    for exp, html in per_exp_html.items():
        body += f"<h2 id='exp-{_esc(exp)}'>{_esc(exp)}</h2>" + html
    out.write_text(page_html(f"{mouse_dir.name} — QC report", body, QC_CSS, QC_JS), encoding="utf-8")
    log(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB) in {(time.time() - t_start) / 60:.1f} min")
    return out
