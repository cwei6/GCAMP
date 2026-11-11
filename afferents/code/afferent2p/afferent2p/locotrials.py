"""Locomotion-evoked activity — the readout of Yu Yong's lab-meeting pipeline
(slides 5-12 of LabMeeting_261012) applied to the afferent recordings.

The afferent activity we test is spontaneous and running is its main
trigger, so the effect of the stimulation (RearChr: post vs pre, acute) or of
the days (Treg / across days, chronic) is read on the locomotion-evoked
responses, not only on the overall event rate.

Trials (slide 5): one trial = 1 min around a running onset — a pre-bout window
of 5 frames (10 s) with no running, a running window of 15 frames (30 s) and a
post-bout window of 10 frames (20 s) with no running (frames are 1.984 s, so
Yu's 5 / 15 / 10 frames are the same seconds as ours). That is the STRICT set.
Afferent sessions often yield few strict trials, so a RELAXED set is scored
too: >= 10 s still before the onset, bout >= 4 s, nothing required after. Both
are reported with their N everywhere.

Responses per trial: the population signal (mean dF/F of the kept CFUs),
the event-count trace (kept events starting per frame) and every CFU's dF/F;
baseline = mean over the pre-bout window, peak = max in the running window.
"Responsive" trial, three rules side by side (the lab should pick one):
  yu   peak >= 1.5 x baseline on F/F0 (= 1 + dF/F), the slide's definition —
       with afferent dF/F baselines near 0 it needs a peak of +0.5 dF/F;
  z    peak - baseline >= 2 x the run's still-period noise (robust SD of the
       population dF/F over still frames) — the default here;
  ev   at least one kept event starts in the running window.

Also: per-CFU locomotion responsiveness against circularly shifted onsets
(L4), rest-only activity >= 20 s away from any running (L5), Loco2Fluor and
LocoPupil2Fluor ridge GLMs with held-out deviance explained and
leave-one-predictor-out (L6), and the behaviour check (L7).

Everything here works on afferent2p.runbundle dicts (recorded-frame grid).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .runbundle import population_dff, robust_sd

PRE_FRAMES, RUN_FRAMES, POST_FRAMES = 5, 15, 10          # Yu's windows (frames of ~2 s)
RUN_FRAC, POST_STILL_FRAC = 0.8, 0.8                      # strict: >= 80 % of the run window running, >= 80 % of the post window still
RELAXED_MIN_BOUT = 2                                      # frames (4 s)
YU_FACTOR, Z_FACTOR = 1.5, 2.0
REST_MARGIN_FRAMES = 10                                   # 20 s away from any running (L5)
RULES = ("z", "yu", "ev")
PALETTE = ["#4f5d75", "#d7263d", "#1b998b", "#f46036", "#7b2cbf", "#2a78d6"]


# ----------------------------------------------------------------------------- bouts and trials

def bouts(b: dict, merge_gap: int = 1, min_len: int = 2) -> list[dict]:
    """Running bouts from the per-frame running flag: gaps <= merge_gap frames
    inside a bout are bridged, bouts shorter than min_len frames dropped."""
    run = np.asarray(b["running"], bool).copy()
    T = run.size
    if merge_gap > 0 and run.any():
        d = np.diff(np.r_[0, run.astype(int), 0])
        starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        for s, e in zip(starts[1:], ends[:-1]):         # gap between consecutive bouts = [e, s)
            if s - e <= merge_gap:
                run[e:s] = True
    d = np.diff(np.r_[0, run.astype(int), 0])
    starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    out = []
    spf = b["spf"]
    speed = np.asarray(b["speed"], float)
    for i, (s, e) in enumerate(zip(starts, ends)):
        if e - s < min_len:
            continue
        prev_end = ends[i - 1] if i > 0 else 0
        next_start = starts[i + 1] if i + 1 < starts.size else T
        sp = speed[s:e]
        out.append({"start": int(s), "end": int(e), "duration_s": float((e - s) * spf),
                    "still_before_frames": int(s - prev_end), "still_after_frames": int(next_start - e),
                    "speed_mean": float(np.nanmean(sp)) if np.isfinite(sp).any() else np.nan,
                    "speed_peak": float(np.nanmax(sp)) if np.isfinite(sp).any() else np.nan})
    b["running_merged"] = run
    return out


def find_trials(b: dict, bts: list[dict], pre: int = PRE_FRAMES, run: int = RUN_FRAMES, post: int = POST_FRAMES,
                run_frac: float = RUN_FRAC, post_still_frac: float = POST_STILL_FRAC,
                relaxed_min_bout: int = RELAXED_MIN_BOUT) -> list[dict]:
    """Trials around bout onsets. kind = 'strict' (Yu) or 'relaxed'; a strict
    trial is also in the relaxed set (one row, kind='strict'). Trials whose
    windows contain an excluded frame, or do not fit in the run, are dropped
    and counted in b['trials_dropped']."""
    T = b["T"]
    running = np.asarray(b.get("running_merged", b["running"]), bool)
    excluded = np.asarray(b["excluded"], bool)
    out, dropped = [], 0
    for i, bt in enumerate(bts):
        t = bt["start"]
        if bt["still_before_frames"] < pre or t - pre < 0:
            continue
        if t + run > T:
            dropped += 1
            continue
        pre_w, run_w = slice(t - pre, t), slice(t, t + run)
        if running[pre_w].any():
            continue
        frac_run = float(running[run_w].mean())
        strict = (bt["end"] - bt["start"] >= relaxed_min_bout and frac_run >= run_frac and t + run + post <= T
                  and float((~running[t + run:t + run + post]).mean()) >= post_still_frac)
        relaxed = bt["end"] - bt["start"] >= relaxed_min_bout
        if not (strict or relaxed):
            continue
        last = t + run + post if strict else t + run
        if excluded[t - pre:last].any():
            dropped += 1
            continue
        out.append({"trial": len(out) + 1, "bout": i + 1, "onset": int(t), "onset_s": float(t * b["spf"]),
                    "kind": "strict" if strict else "relaxed", "bout_duration_s": bt["duration_s"],
                    "run_fraction": frac_run, "speed_mean_run": float(np.nanmean(b["speed"][run_w])) if np.isfinite(b["speed"][run_w]).any() else np.nan,
                    "speed_peak_run": float(np.nanmax(b["speed"][run_w])) if np.isfinite(b["speed"][run_w]).any() else np.nan})
    b["trials_dropped"] = dropped
    return out


# ----------------------------------------------------------------------------- responses

def still_noise(b: dict, pop: np.ndarray) -> float:
    """Robust SD of the population dF/F over still, non-excluded frames (the
    noise the z rule compares a peak with)."""
    ok = ~np.asarray(b["running"], bool) & ~np.asarray(b["excluded"], bool)
    sd = robust_sd(pop[ok])
    if not np.isfinite(sd) or sd <= 0:
        sd = robust_sd(pop)
    return sd


def _window_stats(trace: np.ndarray, t: int, pre: int = PRE_FRAMES, run: int = RUN_FRAMES):
    base = trace[t - pre:t]; rw = trace[t:t + run]
    b = float(np.nanmean(base)) if np.isfinite(base).any() else np.nan
    p = float(np.nanmax(rw)) if np.isfinite(rw).any() else np.nan
    return b, p


def responsive(base: float, peak: float, sd: float, n_events_run: int, rule: str) -> bool:
    if rule == "ev":
        return n_events_run >= 1
    if not (np.isfinite(base) and np.isfinite(peak)):
        return False
    if rule == "yu":
        return (1.0 + peak) >= YU_FACTOR * (1.0 + base)
    return np.isfinite(sd) and (peak - base) >= Z_FACTOR * sd


def trial_responses(b: dict, trials: list[dict], pre: int = PRE_FRAMES, run: int = RUN_FRAMES, post: int = POST_FRAMES) -> dict:
    """Per-trial rows + the onset-aligned traces. Returns {rows, pop_traces (n, pre+run+post),
    ev_traces, cfu_resp (n_trials, nCFU) peak-baseline, cfu_base, cfu_peak, sd_still, t_s}."""
    pop = population_dff(b)
    sd = still_noise(b, pop)
    onsets = np.asarray(b["onsets"], float)
    pupil = np.asarray(b["pupil"], float)
    L = pre + run + post
    rows, pop_tr, ev_tr, pup_tr = [], [], [], []
    n_cfu = b["n_cfu"]
    cfu_base = np.full((len(trials), n_cfu), np.nan); cfu_peak = np.full((len(trials), n_cfu), np.nan)
    cfu_ev = np.zeros((len(trials), n_cfu))
    for i, tr in enumerate(trials):
        t = tr["onset"]
        seg = np.full(L, np.nan); e_seg = np.full(L, np.nan); p_seg = np.full(L, np.nan)
        n = min(L, b["T"] - (t - pre))
        seg[:n] = pop[t - pre:t - pre + n]; e_seg[:n] = onsets[t - pre:t - pre + n]; p_seg[:n] = pupil[t - pre:t - pre + n]
        pop_tr.append(seg); ev_tr.append(e_seg); pup_tr.append(p_seg)
        base, peak = _window_stats(pop, t, pre, run)
        n_pre = int(onsets[t - pre:t].sum()); n_run = int(onsets[t:t + run].sum())
        n_post = int(onsets[t + run:t + run + post].sum()) if t + run + post <= b["T"] else np.nan
        pb, pp = _window_stats(pupil, t, pre, run)
        row = dict(tr)
        row.update({"run": b["run"], "condition": b["condition"], "pop_baseline": base, "pop_peak": peak,
                    "pop_response": peak - base if np.isfinite(base) and np.isfinite(peak) else np.nan,
                    "noise_sd_still": sd, "n_events_pre": n_pre, "n_events_run": n_run, "n_events_post": n_post,
                    "events_per_min_run": n_run / (run * b["spf"] / 60.0),
                    "pupil_baseline": pb, "pupil_peak": pp,
                    "pupil_dilation": pp - pb if np.isfinite(pb) and np.isfinite(pp) else np.nan})
        for rule in RULES:
            row[f"responsive_{rule}"] = int(responsive(base, peak, sd, n_run, rule))
        if n_cfu and b["cfu_dff"] is not None:
            for k in range(n_cfu):
                cb, cp = _window_stats(b["cfu_dff"][k], t, pre, run)
                cfu_base[i, k], cfu_peak[i, k] = cb, cp
                cfu_ev[i, k] = b["cfu_onsets"][k, t:t + run].sum() if b["cfu_onsets"] is not None else 0
            resp = cfu_peak[i] - cfu_base[i]
            row["n_cfu_responding"] = int(np.nansum(resp >= Z_FACTOR * np.nan_to_num(sd, nan=np.inf)))
            row["frac_cfu_responding"] = row["n_cfu_responding"] / n_cfu
        rows.append(row)
    t_s = (np.arange(L) - pre + 0.5) * b["spf"]
    return {"rows": rows, "pop_traces": np.array(pop_traces_or_empty(pop_tr, L)), "ev_traces": np.array(pop_traces_or_empty(ev_tr, L)),
            "pupil_traces": np.array(pop_traces_or_empty(pup_tr, L)), "cfu_base": cfu_base, "cfu_peak": cfu_peak,
            "cfu_events": cfu_ev, "sd_still": sd, "t_s": t_s, "pop": pop}


def pop_traces_or_empty(lst, L):
    return lst if lst else np.zeros((0, L))


# ----------------------------------------------------------------------------- per-CFU responsiveness (L4)

def cfu_responsiveness(b: dict, trials: list[dict], n_shift: int = 500, pre: int = PRE_FRAMES, run: int = RUN_FRAMES,
                       rng=None) -> list[dict]:
    """Each kept CFU's mean dF/F in the running windows minus the pre-bout
    windows (and its event probability per running window), against a null of
    the same windows at circularly shifted onsets. responsive = p < 0.05
    (one-sided, dF/F statistic)."""
    rng = rng or np.random.default_rng(0)
    n_cfu = b["n_cfu"]
    if not trials or n_cfu == 0 or b["cfu_dff"] is None:
        return []
    T = b["T"]
    dff = np.asarray(b["cfu_dff"], float)
    on = np.asarray(b["cfu_onsets"], float) if b["cfu_onsets"] is not None else np.zeros_like(dff)
    excluded = np.asarray(b["excluded"], bool)
    onsets = np.array([tr["onset"] for tr in trials], int)

    def stat(starts):
        d = np.zeros(n_cfu); e = np.zeros(n_cfu); n = 0
        for t in starts:
            if t - pre < 0 or t + run > T or excluded[t - pre:t + run].any():
                continue
            base = np.nanmean(dff[:, t - pre:t], axis=1); rw = np.nanmean(dff[:, t:t + run], axis=1)
            d += np.nan_to_num(rw - base); e += (on[:, t:t + run].sum(axis=1) > 0); n += 1
        if n == 0:
            return np.full(n_cfu, np.nan), np.full(n_cfu, np.nan), 0
        return d / n, e / n, n
    obs_d, obs_e, n_obs = stat(onsets)
    if n_obs == 0:
        return []
    null_d = np.full((n_shift, n_cfu), np.nan); null_e = np.full((n_shift, n_cfu), np.nan)
    for s in range(n_shift):
        shift = int(rng.integers(pre + run, T - run)) if T > 2 * run + pre else int(rng.integers(0, T))
        d, e, n = stat((onsets + shift) % T)
        if n:
            null_d[s], null_e[s] = d, e
    rows = []
    for k in range(n_cfu):
        nd = null_d[:, k]; nd = nd[np.isfinite(nd)]
        ne = null_e[:, k]; ne = ne[np.isfinite(ne)]
        p_d = float((np.sum(nd >= obs_d[k]) + 1) / (nd.size + 1)) if nd.size else np.nan
        p_e = float((np.sum(ne >= obs_e[k]) + 1) / (ne.size + 1)) if ne.size else np.nan
        rows.append({"run": b["run"], "condition": b["condition"], "cfu": int(b["cfu_ids"][k]) if len(b.get("cfu_ids", [])) > k else k + 1,
                     "n_trials": n_obs,
                     "dff_run_minus_pre": float(obs_d[k]), "p_dff": p_d,
                     "event_prob_run": float(obs_e[k]), "p_events": p_e,
                     "null_dff_p95": float(np.percentile(nd, 95)) if nd.size else np.nan,
                     "loco_responsive": int(np.isfinite(p_d) and p_d < 0.05),
                     "n_events_total": int(on[k].sum())})
    return rows


# ----------------------------------------------------------------------------- rest-only activity (L5)

def rest_only(b: dict, margin: int = REST_MARGIN_FRAMES) -> dict:
    """Events during still frames >= margin frames away from any running, before and after."""
    run = np.asarray(b.get("running_merged", b["running"]), bool)
    T = run.size
    near = run.copy()
    idx = np.flatnonzero(run)
    for m in range(1, margin + 1):
        near[np.clip(idx + m, 0, T - 1)] = True; near[np.clip(idx - m, 0, T - 1)] = True
    ok = ~near & ~np.asarray(b["excluded"], bool) & np.isfinite(b["still"])
    minutes = float(ok.sum()) * b["spf"] / 60.0
    on = np.asarray(b["onsets"], float)
    n_ev = int(on[ok].sum())
    ev = b["events"]
    dmax = np.nan
    if ev is not None and len(ev) and "t0_frame" in ev.columns and "dff_max" in ev.columns:
        t0 = np.clip(ev["t0_frame"].to_numpy(int) - 1, 0, T - 1)
        v = ev["dff_max"].to_numpy(float)[ok[t0]]
        v = v[np.isfinite(v)]
        dmax = float(np.median(v)) if v.size else np.nan
    return {"rest_minutes": minutes, "rest_n_events": n_ev, "rest_events_per_min": n_ev / minutes if minutes >= 0.25 else np.nan,
            "rest_dff_max_median": dmax, "rest_frames_fraction": float(ok.mean())}


# ----------------------------------------------------------------------------- GLMs (L6)

def lagged_design(preds: dict, lags: np.ndarray) -> tuple[np.ndarray, list]:
    """X[t, (p, l)] = pred_p[t - lag_l]; lag > 0 = the predictor's past."""
    cols, names = [], []
    for name, x in preds.items():
        x = np.asarray(x, float)
        T = x.size
        for lag in lags:
            c = np.full(T, np.nan)
            if lag >= 0:
                c[lag:] = x[:T - lag] if lag else x
            else:
                c[:T + lag] = x[-lag:]
            cols.append(c); names.append((name, int(lag)))
    return np.column_stack(cols) if cols else np.zeros((0, 0)), names


def _ridge(X, y, lam):
    n, p = X.shape
    A = X.T @ X + lam * np.eye(p)
    return np.linalg.solve(A, X.T @ y)


def ridge_glm(X: np.ndarray, y: np.ndarray, train_frac: float = 0.75, n_cycle: int = 50, n_fold: int = 5,
              lambdas=None, rng=None) -> dict:
    """Gaussian ridge regression with the penalty chosen by K-fold CV inside each
    random train/test split; returns the mean held-out deviance explained
    (1 - SS_res / SS_tot) over the cycles, the coefficients fitted on all
    rows (standardised predictors) and the chosen penalty."""
    rng = rng or np.random.default_rng(0)
    lambdas = np.logspace(-2, 4, 13) if lambdas is None else lambdas
    ok = np.isfinite(X).all(axis=1) & np.isfinite(y)
    X, y = X[ok], y[ok]
    n, p = X.shape
    if n < max(4 * p, 30) or np.nanstd(y) == 0:
        return {"dev_explained": np.nan, "coef": np.full(p, np.nan), "n": int(n), "lambda": np.nan, "why_not": f"{n} usable frames for {p} coefficients"}
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    ym = y.mean()
    ys = y - ym

    def cv_lambda(idx):
        folds = np.array_split(rng.permutation(idx), n_fold)
        err = np.zeros(lambdas.size)
        for f in folds:
            tr = np.setdiff1d(idx, f)
            if tr.size < p + 2 or f.size == 0:
                continue
            for i, lam in enumerate(lambdas):
                beta = _ridge(Xs[tr], ys[tr], lam * tr.size / n)
                err[i] += float(np.sum((ys[f] - Xs[f] @ beta) ** 2))
        return lambdas[int(np.argmin(err))]
    devs = []
    for _ in range(n_cycle):
        perm = rng.permutation(n)
        n_tr = int(round(train_frac * n))
        tr, te = perm[:n_tr], perm[n_tr:]
        if te.size < 5:
            continue
        lam = cv_lambda(tr)
        beta = _ridge(Xs[tr], ys[tr], lam * tr.size / n)
        pred = Xs[te] @ beta
        ss_res = float(np.sum((ys[te] - pred) ** 2)); ss_tot = float(np.sum((ys[te] - ys[te].mean()) ** 2))
        devs.append(1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan)
    lam = cv_lambda(np.arange(n))
    beta = _ridge(Xs, ys, lam)
    return {"dev_explained": float(np.nanmean(devs)) if devs else np.nan, "coef": beta / sd, "coef_std": beta,      # coef_std: per SD of the predictor
            "n": int(n), "lambda": float(lam)}


def glm_unit(preds: dict, y: np.ndarray, lags: np.ndarray, spf: float, rows_mask: np.ndarray | None = None,
             train_frac: float = 0.75, n_cycle: int = 50, n_fold: int = 5, min_dev: float = 0.2, lopo: bool = True,
             rng=None) -> dict:
    """One response (population or CFU dF/F) against lagged predictors: full fit,
    LOPO relative explanatory value per predictor, kernels."""
    X, names = lagged_design(preds, lags)
    if rows_mask is not None:
        X = X.copy(); X[~rows_mask] = np.nan
    full = ridge_glm(X, y, train_frac, n_cycle, n_fold, rng=rng)
    out = {"dev_explained": full["dev_explained"], "n_frames": full["n"], "well_fitted": int(np.isfinite(full["dev_explained"]) and full["dev_explained"] >= min_dev),
           "lambda": full.get("lambda"), "why_not": full.get("why_not", "")}
    kernels = {}
    coef = full.get("coef_std", full["coef"])           # per SD of the predictor, comparable between runs and predictors
    for p in preds:
        idx = [i for i, (nm, _) in enumerate(names) if nm == p]
        k = coef[idx]; lg = np.array([names[i][1] for i in idx]) * spf
        kernels[p] = (lg, k)
        if np.isfinite(k).any():
            j = int(np.nanargmax(np.abs(k)))
            out[f"peak_coef_{p}"] = float(k[j]); out[f"peak_lag_s_{p}"] = float(lg[j])
        else:
            out[f"peak_coef_{p}"] = np.nan; out[f"peak_lag_s_{p}"] = np.nan
    if lopo and len(preds) > 1 and np.isfinite(full["dev_explained"]):
        for p in preds:
            keep = [i for i, (nm, _) in enumerate(names) if nm != p]
            part = ridge_glm(X[:, keep], y, train_frac, n_cycle, n_fold, rng=rng)
            d_full, d_part = full["dev_explained"], part["dev_explained"]
            out[f"dev_without_{p}"] = d_part
            out[f"rel_value_{p}"] = float(np.clip((d_full - d_part) / d_full, 0.0, 1.0)) if (np.isfinite(d_part) and d_full > 0) else np.nan
    out["kernels"] = kernels
    return out


def glms_for_run(b: dict, trials: list[dict], pop: np.ndarray, min_events_cfu: int = 3, lags_loco_s=(-5.0, 10.0),
                 lags_lp_s=(-20.0, 20.0), log=print, rng=None) -> tuple[list[dict], list[dict]]:
    """Loco2Fluor on the trial frames (speed, lags -5..10 s; Ncycle 50, CV 5) and
    LocoPupil2Fluor on the whole trace (speed + pupil 0-1, lags -20..20 s;
    Ncycle 20, CV 10) for the population and every CFU with >= min_events_cfu
    events. Returns (rows, kernel rows)."""
    rng = rng or np.random.default_rng(0)
    spf = b["spf"]
    T = b["T"]
    speed = np.asarray(b["speed"], float).copy()
    speed[np.asarray(b["excluded"], bool)] = np.nan
    units = [("population", pop)]
    if b["n_cfu"] and b["cfu_dff"] is not None and b["cfu_onsets"] is not None:
        for k in range(b["n_cfu"]):
            if b["cfu_onsets"][k].sum() >= min_events_cfu:
                cid = int(b["cfu_ids"][k]) if len(b.get("cfu_ids", [])) > k else k + 1
                units.append((f"cfu{cid}", np.asarray(b["cfu_dff"][k], float)))
    rows, krows = [], []
    # trial frames for Loco2Fluor
    mask = np.zeros(T, bool)
    for tr in trials:
        mask[max(0, tr["onset"] - PRE_FRAMES):min(T, tr["onset"] + RUN_FRAMES + POST_FRAMES)] = True
    l1 = np.arange(int(np.floor(lags_loco_s[0] / spf)), int(np.ceil(lags_loco_s[1] / spf)) + 1)
    l2 = np.arange(int(np.floor(lags_lp_s[0] / spf)), int(np.ceil(lags_lp_s[1] / spf)) + 1)
    pupil = np.asarray(b["pupil"], float)
    pupil01 = None
    if np.isfinite(pupil).sum() > 20:
        lo, hi = np.nanpercentile(pupil, 1), np.nanpercentile(pupil, 99)
        pupil01 = np.clip((pupil - lo) / max(hi - lo, 1e-6), 0, 1) if hi > lo else None
    for name, y in units:
        y = np.asarray(y, float).copy(); y[np.asarray(b["excluded"], bool)] = np.nan
        if mask.sum() >= 30:
            g = glm_unit({"speed": speed}, y, l1, spf, rows_mask=mask, n_cycle=50, n_fold=5, lopo=False, rng=rng)
            rows.append({"run": b["run"], "condition": b["condition"], "unit": name, "model": "loco2fluor",
                         "n_frames": g["n_frames"], "dev_explained": g["dev_explained"], "well_fitted": g["well_fitted"],
                         "peak_coef_speed": g["peak_coef_speed"], "peak_lag_s_speed": g["peak_lag_s_speed"], "note": g["why_not"]})
            for p, (lg, k) in g["kernels"].items():
                krows += [{"run": b["run"], "condition": b["condition"], "unit": name, "model": "loco2fluor", "predictor": p,
                           "lag_s": float(a), "coef": float(c)} for a, c in zip(lg, k)]
        if pupil01 is not None:
            g = glm_unit({"speed": speed, "pupil": pupil01}, y, l2, spf, n_cycle=20, n_fold=10, lopo=True, rng=rng)
            rows.append({"run": b["run"], "condition": b["condition"], "unit": name, "model": "locopupil2fluor",
                         "n_frames": g["n_frames"], "dev_explained": g["dev_explained"], "well_fitted": g["well_fitted"],
                         "peak_coef_speed": g["peak_coef_speed"], "peak_lag_s_speed": g["peak_lag_s_speed"],
                         "peak_coef_pupil": g["peak_coef_pupil"], "peak_lag_s_pupil": g["peak_lag_s_pupil"],
                         "rel_value_speed": g.get("rel_value_speed", np.nan), "rel_value_pupil": g.get("rel_value_pupil", np.nan),
                         "dev_without_speed": g.get("dev_without_speed", np.nan), "dev_without_pupil": g.get("dev_without_pupil", np.nan),
                         "note": g["why_not"]})
            for p, (lg, k) in g["kernels"].items():
                krows += [{"run": b["run"], "condition": b["condition"], "unit": name, "model": "locopupil2fluor", "predictor": p,
                           "lag_s": float(a), "coef": float(c)} for a, c in zip(lg, k)]
    return rows, krows


# ----------------------------------------------------------------------------- behaviour check (L7)

def behaviour_check(b: dict, bts: list[dict], trials: list[dict]) -> dict:
    run = np.asarray(b.get("running_merged", b["running"]), bool)
    ok = ~np.asarray(b["excluded"], bool) & np.isfinite(b["still"])
    sp = np.asarray(b["speed"], float)
    strict = [t for t in trials if t["kind"] == "strict"]
    out = {"pct_running": 100.0 * float(run[ok].mean()) if ok.any() else np.nan,
           "n_bouts": len(bts), "bouts_per_min": len(bts) / max(b["minutes_used"], 1e-9) if b["minutes_used"] else np.nan,
           "bout_duration_s_mean": float(np.mean([x["duration_s"] for x in bts])) if bts else np.nan,
           "speed_running_mean": float(np.nanmean(sp[run & ok])) if (run & ok).any() else np.nan,
           "n_trials_strict": len(strict), "n_trials_relaxed": len(trials), "n_trials_dropped_excluded": int(b.get("trials_dropped", 0))}
    dil = [t.get("pupil_dilation", np.nan) for t in trials]
    dil = [d for d in dil if np.isfinite(d)]
    out["pupil_dilation_mean"] = float(np.mean(dil)) if dil else np.nan
    return out


# ----------------------------------------------------------------------------- one experiment

def _sem(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return float(np.std(x, ddof=1) / np.sqrt(x.size)) if x.size > 1 else np.nan


def _mwu(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return np.nan
    try:
        from scipy.stats import mannwhitneyu
        return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def analyze_experiment(exp: str, bundles: list[dict], conditions: list[str], out_dir: Path, glm: bool = True,
                       n_shift: int = 500, log=print) -> dict:
    """All locomotion readouts of one experiment. Writes loco_trials_<exp>.csv,
    loco_cfu_<exp>.csv, loco_glm_<exp>.csv (+ kernels), <exp>_loco.png,
    <exp>_loco_cfu.png. Returns {run_rows, cond_rows, lines, traces}."""
    import pandas as pd
    out_dir = Path(out_dir)
    trial_rows, cfu_rows, glm_rows, kern_rows, run_rows = [], [], [], [], []
    traces = {}                                     # condition -> dict(pop, ev, pupil, t_s, kinds)
    per_cond_cfu = {}                               # condition -> peri-onset CFU matrices (H4)
    for b in bundles:
        if not b["has_loco"]:
            run_rows.append({"run": b["run"], "experiment": exp, "condition": b["condition"], "note": "no locomotion data"})
            continue
        bts = bouts(b)
        trials = find_trials(b, bts)
        resp = trial_responses(b, trials)
        for r in resp["rows"]:
            r["experiment"] = exp
        trial_rows += resp["rows"]
        cr = cfu_responsiveness(b, trials, n_shift=n_shift)
        for r in cr:
            r["experiment"] = exp
        cfu_rows += cr
        rest = rest_only(b)
        chk = behaviour_check(b, bts, resp["rows"])
        row = {"run": b["run"], "experiment": exp, "condition": b["condition"], "minutes_used": b["minutes_used"],
               "n_cfu": b["n_cfu"], "noise_sd_still": resp["sd_still"], **chk, **rest}
        for kind in ("strict", "relaxed"):
            sel = [r for r in resp["rows"] if (r["kind"] == "strict") or kind == "relaxed"]
            row[f"n_trials_{kind}"] = len(sel)
            for rule in RULES:
                n_resp = int(sum(r[f"responsive_{rule}"] for r in sel))
                row[f"n_responsive_{rule}_{kind}"] = n_resp
                row[f"pct_responsive_{rule}_{kind}"] = 100.0 * n_resp / len(sel) if sel else np.nan
            pk = [r["pop_peak"] for r in sel if r["responsive_z"]]
            row[f"peak_responsive_mean_{kind}"] = float(np.nanmean(pk)) if pk else np.nan
            row[f"response_mean_{kind}"] = float(np.nanmean([r["pop_response"] for r in sel])) if sel else np.nan
            row[f"events_per_trial_{kind}"] = float(np.mean([r["n_events_run"] for r in sel])) if sel else np.nan
            ev_resp = [r["n_events_run"] for r in sel if r["responsive_z"]]
            row[f"events_per_responsive_trial_{kind}"] = float(np.mean(ev_resp)) if ev_resp else np.nan
        if cr:
            row["n_cfu_loco_responsive"] = int(sum(r["loco_responsive"] for r in cr))
            row["pct_cfu_loco_responsive"] = 100.0 * row["n_cfu_loco_responsive"] / len(cr)
        if glm:
            try:
                gr, kr = glms_for_run(b, trials, resp["pop"], log=log)
                for r in gr:
                    r["experiment"] = exp
                glm_rows += gr; kern_rows += kr
                popg = {r["model"]: r for r in gr if r["unit"] == "population"}
                for model, short in (("loco2fluor", "l2f"), ("locopupil2fluor", "lp2f")):
                    g = popg.get(model)
                    if g:
                        row[f"dev_{short}_population"] = g["dev_explained"]
                        row[f"peak_coef_speed_{short}"] = g.get("peak_coef_speed")
                        if model == "locopupil2fluor":
                            row["rel_value_speed"] = g.get("rel_value_speed"); row["rel_value_pupil"] = g.get("rel_value_pupil")
                fitted = [r for r in gr if r["unit"] != "population" and r["model"] == "loco2fluor"]
                row["n_cfu_glm"] = len(fitted)
                row["n_cfu_well_fitted"] = int(sum(r["well_fitted"] for r in fitted))
                dv = [r["dev_explained"] for r in fitted if r["well_fitted"]]
                row["dev_l2f_cfu_median"] = float(np.median(dv)) if dv else np.nan
            except Exception as e:
                log(f"    GLM failed for {b['run']} ({e})")
        run_rows.append(row)
        tr = traces.setdefault(b["condition"], {"pop": [], "ev": [], "pupil": [], "kinds": [], "t_s": resp["t_s"], "runs": []})
        tr["pop"].append(resp["pop_traces"]); tr["ev"].append(resp["ev_traces"]); tr["pupil"].append(resp["pupil_traces"])
        tr["kinds"] += [r["kind"] for r in resp["rows"]]; tr["runs"] += [b["run"]] * len(resp["rows"])
        if trials and b["n_cfu"] and b["cfu_dff"] is not None:
            L = PRE_FRAMES + RUN_FRAMES + POST_FRAMES
            M = np.full((len(trials), b["n_cfu"], L), np.nan)
            for i, t in enumerate(trials):
                t0 = t["onset"] - PRE_FRAMES; n = min(L, b["T"] - t0)
                M[i, :, :n] = b["cfu_dff"][:, t0:t0 + n]
            per_cond_cfu.setdefault(b["condition"], []).append(np.nanmean(M, axis=0))
    for c, tr in traces.items():
        for k in ("pop", "ev", "pupil"):
            tr[k] = np.vstack([x for x in tr[k] if x.size]) if any(x.size for x in tr[k]) else np.zeros((0, tr["t_s"].size))
    # ---- per condition (trials pooled over the runs of the condition)
    cond_rows = []
    tdf = pd.DataFrame(trial_rows) if trial_rows else pd.DataFrame()
    rdf = pd.DataFrame(run_rows)
    for c in conditions:
        sub = rdf[rdf["condition"] == c] if len(rdf) else rdf
        row = {"experiment": exp, "condition": c, "n_runs": int(len(sub))}
        for k in ("pct_running", "n_bouts", "bout_duration_s_mean", "speed_running_mean", "pupil_dilation_mean",
                  "rest_events_per_min", "rest_dff_max_median", "dev_l2f_population", "dev_lp2f_population",
                  "rel_value_speed", "rel_value_pupil", "pct_cfu_loco_responsive", "dev_l2f_cfu_median"):
            if k in sub.columns:
                v = sub[k].to_numpy(float); v = v[np.isfinite(v)]
                row[f"{k}_mean"] = float(v.mean()) if v.size else np.nan
                row[f"{k}_sem"] = _sem(v)
        if "rest_n_events" in sub.columns:
            row["rest_n_events"] = int(sub["rest_n_events"].sum()); row["rest_minutes"] = float(sub["rest_minutes"].sum())
            row["rest_events_per_min"] = row["rest_n_events"] / row["rest_minutes"] if row["rest_minutes"] >= 0.5 else np.nan
        if "n_cfu_loco_responsive" in sub.columns:
            row["n_cfu_loco_responsive"] = int(sub["n_cfu_loco_responsive"].fillna(0).sum())
            row["n_cfu_tested"] = int(sum(1 for r in cfu_rows if r["condition"] == c))
            row["pct_cfu_loco_responsive"] = 100.0 * row["n_cfu_loco_responsive"] / row["n_cfu_tested"] if row["n_cfu_tested"] else np.nan
        for kind in ("strict", "relaxed"):
            t = tdf[(tdf["condition"] == c) & ((tdf["kind"] == "strict") | (kind == "relaxed"))] if len(tdf) else tdf
            row[f"n_trials_{kind}"] = int(len(t))
            for rule in RULES:
                n_resp = int(t[f"responsive_{rule}"].sum()) if len(t) else 0
                row[f"n_responsive_{rule}_{kind}"] = n_resp
                row[f"pct_responsive_{rule}_{kind}"] = 100.0 * n_resp / len(t) if len(t) else np.nan
            if len(t):
                rz = t[t["responsive_z"] == 1]
                row[f"peak_responsive_mean_{kind}"] = float(rz["pop_peak"].mean()) if len(rz) else np.nan
                row[f"peak_responsive_sem_{kind}"] = _sem(rz["pop_peak"]) if len(rz) else np.nan
                row[f"response_mean_{kind}"] = float(t["pop_response"].mean()); row[f"response_sem_{kind}"] = _sem(t["pop_response"])
                row[f"events_per_trial_{kind}"] = float(t["n_events_run"].mean()); row[f"events_per_trial_sem_{kind}"] = _sem(t["n_events_run"])
                row[f"events_per_responsive_trial_{kind}"] = float(rz["n_events_run"].mean()) if len(rz) else np.nan
                row[f"pupil_dilation_{kind}"] = float(t["pupil_dilation"].mean()) if t["pupil_dilation"].notna().any() else np.nan
            else:
                for k in ("peak_responsive_mean", "peak_responsive_sem", "response_mean", "response_sem", "events_per_trial",
                          "events_per_trial_sem", "events_per_responsive_trial", "pupil_dilation"):
                    row[f"{k}_{kind}"] = np.nan
        cond_rows.append(row)
    # ---- tests between the baseline and each condition (trials as the unit, within one animal: descriptive)
    base = conditions[0]
    lines = [f"  locomotion-evoked activity ({exp}):"]
    for c in conditions:
        r = next((x for x in cond_rows if x["condition"] == c), None)
        if r is None:
            continue
        lines.append(f"    {c:>8}: {r['n_trials_strict']} strict / {r['n_trials_relaxed']} relaxed trials in {r['n_runs']} run(s); "
                     f"responsive (z rule) {r['pct_responsive_z_strict']:.0f} % strict, {r['pct_responsive_z_relaxed']:.0f} % relaxed "
                     f"(yu rule {r['pct_responsive_yu_relaxed']:.0f} %, >= 1 event {r['pct_responsive_ev_relaxed']:.0f} %); "
                     f"response {r['response_mean_relaxed']:.3f} dF/F, {r['events_per_trial_relaxed']:.2f} events per trial; "
                     f"rest-only {r['rest_events_per_min']:.3f}/min ({r.get('rest_n_events', 0)} events in {r.get('rest_minutes', 0):.1f} min)"
                     + (f"; {r['pct_cfu_loco_responsive']:.0f} % of {r['n_cfu_tested']} CFUs locomotion-responsive" if r.get("n_cfu_tested") else ""))
    for c in conditions[1:]:
        if not len(tdf):
            break
        tb = tdf[tdf["condition"] == base]; tc = tdf[tdf["condition"] == c]
        p_resp = _mwu(tb["pop_response"], tc["pop_response"]); p_ev = _mwu(tb["n_events_run"], tc["n_events_run"])
        rb = next((x for x in cond_rows if x["condition"] == base), {}); rc = next((x for x in cond_rows if x["condition"] == c), {})
        def ratio(k):
            a, bb = rc.get(k, np.nan), rb.get(k, np.nan)
            return a / bb if (np.isfinite(a) and np.isfinite(bb) and bb > 0) else np.nan
        lines.append(f"    {c} vs {base} (relaxed trials, Mann-Whitney, trials as the unit — descriptive): response p = {p_resp:.2g}, "
                     f"events per trial p = {p_ev:.2g}; acute ratios: responsive-trial peak x{ratio('peak_responsive_mean_relaxed'):.2f}, "
                     f"events per trial x{ratio('events_per_trial_relaxed'):.2f}, rest-only rate x{ratio('rest_events_per_min'):.2f}")
    n_all = sum(r.get("n_trials_relaxed", 0) for r in cond_rows)
    if n_all == 0:
        lines.append("    no locomotion trials: the animal never stood still >= 10 s before a running bout (or no behaviour data)")
    # ---- csvs
    if trial_rows:
        pd.DataFrame(trial_rows).to_csv(out_dir / f"loco_trials_{exp}.csv", index=False)
    if cfu_rows:
        pd.DataFrame(cfu_rows).to_csv(out_dir / f"loco_cfu_{exp}.csv", index=False)
    if glm_rows:
        pd.DataFrame(glm_rows).to_csv(out_dir / f"loco_glm_{exp}.csv", index=False)
        pd.DataFrame(kern_rows).to_csv(out_dir / f"loco_glm_kernels_{exp}.csv", index=False)
    # ---- figures
    try:
        loco_figure(exp, conditions, traces, tdf, rdf, cond_rows, glm_rows, kern_rows, out_dir / f"{exp}_loco.png")
    except Exception as e:
        log(f"    loco figure failed ({e})")
    try:
        if per_cond_cfu:
            peri_onset_heatmap(exp, conditions, per_cond_cfu, bundles[0]["spf"], out_dir / f"{exp}_loco_cfu.png")
    except Exception as e:
        log(f"    peri-onset heatmap failed ({e})")
    return {"run_rows": run_rows, "cond_rows": cond_rows, "lines": lines, "traces": traces, "n_trials": n_all,
            "trial_rows": trial_rows, "cfu_rows": cfu_rows, "glm_rows": glm_rows}


# ----------------------------------------------------------------------------- figures

def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8)


def _bar_dots(ax, groups: list, labels: list, colors: list, ylabel: str = "", title: str = "", fmt=None):
    """Bars = mean ± sem, dots = the values (runs / trials)."""
    rng = np.random.default_rng(1)
    for i, (g, lab, col) in enumerate(zip(groups, labels, colors)):
        g = np.asarray(g, float); g = g[np.isfinite(g)]
        if g.size:
            ax.bar(i, g.mean(), color=col, alpha=0.35, width=0.7)
            if g.size > 1:
                ax.errorbar(i, g.mean(), yerr=_sem(g), color=col, capsize=3, lw=1.2)
            ax.scatter(i + rng.uniform(-0.18, 0.18, g.size), g, s=14, color=col, zorder=3, edgecolor="white", lw=0.4)
        else:
            ax.text(i, 0, "n/a", ha="center", va="bottom", fontsize=7, color="#888")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8); ax.set_title(title, fontsize=9)
    _style(ax)


def loco_figure(exp: str, conditions: list, traces: dict, tdf, rdf, cond_rows: list, glm_rows: list, kern_rows: list, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(conditions)}
    fig, axes = plt.subplots(3, 4, figsize=(17, 11.5))
    # row 1: onset-aligned traces, pies, trials per run
    ax = axes[0, 0]
    kind_used = "strict" if all(sum(1 for k in traces.get(c, {}).get("kinds", []) if k == "strict") >= 3 for c in conditions if c in traces) and traces else "relaxed"
    for c in conditions:
        tr = traces.get(c)
        if not tr or tr["pop"].shape[0] == 0:
            continue
        sel = np.array([k == "strict" for k in tr["kinds"]]) if kind_used == "strict" else np.ones(len(tr["kinds"]), bool)
        M = tr["pop"][sel]
        if M.shape[0] == 0:
            continue
        m = np.nanmean(M, axis=0); s = np.nanstd(M, axis=0, ddof=1) / np.sqrt(max(M.shape[0], 1)) if M.shape[0] > 1 else np.zeros_like(m)
        ax.plot(tr["t_s"], m, color=cols[c], label=f"{c} (n={M.shape[0]})")
        ax.fill_between(tr["t_s"], m - s, m + s, color=cols[c], alpha=0.2, lw=0)
    ax.axvspan(0, RUN_FRAMES * (traces[conditions[0]]["t_s"][1] - traces[conditions[0]]["t_s"][0]) if traces.get(conditions[0]) else 30, color="#f2c94c", alpha=0.15, lw=0)
    ax.axvline(0, color="#999", lw=0.8, ls=":")
    ax.set_xlabel("time from running onset (s)", fontsize=8); ax.set_ylabel("population dF/F (mean of kept CFUs)", fontsize=8)
    ax.set_title(f"onset-aligned population dF/F, mean ± sem ({kind_used} trials)", fontsize=9); ax.legend(fontsize=7, frameon=False); _style(ax)
    ax = axes[0, 1]
    for c in conditions:
        tr = traces.get(c)
        if not tr or tr["ev"].shape[0] == 0:
            continue
        spf = tr["t_s"][1] - tr["t_s"][0]
        m = np.nanmean(tr["ev"], axis=0) / (spf / 60.0)
        ax.plot(tr["t_s"], m, color=cols[c], label=f"{c} (n={tr['ev'].shape[0]})", drawstyle="steps-mid")
    ax.axvline(0, color="#999", lw=0.8, ls=":")
    ax.set_xlabel("time from running onset (s)", fontsize=8); ax.set_ylabel("kept events / min", fontsize=8)
    ax.set_title("event-count trace (relaxed trials)", fontsize=9); ax.legend(fontsize=7, frameon=False); _style(ax)
    ax = axes[0, 2]
    ax.axis("off")
    n = len(conditions)
    for i, c in enumerate(conditions):
        r = next((x for x in cond_rows if x["condition"] == c), None)
        if not r or not r.get("n_trials_relaxed"):
            continue
        sub = ax.inset_axes([i / n, 0.1, 1 / n, 0.8])
        pct = r["pct_responsive_z_relaxed"]
        sub.pie([pct, 100 - pct], colors=[cols[c], "#e6e6e6"], startangle=90, counterclock=False,
                wedgeprops={"lw": 0.5, "edgecolor": "white"})
        sub.set_title(f"{c}\n{pct:.0f} % of {r['n_trials_relaxed']}", fontsize=8)
    ax.set_title("responsive trials (z rule, relaxed)", fontsize=9)
    ax = axes[0, 3]
    if len(rdf) and "n_trials_relaxed" in rdf.columns:
        x = 0; ticks = []; labs = []
        for c in conditions:
            sub = rdf[rdf["condition"] == c]
            xs = np.arange(len(sub)) + x
            ax.bar(xs, sub["n_trials_relaxed"], color="#dddddd", width=0.8, label="relaxed" if c == conditions[0] else None)
            ax.bar(xs, sub["n_trials_strict"], color=cols[c], width=0.8, alpha=0.6, label="strict" if c == conditions[0] else None)
            ax.scatter(xs, sub["n_responsive_z_relaxed"], color="#111", s=16, zorder=3, label="responsive (z)" if c == conditions[0] else None)
            ticks.append(x + (len(sub) - 1) / 2); labs.append(c); x += len(sub) + 1
        ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=8); ax.legend(fontsize=7, frameon=False)
    ax.set_ylabel("trials per run", fontsize=8); ax.set_title("N trials per recording (bar) and responsive (dot)", fontsize=9); _style(ax)
    # row 2: per-trial amplitude and events, per-CFU responsiveness, rest-only
    if len(tdf):
        _bar_dots(axes[1, 0], [tdf.loc[(tdf["condition"] == c) & (tdf["responsive_z"] == 1), "pop_peak"] for c in conditions], conditions,
                  [cols[c] for c in conditions], "peak population dF/F", "responsive-trial peak (dot = trial)")
        _bar_dots(axes[1, 1], [tdf.loc[tdf["condition"] == c, "n_events_run"] for c in conditions], conditions,
                  [cols[c] for c in conditions], "kept events in the running window", "events per trial (dot = trial)")
    else:
        axes[1, 0].axis("off"); axes[1, 1].axis("off")
    ax = axes[1, 2]
    vals = [next((x.get("pct_cfu_loco_responsive", np.nan) for x in cond_rows if x["condition"] == c), np.nan) for c in conditions]
    ns = [next((x.get("n_cfu_tested", 0) for x in cond_rows if x["condition"] == c), 0) for c in conditions]
    ax.bar(range(len(conditions)), [v if np.isfinite(v) else 0 for v in vals], color=[cols[c] for c in conditions], alpha=0.7)
    top = max([v for v in vals if np.isfinite(v)] + [5.0])
    ax.set_ylim(0, top * 1.3)
    for i, (v, m) in enumerate(zip(vals, ns)):
        ax.text(i, (v if np.isfinite(v) else 0) + top * 0.04, f"{m} CFUs", ha="center", fontsize=7)
    ax.set_xticks(range(len(conditions))); ax.set_xticklabels(conditions, fontsize=8); ax.set_ylabel("% CFUs", fontsize=8)
    ax.set_title("locomotion-responsive CFUs (vs shifted onsets, p < 0.05)", fontsize=9); _style(ax)
    if len(rdf) and "rest_events_per_min" in rdf.columns:
        _bar_dots(axes[1, 3], [rdf.loc[rdf["condition"] == c, "rest_events_per_min"] for c in conditions], conditions,
                  [cols[c] for c in conditions], "events / min at rest", "rest-only activity (>= 20 s from running; dot = run)")
    # row 3: behaviour check, pupil dilation, GLM
    if len(rdf) and "pct_running" in rdf.columns:
        _bar_dots(axes[2, 0], [rdf.loc[rdf["condition"] == c, "pct_running"] for c in conditions], conditions,
                  [cols[c] for c in conditions], "% of frames running", "behaviour check: running (dot = run)")
        _bar_dots(axes[2, 1], [rdf.loc[rdf["condition"] == c, "n_bouts"] for c in conditions], conditions,
                  [cols[c] for c in conditions], "running bouts per recording", "behaviour check: bouts (dot = run)")
    if len(tdf) and tdf["pupil_dilation"].notna().any():
        _bar_dots(axes[2, 2], [tdf.loc[tdf["condition"] == c, "pupil_dilation"] for c in conditions], conditions,
                  [cols[c] for c in conditions], "pupil peak − pre-bout (px)", "locomotion-related pupil dilation (dot = trial)")
    else:
        axes[2, 2].axis("off")
    ax = axes[2, 3]
    if glm_rows:
        import pandas as pd
        g = pd.DataFrame(glm_rows)
        pop = g[(g["unit"] == "population")]
        labels, groups, colors = [], [], []
        for model, short in (("loco2fluor", "speed→dF/F"), ("locopupil2fluor", "speed+pupil→dF/F")):
            for c in conditions:
                v = pop.loc[(pop["model"] == model) & (pop["condition"] == c), "dev_explained"]
                labels.append(f"{short}\n{c}"); groups.append(v); colors.append(cols[c])
        _bar_dots(ax, groups, labels, colors, "held-out deviance explained", "GLM on the population dF/F (dot = run)")
        ax.axhline(0.2, color="#999", lw=0.8, ls=":")
        ax.tick_params(axis="x", labelsize=6)
    else:
        ax.axis("off")
    fig.suptitle(f"{exp}: locomotion-evoked activity (Yu's trial pipeline; dots = trials / runs within one animal)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=105)
    plt.close(fig)


def peri_onset_heatmap(exp: str, conditions: list, per_cond: dict, spf: float, out_png: Path):
    """H4: CFU x time from the running onset, trial-averaged per run and
    normalised 0-1 per CFU over the whole figure, one panel per condition
    (rows = that condition's CFUs of each run, sorted by response)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    conds = [c for c in conditions if c in per_cond]
    if not conds:
        return
    mats = {c: np.vstack(per_cond[c]) for c in conds}
    L = next(iter(mats.values())).shape[1]
    t = (np.arange(L) - PRE_FRAMES + 0.5) * spf
    fig, axes = plt.subplots(1, len(conds), figsize=(4.2 * len(conds) + 1, 5.5), squeeze=False)
    for ax, c in zip(axes[0], conds):
        M = mats[c].copy()
        lo = np.nanmin(M, axis=1, keepdims=True); hi = np.nanmax(M, axis=1, keepdims=True)
        N = (M - lo) / np.maximum(hi - lo, 1e-9)
        resp = np.nanmean(N[:, PRE_FRAMES:PRE_FRAMES + RUN_FRAMES], axis=1) - np.nanmean(N[:, :PRE_FRAMES], axis=1)
        order = np.argsort(-np.nan_to_num(resp, nan=-9))
        im = ax.imshow(N[order], aspect="auto", cmap="viridis", vmin=0, vmax=1, extent=[t[0], t[-1], N.shape[0], 0], interpolation="nearest")
        ax.axvline(0, color="white", lw=0.9, ls="--"); ax.axvline(RUN_FRAMES * spf, color="white", lw=0.6, ls=":")
        ax.set_title(f"{c}: {N.shape[0]} CFU rows (per run)", fontsize=9); ax.set_xlabel("time from running onset (s)", fontsize=8)
        ax.set_ylabel("CFU (sorted by response)", fontsize=8); ax.tick_params(labelsize=7)
    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.6, label="dF/F, 0-1 per CFU")
    fig.suptitle(f"{exp}: peri-locomotion-onset CFU activity (trial-averaged)", fontsize=10)
    fig.savefig(out_png, dpi=105, bbox_inches="tight")
    plt.close(fig)
