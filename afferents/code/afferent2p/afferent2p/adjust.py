"""Behaviour-adjusted comparison of event rates.

The animal is not in the same state in the two conditions (it may run through
all of 'post' and sit through all of 'pre'), and afferent events follow the
running state. Matching (still frames only, pupil terciles) throws away most of
the data when the states hardly overlap; a regression keeps all of it:

    events in bin ~ Poisson( minutes_in_bin * exp(b0 + b_cond * [condition]
                             + b_speed * speed_z + b_pupil * pupil_z + b_run * running) )

on bins of `bin_s` seconds (excluded frames leave the exposure), so
exp(b_cond) is the condition's rate ratio AT THE SAME running speed, running
state and pupil size. Its confidence interval comes from a block bootstrap
(blocks of `block_s` seconds resampled within each run: neighbouring bins are
not independent). The share of the rate variation that behaviour explains
(deviance explained, 0..1) says how locomotion-driven the afferents are, and
the locomotion-onset-triggered rate shows it directly.
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------------- Poisson regression (IRLS)

def poisson_fit(X: np.ndarray, y: np.ndarray, offset: np.ndarray, max_iter: int = 50, tol: float = 1e-8,
                ridge: float = 0.05):
    """Maximum-likelihood Poisson regression with a log link and an offset,
    with a small ridge penalty on the non-intercept coefficients (keeps the
    fit finite when a condition has no events in one behavioural state; it is
    negligible with hundreds of bins). Returns (beta, deviance, null
    deviance). X includes the intercept column."""
    X = np.asarray(X, float); y = np.asarray(y, float); off = np.asarray(offset, float)
    n, p = X.shape
    beta = np.zeros(p)
    beta[0] = np.log(max(y.sum(), 0.5) / max(np.exp(off).sum(), 1e-9))
    pen = np.full(p, ridge); pen[0] = 0.0
    for _ in range(max_iter):
        eta = X @ beta + off
        mu = np.exp(np.clip(eta, -30, 30))
        z = eta - off + (y - mu) / np.maximum(mu, 1e-12)
        Wd = mu
        XtW = X.T * Wd
        try:
            new = np.linalg.solve(XtW @ X + np.diag(pen) + 1e-9 * np.eye(p), XtW @ z)
        except np.linalg.LinAlgError:
            break
        new = np.clip(new, -12, 12)
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new
    eta = X @ beta + off
    mu = np.exp(np.clip(eta, -30, 30))
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = 2 * np.sum(np.where(y > 0, y * np.log(y / mu), 0.0) - (y - mu))
        mu0 = np.exp(off) * y.sum() / max(np.exp(off).sum(), 1e-9)
        dev0 = 2 * np.sum(np.where(y > 0, y * np.log(y / mu0), 0.0) - (y - mu0))
    return beta, float(dev), float(dev0)


# ----------------------------------------------------------------------------- bins

def bin_run(onsets: np.ndarray, excluded: np.ndarray, speed: np.ndarray, still: np.ndarray, pupil: np.ndarray,
            spf: float, bin_s: float = 20.0) -> dict:
    """Per-bin counts and covariates for one run (recorded-frame grid)."""
    T = onsets.size
    nb = max(1, int(round(bin_s / spf)))
    n_bins = int(np.ceil(T / nb))
    rows = {k: np.full(n_bins, np.nan) for k in ("events", "minutes", "speed", "running", "pupil")}
    for b in range(n_bins):
        sl = slice(b * nb, min(T, (b + 1) * nb))
        ok = ~excluded[sl]
        rows["minutes"][b] = ok.sum() * spf / 60.0
        rows["events"][b] = onsets[sl][ok].sum()
        if ok.any():
            rows["speed"][b] = np.nanmean(speed[sl][ok]) if np.isfinite(speed[sl][ok]).any() else np.nan
            rows["running"][b] = np.nanmean((still[sl][ok] < 0.5).astype(float)) if np.isfinite(still[sl][ok]).any() else np.nan
            rows["pupil"][b] = np.nanmean(pupil[sl][ok]) if np.isfinite(pupil[sl][ok]).any() else np.nan
    return rows


def adjusted_ratio(runs: list, conditions: list, spf: float, bin_s: float = 20.0, block_s: float = 120.0,
                   n_boot: int = 300, rng=None) -> dict:
    """runs: [{condition, onsets (events per recorded frame), excluded (bool),
    speed, still, pupil}] of one experiment. Returns the behaviour-adjusted
    rate ratios of every condition vs the first, with bootstrap CIs, plus the
    raw ratios on the same bins and the deviance explained by behaviour."""
    rng = rng or np.random.default_rng(0)
    base = conditions[0]
    tabs = []
    for ri, r in enumerate(runs):
        t = bin_run(r["onsets"], r["excluded"], r["speed"], r["still"], r["pupil"], spf, bin_s)
        n = t["events"].size
        t["cond"] = np.full(n, conditions.index(r["condition"]) if r["condition"] in conditions else -1)
        t["run"] = np.full(n, ri)
        t["block"] = np.arange(n) // max(1, int(round(block_s / bin_s)))
        tabs.append(t)
    keys = list(tabs[0])
    d = {k: np.concatenate([t[k] for t in tabs]) for k in keys}
    use = (d["minutes"] > 0.25 * bin_s / 60.0) & (d["cond"] >= 0)
    have_speed = np.isfinite(d["speed"][use]).mean() > 0.8
    have_pupil = np.isfinite(d["pupil"][use]).mean() > 0.8
    for k, have in (("speed", have_speed), ("running", have_speed), ("pupil", have_pupil)):
        if have:
            use &= np.isfinite(d[k])
    idx = np.flatnonzero(use)
    out = {"n_bins": int(idx.size), "bin_s": bin_s, "covariates": []}
    n_ev_cond = [float(d["events"][idx][d["cond"][idx] == ci].sum()) for ci in range(len(conditions))]
    out["events_per_condition"] = n_ev_cond
    if idx.size < 30 or min(n_ev_cond) < 5:
        out["why_not"] = (f"only {idx.size} bins" if idx.size < 30 else
                          f"fewer than 5 events in a condition ({', '.join(f'{c}: {int(n)}' for c, n in zip(conditions, n_ev_cond))})")
        return out

    def design(sel):
        cols = [np.ones(sel.size)]
        names = ["intercept"]
        for ci, c in enumerate(conditions[1:], 1):
            cols.append((d["cond"][sel] == ci).astype(float)); names.append(f"cond:{c}")
        if have_speed:
            sp = np.log1p(np.maximum(d["speed"][sel], 0))
            cols.append((sp - sp_mean) / sp_sd); names.append("speed_z")
            cols.append(d["running"][sel]); names.append("running")
        if have_pupil:
            cols.append((d["pupil"][sel] - pu_mean) / pu_sd); names.append("pupil_z")
        return np.column_stack(cols), names

    sp_all = np.log1p(np.maximum(d["speed"][idx], 0)) if have_speed else np.zeros(1)
    sp_mean, sp_sd = (float(np.nanmean(sp_all)), float(np.nanstd(sp_all) + 1e-9)) if have_speed else (0.0, 1.0)
    pu_mean, pu_sd = (float(np.nanmean(d["pupil"][idx])), float(np.nanstd(d["pupil"][idx]) + 1e-9)) if have_pupil else (0.0, 1.0)
    X, names = design(idx)
    y = d["events"][idx]; off = np.log(d["minutes"][idx])
    beta, dev, dev0 = poisson_fit(X, y, off)
    # behaviour-only and condition-only models for the deviance shares
    beh_cols = [i for i, nm in enumerate(names) if nm in ("speed_z", "running", "pupil_z")]
    cond_cols = [i for i, nm in enumerate(names) if nm.startswith("cond:")]
    _, dev_beh, _ = poisson_fit(X[:, [0] + beh_cols], y, off) if beh_cols else (None, dev0, dev0)
    _, dev_cond, _ = poisson_fit(X[:, [0] + cond_cols], y, off) if cond_cols else (None, dev0, dev0)
    out["covariates"] = names
    out["coef"] = {nm: float(b) for nm, b in zip(names, beta)}
    out["deviance_explained_behaviour"] = float(1 - dev_beh / dev0) if dev0 > 0 else np.nan
    out["deviance_explained_condition"] = float(1 - dev_cond / dev0) if dev0 > 0 else np.nan
    out["deviance_explained_full"] = float(1 - dev / dev0) if dev0 > 0 else np.nan
    # raw ratios on the same bins
    for ci, c in enumerate(conditions[1:], 1):
        mb = d["cond"][idx] == 0; mc = d["cond"][idx] == ci
        rb = y[mb].sum() / max(d["minutes"][idx][mb].sum(), 1e-9)
        rc = y[mc].sum() / max(d["minutes"][idx][mc].sum(), 1e-9)
        out[f"raw_ratio_{c}"] = float(rc / rb) if rb > 0 else np.nan
        out[f"adjusted_ratio_{c}"] = float(np.exp(beta[names.index(f"cond:{c}")]))
    # block bootstrap within runs
    blk = d["run"][idx] * 100000 + d["block"][idx]
    ids = np.unique(blk)
    by_run = {}
    for b in ids:
        by_run.setdefault(int(b // 100000), []).append(b)
    boots = {c: [] for c in conditions[1:]}
    for _ in range(n_boot):
        pick = []
        for r_, bl in by_run.items():
            bl = np.asarray(bl)
            pick.append(rng.choice(bl, bl.size, replace=True))
        pick = np.concatenate(pick)
        sel = np.concatenate([np.flatnonzero(blk == b) for b in pick])
        Xb, _ = design(idx[sel])
        try:
            bb, _, _ = poisson_fit(Xb, y[sel], off[sel], max_iter=30)
        except Exception:
            continue
        for c in conditions[1:]:
            boots[c].append(np.exp(bb[names.index(f"cond:{c}")]))
    for c in conditions[1:]:
        if len(boots[c]) >= 20:
            out[f"adjusted_ratio_ci_{c}"] = [float(np.percentile(boots[c], 2.5)), float(np.percentile(boots[c], 97.5))]
    # what the model says the rate is when still vs running, per condition
    if have_speed:
        for ci, c in enumerate(conditions):
            for state, run_val in (("still", 0.0), ("running", 1.0)):
                sp_val = np.nanmedian(sp_all[d["running"][idx] < 0.5]) if state == "still" else np.nanmedian(sp_all[d["running"][idx] >= 0.5])
                if not np.isfinite(sp_val):
                    continue
                x = {"intercept": 1.0, "speed_z": (sp_val - sp_mean) / sp_sd, "running": run_val, "pupil_z": 0.0}
                for cj, cc in enumerate(conditions[1:], 1):
                    x[f"cond:{cc}"] = 1.0 if cj == ci else 0.0
                eta = sum(out["coef"][nm] * x.get(nm, 0.0) for nm in names)
                out[f"model_rate_{state}_{c}"] = float(np.exp(eta))
    return out


# ----------------------------------------------------------------------------- locomotion-onset triggered rate

def onset_triggered(onsets: np.ndarray, still: np.ndarray, excluded: np.ndarray, spf: float,
                    before_s: float = 40.0, after_s: float = 60.0, min_still_s: float = 20.0) -> dict | None:
    """Mean event rate around running onsets (still for >= min_still_s, then
    running): time axis and events/min, plus the number of onsets."""
    T = onsets.size
    run = still < 0.5
    nb, na = int(round(before_s / spf)), int(round(after_s / spf))
    ms = int(round(min_still_s / spf))
    starts = [t for t in range(ms, T - na) if run[t] and not run[t - ms:t].any()]
    if len(starts) < 3:
        return None
    acc = np.zeros(nb + na); cnt = np.zeros(nb + na)
    for t in starts:
        seg = onsets[t - nb:t + na]; ok = ~excluded[t - nb:t + na]
        acc[ok] += seg[ok]; cnt[ok] += 1
    rate = np.where(cnt > 0, acc / np.maximum(cnt, 1) / (spf / 60.0), np.nan)
    return {"t_s": (np.arange(-nb, na) + 0.5) * spf, "rate": rate, "n_onsets": len(starts)}
