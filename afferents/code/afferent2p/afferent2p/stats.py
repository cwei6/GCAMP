"""Scale-aware statistics for comparing event counts between conditions
(used by the main comparison in aqua_analyze.py and by the trigger
comparison in afferent2p.triggers)."""
from __future__ import annotations

import numpy as np


def poisson_rate_ci(n: float, minutes: float, alpha: float = 0.05):
    """Exact (Garwood) 95 % interval for a rate from n events in `minutes`."""
    from scipy.stats import chi2
    if not np.isfinite(n) or not np.isfinite(minutes) or minutes <= 0:
        return np.nan, np.nan
    n = float(n)
    lo = chi2.ppf(alpha / 2, 2 * n) / 2 if n > 0 else 0.0
    hi = chi2.ppf(1 - alpha / 2, 2 * n + 2) / 2
    return lo / minutes, hi / minutes


def rate_ratio_ci(nb: float, tb: float, nc: float, tc: float, alpha: float = 0.05):
    """Exact conditional 95 % interval for the ratio of two Poisson rates
    (condition / baseline): given nb + nc events, nc is binomial with
    p = rc*tc / (rc*tc + rb*tb); Clopper-Pearson on p, mapped back. An
    infinite upper bound when the baseline has no events, 0 lower bound when
    the condition has none; NaN when both are empty."""
    from scipy.stats import beta
    if any(not np.isfinite(v) for v in (nb, tb, nc, tc)) or tb <= 0 or tc <= 0:
        return np.nan, np.nan
    nb, nc = float(nb), float(nc)
    if nb + nc == 0:
        return np.nan, np.nan
    p_lo = beta.ppf(alpha / 2, nc, nb + 1) if nc > 0 else 0.0
    p_hi = beta.ppf(1 - alpha / 2, nc + 1, nb) if nb > 0 else 1.0
    lo = (p_lo / (1 - p_lo)) * tb / tc if p_lo < 1 else np.inf
    hi = (p_hi / (1 - p_hi)) * tb / tc if p_hi < 1 else np.inf
    return float(lo), float(hi)


def change_measures(nb, tb, nc, tc, pseudo: float = 0.5) -> dict:
    """Scale-aware descriptions of baseline -> condition from the event counts
    and minutes: the rate difference, the log2 ratio with a pseudo-count
    (defined when a condition has no events), the exact interval of the
    ratio, the symmetric modulation index (rc - rb) / (rc + rb) in [-1, 1],
    and a low-count flag (fewer than 10 events in both conditions together)."""
    out = {}
    try:
        nb, tb, nc, tc = float(nb), float(tb), float(nc), float(tc)
    except (TypeError, ValueError):
        return out
    if not all(np.isfinite(v) for v in (nb, tb, nc, tc)) or tb <= 0 or tc <= 0:
        return out
    rb, rc = nb / tb, nc / tc
    out["n_events_baseline"] = int(nb); out["n_events_condition"] = int(nc)
    out["delta_rate"] = rc - rb
    out["log2_ratio_pc"] = float(np.log2(((nc + pseudo) / tc) / ((nb + pseudo) / tb)))
    lo, hi = rate_ratio_ci(nb, tb, nc, tc)
    out["ratio_ci_lo"], out["ratio_ci_hi"] = lo, hi
    out["modulation_index"] = (rc - rb) / (rc + rb) if (rc + rb) > 0 else np.nan
    out["low_counts"] = bool(nb + nc < 10 or min(nb, nc) < 3)      # too few to SIZE the change (direction may still be clear)
    out["rate_baseline_ci_lo"], out["rate_baseline_ci_hi"] = poisson_rate_ci(nb, tb)
    out["rate_condition_ci_lo"], out["rate_condition_ci_hi"] = poisson_rate_ci(nc, tc)
    return out




def fraction_change(kb: float, nb: float, kc: float, nc: float, alpha: float = 0.05) -> dict:
    """Responsive fraction baseline -> condition: the two fractions, their
    difference, the risk ratio with a log (Katz) interval (+0.5 continuity
    when a count is 0), Fisher's exact p (two-sided) and a low-count flag."""
    out = {}
    try:
        kb, nb, kc, nc = float(kb), float(nb), float(kc), float(nc)
    except (TypeError, ValueError):
        return out
    if not all(np.isfinite(v) for v in (kb, nb, kc, nc)) or nb <= 0 or nc <= 0:
        return out
    fb, fc = kb / nb, kc / nc
    out["frac_baseline"], out["frac_condition"] = fb, fc
    out["frac_delta"] = fc - fb
    a, b_, c, d = kc, nc - kc, kb, nb - kb
    cc = 0.5 if min(a, c) == 0 else 0.0
    rr = ((a + cc) / (nc + 2 * cc)) / ((c + cc) / (nb + 2 * cc))
    se = np.sqrt(max(1.0 / (a + cc) - 1.0 / (nc + 2 * cc) + 1.0 / (c + cc) - 1.0 / (nb + 2 * cc), 0.0))
    from scipy.stats import norm, fisher_exact
    z = norm.ppf(1 - alpha / 2)
    out["frac_ratio"] = rr
    out["frac_ratio_ci_lo"], out["frac_ratio_ci_hi"] = float(rr * np.exp(-z * se)), float(rr * np.exp(z * se))
    try:
        out["frac_p_fisher"] = float(fisher_exact([[a, b_], [c, d]], alternative="two-sided")[1])
    except Exception:
        out["frac_p_fisher"] = np.nan
    out["frac_low_counts"] = bool(nb < 5 or nc < 5)
    return out


def amplitude_change(xb, xc, n_boot: int = 2000, seed: int = 0, alpha: float = 0.05) -> dict:
    """A continuous per-trial measure baseline -> condition: means, their
    difference with a bootstrap percentile interval, and the Mann-Whitney p."""
    out = {}
    xb = np.asarray(xb, float); xc = np.asarray(xc, float)
    xb, xc = xb[np.isfinite(xb)], xc[np.isfinite(xc)]
    out["amp_n_baseline"], out["amp_n_condition"] = int(xb.size), int(xc.size)
    out["amp_baseline"] = float(xb.mean()) if xb.size else np.nan
    out["amp_condition"] = float(xc.mean()) if xc.size else np.nan
    out["amp_delta"] = out["amp_condition"] - out["amp_baseline"] if xb.size and xc.size else np.nan
    if xb.size >= 2 and xc.size >= 2:
        rng = np.random.default_rng(seed)
        d = np.array([rng.choice(xc, xc.size).mean() - rng.choice(xb, xb.size).mean() for _ in range(n_boot)])
        out["amp_delta_ci_lo"], out["amp_delta_ci_hi"] = float(np.percentile(d, 100 * alpha / 2)), float(np.percentile(d, 100 * (1 - alpha / 2)))
        try:
            from scipy.stats import mannwhitneyu
            out["amp_p_mwu"] = float(mannwhitneyu(xb, xc, alternative="two-sided").pvalue)
        except Exception:
            out["amp_p_mwu"] = np.nan
    else:
        out["amp_delta_ci_lo"] = out["amp_delta_ci_hi"] = out["amp_p_mwu"] = np.nan
    return out
