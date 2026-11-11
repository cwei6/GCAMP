"""QC: pre/post projections, z-drift score, registration shift plot."""
from __future__ import annotations
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _percentile_image(mov: np.memmap, start: int, n: int, q: float,
                      max_frames: int = 2000) -> np.ndarray:
    step = max(1, n // max_frames)
    sub = np.asarray(mov[start:start + n:step]).astype(np.float32)
    return np.percentile(sub, q, axis=0)


def _norm(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, [1, 99.5])
    return np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1)


def qc_report(
    plane0: Path,
    out_dir: Path,
    run_names: list[str],
    run_conditions: list[str],
    run_nframes: list[int],
    Ly: int,
    Lx: int,
    loco: list[dict] | None = None,
    fs: float | None = None,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    total = sum(run_nframes)
    mov = np.memmap(plane0 / "data.bin", dtype=np.int16, mode="r",
                    shape=(total, Ly, Lx))

    starts = np.r_[0, np.cumsum(run_nframes)[:-1]]
    nrun = len(run_names)
    means, structs = [], []
    for s, n in zip(starts, run_nframes):
        step = max(1, n // 2000)
        sub = np.asarray(mov[s:s + n:step]).astype(np.float32)
        means.append(sub.mean(axis=0))
        structs.append(np.percentile(sub, 10, axis=0))  # ~static structure

    # pairwise correlation of mean images between runs (z-drift / FOV-change score)
    def corr(a, b):
        a, b = a.ravel() - a.mean(), b.ravel() - b.mean()
        return float(a @ b / np.sqrt((a @ a) * (b @ b)))

    # every run vs the first run (baseline): same-field / z-drift score
    drift_scores = {}
    for j in range(1, nrun):
        drift_scores[f"{run_names[j]} [{run_conditions[j]}] vs {run_names[0]} [{run_conditions[0]}]"] = \
            corr(means[0], means[j])

    # figure: mean + structural image per run
    fig, axes = plt.subplots(2, nrun, figsize=(4 * nrun, 8), squeeze=False)
    for i in range(nrun):
        axes[0][i].imshow(_norm(means[i]), cmap="gray")
        axes[0][i].set_title(f"{run_conditions[i]}: mean\n{run_names[i]}", fontsize=9)
        axes[1][i].imshow(_norm(structs[i]), cmap="gray")
        axes[1][i].set_title("10th pct (structure)", fontsize=9)
        for ax in (axes[0][i], axes[1][i]):
            ax.axis("off")
    fig.suptitle("Registered projections per run — check the SAME axons are visible pre & post")
    fig.tight_layout()
    fig.savefig(out_dir / "qc_projections.png", dpi=150)
    plt.close(fig)

    # registration shifts, if available
    try:
        ops = {}
        for name in ("ops.npy", "reg_outputs.npy"):
            p = plane0 / name
            if p.exists():
                ops.update(np.load(p, allow_pickle=True).item())
        if "yoff" in ops:
            fig, ax = plt.subplots(figsize=(10, 3))
            ax.plot(ops["yoff"], lw=0.5, label="y shift")
            ax.plot(ops["xoff"], lw=0.5, label="x shift")
            for b in np.cumsum(run_nframes)[:-1]:
                ax.axvline(b, color="k", ls="--", lw=0.8)
            ax.set_xlabel("frame"); ax.set_ylabel("shift (px)"); ax.legend()
            ax.set_title("Rigid registration shifts (dashed = run boundaries)")
            fig.tight_layout()
            fig.savefig(out_dir / "qc_shifts.png", dpi=150)
            plt.close(fig)
    except Exception as e:
        print(f"  qc: could not plot shifts ({e})")

    # locomotion overview
    if loco is not None and fs:
        fig, axes = plt.subplots(nrun, 1, figsize=(10, 2.2 * nrun), squeeze=False)
        for i, ld in enumerate(loco):
            t = np.arange(len(ld["velocity"])) / fs
            ax = axes[i][0]
            ax.plot(t, ld["velocity"], lw=0.5)
            ax.fill_between(t, ax.get_ylim()[0], ax.get_ylim()[1],
                            where=ld["still"], alpha=0.15, color="green",
                            label="still")
            ax.set_ylabel("vel"); ax.legend(loc="upper left", fontsize=8)
            if ld.get("pupil") is not None:
                ax2 = ax.twinx()
                ax2.plot(t[:len(ld["pupil"])], ld["pupil"], lw=0.6, color="purple", alpha=0.8,
                         label="pupil diam (px)")
                ax2.set_ylabel("pupil (px)", color="purple"); ax2.legend(loc="upper right", fontsize=8)
            ax.set_title(f"{run_conditions[i]}: {run_names[i]}", fontsize=9)
        axes[-1][0].set_xlabel("time (s)")
        fig.tight_layout()
        fig.savefig(out_dir / "qc_locomotion.png", dpi=150)
        plt.close(fig)
        # behaviour summary per run: still fraction, pupil, pupil-vs-running check
        blines = ["behaviour per run", "-" * 40]
        for i, ld in enumerate(loco):
            line = (f"[{run_conditions[i]:>4}] {run_names[i]}: still {100 * ld['still'].mean():.0f}% "
                    f"of frames, median |vel| {np.nanmedian(np.abs(ld['velocity'])):.1f}")
            if ld.get("pupil") is not None:
                pu = ld["pupil"]; ok = np.isfinite(pu)
                line += (f" | pupil median {np.nanmedian(pu):.1f} px, still {np.nanmedian(pu[ok & ld['still']]):.1f} "
                         f"vs moving {np.nanmedian(pu[ok & ~ld['still']]):.1f}; valid {100 * ok.mean():.1f}%")
                if ok.sum() > 100 and np.std(ld["velocity"]) > 0:
                    r = np.corrcoef(pu[ok], np.abs(ld["velocity"])[ok])[0, 1]
                    line += f"; r(pupil,|vel|) = {r:.2f}" + ("" if r > 0.15 else "  <- low: check eye/quad alignment")
            blines.append(line)
        (out_dir / "qc_behavior.txt").write_text("\n".join(blines), encoding="utf-8")
        print("\n".join(blines))

    # text report
    lines = ["afferent2p QC report", "=" * 40]
    for name, cond, n in zip(run_names, run_conditions, run_nframes):
        lines.append(f"[{cond:>4}] {name}: {n} frames")
    if drift_scores:
        lines.append("")
        lines.append("Mean-image correlation of each run vs the first run (same-field / z-drift score):")
        for k, v in drift_scores.items():
            flag = "OK" if v >= 0.8 else ("CHECK" if v >= 0.6 else "BAD — FOV likely changed")
            lines.append(f"  {k}: r = {v:.3f}   [{flag}]")
        lines.append("")
        lines.append("r >= ~0.8: same field, comparison valid. r < ~0.6: inspect "
                     "qc_projections.png — the plane may have drifted in z (or, across days, "
                     "the FOV was not re-found); per-afferent comparison invalid for that run.")
    else:
        lines.append("(single condition — no pre-vs-post drift score)")
    (out_dir / "qc_report.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    del mov
    return drift_scores
