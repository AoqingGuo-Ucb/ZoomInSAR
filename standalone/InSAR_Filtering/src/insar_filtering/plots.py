"""Filtering quality-control figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_filter_result(
    destination: Path,
    pair: str,
    original: np.ndarray,
    filtered: np.ndarray,
    coherence: np.ndarray,
    concentration: np.ndarray,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    original_phase, filtered_phase = np.angle(original), np.angle(filtered)
    residual = np.angle(filtered * np.conj(original))
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    panels = [
        (original_phase, "Original wrapped phase", "twilight", -np.pi, np.pi),
        (filtered_phase, "Filtered wrapped phase", "twilight", -np.pi, np.pi),
        (residual, "Circular phase change", "RdBu_r", -np.pi, np.pi),
        (coherence, "Input coherence", "gray", 0, 1),
        (concentration, "Nonlocal phase concentration", "viridis", 0, 1),
        (concentration - coherence, "Concentration minus coherence", "RdBu_r", -1, 1),
    ]
    for axis, (data, title, cmap, vmin, vmax) in zip(axes.flat, panels):
        image = axis.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin="upper", interpolation="nearest")
        axis.set_title(title); axis.set_xlabel("Range column"); axis.set_ylabel("Azimuth row")
        fig.colorbar(image, ax=axis, shrink=0.82)
    fig.suptitle(f"InSAR phase filtering: {pair}", fontsize=15)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def plot_summary(destination: Path, rows: list[dict]) -> None:
    pairs = [row["pair"] for row in rows]
    x = np.arange(len(rows))
    before = [row["input_circular_variance"] for row in rows]
    after = [row["filtered_circular_variance"] for row in rows]
    concentration = [row["median_phase_concentration"] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    width = 0.38
    axes[0].bar(x - width / 2, before, width, label="Input")
    axes[0].bar(x + width / 2, after, width, label="Filtered")
    axes[0].set_title("Local circular variance (lower is smoother)"); axes[0].legend()
    axes[1].bar(x, concentration, color="seagreen")
    axes[1].set_ylim(0, 1); axes[1].set_title("Median nonlocal phase concentration")
    labels = [pair[2:8] + "\n" + pair[11:17] for pair in pairs]
    for axis in axes:
        axis.set_xticks(x, labels, rotation=45, ha="right", fontsize=8); axis.grid(axis="y", alpha=0.25)
    fig.savefig(destination, dpi=180)
    plt.close(fig)
