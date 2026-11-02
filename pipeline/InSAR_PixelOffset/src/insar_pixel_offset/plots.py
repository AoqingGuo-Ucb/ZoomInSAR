"""Quality-control figures for geocoded pixel-offset tracking."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _limit(*arrays: np.ndarray) -> float:
    values = np.concatenate([array[np.isfinite(array)] for array in arrays])
    return max(float(np.nanpercentile(np.abs(values), 98)), 1e-6) if values.size else 1.0


def plot_pair(
    destination: Path,
    pair_label: str,
    x_offset: np.ndarray,
    y_offset: np.ndarray,
    correlation: np.ndarray,
    unit: str,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    limit = _limit(x_offset, y_offset)
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for axis, values, title in zip(
        axes[:2], (x_offset, y_offset), (f"Map x offset ({unit})", f"Map y offset ({unit})"),
    ):
        image = axis.imshow(values, cmap="RdBu_r", vmin=-limit, vmax=limit, origin="upper")
        axis.set_title(title); axis.set_xlabel("Offset-grid column"); axis.set_ylabel("Offset-grid row")
        figure.colorbar(image, ax=axis, shrink=0.82)
    quality = axes[2].imshow(correlation, cmap="viridis", vmin=0, vmax=1, origin="upper")
    axes[2].set_title("Peak normalized correlation")
    axes[2].set_xlabel("Offset-grid column"); axes[2].set_ylabel("Offset-grid row")
    figure.colorbar(quality, ax=axes[2], shrink=0.82)
    figure.suptitle(f"NISAR GSLC pixel offsets: {pair_label}")
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def plot_summary(destination: Path, records: list[dict], x_stack: np.ndarray,
                 y_stack: np.ndarray, correlation_stack: np.ndarray, unit: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    dates = [datetime.fromisoformat(record["secondary_date"]) for record in records]
    median_x = np.nanmedian(x_stack, axis=(1, 2))
    median_y = np.nanmedian(y_stack, axis=(1, 2))
    median_quality = np.nanmedian(correlation_stack, axis=(1, 2))
    figure, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    axes[0].plot(dates, median_x, "o-", color="#b7374a", label="Map x")
    axes[0].plot(dates, median_y, "o-", color="#286d9a", label="Map y")
    axes[0].axhline(0, color="#555", linewidth=0.8); axes[0].set_ylabel(f"Median offset ({unit})")
    axes[0].legend(); axes[0].grid(alpha=0.25)
    axes[1].plot(dates, median_quality, "o-", color="#38845c")
    axes[1].set_ylim(0, 1); axes[1].set_ylabel("Median correlation"); axes[1].grid(alpha=0.25)
    axes[2].bar(dates, [record["valid_offsets"] for record in records], width=5, color="#6a6485")
    axes[2].set_ylabel("Valid offsets"); axes[2].set_xlabel("Secondary acquisition date"); axes[2].grid(alpha=0.25)
    figure.suptitle("NISAR GSLC pixel-offset tracking summary")
    figure.savefig(destination, dpi=180)
    plt.close(figure)
