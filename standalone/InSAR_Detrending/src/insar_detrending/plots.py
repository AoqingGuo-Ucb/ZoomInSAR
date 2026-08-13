"""Detrending diagnostic figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _limits(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return -1, 1
    maximum = max(float(np.nanpercentile(np.abs(finite), 98)), 1e-6)
    return -maximum, maximum


def plot_case(path: Path, name: str, before: np.ndarray, surface: np.ndarray,
              after: np.ndarray, fit_mask: np.ndarray,
              reference: dict | None = None,
    components: dict[str, np.ndarray] | None = None) -> None:
    components = components or {}
    fig, axes = plt.subplots(3, 4, figsize=(20, 14), constrained_layout=True)
    panels = (
        (axes[0, 0], before, "Before detrending", "RdBu_r"),
        (axes[0, 1], components.get("orbit", np.full_like(before, np.nan)), "Orbit ramp", "RdBu_r"),
        (axes[0, 2], components.get("after_orbit", np.full_like(before, np.nan)), "After orbit correction", "RdBu_r"),
        (axes[0, 3], components.get("terrain", np.full_like(before, np.nan)), "DEM-correlated term", "RdBu_r"),
        (axes[1, 0], components.get("after_dem", np.full_like(before, np.nan)), "After orbit + DEM correction", "RdBu_r"),
        (axes[1, 1], components.get("turbulent_qa", np.full_like(before, np.nan)), "Non-DEM low-frequency residual (QA only)", "RdBu_r"),
        (
            axes[1, 2],
            components.get("after_turbulent_qa", np.full_like(before, np.nan)),
            "After subtracting QA residual (diagnostic)",
            "RdBu_r",
        ),
        (axes[1, 3], surface, "Removed correction", "RdBu_r"),
        (axes[2, 0], after, "Final after accepted correction", "RdBu_r"),
        (axes[2, 1], fit_mask.astype(float), "Final robust fit mask", "gray"),
        (axes[2, 2], np.isfinite(before).astype(float), "Valid pixels", "gray"),
    )
    for axis, values, title, cmap in panels:
        if title == "Final robust fit mask":
            image = axis.imshow(values, cmap=cmap, vmin=0, vmax=1, origin="upper", interpolation="nearest")
        elif title == "Valid pixels":
            image = axis.imshow(values, cmap=cmap, vmin=0, vmax=1, origin="upper", interpolation="nearest")
        else:
            vmin, vmax = _limits(values)
            image = axis.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, origin="upper", interpolation="nearest")
        if reference is not None:
            axis.plot(
                reference["col"], reference["row"], marker="+", color="yellow",
                markersize=14, markeredgewidth=2.2,
            )
            axis.plot(
                reference["col"], reference["row"], marker="o", color="black",
                markersize=12, markerfacecolor="none", markeredgewidth=1.5,
            )
        axis.set_title(title); axis.set_xlabel("Range column"); axis.set_ylabel("Azimuth row")
        fig.colorbar(image, ax=axis, shrink=0.82)
    axes[2, 3].axis("off")
    fig.suptitle(name, fontsize=15); fig.savefig(path, dpi=180); plt.close(fig)


def plot_summary(path: Path, names: list[str], before: list[float], after: list[float], ranges: list[float]) -> None:
    x = np.arange(len(names)); width = 0.38
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    axes[0].bar(x - width / 2, before, width, label="Before")
    axes[0].bar(x + width / 2, after, width, label="After")
    axes[0].set_title("Robust spatial standard deviation (rad)"); axes[0].legend()
    axes[1].bar(x, ranges, color="darkorange"); axes[1].set_title("Removed trend robust range (rad)")
    labels = [name[2:8] + "\n" + name[11:17] for name in names]
    for axis in axes:
        axis.set_xticks(x, labels, rotation=45, ha="right", fontsize=8); axis.grid(axis="y", alpha=0.25)
    fig.savefig(path, dpi=180); plt.close(fig)
