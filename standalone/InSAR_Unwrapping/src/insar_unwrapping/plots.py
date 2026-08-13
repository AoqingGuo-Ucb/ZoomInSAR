"""Non-interactive diagnostic figures for unwrapping outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _phase_limits(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return -np.pi, np.pi
    low, high = np.nanpercentile(finite, [2, 98])
    if high <= low:
        return float(low - 1), float(high + 1)
    return float(low), float(high)


def plot_interferogram_diagnostics(
    output: str | Path,
    pair: str,
    wrapped: np.ndarray,
    coherence: np.ndarray,
    branch_cuts: np.ndarray,
    before_closure: np.ndarray,
    closure_cycles: np.ndarray,
    final: np.ndarray,
) -> Path:
    """Create a six-panel quality-control figure for one interferogram."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    panels = [
        (wrapped, "Wrapped phase", "twilight", (-np.pi, np.pi)),
        (coherence, "Coherence", "gray", (0.0, 1.0)),
        (branch_cuts.astype(float), "Branch cuts", "Reds", (0.0, 1.0)),
        (before_closure, "Unwrapped before closure repair", "turbo", _phase_limits(before_closure)),
        (closure_cycles, "Closure correction (cycles)", "RdBu_r", None),
        (final, "Final unwrapped phase", "turbo", _phase_limits(final)),
    ]
    for axis, (values, title, cmap, limits) in zip(axes.flat, panels):
        kwargs = {"cmap": cmap, "origin": "upper", "interpolation": "nearest"}
        if limits:
            kwargs.update(vmin=limits[0], vmax=limits[1])
        elif np.any(np.isfinite(values)):
            maximum = max(1.0, float(np.nanmax(np.abs(values))))
            kwargs.update(vmin=-maximum, vmax=maximum)
        image = axis.imshow(values, **kwargs)
        axis.set_title(title)
        axis.set_xlabel("Range column")
        axis.set_ylabel("Azimuth row")
        fig.colorbar(image, ax=axis, shrink=0.82)
    fig.suptitle(f"Interferogram {pair}", fontsize=15)
    destination = output / f"{pair}_diagnostics.png"
    fig.savefig(destination, dpi=180)
    plt.close(fig)
    return destination


def plot_dataset_summary(
    output: str | Path,
    pairs: list[str],
    quality_rows: list[dict],
    closure_summary: dict | None,
) -> Path:
    """Plot residue, branch-cut, coverage, and closure-repair summaries."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    x = np.arange(len(pairs))
    positive = np.array([row["positive_residues"] for row in quality_rows])
    negative = np.array([row["negative_residues"] for row in quality_rows])
    cuts = np.array([row["branch_cut_pixels"] for row in quality_rows])
    coverage = np.array([
        100.0 * row["unwrapped_pixels"] / max(row["valid_pixels"], 1)
        for row in quality_rows
    ])
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    width = 0.38
    axes[0, 0].bar(x - width / 2, positive, width, label="Positive")
    axes[0, 0].bar(x + width / 2, negative, width, label="Negative")
    axes[0, 0].set_title("Phase residues")
    axes[0, 0].legend()
    axes[0, 1].bar(x, cuts, color="firebrick")
    axes[0, 1].set_title("Branch-cut pixels")
    axes[1, 0].bar(x, coverage, color="seagreen")
    axes[1, 0].set_ylim(0, 105)
    axes[1, 0].set_ylabel("Coverage (%)")
    axes[1, 0].set_title("Unwrapped valid-pixel coverage")
    if closure_summary:
        values = [closure_summary["bad_pixels_before"], closure_summary["bad_pixels_after"]]
        bars = axes[1, 1].bar(["Before", "After"], values, color=["darkorange", "steelblue"])
        axes[1, 1].bar_label(bars)
        axes[1, 1].set_title("Integer closure-error pixels")
    else:
        axes[1, 1].text(0.5, 0.5, "Closure repair disabled", ha="center", va="center")
        axes[1, 1].set_axis_off()
    short_labels = [pair[2:8] + "\n" + pair[11:17] for pair in pairs]
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xticks(x, short_labels, rotation=45, ha="right", fontsize=8)
        axis.grid(axis="y", alpha=0.25)
    fig.suptitle("Unwrapping quality summary", fontsize=16)
    destination = output / "dataset_quality_summary.png"
    fig.savefig(destination, dpi=180)
    plt.close(fig)
    return destination


def plot_water_mask(
    output: str | Path,
    longitude: np.ndarray,
    latitude: np.ndarray,
    elevation: np.ndarray,
    water_mask: np.ndarray,
    valid_geo: np.ndarray,
    water_max_elevation: float,
) -> Path:
    """Plot geolocation rasters and water classification used by the pipeline."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    classification = np.full(water_mask.shape, np.nan)
    classification[valid_geo & ~water_mask] = 0.0
    classification[water_mask] = 1.0
    longitude_plot = np.where(valid_geo, longitude, np.nan)
    latitude_plot = np.where(valid_geo, latitude, np.nan)
    elevation_plot = np.where(valid_geo, elevation, np.nan)
    fig, axes = plt.subplots(1, 4, figsize=(19, 5), constrained_layout=True)
    for axis, data, title, cmap in (
        (axes[0], longitude_plot, "Longitude", "viridis"),
        (axes[1], latitude_plot, "Latitude", "viridis"),
        (axes[2], elevation_plot, "USGS DEM elevation (m)", "terrain"),
        (axes[3], classification,
         f"Water mask (elevation <= {water_max_elevation:g} m)", "Blues"),
    ):
        image = axis.imshow(data, origin="upper", cmap=cmap, interpolation="nearest")
        axis.set_title(title); axis.set_xlabel("Range column"); axis.set_ylabel("Azimuth row")
        fig.colorbar(image, ax=axis, shrink=0.82)
    destination = output / "water_mask.png"
    fig.savefig(destination, dpi=180)
    plt.close(fig)
    return destination
