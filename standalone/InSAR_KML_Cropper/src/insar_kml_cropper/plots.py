"""Geographic quality-control plots for KML cropping."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from .kml import kml_coordinate_rings


def plot_kml_crop_overlay(
    destination: str | Path,
    kml_file: str | Path,
    longitude: np.ndarray,
    latitude: np.ndarray,
    coherence: np.ndarray | None,
    expanded_bounds: tuple[float, float, float, float],
) -> Path:
    """Overlay KML, expanded bounds, and cropped pixels in geographic space."""
    destination, kml_file = Path(destination), Path(kml_file)
    valid = np.isfinite(longitude) & np.isfinite(latitude) & (longitude != 0) & (latitude != 0)
    if coherence is not None:
        valid &= np.isfinite(coherence)
    fig, axis = plt.subplots(figsize=(10, 8), constrained_layout=True)
    if coherence is not None and np.any(valid):
        points = axis.scatter(
            longitude[valid], latitude[valid], c=coherence[valid], s=5,
            cmap="gray", vmin=0, vmax=1, linewidths=0, rasterized=True,
            label="Cropped coherence pixels",
        )
        fig.colorbar(points, ax=axis, shrink=0.82, label="Coherence")
    elif np.any(valid):
        axis.scatter(longitude[valid], latitude[valid], s=3, color="0.65",
                     linewidths=0, rasterized=True, label="Cropped pixels")
    for index, ring in enumerate(kml_coordinate_rings(kml_file)):
        coordinates = np.asarray(ring)
        axis.plot(coordinates[:, 0], coordinates[:, 1], color="red", linewidth=2.2,
                  label="Original KML" if index == 0 else None)
    west, east, south, north = expanded_bounds
    axis.add_patch(Rectangle(
        (west, south), east - west, north - south, fill=False,
        edgecolor="dodgerblue", linewidth=2, linestyle="--",
        label="Expanded KML bounds",
    ))
    if np.any(valid):
        pixel_west, pixel_east = float(np.nanmin(longitude[valid])), float(np.nanmax(longitude[valid]))
        pixel_south, pixel_north = float(np.nanmin(latitude[valid])), float(np.nanmax(latitude[valid]))
        axis.add_patch(Rectangle(
            (pixel_west, pixel_south), pixel_east - pixel_west, pixel_north - pixel_south,
            fill=False, edgecolor="limegreen", linewidth=1.6, linestyle=":",
            label="Final crop pixel extent",
        ))
    axis.set_title(f"KML crop overlay: {kml_file.stem}")
    axis.set_xlabel("Longitude (degrees)"); axis.set_ylabel("Latitude (degrees)")
    axis.set_aspect("equal", adjustable="box"); axis.grid(alpha=0.25)
    axis.legend(loc="best")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=200)
    plt.close(fig)
    return destination
