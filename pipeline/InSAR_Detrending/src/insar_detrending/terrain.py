"""DEM sampling for elevation-correlated phase correction."""

from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform


def _nanmean_filter_axis(values: np.ndarray, radius: int, axis: int) -> np.ndarray:
    if radius <= 0:
        return values.copy()
    moved = np.moveaxis(values, axis, 0)
    valid = np.isfinite(moved)
    filled = np.where(valid, moved, 0.0)
    weights = valid.astype(float)
    pad_width = [(0, 0)] * filled.ndim
    pad_width[0] = (radius, radius)
    filled_pad = np.pad(filled, pad_width, mode="constant", constant_values=0.0)
    weights_pad = np.pad(weights, pad_width, mode="constant", constant_values=0.0)
    leading = np.zeros((1, *filled_pad.shape[1:]), dtype=float)
    sums = np.concatenate([leading, np.cumsum(filled_pad, axis=0)], axis=0)
    counts = np.concatenate([leading, np.cumsum(weights_pad, axis=0)], axis=0)
    window = 2 * radius + 1
    sums = sums[window:] - sums[:-window]
    counts = counts[window:] - counts[:-window]
    with np.errstate(invalid="ignore", divide="ignore"):
        filtered = sums / counts
    filtered[counts == 0] = np.nan
    return np.moveaxis(filtered, 0, axis)


def lowpass_elevation(elevation: np.ndarray, radius_pixels: int, passes: int = 3) -> np.ndarray:
    """Return a low-frequency DEM field for terrain-correlated trend fitting.

    The repeated box filters suppress local drainage and ridge texture so the
    DEM term represents broad topographic dependence rather than detailed relief.
    """
    radius_pixels = int(radius_pixels)
    if radius_pixels <= 0:
        return np.asarray(elevation, dtype=float).copy()
    smoothed = np.asarray(elevation, dtype=float).copy()
    for _ in range(max(1, int(passes))):
        smoothed = _nanmean_filter_axis(smoothed, radius_pixels, axis=0)
        smoothed = _nanmean_filter_axis(smoothed, radius_pixels, axis=1)
    smoothed[~np.isfinite(elevation)] = np.nan
    return smoothed


def sample_dem(dataset: Path, shape: tuple[int, int], dem_path: Path) -> np.ndarray:
    geo = dataset / "GEO"
    lon_files, lat_files = sorted(geo.glob("*.lon")), sorted(geo.glob("*.lat"))
    if len(lon_files) != 1 or len(lat_files) != 1:
        raise ValueError(f"{geo} must contain exactly one *.lon and one *.lat")
    dtype = np.dtype("f4").newbyteorder(">")
    longitude = np.asarray(np.memmap(lon_files[0], dtype=dtype, mode="r", shape=shape), float)
    latitude = np.asarray(np.memmap(lat_files[0], dtype=dtype, mode="r", shape=shape), float)
    valid = (
        np.isfinite(longitude) & np.isfinite(latitude)
        & (longitude >= -180) & (longitude <= 180)
        & (latitude >= -90) & (latitude <= 90)
        & ~((longitude == 0) & (latitude == 0))
    )
    elevation = np.full(shape, np.nan, dtype=float)
    with rasterio.open(dem_path) as dem:
        xs, ys = transform("EPSG:4326", dem.crs,
                           longitude[valid].tolist(), latitude[valid].tolist())
        samples = np.fromiter(
            (value[0] for value in dem.sample(zip(xs, ys), masked=False)),
            dtype=float, count=int(np.count_nonzero(valid)),
        )
        if dem.nodata is not None:
            samples[np.isclose(samples, dem.nodata)] = np.nan
        elevation[valid] = samples
    return elevation
