"""High-resolution DEM water masking at InSAR longitude/latitude pixels."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform

from .io import read_gamma


def find_lon_lat(geo_dir: Path) -> tuple[Path, Path]:
    longitude = sorted(geo_dir.glob("*.lon"))
    latitude = sorted(geo_dir.glob("*.lat"))
    if len(longitude) != 1 or len(latitude) != 1:
        raise ValueError(f"{geo_dir} must contain exactly one *.lon and one *.lat file")
    return longitude[0], latitude[0]


def water_mask_from_dem(
    dataset: str | Path,
    shape: tuple[int, int],
    dem_path: str | Path,
    water_max_elevation: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample a DEM at every lon/lat and classify elevation <= threshold as water.

    Returns water mask, valid-coordinate mask, longitude, latitude, and sampled
    elevation. The DEM must cover every valid geolocation pixel so missing
    coverage cannot silently be mislabeled as water.
    """
    dataset, dem_path = Path(dataset), Path(dem_path)
    if not dem_path.is_file():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    lon_file, lat_file = find_lon_lat(dataset / "GEO")
    longitude = read_gamma(lon_file, shape, "f4").astype(float)
    latitude = read_gamma(lat_file, shape, "f4").astype(float)
    valid_geo = (
        np.isfinite(longitude) & np.isfinite(latitude)
        & (longitude >= -180.0) & (longitude <= 180.0)
        & (latitude >= -90.0) & (latitude <= 90.0)
        & ~((longitude == 0.0) & (latitude == 0.0))
    )
    elevation = np.full(shape, np.nan, dtype=float)
    with rasterio.open(dem_path) as dem:
        source_crs = "EPSG:4326"
        xs, ys = transform(
            source_crs, dem.crs,
            longitude[valid_geo].tolist(), latitude[valid_geo].tolist(),
        )
        xs, ys = np.asarray(xs), np.asarray(ys)
        inside = (
            (xs >= dem.bounds.left) & (xs <= dem.bounds.right)
            & (ys >= dem.bounds.bottom) & (ys <= dem.bounds.top)
        )
        if not np.all(inside):
            raise ValueError(
                f"DEM {dem_path} does not cover {int(np.count_nonzero(~inside))} valid lon/lat pixels"
            )
        samples = np.fromiter(
            (value[0] for value in dem.sample(zip(xs, ys), masked=False)),
            dtype=float, count=len(xs),
        )
        if dem.nodata is not None:
            samples[np.isclose(samples, dem.nodata)] = np.nan
        elevation[valid_geo] = samples
    missing_dem = valid_geo & ~np.isfinite(elevation)
    if np.any(missing_dem):
        raise ValueError(
            f"DEM has NoData at {int(np.count_nonzero(missing_dem))} valid lon/lat pixels"
        )
    water = valid_geo & (elevation <= float(water_max_elevation))
    return water, valid_geo, longitude, latitude, elevation
