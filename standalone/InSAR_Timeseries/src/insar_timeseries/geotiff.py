"""Geographic raster and representative-point exports for QGIS."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_bounds
from scipy.interpolate import griddata


def _read_lon_lat_from_geo(geo: Path, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    lon_files, lat_files = sorted(geo.glob("*.lon")), sorted(geo.glob("*.lat"))
    if len(lon_files) != 1 or len(lat_files) != 1:
        raise ValueError(f"{geo} must contain exactly one *.lon and one *.lat")
    dtype = np.dtype("f4").newbyteorder(">")
    longitude = np.asarray(np.memmap(lon_files[0], dtype=dtype, mode="r", shape=shape), float)
    latitude = np.asarray(np.memmap(lat_files[0], dtype=dtype, mode="r", shape=shape), float)
    return longitude, latitude


def load_geolocation(dataset: Path, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Locate GEO rasters through the Detrending/Unwrapping provenance chain."""
    summary = json.loads((dataset / "run_summary.json").read_text(encoding="utf-8"))
    candidates = [dataset / "GEO"]
    config = summary.get("config", {})
    if config.get("geo_root"):
        candidates.append(Path(config["geo_root"]) / dataset.name / "GEO")
    source = summary.get("source_unwrapping_summary", summary)
    if source.get("input_directory"):
        candidates.append(Path(source["input_directory"]) / "GEO")
    for geo in candidates:
        if geo.is_dir() and list(geo.glob("*.lon")) and list(geo.glob("*.lat")):
            return _read_lon_lat_from_geo(geo, shape)
    raise FileNotFoundError(f"Cannot locate one *.lon and one *.lat for {dataset}")


def export_velocity_geotiff(
    path: Path, velocity_m_per_year: np.ndarray,
    longitude: np.ndarray, latitude: np.ndarray,
) -> dict:
    """Resample radar geometry to a north-up EPSG:4326 GeoTIFF in mm/year."""
    valid = (
        np.isfinite(velocity_m_per_year) & np.isfinite(longitude) & np.isfinite(latitude)
        & (longitude >= -180) & (longitude <= 180)
        & (latitude >= -90) & (latitude <= 90)
        & ~((longitude == 0) & (latitude == 0))
    )
    if np.count_nonzero(valid) < 4:
        raise ValueError("Too few valid velocity/lon/lat pixels for GeoTIFF export")
    west, east = float(np.nanmin(longitude[valid])), float(np.nanmax(longitude[valid]))
    south, north = float(np.nanmin(latitude[valid])), float(np.nanmax(latitude[valid]))
    height, width = velocity_m_per_year.shape
    transform = from_bounds(west, south, east, north, width, height)
    xs = west + (np.arange(width) + 0.5) * (east - west) / width
    ys = north - (np.arange(height) + 0.5) * (north - south) / height
    grid_x, grid_y = np.meshgrid(xs, ys)
    points = np.column_stack((longitude[valid], latitude[valid]))
    values = velocity_m_per_year[valid] * 1000.0
    regular = griddata(points, values, (grid_x, grid_y), method="linear")
    nodata = np.float32(-9999.0)
    output = np.where(np.isfinite(regular), regular, nodata).astype(np.float32)
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=1,
        dtype="float32", crs="EPSG:4326", transform=transform, nodata=nodata,
        compress="deflate", predictor=3, tiled=True,
    ) as dst:
        dst.write(output, 1)
        dst.set_band_description(1, "Mean LOS velocity (mm/year)")
        dst.update_tags(
            AREA_OR_POINT="Point", units="mm/year",
            quantity="mean_line_of_sight_velocity",
        )
        dst.update_tags(1, units="mm/year", quantity="mean_line_of_sight_velocity")
    return {
        "path": str(path), "crs": "EPSG:4326", "units": "mm/year",
        "bounds": [west, south, east, north], "width": width, "height": height,
        "valid_output_pixels": int(np.count_nonzero(output != nodata)),
    }
