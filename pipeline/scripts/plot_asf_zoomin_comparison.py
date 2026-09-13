#!/usr/bin/env python3
"""Plot ASF/OPERA and ZoomInSAR LOS-velocity and point time-series comparisons.

The map region, pixel spacing, and geographic reference are taken from the
ZoomInSAR velocity GeoTIFF.  ASF and DEM rasters are reprojected onto that
grid, so the two map panels are directly comparable even if their original
projections or postings differ.

ASF NetCDFs may be supplied locally or downloaded automatically with an
Earthdata login. ``--auto-asf-overlap`` uses the overlap between the ZoomInSAR
date series and available compatible OPERA DISP products. ``--list-asf-subdatasets``
is useful when a NetCDF contains more than one raster layer.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pygmt
import rasterio
import xarray as xr
from rasterio.transform import array_bounds
from rasterio.warp import Resampling, reproject


DATE_PAIR = re.compile(r"(20\d{6})[_-](20\d{6})")


@dataclass(frozen=True)
class GridReference:
    """A geographic target grid based on ZoomInSAR's exported velocity TIFF."""

    data: np.ndarray
    transform: rasterio.Affine
    crs: rasterio.crs.CRS
    region: tuple[float, float, float, float]


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description=("Compare ASF/OPERA and ZoomInSAR LOS velocity over DEM "
                     "hillshade, then compare a point displacement time series.")
    )
    source = command.add_mutually_exclusive_group(required=True)
    source.add_argument("--asf-nc", type=Path,
                        help="Downloaded ASF/OPERA NetCDF used for panel (a).")
    source.add_argument("--download-asf-granule", metavar="GRANULE",
                        help=("ASF granule name to download automatically, for example "
                              "20221107_20221213.unw.nc. Requires an Earthdata login."))
    source.add_argument("--auto-asf-overlap", action="store_true",
                        help=("Find/download OPERA DISP products on the linked ASF track during the "
                              "ZoomInSAR/ASF overlapping period, then derive ASF velocity from that stack."))
    command.add_argument("--asf-download-dir", type=Path,
                         help="Directory for downloaded ASF files (inferred from --project-root when omitted).")
    command.add_argument("--asf-reference-granule", default="20221107_20221213.unw.nc",
                         help=("ASF granule used to identify the correct orbit/direction for "
                               "--auto-asf-overlap (default: 20221107_20221213.unw.nc)."))
    command.add_argument("--asf-max-products", type=int, default=600,
                         help="Safety limit for automatic overlap search results (default: 600).")
    command.add_argument("--project-root", type=Path,
                         help=("One-folder shortcut: project folder containing work/ and ZoomInSAR_Results/. "
                               "It fills the DEM, result, download, and output paths."))
    command.add_argument("--dataset",
                         help="Dataset name after Dataset_ (defaults to the --project-root folder name).")
    command.add_argument("--zoom-velocity", type=Path,
                         help="ZoomInSAR mean_los_velocity_mm_per_year.tif for panel (b).")
    command.add_argument("--dem", type=Path,
                         help="DEM GeoTIFF; it is resampled to the ZoomInSAR grid.")
    command.add_argument("--output", type=Path,
                         help="Output PNG or PDF path.")
    command.add_argument("--asf-velocity-subdataset", default=None,
                         help="Case-insensitive substring selecting the ASF velocity layer.")
    command.add_argument("--list-asf-subdatasets", action="store_true",
                         help="Print all raster layers in --asf-nc and exit.")
    command.add_argument("--asf-velocity-scale", type=float, default=1000.0,
                         help="Multiply ASF velocity by this value (default: m/year to mm/year).")
    command.add_argument("--vmax-mm-year", type=float, default=30.0,
                         help="Shared symmetric velocity limit in mm/year (default: 30).")
    command.add_argument("--transparency", type=float, default=35.0,
                         help="Percent transparency for velocity over the DEM (default: 35).")
    command.add_argument("--point-lon", type=float, help="Representative point longitude.")
    command.add_argument("--point-lat", type=float, help="Representative point latitude.")
    command.add_argument("--point-file", type=Path,
                         help="CSV or GeoJSON point file; defaults to representative_points.csv beside time series.")
    command.add_argument("--zoom-timeseries-dir", type=Path,
                         help=("ZoomInSAR InSAR_Timeseries/OUTPUT/Dataset_* directory. "
                               "Required for panel (c)."))
    command.add_argument("--zoom-data-dir", type=Path,
                         help=("ZoomInSAR Data/ROI/Dataset_* directory containing GEO longitude/latitude files. "
                               "Required for panel (c)."))
    command.add_argument("--asf-timeseries-nc", type=Path, nargs="*", default=[],
                         help=("Chronological ASF displacement NetCDF files for panel (c). "
                               "Their second YYYYMMDD filename date is used as the acquisition date."))
    command.add_argument("--asf-displacement-subdataset", default=None,
                         help="Case-insensitive substring selecting an ASF displacement layer for panel (c).")
    command.add_argument("--asf-displacement-scale", type=float, default=1.0,
                         help="Multiply ASF displacement by this value (default: metres remain metres).")
    command.add_argument("--no-timeseries", action="store_true",
                         help="Create panels (a) and (b) only.")
    return command


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {path}")


def apply_project_defaults(args: argparse.Namespace) -> None:
    """Fill conventional ZoomInSAR paths from one project folder when requested."""
    if not args.project_root:
        args.asf_download_dir = args.asf_download_dir or Path("ASF_downloads")
        args.output = args.output or Path("asf_zoomin_comparison.png")
        if not args.zoom_velocity or not args.dem:
            raise ValueError("Provide --project-root, or provide both --zoom-velocity and --dem.")
        return
    root = args.project_root.expanduser()
    dataset = args.dataset or root.name
    dataset_name = dataset if dataset.startswith("Dataset_") else f"Dataset_{dataset}"
    results = root / "ZoomInSAR_Results"
    timeseries = results / "InSAR_Timeseries" / "OUTPUT" / dataset_name
    data = results / "Data" / "ROI" / dataset_name
    args.zoom_velocity = args.zoom_velocity or timeseries / "mean_los_velocity_mm_per_year.tif"
    args.dem = args.dem or root / "work" / "dem.tif"
    args.zoom_timeseries_dir = args.zoom_timeseries_dir or timeseries
    args.zoom_data_dir = args.zoom_data_dir or data
    args.asf_download_dir = args.asf_download_dir or root / "ASF_overlap_downloads"
    args.output = args.output or root / "asf_vs_zoomin_overlap.png"


def download_asf_granule(granule: str, destination: Path) -> Path:
    """Download one ASF catalog granule using the recipient's Earthdata login.

    ASF recommends ``.netrc`` authentication. An ``EARTHDATA_TOKEN`` environment
    variable is also accepted so neither a password nor token needs to appear in
    a command history or a shared script.
    """
    asf_search, session = asf_session()
    destination.mkdir(parents=True, exist_ok=True)
    existing = destination / Path(granule).name
    if existing.is_file() and existing.stat().st_size:
        print(f"Using existing ASF download: {existing}")
        return existing
    results = asf_search.granule_search([granule])
    if not results:
        raise ValueError(f"ASF did not find a granule named {granule!r}.")
    # With no token, asf_search/requests uses the recipient's ~/.netrc credentials.
    results.download(path=str(destination), session=session)
    candidates = sorted(destination.glob(f"*{Path(granule).stem}*"))
    netcdfs = [path for path in candidates if path.is_file() and path.suffix.lower() == ".nc"]
    if len(netcdfs) == 1:
        return netcdfs[0]
    if existing.is_file():
        return existing
    raise FileNotFoundError(
        f"ASF download completed but the expected NetCDF was not found in {destination}. "
        f"Found: {[path.name for path in candidates]}"
    )


def asf_session():
    """Create an ASF session using .netrc or an optional Earthdata token."""
    try:
        import asf_search
    except ImportError as error:
        raise RuntimeError(
            "Automatic ASF download needs asf_search. Install it with: "
            "conda install -c conda-forge asf-search"
        ) from error
    session = asf_search.ASFSession()
    token = os.environ.get("EARTHDATA_TOKEN")
    if token:
        session.auth_with_token(token)
    return asf_search, session


def product_property(product, *names: str) -> str | None:
    properties = getattr(product, "properties", {})
    for name in names:
        value = properties.get(name)
        if value not in (None, ""):
            return str(value)
    return None


def product_date(product) -> datetime:
    value = product_property(product, "stopTime", "sceneDate", "startTime")
    if value:
        try:
            return datetime.fromisoformat(value.replace("Z", "+00:00")).replace(tzinfo=None)
        except ValueError:
            pass
    name = product_property(product, "fileName", "fileID", "sceneName")
    if name:
        return date_from_asf_filename(Path(name))
    raise ValueError(f"ASF product does not expose an acquisition date: {product}")


def product_filename(product) -> str:
    return product_property(product, "fileName", "fileID", "sceneName", "url") or ""


def overlap_wkt(region: tuple[float, float, float, float]) -> str:
    west, east, south, north = region
    return (f"POLYGON(({west} {south},{east} {south},{east} {north},"
            f"{west} {north},{west} {south}))")


def derive_asf_velocity(displacements: np.ndarray, dates: list[datetime]) -> np.ndarray:
    """Fit mm/year velocity independently at every valid raster cell."""
    if len(dates) < 2:
        raise ValueError("At least two ASF displacement dates are required to derive velocity.")
    years = np.array([(date - dates[0]).total_seconds() / (365.25 * 86400) for date in dates])
    valid = np.isfinite(displacements)
    count = valid.sum(axis=0)
    safe_count = np.where(count, count, 1)
    mean_time = (years[:, None, None] * valid).sum(axis=0) / safe_count
    mean_value = np.nansum(displacements, axis=0) / safe_count
    centered_time = years[:, None, None] - mean_time
    numerator = np.nansum(centered_time * (displacements - mean_value) * valid, axis=0)
    denominator = np.sum((centered_time ** 2) * valid, axis=0)
    velocity = np.full(displacements.shape[1:], np.nan, dtype=np.float32)
    good = (count >= 2) & (denominator > 0)
    velocity[good] = numerator[good] / denominator[good]
    return velocity


def download_asf_overlap(zoom_dates: list[datetime], target: GridReference, reference_granule: str,
                         destination: Path, max_products: int) -> list[tuple[Path, datetime]]:
    """Download compatible OPERA DISP rasters for the ZoomInSAR/ASF time overlap.

    The reference granule determines orbit and look direction, avoiding an
    invalid mixture of ascending and descending LOS observations.
    """
    if max_products < 2:
        raise ValueError("--asf-max-products must be at least 2.")
    asf_search, session = asf_session()
    reference = asf_search.granule_search([reference_granule])
    if not reference:
        raise ValueError(f"ASF did not find reference granule {reference_granule!r}.")
    reference_product = reference[0]
    orbit = product_property(reference_product, "relativeOrbit", "pathNumber")
    direction = product_property(reference_product, "flightDirection")
    if not orbit or not direction:
        raise ValueError(
            "ASF reference metadata did not include orbit/direction. Use an ASF granule from the "
            "same track as ZoomInSAR, or download the matching DISP stack manually."
        )
    start, end = min(zoom_dates), max(zoom_dates)
    results = asf_search.search(
        dataset=asf_search.DATASET.OPERA_S1,
        intersectsWith=overlap_wkt(target.region),
        start=start.strftime("%Y-%m-%dT00:00:00Z"),
        end=end.strftime("%Y-%m-%dT23:59:59Z"),
        maxResults=max_products + 1,
    )
    selected = []
    for product in results:
        name = product_filename(product).lower()
        product_orbit = product_property(product, "relativeOrbit", "pathNumber")
        product_direction = product_property(product, "flightDirection")
        if ("disp" in name and name.endswith(".nc") and product_orbit == orbit
                and product_direction and product_direction.lower() == direction.lower()):
            selected.append(product)
    if len(results) > max_products:
        raise ValueError(
            f"ASF search exceeded --asf-max-products={max_products}; narrow the map region or raise the limit."
        )
    if len(selected) < 2:
        raise ValueError(
            "ASF found fewer than two matching OPERA DISP NetCDF products in the overlap. "
            "A velocity/time-series comparison cannot be made for this interval."
        )
    destination.mkdir(parents=True, exist_ok=True)
    records: list[tuple[Path, datetime]] = []
    for product in sorted(selected, key=product_date):
        name = Path(product_filename(product)).name
        target_file = destination / name
        if not target_file.is_file() or not target_file.stat().st_size:
            product.download(path=str(destination), session=session)
        if not target_file.is_file():
            # ASF occasionally normalizes the local filename; find the newest matching NetCDF.
            candidates = sorted(destination.glob("*.nc"), key=lambda path: path.stat().st_mtime)
            if not candidates:
                raise FileNotFoundError(f"ASF did not create a NetCDF in {destination} for {name}.")
            target_file = candidates[-1]
        records.append((target_file, product_date(product)))
    manifest = {
        "reference_granule": reference_granule,
        "relative_orbit": orbit,
        "flight_direction": direction,
        "zoominsar_dates": [date.strftime("%Y-%m-%d") for date in zoom_dates],
        "asf_overlap_dates": [date.strftime("%Y-%m-%d") for _, date in records],
        "files": [str(path) for path, _ in records],
    }
    (destination / "asf_overlap_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return records


def valid_float_data(dataset: rasterio.DatasetReader, band: int = 1) -> np.ndarray:
    """Read a raster as float32 and normalize all declared mask values to NaN."""
    values = dataset.read(band).astype(np.float32)
    nodata = dataset.nodata
    if nodata is not None:
        values[np.isclose(values, nodata)] = np.nan
    values[~np.isfinite(values)] = np.nan
    return values


def read_zoom_grid(path: Path) -> GridReference:
    require_file(path, "ZoomInSAR velocity GeoTIFF")
    with rasterio.open(path) as dataset:
        if dataset.crs is None or not dataset.crs.is_geographic:
            raise ValueError(
                f"ZoomInSAR velocity GeoTIFF must use longitude/latitude coordinates; got {dataset.crs}."
            )
        west, south, east, north = array_bounds(dataset.height, dataset.width, dataset.transform)
        return GridReference(valid_float_data(dataset), dataset.transform, dataset.crs,
                             (west, east, south, north))


def to_dataarray(values: np.ndarray, grid: GridReference, name: str) -> xr.DataArray:
    """Make an ascending-latitude xarray grid, which GMT handles consistently."""
    height, width = values.shape
    xs = grid.transform.c + (np.arange(width) + 0.5) * grid.transform.a
    ys = grid.transform.f + (np.arange(height) + 0.5) * grid.transform.e
    if ys[0] > ys[-1]:
        values, ys = values[::-1, :], ys[::-1]
    return xr.DataArray(values, coords={"y": ys, "x": xs}, dims=("y", "x"), name=name)


def available_subdatasets(path: Path) -> list[str]:
    """Return GDAL-visible NetCDF rasters, including the single-raster fallback."""
    require_file(path, "ASF NetCDF")
    with rasterio.open(path) as container:
        layers = list(container.subdatasets)
        if not layers and container.count:
            layers = [str(path)]
    return layers


def choose_subdataset(path: Path, requested: str | None, kind: str) -> str:
    layers = available_subdatasets(path)
    if not layers:
        raise ValueError(f"No raster layers could be read from {path}.")
    lowered = [(layer, layer.lower()) for layer in layers]
    if requested:
        matches = [layer for layer, text in lowered if requested.lower() in text]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise ValueError(
                f"No ASF layer contains {requested!r}. Run with --list-asf-subdatasets to inspect choices."
            )
        raise ValueError(f"More than one ASF layer contains {requested!r}: {matches}")

    preferred = ("velocity",) if kind == "velocity" else (
        "short_wavelength_displacement", "displacement", "unwrapped",
    )
    for token in preferred:
        matches = [layer for layer, text in lowered if token in text and "mask" not in text]
        if len(matches) == 1:
            return matches[0]
    layer_text = "\n  ".join(layers)
    raise ValueError(
        f"Could not unambiguously select an ASF {kind} layer. Use the relevant --asf-*-subdataset option.\n"
        f"Available layers:\n  {layer_text}"
    )


def reproject_to_zoom(source_name: str, target: GridReference) -> np.ndarray:
    """Resample any source raster to the exact ZoomInSAR grid."""
    with rasterio.open(source_name) as source:
        if source.crs is None:
            raise ValueError(f"Raster has no CRS: {source_name}")
        values = valid_float_data(source)
        output = np.full(target.data.shape, np.nan, dtype=np.float32)
        reproject(
            source=values,
            destination=output,
            src_transform=source.transform,
            src_crs=source.crs,
            src_nodata=np.nan,
            dst_transform=target.transform,
            dst_crs=target.crs,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    return output


def write_velocity_cpt(path: Path, vmax: float) -> None:
    """Blue for negative LOS motion, white at zero, and red for positive motion."""
    path.write_text(
        f"{-vmax:g} 30/70/180 0 255/255/255\n"
        f"0 255/255/255 {vmax:g} 180/25/25\n"
        "B 30/70/180\nF 180/25/25\nN 128\n",
        encoding="utf-8",
    )


def point_from_file(path: Path) -> tuple[float, float]:
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8", newline="") as stream:
            row = next(csv.DictReader(stream), None)
        if not row:
            raise ValueError(f"Point CSV contains no rows: {path}")
        return float(row["longitude"]), float(row["latitude"])
    payload = json.loads(path.read_text(encoding="utf-8"))
    try:
        coordinates = payload["features"][0]["geometry"]["coordinates"]
        return float(coordinates[0]), float(coordinates[1])
    except (KeyError, IndexError, TypeError) as error:
        raise ValueError(f"Cannot read the first point from GeoJSON: {path}") from error


def choose_point(args: argparse.Namespace) -> tuple[float, float]:
    if (args.point_lon is None) != (args.point_lat is None):
        raise ValueError("Provide both --point-lon and --point-lat.")
    if args.point_lon is not None:
        return args.point_lon, args.point_lat
    candidate = args.point_file
    if candidate is None and args.zoom_timeseries_dir:
        candidate = args.zoom_timeseries_dir / "representative_points.csv"
    if candidate is None:
        raise ValueError("Set --point-lon/--point-lat or provide --point-file for panel (c).")
    require_file(candidate, "Representative point file")
    return point_from_file(candidate)


def load_shape(data_dir: Path) -> tuple[int, int]:
    metadata = data_dir / "crop_metadata.json"
    if metadata.is_file():
        payload = json.loads(metadata.read_text(encoding="utf-8"))
        shape = (payload.get("output_shape_lines_width") or payload.get("output_shape")
                 or payload.get("shape"))
        if shape and len(shape) == 2:
            return int(shape[0]), int(shape[1])
    shape_file = data_dir / "GEO" / "dem_seg.par"
    if shape_file.is_file():
        text = shape_file.read_text(encoding="utf-8", errors="ignore")
        width = re.search(r"width\s*[:=]\s*(\d+)", text)
        lines = re.search(r"(?:nlines|azimuth_lines)\s*[:=]\s*(\d+)", text)
        if width and lines:
            return int(lines.group(1)), int(width.group(1))
    raise FileNotFoundError(
        f"Cannot determine raster shape in {data_dir}. Expected crop_metadata.json or GEO/dem_seg.par."
    )


def find_geo_file(geo_dir: Path, suffix: str) -> Path:
    choices = sorted(geo_dir.glob(f"*.{suffix}"))
    if not choices:
        raise FileNotFoundError(f"No *.{suffix} geolocation file found in {geo_dir}")
    preferred = [path for path in choices if path.stem.lower() in {"sweets", "dem_seg"}]
    return preferred[0] if preferred else choices[0]


def nearest_radar_pixel(data_dir: Path, longitude: float, latitude: float) -> tuple[int, int]:
    shape = load_shape(data_dir)
    geo_dir = data_dir / "GEO"
    lon_path, lat_path = find_geo_file(geo_dir, "lon"), find_geo_file(geo_dir, "lat")
    # ZoomInSAR's importer/GAMMA products store these as big-endian float32 grids.
    lon = np.fromfile(lon_path, dtype=">f4").reshape(shape)
    lat = np.fromfile(lat_path, dtype=">f4").reshape(shape)
    distance = (lon - longitude) ** 2 + (lat - latitude) ** 2
    distance[~np.isfinite(distance)] = np.inf
    row, col = np.unravel_index(np.argmin(distance), shape)
    if not np.isfinite(distance[row, col]):
        raise ValueError("No finite longitude/latitude pixel is available for the selected point.")
    return int(row), int(col)


def load_zoom_timeseries(directory: Path, data_dir: Path, longitude: float, latitude: float) -> tuple[list[datetime], np.ndarray]:
    require_file(directory / "cumulative_los_displacement_m.npy", "ZoomInSAR displacement stack")
    require_file(directory / "dates.json", "ZoomInSAR dates")
    stack = np.load(directory / "cumulative_los_displacement_m.npy", mmap_mode="r")
    dates = [datetime.fromisoformat(value) for value in json.loads((directory / "dates.json").read_text())]
    if stack.ndim != 3 or len(dates) != stack.shape[0]:
        raise ValueError("ZoomInSAR cumulative displacement stack and dates.json are inconsistent.")
    row, col = nearest_radar_pixel(data_dir, longitude, latitude)
    values = np.asarray(stack[:, row, col], dtype=float)
    values[~np.isfinite(values)] = np.nan
    return dates, values - values[np.flatnonzero(np.isfinite(values))[0]]


def date_from_asf_filename(path: Path) -> datetime:
    match = DATE_PAIR.search(path.name)
    if not match:
        raise ValueError(f"Cannot find YYYYMMDD_YYYYMMDD in ASF filename: {path.name}")
    return datetime.strptime(match.group(2), "%Y%m%d")


def sample_asf_displacement(paths: Iterable[Path], target: GridReference, longitude: float,
                            latitude: float, requested: str | None) -> tuple[list[datetime], np.ndarray]:
    row, col = rasterio.transform.rowcol(target.transform, longitude, latitude)
    values, dates = [], []
    for path in paths:
        require_file(path, "ASF displacement NetCDF")
        layer = choose_subdataset(path, requested, "displacement")
        values.append(float(reproject_to_zoom(layer, target)[row, col]))
        dates.append(date_from_asf_filename(path))
    order = np.argsort(dates)
    dates = [dates[index] for index in order]
    series = np.asarray(values, dtype=float)[order]
    valid = np.flatnonzero(np.isfinite(series))
    if not len(valid):
        raise ValueError("The selected point is outside all supplied ASF displacement rasters.")
    return dates, series - series[valid[0]]


def load_auto_asf_stack(records: list[tuple[Path, datetime]], target: GridReference,
                        requested: str | None, scale: float) -> tuple[list[datetime], np.ndarray]:
    """Load the automatically downloaded common-reference DISP displacement stack."""
    dates, grids = [], []
    for path, date in records:
        layer = choose_subdataset(path, requested, "displacement")
        grids.append(reproject_to_zoom(layer, target) * scale)
        dates.append(date)
    order = np.argsort(dates)
    return [dates[index] for index in order], np.stack(grids, axis=0)[order]


def plot_map(fig: pygmt.Figure, dem: xr.DataArray, velocity: xr.DataArray, region: tuple[float, ...],
             projection: str, cpt: str, title: str, point: tuple[float, float] | None,
             transparency: float) -> None:
    shade = pygmt.grdgradient(grid=dem, azimuth=315, normalize="t1")
    fig.grdimage(grid=dem, region=region, projection=projection, cmap="gray", shading=shade,
                 frame=["af", f'+t"{title}"'])
    fig.grdimage(grid=velocity, cmap=cpt, transparency=transparency)
    if point:
        fig.plot(x=point[0], y=point[1], style="c0.22c", fill="yellow", pen="0.6p,black")
    fig.basemap(map_scale="jBL+w2k+f+lkm", rose="jTR+w1.5c+f2+l")


def plot_timeseries(fig: pygmt.Figure, zoom_dates: list[datetime], zoom_values: np.ndarray,
                    asf_dates: list[datetime], asf_values: np.ndarray) -> None:
    all_dates = zoom_dates + asf_dates
    all_values = np.concatenate([zoom_values, asf_values])
    finite = all_values[np.isfinite(all_values)]
    span = max(0.01, float(np.nanmax(finite) - np.nanmin(finite)))
    ymin, ymax = float(np.nanmin(finite) - span * 0.15), float(np.nanmax(finite) + span * 0.15)
    region = [min(all_dates).strftime("%Y-%m-%d"), max(all_dates).strftime("%Y-%m-%d"), ymin, ymax]
    fig.basemap(region=region, projection="X25c/7c",
                frame=["pxa1Yf3o", 'ya+l"Relative LOS displacement (m)"', '+t"(c) Representative-point time series"'])
    fig.plot(x=zoom_dates, y=zoom_values, pen="1.3p,black", style="c0.10c", fill="black", label="ZoomInSAR")
    fig.plot(x=asf_dates, y=asf_values, pen="1.3p,180/25/25", style="t0.16c", fill="180/25/25", label="ASF/OPERA")
    fig.legend(position="jTR+o0.2c", box="+gwhite+p0.25p")


def main() -> None:
    args = parser().parse_args()
    apply_project_defaults(args)
    if args.download_asf_granule:
        args.asf_nc = download_asf_granule(args.download_asf_granule, args.asf_download_dir)
    if args.asf_nc:
        require_file(args.asf_nc, "ASF NetCDF")
    if args.list_asf_subdatasets:
        if not args.asf_nc:
            raise ValueError("--list-asf-subdatasets requires --asf-nc or --download-asf-granule.")
        print("\n".join(available_subdatasets(args.asf_nc)))
        return
    if args.vmax_mm_year <= 0:
        raise ValueError("--vmax-mm-year must be positive.")
    if not 0 <= args.transparency <= 100:
        raise ValueError("--transparency must be between 0 and 100.")
    if args.no_timeseries and (args.asf_timeseries_nc or (args.zoom_timeseries_dir and not args.auto_asf_overlap)):
        raise ValueError("--no-timeseries cannot be combined with time-series inputs.")
    if args.auto_asf_overlap and not args.zoom_timeseries_dir:
        raise ValueError("--auto-asf-overlap needs --zoom-timeseries-dir to determine the overlap dates.")
    if args.auto_asf_overlap and not args.no_timeseries and not args.zoom_data_dir:
        raise ValueError("Panel (c) with --auto-asf-overlap needs --zoom-data-dir.")
    if (not args.auto_asf_overlap and not args.no_timeseries
            and (not args.zoom_timeseries_dir or not args.zoom_data_dir or not args.asf_timeseries_nc)):
        raise ValueError(
            "Panel (c) requires --zoom-timeseries-dir, --zoom-data-dir, and at least one "
            "--asf-timeseries-nc file. Use --no-timeseries for a two-map comparison."
        )

    zoom = read_zoom_grid(args.zoom_velocity)
    require_file(args.dem, "DEM GeoTIFF")
    dem_name = str(args.dem)
    dem = to_dataarray(reproject_to_zoom(dem_name, zoom), zoom, "dem")
    zoom_velocity = to_dataarray(zoom.data, zoom, "zoomin_velocity_mm_year")
    point = None if args.no_timeseries else choose_point(args)

    zoom_dates: list[datetime] = []
    zoom_series = np.asarray([])
    asf_dates: list[datetime] = []
    asf_series = np.asarray([])
    if args.zoom_timeseries_dir:
        require_file(args.zoom_timeseries_dir / "dates.json", "ZoomInSAR dates")
        zoom_dates = [datetime.fromisoformat(value) for value in json.loads(
            (args.zoom_timeseries_dir / "dates.json").read_text(encoding="utf-8")
        )]
    if point:
        zoom_dates, zoom_series = load_zoom_timeseries(
            args.zoom_timeseries_dir, args.zoom_data_dir, *point
        )
    if args.auto_asf_overlap:
        records = download_asf_overlap(
            zoom_dates, zoom, args.asf_reference_granule, args.asf_download_dir, args.asf_max_products
        )
        asf_dates, asf_stack = load_auto_asf_stack(
            records, zoom, args.asf_displacement_subdataset, args.asf_displacement_scale
        )
        asf_velocity = to_dataarray(
            derive_asf_velocity(asf_stack, asf_dates) * 1000.0, zoom, "asf_velocity_mm_year"
        )
        if point:
            row, col = rasterio.transform.rowcol(zoom.transform, point[0], point[1])
            asf_series = np.asarray(asf_stack[:, row, col], dtype=float)
            valid = np.flatnonzero(np.isfinite(asf_series))
            if not len(valid):
                raise ValueError("The selected point is outside the ASF overlap stack.")
            asf_series -= asf_series[valid[0]]
    else:
        asf_velocity_layer = choose_subdataset(args.asf_nc, args.asf_velocity_subdataset, "velocity")
        asf_velocity = to_dataarray(reproject_to_zoom(asf_velocity_layer, zoom) * args.asf_velocity_scale,
                                    zoom, "asf_velocity_mm_year")
        if point:
            asf_dates, asf_series = sample_asf_displacement(
                args.asf_timeseries_nc, zoom, *point, args.asf_displacement_subdataset
            )

    with tempfile.TemporaryDirectory(prefix="zoomin_asf_") as temporary:
        cpt = Path(temporary) / "velocity.cpt"
        write_velocity_cpt(cpt, args.vmax_mm_year)
        fig = pygmt.Figure()
        projection = "M12c"
        plot_map(fig, dem, asf_velocity, zoom.region, projection, str(cpt),
                 "(a) ASF/OPERA LOS velocity", point, args.transparency)
        fig.shift_origin(xshift="13.5c")
        plot_map(fig, dem, zoom_velocity, zoom.region, projection, str(cpt),
                 "(b) ZoomInSAR LOS velocity", point, args.transparency)
        fig.shift_origin(xshift="-13.5c", yshift="-10.5c")
        if point:
            plot_timeseries(fig, zoom_dates, zoom_series, asf_dates, asf_series)
        fig.colorbar(position="JBC+w11c/0.4c+o0c/-8.4c+h", cmap=str(cpt),
                     frame=['xaf+l"LOS velocity (mm/year): blue negative, red positive"'])
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=300)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # Provide a short actionable CLI error rather than a GMT traceback.
        raise SystemExit(f"ERROR: {error}") from error
