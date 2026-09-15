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
import time
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
from scipy.interpolate import griddata


DATE_PAIR = re.compile(r"(20\d{6})[_-](20\d{6})")
OPERA_DISP_NAME = re.compile(
    r"OPERA_L3_DISP-S1_IW_F(?P<frame>\d+)_VV_"
    r"(?P<reference>20\d{6})T\d{6}Z_(?P<secondary>20\d{6})T\d{6}Z"
)


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
    command.add_argument("--asf-reference-granule",
                         help=("Optional ASF granule used to force its orbit/direction. Usually omit this; "
                               "automatic mode otherwise selects the best-covered compatible DISP track."))
    command.add_argument("--asf-flight-direction", choices=("ascending", "descending"),
                         help=("ASF flight direction for --auto-asf-overlap. Set this to descending when "
                               "the ZoomInSAR input is descending; it prevents an ascending-track comparison."))
    command.add_argument("--asf-max-products", type=int, default=5000,
                         help=("Safety limit for ASF catalogue search results before track filtering "
                               "(default: 5000; products are filtered before download)."))
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


def readable_netcdf(path: Path) -> bool:
    """Check enough of a NetCDF to reject an interrupted/partial download."""
    if not path.is_file() or not path.stat().st_size:
        return False
    try:
        with rasterio.open(path) as dataset:
            layers = dataset.subdatasets
            if layers:
                with rasterio.open(layers[0]) as layer:
                    layer.read(1, window=((layer.height - 1, layer.height),
                                          (layer.width - 1, layer.width)))
            elif dataset.count:
                dataset.read(1, window=((dataset.height - 1, dataset.height),
                                        (dataset.width - 1, dataset.width)))
            else:
                return False
        return True
    except Exception:
        return False


def download_product_with_retries(product, destination: Path, filename: str, session,
                                  attempts: int = 3) -> Path:
    """Download one ASF NetCDF, retrying network interruptions safely."""
    target = destination / filename
    if readable_netcdf(target):
        return target
    # This target is an ASF product file in the selected download directory.
    # A non-readable file here can only be an incomplete earlier download.
    target.unlink(missing_ok=True)
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            product.download(path=str(destination), session=session)
            if readable_netcdf(target):
                return target
            raise IOError(f"ASF saved an incomplete NetCDF: {target.name}")
        except Exception as error:
            last_error = error
            target.unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(5 * attempt)
    raise RuntimeError(
        f"ASF could not download {filename} after {attempts} attempts. "
        f"Please rerun the same command; completed files will be reused. Last error: {last_error}"
    )


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
    default_downloads = root / "ASF_overlap_downloads"
    if args.asf_flight_direction:
        default_downloads /= args.asf_flight_direction
    args.asf_download_dir = args.asf_download_dir or default_downloads
    args.output = args.output or results / "asf_vs_zoomin.png"


def download_asf_granule(granule: str, destination: Path) -> Path:
    """Download one ASF catalog granule using the recipient's Earthdata login.

    ASF recommends ``.netrc`` authentication. An ``EARTHDATA_TOKEN`` environment
    variable is also accepted so neither a password nor token needs to appear in
    a command history or a shared script.
    """
    asf_search, session = asf_session()
    destination.mkdir(parents=True, exist_ok=True)
    existing = destination / Path(granule).name
    if readable_netcdf(existing):
        print(f"Using existing ASF download: {existing}")
        return existing
    results = asf_search.granule_search([granule])
    if not results:
        raise ValueError(f"ASF did not find a granule named {granule!r}.")
    # With no token, asf_search/requests uses the recipient's ~/.netrc credentials.
    return download_product_with_retries(results[0], destination, existing.name, session)


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


def opera_disp_identity(name: str) -> tuple[str, datetime, datetime]:
    """Get the OPERA frame, reference acquisition, and secondary acquisition."""
    match = OPERA_DISP_NAME.search(Path(name).name)
    if not match:
        raise ValueError(f"Cannot parse OPERA DISP frame/reference dates from {name!r}")
    return (
        match.group("frame"),
        datetime.strptime(match.group("reference"), "%Y%m%d"),
        datetime.strptime(match.group("secondary"), "%Y%m%d"),
    )


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


def download_asf_overlap(zoom_dates: list[datetime], target: GridReference, reference_granule: str | None,
                         flight_direction: str | None,
                         destination: Path, max_products: int) -> list[tuple[Path, datetime]]:
    """Download compatible OPERA DISP rasters for the ZoomInSAR/ASF time overlap.

    Products are grouped by orbit/direction so ascending and descending LOS
    observations are never mixed. A requested direction is honored; otherwise
    the best-covered group is selected. An optional reference granule can force
    a known matching group.
    """
    if max_products < 2:
        raise ValueError("--asf-max-products must be at least 2.")
    asf_search, session = asf_session()
    reference_orbit = reference_direction = None
    if reference_granule:
        reference = asf_search.granule_search([reference_granule])
        if not reference:
            raise ValueError(f"ASF did not find reference granule {reference_granule!r}.")
        reference_orbit = product_property(reference[0], "relativeOrbit", "pathNumber")
        reference_direction = product_property(reference[0], "flightDirection")
        if not reference_orbit or not reference_direction:
            raise ValueError(f"ASF reference granule {reference_granule!r} has no orbit/direction metadata.")
    start, end = min(zoom_dates), max(zoom_dates)
    results = asf_search.search(
        dataset=asf_search.DATASET.OPERA_S1,
        intersectsWith=overlap_wkt(target.region),
        start=start.strftime("%Y-%m-%dT00:00:00Z"),
        end=end.strftime("%Y-%m-%dT23:59:59Z"),
        maxResults=max_products + 1,
    )
    groups: dict[tuple[str, str, str], list] = {}
    for product in results:
        name = product_filename(product).lower()
        product_orbit = product_property(product, "relativeOrbit", "pathNumber")
        product_direction = product_property(product, "flightDirection")
        if "disp" in name and name.endswith(".nc") and product_orbit and product_direction:
            try:
                frame, _, _ = opera_disp_identity(product_filename(product))
            except ValueError:
                continue
            key = (product_orbit, product_direction.lower(), frame)
            groups.setdefault(key, []).append(product)
    if len(results) > max_products:
        raise ValueError(
            f"ASF search exceeded --asf-max-products={max_products}; narrow the map region or raise the limit."
        )
    if reference_orbit:
        matching = [item for item in groups.items()
                    if item[0][0] == reference_orbit and item[0][1] == reference_direction.lower()]
        _, selected = max(matching, key=lambda item: (len(item[1]), item[0])) if matching else (None, [])
        selection_method = f"reference granule {reference_granule}"
    elif flight_direction:
        direction_groups = [item for item in groups.items() if item[0][1] == flight_direction]
        if direction_groups:
            track, selected = max(direction_groups, key=lambda item: (len(item[1]), item[0]))
            selection_method = f"largest {flight_direction} group ({track[0]}, {track[1]})"
        else:
            selected, selection_method = [], f"no {flight_direction} group"
    elif groups:
        track, selected = max(groups.items(), key=lambda item: (len(item[1]), item[0]))
        selection_method = f"largest compatible group ({track[0]}, {track[1]})"
    else:
        selected = []
        selection_method = "none"
    # ASF may return multiple processing versions for one reference-secondary
    # pair. Keep one version, but never collapse different reference groups:
    # their displacement layers have different zero/reference epochs.
    one_per_date = {}
    for product in sorted(selected, key=product_date):
        _, reference, secondary = opera_disp_identity(product_filename(product))
        one_per_date.setdefault((reference.date(), secondary.date()), product)
    selected = list(one_per_date.values())
    if len(selected) < 2:
        raise ValueError(
            "ASF found fewer than two matching OPERA DISP NetCDF products in the overlap. "
            "A velocity/time-series comparison cannot be made for this interval."
        )
    orbit = product_property(selected[0], "relativeOrbit", "pathNumber")
    direction = product_property(selected[0], "flightDirection")
    frame, _, _ = opera_disp_identity(product_filename(selected[0]))
    destination.mkdir(parents=True, exist_ok=True)
    records: list[tuple[Path, datetime]] = []
    for product in sorted(selected, key=product_date):
        name = Path(product_filename(product)).name
        target_file = download_product_with_retries(product, destination, name, session)
        records.append((target_file, product_date(product)))
    reference_groups: dict[datetime, list[datetime]] = {}
    for path, _ in records:
        _, reference, secondary = opera_disp_identity(path.name)
        reference_groups.setdefault(reference, []).append(secondary)
    manifest = {
        "reference_granule": reference_granule,
        "track_selection": selection_method,
        "relative_orbit": orbit,
        "flight_direction": direction,
        "frame": frame,
        "zoominsar_dates": [date.strftime("%Y-%m-%d") for date in zoom_dates],
        "asf_overlap_dates": [date.strftime("%Y-%m-%d") for _, date in records],
        "reference_groups": [
            {
                "reference_date": reference.strftime("%Y-%m-%d"),
                "secondary_dates": [date.strftime("%Y-%m-%d") for date in sorted(secondaries)],
            }
            for reference, secondaries in sorted(reference_groups.items())
        ],
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


def suitable_scale_km(region: tuple[float, float, float, float]) -> float:
    """Choose a familiar scale-bar length that occupies less than half the map width."""
    west, east, south, north = region
    mid_latitude = np.deg2rad((south + north) / 2)
    width_km = abs(east - west) * 111.32 * np.cos(mid_latitude)
    maximum = width_km * 0.45
    candidates = np.array([
        0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50, 100,
    ])
    valid = candidates[candidates <= maximum]
    return float(valid[-1]) if len(valid) else float(candidates[0])


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


def representative_point_candidates(args: argparse.Namespace) -> list[tuple[float, float]]:
    """Return explicit point(s), or all ZoomInSAR representative points in rank order."""
    if (args.point_lon is None) != (args.point_lat is None):
        raise ValueError("Provide both --point-lon and --point-lat.")
    if args.point_lon is not None:
        return [(args.point_lon, args.point_lat)]
    if args.point_file:
        require_file(args.point_file, "Representative point file")
        if args.point_file.suffix.lower() != ".csv":
            return [point_from_file(args.point_file)]
        with args.point_file.open(encoding="utf-8", newline="") as stream:
            return [(float(row["longitude"]), float(row["latitude"])) for row in csv.DictReader(stream)]
    source = args.zoom_timeseries_dir / "representative_points.csv"
    require_file(source, "ZoomInSAR representative-point file")
    with source.open(encoding="utf-8", newline="") as stream:
        return [(float(row["longitude"]), float(row["latitude"])) for row in csv.DictReader(stream)]


def asf_series_at_point(stack: np.ndarray, target: GridReference,
                        point: tuple[float, float]) -> np.ndarray | None:
    row, col = rasterio.transform.rowcol(target.transform, point[0], point[1])
    if not (0 <= row < stack.shape[1] and 0 <= col < stack.shape[2]):
        return None
    values = np.asarray(stack[:, row, col], dtype=float)
    if np.count_nonzero(np.isfinite(values)) < 2:
        return None
    values -= values[np.flatnonzero(np.isfinite(values))[0]]
    return values


def select_common_representative_point(args: argparse.Namespace, asf_stack: np.ndarray,
                                       target: GridReference) -> tuple[tuple[float, float] | None, np.ndarray]:
    """Choose a ZoomInSAR representative point that also has an ASF time series.

    Failure to find a common ASF point is non-fatal: the caller can still plot
    the ZoomInSAR-only time series.
    """
    try:
        candidates = representative_point_candidates(args)
    except (FileNotFoundError, ValueError):
        candidates = []
    for point in candidates:
        series = asf_series_at_point(asf_stack, target, point)
        if series is not None:
            return point, series
    print(
        "WARNING: No ZoomInSAR representative point has at least two valid ASF observations; "
        "continuing with ZoomInSAR-only time series."
    )
    return None, np.asarray([], dtype=float)


def pixel_center(target: GridReference, row: int, col: int) -> tuple[float, float]:
    """Return longitude/latitude of a target-grid pixel center."""
    x, y = rasterio.transform.xy(target.transform, row, col, offset="center")
    return float(x), float(y)


def interactive_select_two_points(
    args: argparse.Namespace, target: GridReference
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Interactively inspect and confirm deformation and reference points.

    Stage 1: left-click as many candidates as desired, then press Enter to
    confirm the deformation point.  Stage 2: continue testing candidates and
    press Enter again to confirm the no-deformation/reference point.  Escape
    cancels the selector.  The time-series panel updates after every click.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "Interactive point selection requires matplotlib. Install it or use "
            "--point-lon/--point-lat for a non-interactive deformation point."
        ) from error

    west, east, south, north = target.region
    fig, (map_ax, ts_ax) = plt.subplots(1, 2, figsize=(13, 5.5))
    image = map_ax.imshow(
        target.data, extent=[west, east, south, north], origin="upper",
        cmap="RdBu_r", vmin=-args.vmax_mm_year, vmax=args.vmax_mm_year,
        aspect="auto",
    )
    fig.colorbar(image, ax=map_ax, label="LOS velocity (mm/year)")
    map_ax.set_xlabel("Longitude")
    map_ax.set_ylabel("Latitude")
    ts_ax.set_xlabel("Time (YYYY)")
    ts_ax.set_ylabel("Relative LOS displacement (m)")
    ts_ax.grid(alpha=0.25)

    state: dict[str, object] = {
        "stage": "deformation", "point": None,
        "deformation": None, "reference": None, "current_marker": None,
    }
    tested_x: list[float] = []
    tested_y: list[float] = []
    tested_artist = map_ax.scatter([], [], s=20, facecolors="none", edgecolors="0.35")
    deformation_artist = None

    def update_instructions() -> None:
        if state["stage"] == "deformation":
            map_ax.set_title("1/2: Test deformation points; Enter = confirm")
            ts_ax.set_title("Deformation-point candidate")
        else:
            map_ax.set_title("2/2: Test reference points; Enter = confirm")
            ts_ax.set_title("No-deformation/reference candidate")

    update_instructions()

    def on_click(event) -> None:
        if event.inaxes is not map_ax or event.xdata is None or event.ydata is None:
            return
        point = (float(event.xdata), float(event.ydata))
        try:
            dates, values = load_zoom_timeseries(
                args.zoom_timeseries_dir, args.zoom_data_dir, *point
            )
        except Exception as error:
            print(f"Candidate could not be sampled: {error}")
            return
        finite = np.isfinite(values)
        if np.count_nonzero(finite) < 2:
            print(f"Candidate lon={point[0]:.6f}, lat={point[1]:.6f} has fewer than two valid observations.")
            return

        tested_x.append(point[0]); tested_y.append(point[1])
        tested_artist.set_offsets(np.column_stack([tested_x, tested_y]))
        if state["current_marker"] is not None:
            state["current_marker"].remove()
        marker = "*" if state["stage"] == "deformation" else "o"
        face = "yellow" if state["stage"] == "deformation" else "white"
        state["current_marker"] = map_ax.scatter(
            [point[0]], [point[1]], marker=marker, s=180 if marker == "*" else 110,
            c=face, edgecolors="black", linewidths=0.9, zorder=6
        )
        state["point"] = point

        ts_ax.clear()
        ts_ax.plot(dates, values, "o-", markersize=3, linewidth=1.2)
        ts_ax.axhline(0, color="0.65", linewidth=0.8)
        role = "Deformation" if state["stage"] == "deformation" else "Reference"
        ts_ax.set_title(f"{role} candidate: {point[0]:.5f}, {point[1]:.5f}")
        ts_ax.set_xlabel("Time (YYYY)")
        ts_ax.set_ylabel("Relative LOS displacement (m)")
        ts_ax.grid(alpha=0.25)
        fig.canvas.draw_idle()
        print(f"Testing {role.lower()} candidate: lon={point[0]:.6f}, lat={point[1]:.6f}. "
              "Click another point to test it, or press Enter to confirm.")

    def on_key(event) -> None:
        nonlocal deformation_artist
        if event.key in ("enter", "return") and state["point"] is not None:
            if state["stage"] == "deformation":
                state["deformation"] = state["point"]
                point = state["deformation"]
                if state["current_marker"] is not None:
                    state["current_marker"].remove()
                    state["current_marker"] = None
                deformation_artist = map_ax.scatter(
                    [point[0]], [point[1]], marker="*", s=180, c="yellow",
                    edgecolors="black", linewidths=0.9, zorder=6,
                    label="Confirmed deformation point"
                )
                state["stage"] = "reference"
                state["point"] = None
                update_instructions()
                ts_ax.clear()
                ts_ax.set_title("Click candidates for the no-deformation/reference point")
                ts_ax.set_xlabel("Time (YYYY)")
                ts_ax.set_ylabel("Relative LOS displacement (m)")
                ts_ax.grid(alpha=0.25)
                fig.canvas.draw_idle()
                print("Deformation point confirmed. Now choose the no-deformation/reference point; "
                      "press Enter a second time to confirm it.")
            else:
                state["reference"] = state["point"]
                print("Reference point confirmed. Both selections are complete.")
                plt.close(fig)
        elif event.key == "escape":
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.tight_layout()
    plt.show()

    if state["deformation"] is None or state["reference"] is None:
        raise ValueError(
            "Both points must be confirmed: click/inspect a deformation point and press Enter, "
            "then click/inspect a no-deformation/reference point and press Enter again."
        )
    return state["deformation"], state["reference"]  # type: ignore[return-value]



def choose_stable_point(
    args: argparse.Namespace, target: GridReference, deforming: tuple[float, float]
) -> tuple[float, float]:
    """Automatically choose a stable point using low velocity and low time variability."""
    values = np.asarray(target.data, dtype=float)
    valid = np.isfinite(values)
    if np.count_nonzero(valid) < 2:
        raise ValueError("ZoomInSAR velocity map has fewer than two valid pixels.")

    drow, dcol = rasterio.transform.rowcol(target.transform, *deforming)
    # Start from the lowest-|velocity| pixels, then use temporal variability to
    # distinguish genuinely quiet locations from pixels that merely have a near-zero trend.
    flat = np.argsort(np.where(valid, np.abs(values), np.inf), axis=None)
    candidates: list[tuple[float, float, float, float]] = []
    min_sep_pixels = max(5.0, 0.05 * np.hypot(*values.shape))
    for index in flat[: min(len(flat), 1000)]:
        row, col = np.unravel_index(index, values.shape)
        if not np.isfinite(values[row, col]):
            continue
        if np.hypot(row - drow, col - dcol) < min_sep_pixels:
            continue
        point = pixel_center(target, int(row), int(col))
        try:
            _, series = load_zoom_timeseries(args.zoom_timeseries_dir, args.zoom_data_dir, *point)
        except Exception:
            continue
        finite = series[np.isfinite(series)]
        if finite.size < 2:
            continue
        variability = float(np.nanstd(finite))
        candidates.append((variability, abs(float(values[row, col])), point[0], point[1]))
        if len(candidates) >= 80:
            break
    if not candidates:
        # Conservative fallback: nearest-to-zero valid pixel distinct from deformation point.
        for index in flat:
            row, col = np.unravel_index(index, values.shape)
            if (int(row), int(col)) != (int(drow), int(dcol)):
                return pixel_center(target, int(row), int(col))
        raise ValueError("Could not find a distinct stable ZoomInSAR point.")

    # Normalize the two criteria so neither velocity nor variability dominates by units.
    var = np.asarray([item[0] for item in candidates])
    vel = np.asarray([item[1] for item in candidates])
    var_scale = max(float(np.nanpercentile(var, 90)), 1e-9)
    vel_scale = max(float(np.nanpercentile(vel, 90)), 1e-9)
    score = var / var_scale + vel / vel_scale
    best = candidates[int(np.nanargmin(score))]
    return float(best[2]), float(best[3])


def choose_zoom_deformation_and_stable_points(
    args: argparse.Namespace, target: GridReference
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Choose both ZoomInSAR points interactively, with two confirmations.

    If --point-lon/--point-lat are supplied, they remain a non-interactive
    deformation-point override and the stable point is selected automatically
    for backward-compatible batch runs.  Otherwise both points are chosen by
    the two-stage interactive selector.
    """
    if (args.point_lon is None) != (args.point_lat is None):
        raise ValueError("Provide both --point-lon and --point-lat.")
    if args.point_lon is not None:
        deforming = (float(args.point_lon), float(args.point_lat))
        print(f"Using command-line deformation point: lon={deforming[0]:.6f}, lat={deforming[1]:.6f}")
        # Preserve batch-mode behavior when an explicit deformation point is supplied.
        stable = choose_stable_point(args, target, deforming)
        return deforming, stable
    return interactive_select_two_points(args, target)

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


def load_radar_geolocation(data_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the native ZoomInSAR longitude/latitude rasters."""
    shape = load_shape(data_dir)
    geo_dir = data_dir / "GEO"
    longitude = np.fromfile(find_geo_file(geo_dir, "lon"), dtype=">f4").reshape(shape)
    latitude = np.fromfile(find_geo_file(geo_dir, "lat"), dtype=">f4").reshape(shape)
    return longitude, latitude


def geocode_radar_velocity(velocity: np.ndarray, data_dir: Path, target: GridReference) -> np.ndarray:
    """Interpolate a native-grid ZoomInSAR velocity map to the comparison GeoTIFF grid."""
    longitude, latitude = load_radar_geolocation(data_dir)
    valid = (np.isfinite(velocity) & np.isfinite(longitude) & np.isfinite(latitude)
             & (longitude != 0) & (latitude != 0))
    if np.count_nonzero(valid) < 4:
        raise ValueError("Too few valid ZoomInSAR pixels to calculate overlap-period velocity.")
    height, width = target.data.shape
    xs = target.transform.c + (np.arange(width) + 0.5) * target.transform.a
    ys = target.transform.f + (np.arange(height) + 0.5) * target.transform.e
    grid_x, grid_y = np.meshgrid(xs, ys)
    points = np.column_stack((longitude[valid], latitude[valid]))
    regular = griddata(points, velocity[valid], (grid_x, grid_y), method="linear")
    return regular.astype(np.float32)


def overlap_zoom_velocity(timeseries_dir: Path, data_dir: Path, start: datetime, end: datetime,
                           target: GridReference) -> np.ndarray:
    """Fit ZoomInSAR velocity from precisely the same date interval as ASF."""
    require_file(timeseries_dir / "cumulative_los_displacement_m.npy", "ZoomInSAR displacement stack")
    require_file(timeseries_dir / "dates.json", "ZoomInSAR dates")
    dates = [datetime.fromisoformat(value) for value in json.loads(
        (timeseries_dir / "dates.json").read_text(encoding="utf-8")
    )]
    stack = np.load(timeseries_dir / "cumulative_los_displacement_m.npy", mmap_mode="r")
    selected = [index for index, date in enumerate(dates) if start <= date <= end]
    if len(selected) < 2:
        raise ValueError(
            "Fewer than two ZoomInSAR acquisitions fall within the ASF overlap dates; "
            "cannot make a matched-period velocity comparison."
        )
    velocity = derive_asf_velocity(np.asarray(stack[selected], dtype=np.float32),
                                   [dates[index] for index in selected])
    return geocode_radar_velocity(velocity * 1000.0, data_dir, target)


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
    """Stitch OPERA DISP reference-date groups into one relative time series.

    Each ``short_wavelength_displacement`` layer is relative to the reference
    acquisition embedded in its own filename. OPERA changes that reference
    periodically. Regressing raw layers from several groups creates artificial
    jumps and an invalid velocity, so each group is added to the already
    stitched grid at its reference date.
    """
    grouped: dict[datetime, list[tuple[datetime, Path]]] = {}
    for path, _ in records:
        _, reference, secondary = opera_disp_identity(path.name)
        grouped.setdefault(reference, []).append((secondary, path))
    if not grouped:
        raise ValueError("No OPERA DISP products were available to build an ASF time series.")

    stitched: dict[datetime, np.ndarray] = {}
    first_reference = min(grouped)
    for reference in sorted(grouped):
        if reference == first_reference:
            base = np.zeros(target.data.shape, dtype=np.float32)
            stitched[reference] = base
        elif reference in stitched:
            base = stitched[reference]
        else:
            raise ValueError(
                "ASF OPERA reference-date chain is discontinuous at "
                f"{reference:%Y-%m-%d}; cannot safely combine these displacement groups. "
                "Delete the incomplete ASF_overlap_downloads folder and rerun the command."
            )
        for secondary, path in sorted(grouped[reference]):
            if secondary in stitched:
                raise ValueError(
                    f"ASF has more than one stitched displacement for {secondary:%Y-%m-%d}; "
                    "choose one OPERA frame only."
                )
            layer = choose_subdataset(path, requested, "displacement")
            relative = reproject_to_zoom(layer, target) * scale
            stitched[secondary] = (base + relative).astype(np.float32)
    dates = sorted(stitched)
    return dates, np.stack([stitched[date] for date in dates], axis=0)


def restrict_to_date_range(dates: list[datetime], values: np.ndarray, start: datetime,
                           end: datetime) -> tuple[list[datetime], np.ndarray]:
    """Keep a point series within the same time interval used for both maps."""
    selected = [index for index, date in enumerate(dates) if start <= date <= end]
    if len(selected) < 2:
        raise ValueError("Fewer than two ZoomInSAR dates fall within the ASF overlap interval.")
    clipped_dates = [dates[index] for index in selected]
    clipped_values = np.asarray(values, dtype=float)[selected]
    finite = np.flatnonzero(np.isfinite(clipped_values))
    if len(finite):
        clipped_values = clipped_values - clipped_values[finite[0]]
    return clipped_dates, clipped_values


def plot_map(fig: pygmt.Figure, dem: xr.DataArray, velocity: xr.DataArray, region: tuple[float, ...],
             projection: str, cpt: str, title: str, points: list[tuple[float, float, str]] | None,
             transparency: float, scale_km: float) -> None:
    shade = pygmt.grdgradient(grid=dem, azimuth=315, normalize="t1")
    fig.grdimage(grid=dem, region=region, projection=projection, cmap="gray", shading=shade,
                 frame=["af", f"+t{title}"])
    fig.grdimage(grid=velocity, cmap=cpt, transparency=transparency)
    if points:
        for longitude, latitude, role in points:
            if role == "Deformation":
                fig.plot(x=longitude, y=latitude, style="a0.30c", fill="yellow", pen="0.75p,black")
            else:
                fig.plot(x=longitude, y=latitude, style="c0.25c", fill="white", pen="0.75p,black")
    # Put the scale bar in the upper-left blank map area, away from the colorbar.
    fig.basemap(map_scale=f"jTL+o0.25c/0.75c+w{scale_km:g}k+f+lkm", rose="jTR+w1.2c+f2+l")


def plot_zoom_timeseries(fig: pygmt.Figure, zoom_dates: list[datetime],
                         deformation_values: np.ndarray, stable_values: np.ndarray,
                         asf_dates: list[datetime] | None = None,
                         asf_values: np.ndarray | None = None) -> None:
    """Plot matched ASF/OPERA and ZoomInSAR displacement at one common point."""
    pieces = [deformation_values, stable_values]
    if asf_values is not None:
        pieces.append(asf_values)
    all_values = np.concatenate(pieces)
    finite = all_values[np.isfinite(all_values)]
    if not len(finite):
        raise ValueError("Both selected ZoomInSAR time series contain only NaN values.")
    span = max(0.01, float(np.nanmax(finite) - np.nanmin(finite)))
    ymin = float(np.nanmin(finite) - span * 0.15)
    ymax = float(np.nanmax(finite) + span * 0.15)
    region = [min(zoom_dates).strftime("%Y-%m-%d"), max(zoom_dates).strftime("%Y-%m-%d"), ymin, ymax]
    fig.basemap(
        region=region, projection="X25c/7c",
        frame=["pxa1Yf3o+lTime (YYYY)", "ya+lRelative LOS displacement (m)",
               "+t(c) Matched representative-point LOS displacement time series"],
    )
    if asf_dates and asf_values is not None:
        fig.plot(x=asf_dates, y=asf_values, pen="1.3p,black",
                 style="c0.10c", fill="black", label="ASF/OPERA")
    fig.plot(x=zoom_dates, y=deformation_values, pen="1.3p,180/25/25",
             style="c0.11c", fill="180/25/25", label="ZoomInSAR")
    fig.plot(x=zoom_dates, y=stable_values, pen="0.9p,90/90/90",
             style="c0.09c", fill="white", label="ZoomInSAR stable point")
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
            and (not args.zoom_timeseries_dir or not args.zoom_data_dir)):
        raise ValueError(
            "Panel (c) requires --zoom-timeseries-dir and --zoom-data-dir. "
            "ASF time-series files are optional."
        )

    zoom = read_zoom_grid(args.zoom_velocity)
    require_file(args.dem, "DEM GeoTIFF")
    dem_name = str(args.dem)
    dem = to_dataarray(reproject_to_zoom(dem_name, zoom), zoom, "dem")
    zoom_velocity = to_dataarray(zoom.data, zoom, "zoomin_velocity_mm_year")
    deformation_point = None
    stable_point = None

    zoom_dates: list[datetime] = []
    zoom_deformation_series = np.asarray([])
    zoom_stable_series = np.asarray([])
    asf_dates: list[datetime] = []
    asf_series = np.asarray([])
    if args.zoom_timeseries_dir:
        require_file(args.zoom_timeseries_dir / "dates.json", "ZoomInSAR dates")
        zoom_dates = [datetime.fromisoformat(value) for value in json.loads(
            (args.zoom_timeseries_dir / "dates.json").read_text(encoding="utf-8")
        )]
    if args.auto_asf_overlap:
        records = download_asf_overlap(
            zoom_dates, zoom, args.asf_reference_granule, args.asf_flight_direction,
            args.asf_download_dir, args.asf_max_products
        )
        asf_dates, asf_stack = load_auto_asf_stack(
            records, zoom, args.asf_displacement_subdataset, args.asf_displacement_scale
        )
        asf_velocity = to_dataarray(
            derive_asf_velocity(asf_stack, asf_dates) * 1000.0, zoom, "asf_velocity_mm_year"
        )
        zoom_velocity = to_dataarray(
            overlap_zoom_velocity(args.zoom_timeseries_dir, args.zoom_data_dir,
                                  min(asf_dates), max(asf_dates), zoom),
            zoom, "zoomin_velocity_overlap_mm_year"
        )
        if not args.no_timeseries:
            common_point, asf_series = select_common_representative_point(args, asf_stack, zoom)
            if common_point is not None:
                # This guarantees that the star in both maps and both curves in
                # panel (c) refer to exactly the same physical pixel/location.
                deformation_point = common_point
                stable_point = choose_stable_point(args, zoom, deformation_point)
                print("Using a common ZoomInSAR/ASF representative point for the comparison.")
            else:
                deformation_point, stable_point = choose_zoom_deformation_and_stable_points(args, zoom)
            zoom_dates, zoom_deformation_series = load_zoom_timeseries(
                args.zoom_timeseries_dir, args.zoom_data_dir, *deformation_point
            )
            stable_dates, zoom_stable_series = load_zoom_timeseries(
                args.zoom_timeseries_dir, args.zoom_data_dir, *stable_point
            )
            if stable_dates != zoom_dates:
                raise ValueError("ZoomInSAR time-series dates differ between selected points.")
            overlap_start, overlap_end = min(asf_dates), max(asf_dates)
            zoom_dates, zoom_deformation_series = restrict_to_date_range(
                zoom_dates, zoom_deformation_series, overlap_start, overlap_end
            )
            stable_dates, zoom_stable_series = restrict_to_date_range(
                stable_dates, zoom_stable_series, overlap_start, overlap_end
            )
            if stable_dates != zoom_dates:
                raise ValueError("ZoomInSAR clipped time-series dates differ between selected points.")
    else:
        asf_velocity_layer = choose_subdataset(args.asf_nc, args.asf_velocity_subdataset, "velocity")
        asf_velocity = to_dataarray(reproject_to_zoom(asf_velocity_layer, zoom) * args.asf_velocity_scale,
                                    zoom, "asf_velocity_mm_year")
        if not args.no_timeseries:
            deformation_point, stable_point = choose_zoom_deformation_and_stable_points(args, zoom)
            zoom_dates, zoom_deformation_series = load_zoom_timeseries(
                args.zoom_timeseries_dir, args.zoom_data_dir, *deformation_point
            )
            stable_dates, zoom_stable_series = load_zoom_timeseries(
                args.zoom_timeseries_dir, args.zoom_data_dir, *stable_point
            )
            if stable_dates != zoom_dates:
                raise ValueError("ZoomInSAR time-series dates differ between selected points.")
            # Keep optional ASF sampling for compatibility/diagnostics, but never require it
            # for the ZoomInSAR time-series figure.
            if args.asf_timeseries_nc:
                try:
                    asf_dates, asf_series = sample_asf_displacement(
                        args.asf_timeseries_nc, zoom, *deformation_point,
                        args.asf_displacement_subdataset
                    )
                except ValueError as error:
                    print(f"WARNING: ASF point time series unavailable ({error}); continuing with ZoomInSAR only.")

    with tempfile.TemporaryDirectory(prefix="zoomin_asf_") as temporary:
        cpt = Path(temporary) / "velocity.cpt"
        write_velocity_cpt(cpt, args.vmax_mm_year)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        asf_map_output = args.output.with_name(f"{args.output.stem}_asf_velocity{args.output.suffix}")
        zoom_map_output = args.output.with_name(f"{args.output.stem}_zoomin_velocity{args.output.suffix}")
        scale_km = suitable_scale_km(zoom.region)
        with pygmt.config(
            FONT_TITLE="18p,Helvetica-Bold",
            FONT_LABEL="14p,Helvetica", FONT_ANNOT_PRIMARY="12p,Helvetica",
        ):
            projection = "M13c"
            asf_figure = pygmt.Figure()
            asf_map_points = None
            if deformation_point:
                asf_map_points = [(*deformation_point, "Deformation")]
            plot_map(asf_figure, dem, asf_velocity, zoom.region, projection, str(cpt),
                     "(a) ASF/OPERA LOS velocity", asf_map_points, args.transparency, scale_km)
            # Place the colorbar fully below the map frame in dedicated blank space.
            asf_figure.colorbar(position="JBC+w9c/0.45c+o0c/-2.0c+h", cmap=str(cpt),
                                frame=["xaf+lLOS velocity (mm/year)"])
            asf_figure.savefig(asf_map_output, dpi=300)

            zoom_figure = pygmt.Figure()
            # Show both time-series locations only on the ZoomInSAR map.
            zoom_map_points = None
            if deformation_point and stable_point:
                zoom_map_points = [(*deformation_point, "Deformation"), (*stable_point, "Stable")]
            plot_map(zoom_figure, dem, zoom_velocity, zoom.region, projection, str(cpt),
                     "(b) ZoomInSAR LOS velocity", zoom_map_points, args.transparency, scale_km)
            # Place the colorbar fully below the map frame in dedicated blank space.
            zoom_figure.colorbar(position="JBC+w9c/0.45c+o0c/-2.0c+h", cmap=str(cpt),
                                 frame=["xaf+lLOS velocity (mm/year)"])
            zoom_figure.savefig(zoom_map_output, dpi=300)
            if deformation_point and stable_point:
                timeseries_output = args.output.with_name(
                    f"{args.output.stem}_timeseries{args.output.suffix}"
                )
                timeseries_figure = pygmt.Figure()
                plot_zoom_timeseries(
                    timeseries_figure, zoom_dates,
                    zoom_deformation_series, zoom_stable_series,
                    asf_dates if len(asf_series) else None,
                    asf_series if len(asf_series) else None,
                )
                timeseries_figure.savefig(timeseries_output, dpi=300)
    print(f"Wrote {asf_map_output}")
    print(f"Wrote {zoom_map_output}")
    if deformation_point and stable_point:
        print(f"ZoomInSAR deformation point: lon={deformation_point[0]:.6f}, lat={deformation_point[1]:.6f}")
        print(f"ZoomInSAR stable point:      lon={stable_point[0]:.6f}, lat={stable_point[1]:.6f}")
        print(f"Wrote {timeseries_output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # Provide a short actionable CLI error rather than a GMT traceback.
        raise SystemExit(f"ERROR: {error}") from error
