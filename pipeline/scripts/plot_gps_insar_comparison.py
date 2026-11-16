#!/usr/bin/env python3
"""Compare Portuguese Bend survey GPS displacements with ZoomInSAR LOS.

The Portuguese Bend movement posting is an Excel report, rather than a simple
time-series table.  This utility reads every dated ``... Positions`` block,
extracts the North/East/Elevation coordinates for one monument, and converts
coordinate changes from US survey feet to metres.  It then projects the GPS
East/North/Up (ENU) displacement onto a user-supplied satellite look vector.

The look vector must point from the ground *toward the satellite*.  With that
definition, positive projected GPS displacement is positive LOS displacement
(motion toward the satellite), matching the ZoomInSAR convention used by this
repository.  Do not compare a GPS vertical or horizontal displacement directly
with an InSAR LOS time series.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from rasterio.warp import transform


US_SURVEY_FOOT_TO_METRE = 1200.0 / 3937.0
XLSX_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


@dataclass(frozen=True)
class GpsObservation:
    """One survey coordinate solution from the movement-posting workbook."""

    date: datetime
    station: str
    north_ft: float
    east_ft: float
    elevation_ft: float


def excel_column_number(reference: str) -> int:
    """Return a one-based Excel column number from a cell reference."""
    result = 0
    for character in "".join(char for char in reference if char.isalpha()):
        result = result * 26 + ord(character.upper()) - ord("A") + 1
    return result


def read_xlsx_cells(path: Path) -> dict[int, dict[int, str]]:
    """Read plain values from the populated report worksheet without openpyxl."""
    with zipfile.ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", XLSX_NS):
                shared_strings.append("".join(text.text or "" for text in item.findall(".//m:t", XLSX_NS)))

        # The source report's populated sheet is normally sheet2.  Select the
        # sheet with the largest XML payload so this also works if Excel names
        # the sheets differently in a later posting.
        sheets = [name for name in archive.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", name)]
        if not sheets:
            raise ValueError("The workbook contains no worksheets")
        sheet_name = max(sheets, key=lambda name: archive.getinfo(name).file_size)
        root = ET.fromstring(archive.read(sheet_name))

    rows: dict[int, dict[int, str]] = {}
    for row in root.findall(".//m:sheetData/m:row", XLSX_NS):
        values: dict[int, str] = {}
        for cell in row.findall("m:c", XLSX_NS):
            raw = cell.find("m:v", XLSX_NS)
            if raw is None or raw.text is None:
                continue
            value = raw.text
            if cell.attrib.get("t") == "s":
                value = shared_strings[int(value)]
            values[excel_column_number(cell.attrib["r"])] = value.strip()
        if values:
            rows[int(row.attrib["r"])] = values
    return rows


def parse_report_date(label: str) -> datetime | None:
    """Parse dates in labels such as ``September 3, 2025 Positions``."""
    cleaned = label.replace("Septemebr", "September").replace("Sept.", "Sep")
    match = re.search(r"([A-Za-z]+)\s+(\d{1,2}),?\s+(\d{4})", cleaned)
    if not match:
        return None
    candidate = " ".join(match.groups())
    for format_string in ("%B %d %Y", "%b %d %Y"):
        try:
            return datetime.strptime(candidate, format_string)
        except ValueError:
            pass
    return None


def numeric(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def parse_portuguese_bend_posting(path: Path) -> dict[str, list[GpsObservation]]:
    """Extract all stations' dated North/East/Elevation coordinate solutions."""
    rows = read_xlsx_cells(path)
    header = rows.get(11, {})
    position_columns: list[tuple[int, datetime]] = []
    for column, value in header.items():
        if not value.lower().endswith("positions"):
            continue
        date = parse_report_date(value)
        # The 2007 sheet says ``Sept. 24, 2007 Positions`` and later pages use
        # the same pattern.  The following columns are North, East, Elevation.
        if date is not None:
            position_columns.append((column, date))
    if not position_columns:
        raise ValueError("Could not find dated coordinate blocks in workbook row 11")

    observations_by_station: dict[str, list[GpsObservation]] = {}
    for north_column, date in position_columns:
        station_column = north_column - 1
        for row_number, row in rows.items():
            if row_number < 14:
                continue
            station = row.get(station_column, "").strip().upper()
            if not station:
                continue
            north_ft = numeric(row.get(north_column))
            east_ft = numeric(row.get(north_column + 1))
            elevation_ft = numeric(row.get(north_column + 2))
            if None not in (north_ft, east_ft, elevation_ft):
                observations_by_station.setdefault(station, []).append(
                    GpsObservation(date, station, north_ft, east_ft, elevation_ft)
                )

    # Later report pages occasionally repeat a survey date.  Keep the last
    # complete coordinate solution for that date, then sort chronologically.
    result: dict[str, list[GpsObservation]] = {}
    for station, observations in observations_by_station.items():
        by_date = {observation.date: observation for observation in observations}
        if len(by_date) >= 2:
            result[station] = sorted(by_date.values(), key=lambda observation: observation.date)
    if not result:
        raise ValueError("No station has at least two usable survey observations")
    return result


def find_single(directory: Path, suffix: str) -> Path:
    matches = sorted(directory.glob(f"*.{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one *.{suffix} file in {directory}; found {len(matches)}")
    return matches[0]


def resolve_dataset(project_root: Path, dataset: str | None) -> tuple[Path, Path, str]:
    """Return results root, time-series directory, and Dataset_* name."""
    results_root = project_root / "ZoomInSAR_Results"
    output_root = results_root / "InSAR_Timeseries" / "OUTPUT"
    if not output_root.is_dir():
        raise FileNotFoundError(f"ZoomInSAR time-series output directory not found: {output_root}")
    if dataset:
        names = [dataset if dataset.startswith("Dataset_") else f"Dataset_{dataset}"]
    else:
        names = [path.name for path in sorted(output_root.glob("Dataset_*")) if path.is_dir()]
    existing = [name for name in names if (output_root / name).is_dir()]
    if len(existing) != 1:
        available = ", ".join(path.name.removeprefix("Dataset_") for path in sorted(output_root.glob("Dataset_*")))
        raise ValueError(f"Select exactly one dataset with --dataset. Available: {available or 'none'}")
    dataset_name = existing[0]
    return results_root, output_root / dataset_name, dataset_name


def valid_geo(longitude: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    return (
        np.isfinite(longitude) & np.isfinite(latitude)
        & (longitude >= -180) & (longitude <= 180)
        & (latitude >= -90) & (latitude <= 90)
        & ~((longitude == 0) & (latitude == 0))
    )


def load_insar_grid(project_root: Path, dataset: str | None) -> dict:
    """Open the InSAR grid once so multiple GPS stations share it."""
    results_root, output_dir, dataset_name = resolve_dataset(project_root, dataset)
    velocity = np.load(output_dir / "mean_los_velocity_m_per_year.npy", mmap_mode="r")
    shape = tuple(velocity.shape)
    dates = [datetime.fromisoformat(value) for value in json.loads((output_dir / "dates.json").read_text())]
    stack = np.load(output_dir / "cumulative_los_displacement_m.npy", mmap_mode="r")
    if stack.shape != (len(dates), *shape):
        raise ValueError("InSAR cumulative stack dimensions do not match dates or velocity grid")

    geo_dir = results_root / "Data" / "ROI" / dataset_name / "GEO"
    lon = np.memmap(find_single(geo_dir, "lon"), dtype=np.dtype(">f4"), mode="r", shape=shape)
    lat = np.memmap(find_single(geo_dir, "lat"), dtype=np.dtype(">f4"), mode="r", shape=shape)
    valid = valid_geo(lon, lat)
    if not np.any(valid):
        raise ValueError("No valid longitude/latitude values in the InSAR GEO directory")
    return {
        "dataset": dataset_name,
        "velocity_m_per_year": velocity,
        "longitude_grid": lon,
        "latitude_grid": lat,
        "valid_geo": valid,
        "shape": shape,
        "dates": np.array(dates, dtype="datetime64[ns]"),
        "stack": stack,
    }


def load_insar_point(insar_grid: dict, longitude: float, latitude: float) -> dict:
    """Extract the nearest valid InSAR time series from an opened grid."""
    lon = insar_grid["longitude_grid"]
    lat = insar_grid["latitude_grid"]
    valid = insar_grid["valid_geo"]
    shape = insar_grid["shape"]
    distance_squared = np.where(valid, (lon - longitude) ** 2 + (lat - latitude) ** 2, np.inf)
    row, column = np.unravel_index(np.argmin(distance_squared), shape)
    distance_m = math.sqrt(float(distance_squared[row, column])) * 111_000.0
    series_m = np.asarray(insar_grid["stack"][:, row, column], dtype=float)
    if not np.any(np.isfinite(series_m)):
        raise ValueError("The nearest valid geographic InSAR pixel has no finite time-series values")
    return {
        "dataset": insar_grid["dataset"],
        "row": int(row),
        "column": int(column),
        "longitude": float(lon[row, column]),
        "latitude": float(lat[row, column]),
        "distance_m": distance_m,
        "dates": insar_grid["dates"],
        "displacement_m": series_m,
    }


def dates_as_days(dates: Iterable[datetime] | np.ndarray) -> np.ndarray:
    return np.asarray([(np.datetime64(date, "ns") - np.datetime64("1970-01-01", "ns")) / np.timedelta64(1, "D") for date in dates], dtype=float)


def interpolation_values(source_days: np.ndarray, source_values: np.ndarray, target_days: np.ndarray) -> np.ndarray:
    """Linearly interpolate only within the observed time span; no extrapolation."""
    result = np.full(target_days.shape, np.nan, dtype=float)
    within = (target_days >= source_days[0]) & (target_days <= source_days[-1])
    result[within] = np.interp(target_days[within], source_days, source_values)
    return result


def plot_comparison(
    output_path: Path,
    station: str,
    gps_dates: np.ndarray,
    gps_enu_m: np.ndarray,
    gps_los_m: np.ndarray,
    insar_dates: np.ndarray,
    insar_los_m: np.ndarray,
    insar_longitude: float,
    insar_latitude: float,
    nearest_distance_m: float,
) -> None:
    """Draw GPS components and the like-for-like LOS comparison."""
    fig, (components, comparison) = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True, layout="constrained")
    gps_dates_plot = gps_dates.astype("datetime64[ms]").astype(object)
    insar_dates_plot = insar_dates.astype("datetime64[ms]").astype(object)
    components.plot(gps_dates_plot, gps_enu_m[:, 0] * 1000, "o-", label="GPS east", color="#b33c3c")
    components.plot(gps_dates_plot, gps_enu_m[:, 1] * 1000, "o-", label="GPS north", color="#2275a5")
    components.plot(gps_dates_plot, gps_enu_m[:, 2] * 1000, "o-", label="GPS up", color="#317a45")
    components.axhline(0, color="0.45", linewidth=0.8)
    components.set_ylabel("GPS displacement (mm)")
    components.set_title(f"GPS ENU displacement: {station}")
    components.legend(ncol=3, frameon=False)
    components.grid(alpha=0.25)

    comparison.plot(gps_dates_plot, gps_los_m * 1000, "o-", color="#303030", label="GPS projected to LOS")
    comparison.plot(insar_dates_plot, insar_los_m * 1000, ".-", color="#c62828", label="ZoomInSAR LOS")
    comparison.axhline(0, color="0.45", linewidth=0.8)
    comparison.set_ylabel("Relative LOS displacement (mm)")
    comparison.set_title(
        "LOS comparison (+ = toward satellite) | "
        f"nearest InSAR pixel: {insar_longitude:.5f}, {insar_latitude:.5f} ({nearest_distance_m:.0f} m)"
    )
    comparison.legend(frameon=False)
    comparison.grid(alpha=0.25)
    comparison.xaxis.set_major_locator(mdates.AutoDateLocator())
    comparison.xaxis.set_major_formatter(mdates.ConciseDateFormatter(comparison.xaxis.get_major_locator()))
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_station_velocity_map(output_path: Path, insar_grid: dict, stations: dict[str, tuple[float, float]]) -> None:
    """Plot mean LOS velocity and all GPS station locations in geographic space."""
    longitude = insar_grid["longitude_grid"]
    latitude = insar_grid["latitude_grid"]
    velocity = np.asarray(insar_grid["velocity_m_per_year"], dtype=float) * 1000.0
    valid = insar_grid["valid_geo"] & np.isfinite(velocity)
    if not np.any(valid):
        raise ValueError("No finite mean LOS velocity values for the station map")
    lon_values, lat_values, velocity_values = longitude[valid], latitude[valid], velocity[valid]
    # Keep the figure responsive for very large results while retaining a
    # regular geographic sample; this does not alter any station comparison.
    if lon_values.size > 180_000:
        sample = np.linspace(0, lon_values.size - 1, 180_000).astype(int)
        lon_values, lat_values, velocity_values = lon_values[sample], lat_values[sample], velocity_values[sample]
    limit = max(float(np.nanpercentile(np.abs(velocity_values), 98)), 1.0)
    fig, axis = plt.subplots(figsize=(10, 8), layout="constrained")
    points = axis.scatter(
        lon_values, lat_values, c=velocity_values, s=5, cmap="RdBu_r", vmin=-limit, vmax=limit,
        marker="s", linewidths=0, rasterized=True,
    )
    for station, (longitude_value, latitude_value) in sorted(stations.items()):
        axis.scatter(longitude_value, latitude_value, s=62, facecolor="#e31a1c", edgecolor="white", linewidth=1.3, zorder=3)
        axis.annotate(station, (longitude_value, latitude_value), xytext=(5, 5), textcoords="offset points",
                      fontsize=7, color="black", weight="bold", zorder=4,
                      bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 0.5})
    colorbar = fig.colorbar(points, ax=axis, shrink=0.84, pad=0.02)
    colorbar.set_label("Mean LOS velocity (mm/year; + toward satellite)")
    axis.set_title(f"ZoomInSAR mean LOS velocity with GPS stations ({insar_grid['dataset'].removeprefix('Dataset_')})")
    axis.set_xlabel("Longitude")
    axis.set_ylabel("Latitude")
    axis.set_aspect("equal", adjustable="box")
    axis.grid(alpha=0.2)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gps-xlsx", required=True, type=Path, help="Portuguese Bend movement-posting workbook")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--station", help="One GPS monitoring-point identifier, e.g. AB50")
    selection.add_argument("--all-stations", action="store_true", help="Process every station with at least two usable survey observations")
    parser.add_argument("--project-root", required=True, type=Path, help="Project containing ZoomInSAR_Results")
    parser.add_argument("--dataset", help="Dataset name without Dataset_; auto-detected only if one is present")
    parser.add_argument("--output-dir", type=Path, help="Default: <project-root>/ZoomInSAR_Results/GPS_InSAR_Comparison")
    parser.add_argument("--gps-crs", default="EPSG:2229", help="GPS coordinate CRS; the supplied posting states NAD83 CA Zone 5 (ft), EPSG:2229")
    parser.add_argument("--look-east", type=float, help="Ground-to-satellite LOS unit-vector east component")
    parser.add_argument("--look-north", type=float, help="Ground-to-satellite LOS unit-vector north component")
    parser.add_argument("--look-up", type=float, help="Ground-to-satellite LOS unit-vector up component")
    parser.add_argument("--incidence-deg", type=float, help="Alternative to components: local incidence angle from vertical")
    parser.add_argument("--look-azimuth-deg", type=float, help="Alternative to components: ground-to-satellite look azimuth, clockwise from north")
    parser.add_argument("--reverse-insar-sign", action="store_true", help="Reverse InSAR LOS only if that product defines positive as away from satellite")
    return parser.parse_args()


def look_vector(args: argparse.Namespace) -> np.ndarray:
    components = (args.look_east, args.look_north, args.look_up)
    angle_values = (args.incidence_deg, args.look_azimuth_deg)
    component_mode = any(value is not None for value in components)
    angle_mode = any(value is not None for value in angle_values)
    if component_mode and angle_mode:
        raise ValueError("Use either --look-east/--look-north/--look-up or --incidence-deg with --look-azimuth-deg, not both")
    if component_mode:
        if any(value is None for value in components):
            raise ValueError("Provide all three LOS components: --look-east, --look-north, --look-up")
        vector = np.asarray(components, dtype=float)
    elif angle_mode:
        if any(value is None for value in angle_values):
            raise ValueError("Provide both --incidence-deg and --look-azimuth-deg")
        incidence = math.radians(args.incidence_deg)
        azimuth = math.radians(args.look_azimuth_deg)
        vector = np.array([
            math.sin(incidence) * math.sin(azimuth),
            math.sin(incidence) * math.cos(azimuth),
            math.cos(incidence),
        ])
    else:
        raise ValueError(
            "A real LOS geometry is required. Supply --look-east/--look-north/--look-up "
            "from the product metadata, or --incidence-deg and --look-azimuth-deg."
        )
    magnitude = float(np.linalg.norm(vector))
    if not math.isfinite(magnitude) or magnitude == 0:
        raise ValueError("LOS vector must be finite and non-zero")
    return vector / magnitude


def process_station(
    station: str,
    observations: list[GpsObservation],
    args: argparse.Namespace,
    vector: np.ndarray,
    insar_grid: dict,
    output_dir: Path,
) -> tuple[float, float]:
    """Write one station's comparison products and return its longitude/latitude."""
    north_ft = np.array([item.north_ft for item in observations])
    east_ft = np.array([item.east_ft for item in observations])
    up_ft = np.array([item.elevation_ft for item in observations])
    gps_dates = np.array([item.date for item in observations], dtype="datetime64[ns]")
    # EPSG:2229 is a US-survey-foot CRS. Rasterio performs its native unit
    # conversion before returning WGS84 longitude/latitude.
    station_longitude, station_latitude = transform(args.gps_crs, "EPSG:4326", [east_ft[0]], [north_ft[0]])
    insar = load_insar_point(insar_grid, station_longitude[0], station_latitude[0])

    gps_enu_m = np.column_stack((
        (east_ft - east_ft[0]) * US_SURVEY_FOOT_TO_METRE,
        (north_ft - north_ft[0]) * US_SURVEY_FOOT_TO_METRE,
        (up_ft - up_ft[0]) * US_SURVEY_FOOT_TO_METRE,
    ))
    gps_los_m = gps_enu_m @ vector
    insar_days = dates_as_days(insar["dates"])
    gps_days = dates_as_days(gps_dates)
    insar_values = np.asarray(insar["displacement_m"], dtype=float)
    finite_insar = np.isfinite(insar_values)
    if finite_insar.sum() < 2:
        raise ValueError("Selected InSAR pixel has fewer than two finite values")

    overlap_start = max(gps_days[0], insar_days[finite_insar][0])
    overlap_end = min(gps_days[-1], insar_days[finite_insar][-1])
    if overlap_end <= overlap_start:
        raise ValueError("GPS and InSAR time spans do not overlap")
    gps_zero = interpolation_values(gps_days, gps_los_m, np.array([overlap_start]))[0]
    insar_zero = interpolation_values(insar_days[finite_insar], insar_values[finite_insar], np.array([overlap_start]))[0]
    if not (math.isfinite(gps_zero) and math.isfinite(insar_zero)):
        raise ValueError("Cannot establish a common reference date within the data overlap")
    gps_los_relative = gps_los_m - gps_zero
    # Re-reference all three ENU components at the same comparison epoch.
    gps_enu_relative = gps_enu_m - np.array([
        interpolation_values(gps_days, gps_enu_m[:, index], np.array([overlap_start]))[0]
        for index in range(3)
    ])
    selected_insar = finite_insar & (insar_days >= overlap_start) & (insar_days <= overlap_end)
    insar_dates = insar["dates"][selected_insar]
    insar_los_relative = insar_values[selected_insar] - insar_zero
    if args.reverse_insar_sign:
        insar_los_relative *= -1.0

    insar_sample_days = insar_days[selected_insar]
    gps_at_insar = interpolation_values(gps_days, gps_los_relative, insar_sample_days)
    residual = insar_los_relative - gps_at_insar
    matched = np.isfinite(residual)
    rmse_mm = float(np.sqrt(np.mean(residual[matched] ** 2)) * 1000.0) if np.any(matched) else None

    safe_station = re.sub(r"[^A-Za-z0-9_.-]+", "_", station)
    figure_path = output_dir / f"gps_insar_{safe_station}.png"
    csv_path = output_dir / f"gps_projected_los_{safe_station}.csv"
    summary_path = output_dir / f"gps_insar_{safe_station}_summary.json"

    plot_comparison(
        figure_path, args.station.strip().upper(), gps_dates, gps_enu_relative, gps_los_relative,
        insar_dates, insar_los_relative, insar["longitude"], insar["latitude"], insar["distance_m"],
    )
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["date", "east_m", "north_m", "up_m", "gps_los_m"])
        writer.writeheader()
        for date, enu, los in zip(gps_dates, gps_enu_relative, gps_los_relative):
            writer.writerow({
                "date": str(date)[:10], "east_m": f"{enu[0]:.6f}", "north_m": f"{enu[1]:.6f}",
                "up_m": f"{enu[2]:.6f}", "gps_los_m": f"{los:.6f}",
            })
    summary = {
        "station": station, "gps_workbook": str(args.gps_xlsx.resolve()),
        "gps_crs": args.gps_crs, "gps_station_longitude": station_longitude[0], "gps_station_latitude": station_latitude[0],
        "los_unit_vector_ground_to_satellite_ENU": {"east": vector[0], "north": vector[1], "up": vector[2]},
        "positive_los_direction": "toward satellite", "comparison_start": str(np.datetime64(int(overlap_start), "D")),
        "comparison_end": str(np.datetime64(int(overlap_end), "D")), "dataset": insar["dataset"],
        "insar_pixel": {key: insar[key] for key in ("row", "column", "longitude", "latitude", "distance_m")},
        "matched_insar_dates": int(np.sum(matched)), "gps_insar_rmse_mm_at_interpolated_gps_dates": rmse_mm,
        "reverse_insar_sign": bool(args.reverse_insar_sign),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"GPS observations: {len(observations)} ({gps_dates[0]} to {gps_dates[-1]})")
    print(f"InSAR dataset/pixel: {insar['dataset']} row {insar['row']}, column {insar['column']}; {insar['distance_m']:.0f} m from GPS station")
    print(f"Wrote figure: {figure_path}")
    print(f"Wrote GPS LOS table: {csv_path}")
    print(f"Wrote summary: {summary_path}")
    return station_longitude[0], station_latitude[0]


def main() -> None:
    args = parse_args()
    vector = look_vector(args)
    observations_by_station = parse_portuguese_bend_posting(args.gps_xlsx)
    if args.station:
        station = args.station.strip().upper()
        if station not in observations_by_station:
            available = ", ".join(sorted(observations_by_station))
            raise ValueError(f"No usable observations found for {station}. Available stations: {available}")
        selected = {station: observations_by_station[station]}
    else:
        selected = observations_by_station

    results_root, _, _ = resolve_dataset(args.project_root, args.dataset)
    output_dir = args.output_dir or results_root / "GPS_InSAR_Comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    insar_grid = load_insar_grid(args.project_root, args.dataset)
    station_locations: dict[str, tuple[float, float]] = {}
    failures: dict[str, str] = {}
    for station, observations in sorted(selected.items()):
        try:
            station_locations[station] = process_station(station, observations, args, vector, insar_grid, output_dir)
        except ValueError as error:
            # A station can legitimately fall outside the InSAR time span or
            # on an invalid pixel. Continue through the remaining stations.
            failures[station] = str(error)
            print(f"Skipped {station}: {error}")

    if not station_locations:
        raise ValueError("No station produced a GPS–InSAR comparison")
    station_map = output_dir / "gps_stations_on_mean_los_velocity.png"
    plot_station_velocity_map(station_map, insar_grid, station_locations)
    print(f"Wrote station map: {station_map}")
    if args.all_stations:
        (output_dir / "all_stations_summary.json").write_text(json.dumps({
            "processed_stations": sorted(station_locations), "skipped_stations": failures,
            "station_locations": {name: {"longitude": point[0], "latitude": point[1]} for name, point in station_locations.items()},
        }, indent=2) + "\n", encoding="utf-8")
        print(f"Completed {len(station_locations)} stations; skipped {len(failures)} stations")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, KeyError, zipfile.BadZipFile) as error:
        sys.exit(f"ERROR: {error}")
