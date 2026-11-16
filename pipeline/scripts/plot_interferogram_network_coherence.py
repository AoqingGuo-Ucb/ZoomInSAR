#!/usr/bin/env python3
"""Plot an interferogram network with links coloured by mean coherence.

The script reads pair dates from SWEETS/Dolphin interferogram names such as
``YYYYMMDD_YYYYMMDD.int.cor.tif`` and calculates each pair's mean valid
coherence.  A true perpendicular-baseline plot additionally needs one
perpendicular baseline value for every acquisition; pass those values in a CSV
with ``date,perpendicular_baseline_m`` columns.  Without that CSV, the vertical
axis deliberately shows acquisition order rather than a falsely labelled
perpendicular baseline.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import rasterio
import h5py
from rasterio.crs import CRS
from rasterio.warp import transform


PAIR_PATTERN = re.compile(r"(20\d{6})[_-](20\d{6})")


def mean_coherence(path: Path) -> float:
    """Calculate a nodata-aware mean without loading a large raster at once."""
    total, count = 0.0, 0
    with rasterio.open(path) as dataset:
        nodata = dataset.nodata
        for _, window in dataset.block_windows(1):
            values = dataset.read(1, window=window, masked=True).astype(float)
            values = np.asarray(values.filled(np.nan))
            valid = np.isfinite(values) & (values >= 0.0) & (values <= 1.0)
            if nodata is not None:
                valid &= ~np.isclose(values, nodata)
            total += float(np.sum(values[valid]))
            count += int(np.sum(valid))
    return total / count if count else float("nan")


def find_coherence_files(project_root: Path) -> list[Path]:
    """Choose the directory containing the largest coherent interferogram set."""
    grouped: dict[Path, list[Path]] = defaultdict(list)
    for path in project_root.rglob("*.tif"):
        if path.name.lower().endswith(".int.cor.tif") and PAIR_PATTERN.search(path.name):
            grouped[path.parent].append(path)
    if not grouped:
        raise FileNotFoundError(f"No YYYYMMDD_YYYYMMDD.int.cor.tif files found below {project_root}")
    directory, files = max(grouped.items(), key=lambda item: len(item[1]))
    print(f"Using interferogram directory: {directory}")
    return sorted(files)


def read_baselines(path: Path) -> dict[datetime, float]:
    result: dict[datetime, float] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            try:
                date = datetime.fromisoformat(row["date"])
                result[date] = float(row["perpendicular_baseline_m"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Invalid baseline CSV row: {row}") from error
    if not result:
        raise ValueError("Baseline CSV has no usable rows")
    return result


def _dataset_path(handle: h5py.File, ending: str, required_text: str = "") -> str | None:
    """Find one HDF5 dataset by final name and a distinguishing parent path."""
    matches: list[str] = []
    def visit(name: str, item: object) -> None:
        if isinstance(item, h5py.Dataset) and name.lower().endswith(ending.lower()) and required_text.lower() in name.lower():
            matches.append(name)
    handle.visititems(visit)
    return sorted(matches)[0] if matches else None


def _array(handle: h5py.File, ending: str, parent: str = "") -> np.ndarray:
    name = _dataset_path(handle, ending, parent)
    if not name:
        raise ValueError(f"Missing HDF5 dataset ending in {ending} under {parent or '/'}")
    return np.asarray(handle[name][...], dtype=float)


def _ecef_los_from_enu(east: float, north: float, up: float, longitude: float, latitude: float) -> np.ndarray:
    """Rotate a local ENU look vector into Earth-centred Earth-fixed coordinates."""
    lon, lat = np.deg2rad([longitude, latitude])
    east_axis = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north_axis = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    up_axis = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    vector = east * east_axis + north * north_axis + up * up_axis
    return vector / np.linalg.norm(vector)


def _scene_baseline_geometry(path: Path) -> tuple[datetime, np.ndarray, np.ndarray, np.ndarray]:
    """Return date, satellite ECEF position, LOS ECEF, and along-track vector."""
    date = datetime.strptime(PAIR_PATTERN.search(path.name).group(1), "%Y%m%d") if PAIR_PATTERN.search(path.name) else datetime.strptime(re.search(r"(20\d{6})", path.name).group(1), "%Y%m%d")
    with h5py.File(path, "r") as handle:
        orbit_time = _array(handle, "/time", "/metadata/orbit")
        position = _array(handle, "/position", "/metadata/orbit")
        velocity = _array(handle, "/velocity", "/metadata/orbit")
        zero_doppler = _array(handle, "/zeroDopplerAzimuthTime", "/metadata/radarGrid")
        los_east = _array(handle, "/losUnitVectorX", "/metadata/radarGrid")
        los_north = _array(handle, "/losUnitVectorY", "/metadata/radarGrid")
        x_coordinates = _array(handle, "/xCoordinates", "/metadata/radarGrid")
        y_coordinates = _array(handle, "/yCoordinates", "/metadata/radarGrid")
        projection_name = _dataset_path(handle, "/projection", "/metadata/radarGrid")
        if not projection_name:
            raise ValueError("Missing radarGrid projection")
        projection_raw = handle[projection_name][()]
        projection_text = projection_raw.decode() if isinstance(projection_raw, bytes) else str(projection_raw)
    height_index, row, col = (size // 2 for size in zero_doppler.shape)
    time_value = float(zero_doppler[height_index, row, col])
    satellite = np.array([np.interp(time_value, orbit_time, position[:, component]) for component in range(3)])
    along_track = np.array([np.interp(time_value, orbit_time, velocity[:, component]) for component in range(3)])
    # Metadata-cube x/y postings align with its column/row dimensions.
    x_value, y_value = float(x_coordinates[col]), float(y_coordinates[row])
    longitude, latitude = transform(CRS.from_user_input(projection_text), "EPSG:4326", [x_value], [y_value])
    east, north = float(los_east[height_index, row, col]), float(los_north[height_index, row, col])
    up = np.sqrt(max(0.0, 1.0 - east * east - north * north))
    return date, satellite, _ecef_los_from_enu(east, north, up, longitude[0], latitude[0]), along_track / np.linalg.norm(along_track)


def auto_nisar_baselines(project_root: Path, dates: list[datetime]) -> dict[datetime, float]:
    """Compute signed effective perpendicular baselines from NISAR GSLC metadata."""
    h5_files = sorted({*project_root.rglob("*.h5"), *project_root.rglob("*.H5")})
    by_date: dict[datetime, Path] = {}
    for path in h5_files:
        match = re.search(r"NISAR_L2_.*?(20\d{6})T", path.name)
        if match:
            by_date.setdefault(datetime.strptime(match.group(1), "%Y%m%d"), path)
    missing = [date for date in dates if date not in by_date]
    if missing:
        raise ValueError(f"No NISAR GSLC file was found for {len(missing)} interferogram acquisition date(s)")
    geometry = {date: _scene_baseline_geometry(by_date[date]) for date in dates}
    reference_date = dates[0]
    _, reference_position, reference_los, reference_track = geometry[reference_date]
    cross_track = np.cross(reference_los, reference_track)
    cross_track /= np.linalg.norm(cross_track)
    baselines = {reference_date: 0.0}
    for date in dates[1:]:
        _, satellite, _, _ = geometry[date]
        baselines[date] = float(np.dot(satellite - reference_position, cross_track))
    return baselines


def parse_pair(path: Path) -> tuple[datetime, datetime]:
    match = PAIR_PATTERN.search(path.name)
    if not match:
        raise ValueError(f"Could not read two dates from {path.name}")
    return tuple(datetime.strptime(value, "%Y%m%d") for value in match.groups())  # type: ignore[return-value]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--interferogram-dir", type=Path, help="Directory containing *.int.cor.tif files")
    source.add_argument("--project-root", type=Path, help="Search this project recursively for coherence rasters")
    parser.add_argument("--baseline-csv", type=Path, help="Optional CSV with date,perpendicular_baseline_m")
    parser.add_argument("--output-dir", type=Path, help="Default: <project-root>/ZoomInSAR_Results/Network_Quality")
    parser.add_argument("--title", default="Interferogram network coloured by mean coherence")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.interferogram_dir:
        coherence_files = sorted(args.interferogram_dir.expanduser().glob("*.int.cor.tif"))
        root = args.interferogram_dir.expanduser().resolve().parent
    else:
        root = args.project_root.expanduser().resolve()
        coherence_files = find_coherence_files(root)
    if not coherence_files:
        raise FileNotFoundError("No *.int.cor.tif coherence files found")

    records = []
    for path in coherence_files:
        first, second = parse_pair(path)
        coherence = mean_coherence(path)
        if np.isfinite(coherence):
            records.append((first, second, coherence, path))
    if not records:
        raise ValueError("No coherence raster contains valid values in [0, 1]")
    dates = sorted({date for first, second, _, _ in records for date in (first, second)})
    baselines = read_baselines(args.baseline_csv) if args.baseline_csv else None
    baseline_source = str(args.baseline_csv) if args.baseline_csv else None
    if baselines is None and args.project_root:
        try:
            baselines = auto_nisar_baselines(root, dates)
            baseline_source = "calculated from NISAR GSLC orbit and radar-grid LOS metadata"
            print("Calculated effective perpendicular baselines from NISAR GSLC metadata")
        except (OSError, ValueError, KeyError, IndexError) as error:
            print(f"NISAR perpendicular-baseline calculation unavailable: {error}")
    if baselines:
        missing = [date.strftime("%Y-%m-%d") for date in dates if date not in baselines]
        if missing:
            raise ValueError(f"Baseline CSV lacks {len(missing)} acquisition dates, e.g. {', '.join(missing[:5])}")
        y = baselines
        ylabel = "Perpendicular baseline (m)"
    else:
        y = {date: index for index, date in enumerate(dates)}
        ylabel = "Acquisition order (perpendicular baseline not supplied)"

    figure, axis = plt.subplots(figsize=(12, 7), layout="constrained")
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(0.0, 1.0)
    segments = [[(mdates.date2num(first), y[first]), (mdates.date2num(second), y[second])] for first, second, _, _ in records]
    lines = LineCollection(segments, colors=[cmap(norm(coherence)) for _, _, coherence, _ in records], linewidths=1.8, alpha=0.9)
    axis.add_collection(lines)
    axis.scatter(dates, [y[date] for date in dates], s=30, c="black", zorder=3)
    axis.set_xlim(min(dates), max(dates))
    values = np.array([y[date] for date in dates], dtype=float)
    padding = max((values.max() - values.min()) * 0.08, 1.0)
    axis.set_ylim(values.min() - padding, values.max() + padding)
    axis.set_title(args.title)
    axis.set_xlabel("Acquisition date")
    axis.set_ylabel(ylabel)
    axis.grid(alpha=0.25)
    axis.xaxis.set_major_locator(mdates.AutoDateLocator())
    axis.xaxis.set_major_formatter(mdates.ConciseDateFormatter(axis.xaxis.get_major_locator()))
    colorbar = figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axis, pad=0.02)
    colorbar.set_label("Mean interferogram coherence")

    output_dir = args.output_dir or root / "ZoomInSAR_Results" / "Network_Quality"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = output_dir / "interferogram_network_mean_coherence.png"
    figure.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    csv_path = output_dir / "interferogram_mean_coherence.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["reference_date", "secondary_date", "mean_coherence", "coherence_file"])
        writer.writeheader()
        for first, second, coherence, path in records:
            writer.writerow({"reference_date": first.date(), "secondary_date": second.date(), "mean_coherence": f"{coherence:.6f}", "coherence_file": path})
    (output_dir / "network_plot_summary.json").write_text(json.dumps({
        "interferograms": len(records), "acquisitions": len(dates), "vertical_axis": ylabel,
        "coherence_directory": str(coherence_files[0].parent), "baseline_source": baseline_source,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote network figure: {figure_path}")
    print(f"Wrote coherence table: {csv_path}")


if __name__ == "__main__":
    main()
