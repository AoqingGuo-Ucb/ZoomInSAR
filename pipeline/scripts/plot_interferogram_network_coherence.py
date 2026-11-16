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
        "coherence_directory": str(coherence_files[0].parent), "baseline_csv": str(args.baseline_csv) if args.baseline_csv else None,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote network figure: {figure_path}")
    print(f"Wrote coherence table: {csv_path}")


if __name__ == "__main__":
    main()
