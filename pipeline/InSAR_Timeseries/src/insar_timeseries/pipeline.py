"""End-to-end SBAS time-series workflow using InSAR_Unwrapping OUTPUT only."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from tqdm.auto import tqdm

from .io import discover_unwrapped, read_shape, read_unwrapped
from .network import select_connected_network
from .geotiff import export_velocity_geotiff, load_geolocation
from .plots import (plot_displacement_maps, plot_network, plot_representative_timeseries,
                    plot_velocity_and_residual, select_representative_points)
from .quality import EdgeQuality, load_edge_quality
from .solver import invert_timeseries


@dataclass
class TimeseriesConfig:
    max_baseline_days: int = 12
    wavelength_m: float = 0.056
    phase_sign: float = 1.0
    make_plots: bool = True
    fail_on_low_quality_bridge: bool = True


def run_dataset(unwrapping_dataset: str | Path, output: str | Path,
                config: TimeseriesConfig | None = None) -> dict:
    unwrapping_dataset, output = Path(unwrapping_dataset), Path(output)
    config = config or TimeseriesConfig()
    shape = read_shape(unwrapping_dataset)
    available = discover_unwrapped(unwrapping_dataset)
    all_pairs = [(first, second) for _, first, second in available]
    quality_by_pair = load_edge_quality(unwrapping_dataset)
    qualities = [quality_by_pair.get(f"{a:%Y%m%d}-{b:%Y%m%d}", EdgeQuality()) for a, b in all_pairs]
    selected, edge_details, dates = select_connected_network(
        all_pairs, config.max_baseline_days, qualities
    )
    bridges = [edge for edge in edge_details if edge.selected and "connectivity_bridge" in edge.reason]
    low_quality = [edge for edge in edge_details if edge.selected and edge.quality_tier >= 2]
    if config.fail_on_low_quality_bridge and low_quality:
        pairs_text = ", ".join(f"{e.first:%Y%m%d}-{e.second:%Y%m%d}" for e in low_quality)
        raise ValueError(
            f"{unwrapping_dataset.name}: connected inversion requires low-quality interferogram bridge(s): "
            f"{pairs_text}. Use --allow-low-quality-bridges only if you accept low-weight but still structurally necessary edges."
        )
    selected_files = [available[index][0] for index in selected]
    selected_pairs = [all_pairs[index] for index in selected]
    selected_weights = np.array([edge_details[index].inversion_weight for index in selected])
    phase_stack = np.stack([read_unwrapped(path, shape) for path in tqdm(
        selected_files, total=len(selected_files), unit="interferogram",
        desc=f"Loading time series {unwrapping_dataset.name}", dynamic_ncols=True,
    )])
    displacement, velocity, residual, stats = invert_timeseries(
        phase_stack, selected_pairs, dates, config.wavelength_m, config.phase_sign,
        edge_weights=selected_weights,
    )

    output.mkdir(parents=True, exist_ok=True)
    (output / "diagnostics").mkdir(exist_ok=True)
    (output / "figures").mkdir(exist_ok=True)
    np.save(output / "cumulative_los_displacement_m.npy", displacement.astype(np.float32))
    np.save(output / "mean_los_velocity_m_per_year.npy", velocity.astype(np.float32))
    np.save(output / "diagnostics" / "phase_residual_rms_rad.npy", residual.astype(np.float32))
    longitude, latitude = load_geolocation(unwrapping_dataset, shape)
    geotiff = export_velocity_geotiff(
        output / "mean_los_velocity_mm_per_year.tif", velocity, longitude, latitude
    )
    points = select_representative_points(velocity, longitude, latitude)
    with (output / "representative_points.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["id", "row", "col", "longitude", "latitude", "velocity_m_per_year",
                  "velocity_mm_per_year", "selection"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(points)
    point_geojson = {
        "type": "FeatureCollection",
        "features": [{"type": "Feature", "geometry": {"type": "Point",
            "coordinates": [p["longitude"], p["latitude"]]},
            "properties": {k: v for k, v in p.items() if k not in ("longitude", "latitude")}}
            for p in points],
    }
    (output / "representative_points.geojson").write_text(
        json.dumps(point_geojson, indent=2), encoding="utf-8"
    )
    (output / "dates.json").write_text(
        json.dumps([date.strftime("%Y-%m-%d") for date in dates], indent=2), encoding="utf-8"
    )
    with (output / "interferogram_network.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "index", "first", "second", "baseline_days", "selected", "reason",
            "quality_tier", "inversion_weight", "qc_status", "detrending_accepted", "filename",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for edge in edge_details:
            row = asdict(edge); row["first"] = edge.first.strftime("%Y-%m-%d")
            row["second"] = edge.second.strftime("%Y-%m-%d")
            row["filename"] = available[edge.index][0].name
            writer.writerow(row)
    if config.make_plots:
        plot_network(output / "figures" / "selected_network.png", dates, edge_details,
                     config.max_baseline_days)
        plot_velocity_and_residual(output / "figures" / "velocity_and_residual.png", velocity, residual, points)
        plot_displacement_maps(output / "figures" / "cumulative_displacement_maps.png", displacement, dates, points)
        plot_representative_timeseries(output / "figures" / "representative_timeseries.png",
                                       displacement, dates, points,
                                       [(edge.first, edge.second) for edge in low_quality])
    summary = {
        "input_directory": str(unwrapping_dataset),
        "input_contract": "InSAR_Unwrapping or InSAR_Detrending Dataset_*/unwrapped/*.unw",
        "output_directory": str(output), "shape_lines_width": list(shape),
        "config": asdict(config), "available_interferograms": len(available),
        "selected_interferograms": len(selected), "acquisition_dates": len(dates),
        "network_rank": len(dates) - 1,
        "independent_cycles": len(selected) - len(dates) + 1,
        "connectivity_bridges": len(bridges),
        "bridge_pairs": [f"{e.first:%Y-%m-%d}_{e.second:%Y-%m-%d}" for e in bridges],
        "low_quality_selected": len(low_quality),
        "low_quality_pairs": [f"{e.first:%Y-%m-%d}_{e.second:%Y-%m-%d}" for e in low_quality],
        "velocity_geotiff": geotiff,
        "representative_points": points,
        "inversion": asdict(stats),
    }
    (output / "timeseries_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def run_all_datasets(unwrapping_root: str | Path, output_root: str | Path,
                     config: TimeseriesConfig | None = None) -> list[dict]:
    unwrapping_root, output_root = Path(unwrapping_root), Path(output_root)
    datasets = sorted(path for path in unwrapping_root.glob("Dataset_*") if path.is_dir())
    if not datasets:
        raise FileNotFoundError(f"No Dataset_* directories in {unwrapping_root}")
    summaries = [run_dataset(path, output_root / path.name, config) for path in tqdm(
        datasets, total=len(datasets), unit="dataset", desc="Time-series datasets",
        dynamic_ncols=True,
    )]
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "batch_summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    return summaries
