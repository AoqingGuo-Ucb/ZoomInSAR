"""Complete batch unwrapping and correction workflow."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from tqdm.auto import tqdm

from .closure import repair_phase_closure
from .corrections import (
    repair_edge_connected_cycles,
    repair_integer_cycle_regions,
    repair_region_graph_cycles,
)
from .io import discover_shape, pair_name, read_coherence, read_wrapped, write_unwrapped
from .mcf import unwrap_mcf
from .plots import plot_dataset_summary, plot_interferogram_diagnostics, plot_water_mask
from .quality_control import select_network_safe_interferograms
from .water import water_mask_from_dem


@dataclass
class PipelineConfig:
    coherence_threshold: float = 0.05
    max_branch_cut_length: float = 15.0
    spatial_repair: bool = True
    spatial_sigma: float = 12.0
    spatial_min_pixels: int = 40
    region_graph_repair: bool = True
    region_edge_sigma: float = 1.5
    region_edge_threshold_cycles: float = 0.25
    region_boundary_width: int = 5
    region_max_iterations: int = 3
    region_open_boundary_max_gap: int = 12
    region_open_boundary_snap: int = 10
    edge_strip_repair: bool = True
    edge_strip_min_pixels: int = 40
    edge_strip_boundary_width: int = 5
    closure_repair: bool = True
    closure_iterations: int = 8
    closure_min_pixels: int = 20
    closure_boundary_width: int = 3
    auto_exclude_bad: bool = True
    qc_max_residual_fraction: float = 0.04
    qc_max_phase_range_cycles: float = 2.0
    make_plots: bool = True
    mask_water: bool = True
    dem_path: str = "../Data/DEM/rasters_USGS10m/output_USGS10m.tif"
    water_max_elevation: float = 0.0
    selected_pair: str | None = None


def _match_files(int_dir: Path) -> list[tuple[Path, Path, tuple[str, str]]]:
    coherence_by_pair = {pair_name(path): path for path in int_dir.glob("*.coh")}
    result = []
    for phase_file in sorted(int_dir.glob("*.filt")):
        pair = pair_name(phase_file)
        if pair not in coherence_by_pair:
            raise FileNotFoundError(f"No coherence file matches {phase_file.name}")
        result.append((phase_file, coherence_by_pair[pair], pair))
    if not result:
        raise FileNotFoundError(f"No *.filt files found in {int_dir}")
    return result


def run_pipeline(
    data_dir: str | Path,
    output_dir: str | Path,
    shape: tuple[int, int] | None = None,
    config: PipelineConfig | None = None,
) -> dict:
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    config = config or PipelineConfig()
    shape = discover_shape(data_dir, shape)
    int_dir = data_dir / "INT" if (data_dir / "INT").is_dir() else data_dir
    files = _match_files(int_dir)
    if config.selected_pair:
        files = [item for item in files if f"{item[2][0]}-{item[2][1]}" == config.selected_pair]
        if not files:
            raise FileNotFoundError(f"No interferogram found for pair {config.selected_pair} in {data_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "unwrapped").mkdir(exist_ok=True)
    (output_dir / "diagnostics").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    if config.mask_water:
        water_mask, valid_geo, longitude, latitude, elevation = water_mask_from_dem(
            data_dir, shape, config.dem_path, config.water_max_elevation
        )
    else:
        water_mask = np.zeros(shape, dtype=bool)
        valid_geo = np.ones(shape, dtype=bool)
        longitude = latitude = elevation = np.full(shape, np.nan)
    np.save(output_dir / "diagnostics" / "water_mask.npy", water_mask)
    np.save(output_dir / "diagnostics" / "valid_geolocation_mask.npy", valid_geo)
    np.save(output_dir / "diagnostics" / "dem_elevation.npy", elevation.astype(np.float32))
    if config.make_plots:
        plot_water_mask(
            output_dir / "figures", longitude, latitude, elevation, water_mask,
            valid_geo, config.water_max_elevation,
        )

    stack, wrapped_stack, coherence_stack, cut_stack, pairs, rows = [], [], [], [], [], []
    for phase_file, coherence_file, pair in tqdm(
        files, total=len(files), unit="interferogram",
        desc=f"Unwrapping {data_dir.name}", dynamic_ncols=True,
    ):
        wrapped = read_wrapped(phase_file, shape)
        coherence = read_coherence(coherence_file, shape)
        # Exclude water before every spatial and temporal unwrapping operation.
        wrapped = wrapped.copy(); coherence = coherence.copy()
        wrapped[water_mask] = np.nan
        coherence[water_mask] = np.nan
        unwrapped, offsets, cuts, diagnostics = unwrap_mcf(
            wrapped, coherence, config.coherence_threshold, config.max_branch_cut_length
        )
        spatial_map = np.zeros(shape, dtype=np.int16)
        spatial_stats = None
        region_stats = None
        edge_stats = None
        if config.spatial_repair:
            unwrapped, spatial_map, spatial_stats = repair_integer_cycle_regions(
                unwrapped, sigma=config.spatial_sigma, min_pixels=config.spatial_min_pixels
            )
        if config.region_graph_repair:
            unwrapped, region_map, region_stats = repair_region_graph_cycles(
                unwrapped,
                coherence,
                cuts,
                edge_sigma=config.region_edge_sigma,
                edge_threshold_cycles=config.region_edge_threshold_cycles,
                boundary_width=config.region_boundary_width,
                min_pixels=config.spatial_min_pixels,
                max_iterations=config.region_max_iterations,
                open_boundary_max_gap=config.region_open_boundary_max_gap,
                open_boundary_snap=config.region_open_boundary_snap,
            )
            spatial_map += region_map
        edge_map = np.zeros(shape,dtype=np.int16)
        if config.edge_strip_repair:
            unwrapped,edge_map,edge_stats=repair_edge_connected_cycles(
                unwrapped,coherence,
                boundary_width=config.edge_strip_boundary_width,
                min_pixels=config.edge_strip_min_pixels,
            )
            spatial_map += edge_map
        stack.append(unwrapped); wrapped_stack.append(wrapped); coherence_stack.append(coherence)
        cut_stack.append(cuts); pairs.append(pair)
        np.save(output_dir / "diagnostics" / f"{pair[0]}-{pair[1]}_offsets.npy", offsets)
        np.save(output_dir / "diagnostics" / f"{pair[0]}-{pair[1]}_branch_cuts.npy", cuts)
        np.save(output_dir / "diagnostics" / f"{pair[0]}-{pair[1]}_spatial_cycles.npy", spatial_map)
        if config.region_graph_repair:
            np.save(output_dir / "diagnostics" / f"{pair[0]}-{pair[1]}_region_graph_cycles.npy", region_map)
        if config.edge_strip_repair:
            np.save(output_dir / "diagnostics" / f"{pair[0]}-{pair[1]}_edge_cycles.npy", edge_map)
        row = {"pair": f"{pair[0]}-{pair[1]}", **asdict(diagnostics)}
        if spatial_stats:
            row.update({f"spatial_{key}": value for key, value in asdict(spatial_stats).items()})
        if region_stats:
            row.update({f"region_graph_{key}": value for key, value in asdict(region_stats).items()})
        if edge_stats:
            row.update({f"edge_strip_{key}": value for key,value in asdict(edge_stats).items()})
        rows.append(row)

    unwrapped_stack = np.stack(stack)
    coherence_array = np.stack(coherence_stack)
    np.save(output_dir / "unwrapped_before_closure.npy", unwrapped_stack.astype(np.float32))

    # Repair closure on the complete network before quality-control exclusion.
    # Excluding a bad edge first can destroy the only triangle capable of
    # diagnosing and repairing that edge. QC must assess the fully corrected
    # result, not the preliminary unwrapped stack.
    closure_stats = None
    closure_map = np.zeros(unwrapped_stack.shape, dtype=np.int16)
    if config.closure_repair and len(pairs) >= 3:
        unwrapped_stack, closure_map, closure_stats = repair_phase_closure(
            unwrapped_stack, pairs, coherence_array,
            max_iterations=config.closure_iterations,
            min_component_pixels=config.closure_min_pixels,
            boundary_width=config.closure_boundary_width,
        )
    np.save(output_dir / "closure_cycle_corrections.npy", closure_map)

    if config.auto_exclude_bad:
        kept_indices, qc_rows = select_network_safe_interferograms(
            unwrapped_stack, pairs, sigma=config.spatial_sigma,
            max_residual_fraction=config.qc_max_residual_fraction,
            max_phase_range_cycles=config.qc_max_phase_range_cycles,
        )
    else:
        kept_indices = list(range(len(pairs)))
        qc_rows = [{
            "pair": f"{a}-{b}", "spatial_residual_fraction": None,
            "robust_phase_range_cycles": None, "bad_candidate": False,
            "included": True, "status": "included_qc_disabled",
        } for a, b in pairs]
    kept_set = set(kept_indices)
    for index, row in enumerate(qc_rows):
        rows[index].update({f"qc_{key}": value for key, value in row.items() if key != "pair"})
    with (output_dir / "interferogram_exclusion_report.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(qc_rows[0])); writer.writeheader(); writer.writerows(qc_rows)
    np.save(output_dir / "unwrapped_final.npy", unwrapped_stack.astype(np.float32))
    for index, pair in enumerate(pairs):
        destination = output_dir / "unwrapped" / f"{pair[0]}-{pair[1]}.unw"
        if index in kept_set:
            write_unwrapped(destination, unwrapped_stack[index])
        elif destination.exists():
            destination.unlink()

    if config.make_plots:
        for index, pair in enumerate(pairs):
            pair_text = f"{pair[0]}-{pair[1]}"
            plot_interferogram_diagnostics(
                output_dir / "figures", pair_text, wrapped_stack[index], coherence_array[index],
                cut_stack[index], stack[index], closure_map[index], unwrapped_stack[index],
            )
        plot_dataset_summary(
            output_dir / "figures", [f"{a}-{b}" for a, b in pairs], rows,
            asdict(closure_stats) if closure_stats else None,
        )

    fieldnames = sorted({key for row in rows for key in row})
    with (output_dir / "unwrap_quality.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames); writer.writeheader(); writer.writerows(rows)
    summary = {
        "input_directory": str(data_dir), "output_directory": str(output_dir),
        "shape_lines_width": list(shape), "interferograms": len(files),
        "processing_order": [
            "mcf_unwrap", "spatial_repair", "region_graph_repair", "edge_strip_repair",
            "full_network_closure_repair", "post_repair_quality_control",
            "network_safe_exclusion",
        ],
        "pairs": [f"{pairs[i][0]}-{pairs[i][1]}" for i in kept_indices],
        "excluded_pairs": [f"{a}-{b}" for i, (a, b) in enumerate(pairs) if i not in kept_set],
        "all_input_pairs": [f"{a}-{b}" for a, b in pairs], "config": asdict(config),
        "water_mask": {
            "enabled": bool(config.mask_water),
            "water_pixels": int(np.count_nonzero(water_mask)),
            "valid_geolocation_pixels": int(np.count_nonzero(valid_geo)),
            "method": "dem_elevation_from_lon_lat" if config.mask_water else None,
            "dem_path": config.dem_path if config.mask_water else None,
            "water_max_elevation": float(config.water_max_elevation),
        },
        "closure": asdict(closure_stats) if closure_stats else None,
    }
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def run_all_datasets(
    roi_dir: str | Path,
    output_root: str | Path,
    config: PipelineConfig | None = None,
) -> list[dict]:
    """Process every Dataset_* directory under an ROI directory."""
    roi_dir, output_root = Path(roi_dir), Path(output_root)
    datasets = sorted(path for path in roi_dir.glob("Dataset_*") if path.is_dir())
    if not datasets:
        raise FileNotFoundError(f"No Dataset_* directories found in {roi_dir}")
    output_root.mkdir(parents=True, exist_ok=True)
    summaries = [run_pipeline(dataset, output_root / dataset.name, config=config) for dataset in tqdm(
        datasets, total=len(datasets), unit="dataset", desc="Unwrapping datasets",
        dynamic_ncols=True,
    )]
    (output_root / "batch_summary.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8"
    )
    return summaries
