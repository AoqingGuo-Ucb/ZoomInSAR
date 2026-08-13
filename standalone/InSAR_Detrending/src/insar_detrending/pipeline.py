"""Batch detrending with InSAR_Unwrapping-compatible output."""

from __future__ import annotations

import csv
import json
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
import xml.etree.ElementTree as ET

from matplotlib.path import Path as MplPath
import numpy as np
from tqdm.auto import tqdm

from .detrend import build_consensus_stable_mask, robust_polynomial_detrend
from .io import discover_unwrapped, read_shape, read_unwrapped, write_unwrapped
from .plots import plot_case, plot_summary
from .terrain import lowpass_elevation, sample_dem


@dataclass
class DetrendConfig:
    degree: int = 1
    orbit_model: str = "plane"
    mad_scale: float = 3.5
    max_iterations: int = 5
    stable_valid_fraction: float = 0.8
    stable_dispersion_percentile: float = 40.0
    fit_edge_guard_pixels: int = 0
    output_edge_guard_pixels: int = 0
    reference_border_fraction: float = 0.0
    auto_reference_point: bool = True
    reference_window_pixels: int = 9
    reference_pixel_count: int = 1000
    reference_kml: str | None = None
    preserve_median: bool = False
    exclude_mask: str | None = None
    make_plots: bool = True
    terrain_correction: bool = True
    terrain_degree: int = 2
    terrain_fit_method: str = "spatial"
    terrain_bins: int = 30
    terrain_min_bin_pixels: int = 100
    terrain_local_radius_pixels: int = 80
    terrain_strength: float = 0.3
    terrain_max_range_fraction: float = 0.6
    terrain_local_guard_pixels: int = 96
    terrain_local_guard_tolerance: float = 0.15
    terrain_local_guard_min_pixels: int = 200
    terrain_smoothing_pixels: int = 30
    terrain_smoothing_passes: int = 1
    turbulent_qa: bool = True
    turbulent_smoothing_pixels: int = 60
    geo_root: str = "../InSAR_Filtering/OUTPUT"
    dem_path: str = "../Data/DEM/rasters_USGS10m/output_USGS10m.tif"
    selected_pair: str | None = None


def _load_lon_lat(filtering_dataset: Path, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    geo = filtering_dataset / "GEO"
    lon_files, lat_files = sorted(geo.glob("*.lon")), sorted(geo.glob("*.lat"))
    if len(lon_files) != 1 or len(lat_files) != 1:
        raise ValueError(f"{geo} must contain exactly one *.lon and one *.lat")
    dtype = np.dtype("f4").newbyteorder(">")
    longitude = np.asarray(np.memmap(lon_files[0], dtype=dtype, mode="r", shape=shape), float)
    latitude = np.asarray(np.memmap(lat_files[0], dtype=dtype, mode="r", shape=shape), float)
    return longitude, latitude


def _kml_coordinate_rings(path: Path) -> list[np.ndarray]:
    root = ET.parse(path).getroot()
    rings = []
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] != "coordinates" or not element.text:
            continue
        coordinates = []
        for token in element.text.split():
            parts = token.split(",")
            if len(parts) < 2:
                continue
            try:
                coordinates.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
        if len(coordinates) >= 3:
            rings.append(np.asarray(coordinates, dtype=float))
    if not rings:
        raise ValueError(f"No polygon coordinates found in {path}")
    return rings


def _reference_region_mask(reference_kml: str | None, filtering_dataset: Path,
                           shape: tuple[int, int]) -> np.ndarray | None:
    if not reference_kml:
        return None
    longitude, latitude = _load_lon_lat(filtering_dataset, shape)
    valid = (
        np.isfinite(longitude) & np.isfinite(latitude)
        & (longitude >= -180) & (longitude <= 180)
        & (latitude >= -90) & (latitude <= 90)
        & ~((longitude == 0) & (latitude == 0))
    )
    points = np.column_stack([longitude[valid], latitude[valid]])
    inside_valid = np.zeros(points.shape[0], dtype=bool)
    for ring in _kml_coordinate_rings(Path(reference_kml)):
        inside_valid |= MplPath(ring[:, :2]).contains_points(points)
    mask = np.zeros(shape, dtype=bool)
    mask[valid] = inside_valid
    if not np.any(mask):
        raise ValueError(f"Reference KML does not overlap valid pixels: {reference_kml}")
    return mask


def _load_mean_coherence(files: list[Path], filtering_dataset: Path,
                         shape: tuple[int, int]) -> np.ndarray | None:
    coherence_values = []
    diagnostics = filtering_dataset / "diagnostics"
    for path in files:
        coherence_path = diagnostics / f"{path.stem}_input_coherence.npy"
        if not coherence_path.exists():
            continue
        coherence = np.load(coherence_path).astype(float)
        if coherence.shape != shape:
            continue
        coherence_values.append(coherence)
    if not coherence_values:
        return None
    with np.errstate(invalid="ignore"):
        return np.nanmean(np.stack(coherence_values), axis=0)


def _select_reference_pixels(stack: np.ndarray, stable_mask: np.ndarray,
                             mean_coherence: np.ndarray | None,
                             window_pixels: int,
                             pixel_count: int) -> tuple[dict | None, np.ndarray | None]:
    radius = max(1, int(window_pixels) // 2)
    valid_fraction = np.mean(np.isfinite(stack), axis=0)
    centered = stack - np.nanmedian(stack, axis=(1, 2), keepdims=True)
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        temporal_center = np.nanmedian(centered, axis=0)
        temporal_mad = np.nanmedian(np.abs(centered - temporal_center), axis=0)
    candidate = stable_mask & np.isfinite(temporal_mad) & (valid_fraction >= 0.8)
    candidate[:radius, :] = False
    candidate[-radius:, :] = False
    candidate[:, :radius] = False
    candidate[:, -radius:] = False
    if mean_coherence is not None:
        candidate &= np.isfinite(mean_coherence) & (mean_coherence > 0)
    if not np.any(candidate):
        return None, None
    dispersion = temporal_mad.copy()
    dispersion[~np.isfinite(dispersion)] = np.nan
    disp_scale = np.nanpercentile(dispersion[candidate], 95)
    if not np.isfinite(disp_scale) or disp_scale <= 0:
        disp_scale = 1.0
    if mean_coherence is None:
        score = -dispersion / disp_scale
    else:
        score = mean_coherence - 0.25 * dispersion / disp_scale
    score[~candidate] = -np.inf
    flat_order = np.argsort(score.ravel())[::-1]
    finite_order = flat_order[np.isfinite(score.ravel()[flat_order])]
    if finite_order.size < 5:
        return None, None
    selected_count = min(max(5, int(pixel_count)), int(finite_order.size))
    selected = finite_order[:selected_count]
    reference_mask = np.zeros(score.shape, dtype=bool)
    reference_mask.ravel()[selected] = True
    for flat_index in flat_order[:2000]:
        if not np.isfinite(score.ravel()[flat_index]):
            break
        row, col = np.unravel_index(flat_index, score.shape)
        patch = stable_mask[row - radius:row + radius + 1, col - radius:col + radius + 1]
        if np.count_nonzero(patch) >= max(5, patch.size // 2):
            return {
                "row": int(row),
                "col": int(col),
                "window_pixels": int(radius * 2 + 1),
                "reference_pixels": int(np.count_nonzero(reference_mask)),
                "mean_coherence": (
                    None if mean_coherence is None
                    else float(mean_coherence[row, col])
                ),
                "median_reference_coherence": (
                    None if mean_coherence is None
                    else float(np.nanmedian(mean_coherence[reference_mask]))
                ),
                "temporal_mad_rad": float(temporal_mad[row, col]),
                "median_reference_temporal_mad_rad": float(np.nanmedian(temporal_mad[reference_mask])),
                "valid_fraction": float(valid_fraction[row, col]),
            }, reference_mask
    return None, None


def _reference_value(image: np.ndarray, reference: dict | None,
                     reference_mask: np.ndarray | None) -> float:
    if not reference:
        return 0.0
    if reference_mask is not None and np.any(reference_mask):
        value = float(np.nanmedian(image[reference_mask]))
        return value if np.isfinite(value) else 0.0
    radius = int(reference["window_pixels"]) // 2
    row, col = int(reference["row"]), int(reference["col"])
    patch = image[row - radius:row + radius + 1, col - radius:col + radius + 1]
    value = float(np.nanmedian(patch))
    return value if np.isfinite(value) else 0.0


def _remove_stale_case_outputs(output: Path, current_stems: set[str]) -> list[str]:
    """Remove generated per-interferogram products absent from the current input."""
    removed = []
    patterns = (
        (output / "unwrapped", "*.unw", ""),
        (output / "diagnostics", "*_removed_trend.npy", "_removed_trend"),
        (output / "diagnostics", "*_orbit_ramp.npy", "_orbit_ramp"),
        (output / "diagnostics", "*_after_orbit.npy", "_after_orbit"),
        (output / "diagnostics", "*_dem_correlated.npy", "_dem_correlated"),
        (output / "diagnostics", "*_after_dem.npy", "_after_dem"),
        (output / "diagnostics", "*_turbulent_qa.npy", "_turbulent_qa"),
        (output / "diagnostics", "*_after_turbulent_qa.npy", "_after_turbulent_qa"),
        (output / "diagnostics", "*_fit_mask.npy", "_fit_mask"),
        (output / "figures", "*_detrending.png", "_detrending"),
    )
    for directory, pattern, suffix in patterns:
        if not directory.exists():
            continue
        for path in directory.glob(pattern):
            stem = path.stem
            if suffix and stem.endswith(suffix):
                stem = stem[:-len(suffix)]
            if stem not in current_stems:
                path.unlink()
                removed.append(str(path))
    return removed


def _edge_guard_mask(shape: tuple[int, int], pixels: int) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    if pixels <= 0:
        return mask
    guard = int(pixels)
    if guard * 2 >= min(shape):
        raise ValueError("edge guard is too large for the image shape")
    mask[:guard, :] = True
    mask[-guard:, :] = True
    mask[:, :guard] = True
    mask[:, -guard:] = True
    return mask


def run_dataset(unwrapping_dataset: str | Path, output: str | Path,
                config: DetrendConfig | None = None) -> dict:
    unwrapping_dataset, output = Path(unwrapping_dataset), Path(output)
    config = config or DetrendConfig()
    shape = read_shape(unwrapping_dataset)
    files = discover_unwrapped(unwrapping_dataset)
    stack = np.stack([read_unwrapped(path, shape) for path in files])
    processing_indices = [
        index for index, path in enumerate(files)
        if config.selected_pair is None or config.selected_pair in path.name
    ]
    if not processing_indices:
        raise FileNotFoundError(f"No interferogram found for pair {config.selected_pair}")
    elevation = None
    if config.terrain_correction:
        cached_elevation = unwrapping_dataset / "diagnostics" / "dem_elevation.npy"
        if cached_elevation.exists():
            elevation = np.load(cached_elevation)
        else:
            elevation = sample_dem(
                Path(config.geo_root) / unwrapping_dataset.name,
                shape, Path(config.dem_path),
            )
            np.save(cached_elevation, elevation.astype(np.float32))
        if config.terrain_smoothing_pixels > 0:
            elevation = lowpass_elevation(
                elevation, config.terrain_smoothing_pixels,
                passes=config.terrain_smoothing_passes,
            )
            smooth_name = (
                f"dem_elevation_lowpass_r{config.terrain_smoothing_pixels}"
                f"_p{config.terrain_smoothing_passes}.npy"
            )
            np.save(
                unwrapping_dataset / "diagnostics" / smooth_name,
                elevation.astype(np.float32),
            )
    stable_mask = build_consensus_stable_mask(
        stack, config.stable_valid_fraction, config.stable_dispersion_percentile
    )
    if config.fit_edge_guard_pixels > 0:
        stable_mask &= ~_edge_guard_mask(shape, config.fit_edge_guard_pixels)
    output_edge_mask = _edge_guard_mask(shape, config.output_edge_guard_pixels)
    if np.any(output_edge_mask):
        stack = stack.copy()
        stack[:, output_edge_mask] = np.nan
    if config.reference_border_fraction > 0:
        if not 0 < config.reference_border_fraction < 0.5:
            raise ValueError("reference_border_fraction must be between 0 and 0.5")
        border_rows = max(1, int(round(shape[0] * config.reference_border_fraction)))
        border_cols = max(1, int(round(shape[1] * config.reference_border_fraction)))
        border_mask = np.zeros(shape, dtype=bool)
        border_mask[:border_rows, :] = True
        border_mask[-border_rows:, :] = True
        border_mask[:, :border_cols] = True
        border_mask[:, -border_cols:] = True
        stable_mask &= border_mask
    exclude = None
    if config.exclude_mask:
        exclude = np.load(config.exclude_mask).astype(bool)
        if exclude.shape != shape:
            raise ValueError(f"Exclude mask shape {exclude.shape} does not match {shape}")
    reference_region = _reference_region_mask(
        config.reference_kml, Path(config.geo_root) / unwrapping_dataset.name, shape
    )
    output.mkdir(parents=True, exist_ok=True)
    (output / "unwrapped").mkdir(exist_ok=True)
    (output / "diagnostics").mkdir(exist_ok=True)
    (output / "figures").mkdir(exist_ok=True)
    removed_stale_outputs = _remove_stale_case_outputs(
        output, {path.stem for path in files}
    )
    np.save(output / "diagnostics" / "consensus_stable_mask.npy", stable_mask)
    if reference_region is not None:
        np.save(output / "diagnostics" / "reference_region_mask.npy", reference_region)
    reference = None
    reference_mask = None
    if config.auto_reference_point:
        reference_candidate_mask = stable_mask.copy()
        if exclude is not None:
            reference_candidate_mask &= ~exclude
        if reference_region is not None:
            reference_candidate_mask &= reference_region
        mean_coherence = _load_mean_coherence(
            files, Path(config.geo_root) / unwrapping_dataset.name, shape
        )
        reference, reference_mask = _select_reference_pixels(
            stack, reference_candidate_mask,
            mean_coherence, config.reference_window_pixels,
            config.reference_pixel_count,
        )
    if reference_mask is not None:
        np.save(output / "diagnostics" / "reference_mask.npy", reference_mask)
    (output / "diagnostics" / "reference_point.json").write_text(
        json.dumps(reference, indent=2), encoding="utf-8"
    )
    corrected_stack, rows = [], []
    cases = [(files[index], stack[index]) for index in processing_indices]
    for index, (path, image) in enumerate(tqdm(
        cases, total=len(cases), unit="interferogram",
        desc=f"Detrending {unwrapping_dataset.name}", dynamic_ncols=True,
    )):
        corrected, surface, final_mask, stats, components = robust_polynomial_detrend(
            image, config.degree, stable_mask, exclude, config.mad_scale,
            config.max_iterations, config.preserve_median,
            elevation=elevation, terrain_degree=config.terrain_degree,
            terrain_fit_method=config.terrain_fit_method,
            terrain_bins=config.terrain_bins,
            terrain_min_bin_pixels=config.terrain_min_bin_pixels,
            terrain_local_radius_pixels=config.terrain_local_radius_pixels,
            terrain_strength=config.terrain_strength,
            terrain_max_range_fraction=config.terrain_max_range_fraction,
            terrain_local_guard_pixels=config.terrain_local_guard_pixels,
            terrain_local_guard_tolerance=config.terrain_local_guard_tolerance,
            terrain_local_guard_min_pixels=config.terrain_local_guard_min_pixels,
            orbit_model=config.orbit_model,
            turbulent_qa=config.turbulent_qa,
            turbulent_smoothing_pixels=config.turbulent_smoothing_pixels,
        )
        reference_offset = _reference_value(corrected, reference, reference_mask)
        corrected = corrected - reference_offset
        surface = surface + reference_offset
        components["after_orbit"] = components["after_orbit"] - reference_offset
        components["after_dem"] = components["after_dem"] - reference_offset
        components["after_turbulent_qa"] = components["after_turbulent_qa"] - reference_offset
        if np.any(output_edge_mask):
            corrected[output_edge_mask] = np.nan
            surface[output_edge_mask] = np.nan
            for component in components.values():
                component[output_edge_mask] = np.nan
        corrected_stack.append(corrected)
        write_unwrapped(output / "unwrapped" / path.name, corrected)
        np.save(output / "diagnostics" / f"{path.stem}_removed_trend.npy", surface.astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_orbit_ramp.npy", components["orbit"].astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_after_orbit.npy", components["after_orbit"].astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_dem_correlated.npy", components["terrain"].astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_after_dem.npy", components["after_dem"].astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_turbulent_qa.npy", components["turbulent_qa"].astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_after_turbulent_qa.npy", components["after_turbulent_qa"].astype(np.float32))
        np.save(output / "diagnostics" / f"{path.stem}_fit_mask.npy", final_mask)
        row = {
            "interferogram": path.stem,
            "reference_offset_rad": reference_offset,
            **asdict(stats),
        }; rows.append(row)
        if config.make_plots:
            plot_case(output / "figures" / f"{path.stem}_detrending.png", path.stem,
                      image, surface, corrected, final_mask, reference, components)
    corrected_array = np.stack(corrected_stack).astype(np.float32)
    np.save(output / "detrended_unwrapped_stack.npy", corrected_array)
    with (output / "detrending_quality.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    if config.make_plots:
        plot_summary(
            output / "figures" / "detrending_summary.png",
            [row["interferogram"] for row in rows],
            [row["before_robust_std_rad"] for row in rows],
            [row["after_robust_std_rad"] for row in rows],
            [row["correction_robust_range_rad"] for row in rows],
        )
    source_summary = json.loads((unwrapping_dataset / "run_summary.json").read_text(encoding="utf-8"))
    summary = {
        "input_directory": str(unwrapping_dataset),
        "input_contract": "InSAR_Unwrapping/OUTPUT/Dataset_*/unwrapped/*.unw",
        "output_directory": str(output), "shape_lines_width": list(shape),
        "interferograms": len(files), "config": asdict(config),
        "stable_pixels": int(np.count_nonzero(stable_mask)),
        "reference_point": reference,
        "removed_stale_outputs": removed_stale_outputs,
        "source_unwrapping_summary": source_summary,
    }
    # Keep this conventional filename so InSAR_Timeseries can consume either source.
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def run_all_datasets(unwrapping_root: str | Path, output_root: str | Path,
                     config: DetrendConfig | None = None) -> list[dict]:
    unwrapping_root, output_root = Path(unwrapping_root), Path(output_root)
    datasets = sorted(path for path in unwrapping_root.glob("Dataset_*") if path.is_dir())
    if not datasets:
        raise FileNotFoundError(f"No Dataset_* directories in {unwrapping_root}")
    summaries = [run_dataset(path, output_root / path.name, config) for path in tqdm(
        datasets, total=len(datasets), unit="dataset", desc="Detrending datasets",
        dynamic_ncols=True,
    )]
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "batch_summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    return summaries
