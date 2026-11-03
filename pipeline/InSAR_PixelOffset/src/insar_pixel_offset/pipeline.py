"""Run pairwise geocoded SLC/GSLC pixel-offset tracking."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import h5py
import rasterio
from tqdm.auto import tqdm

from .gslc import GslcScene, coordinate_steps, discover_scenes, locate_complex_dataset
from .plots import plot_pair, plot_summary
from .tracking import grid_centers, normalized_offset


@dataclass
class PixelOffsetConfig:
    frequency: str = "A"
    polarization: str = "HH"
    window_pixels: int = 64
    search_pixels: int = 12
    stride_pixels: int = 32
    min_correlation: float = 0.25
    min_baseline_days: int = 0
    max_baseline_days: int = 24


def temporal_pairs(scenes: list[GslcScene], min_baseline_days: int,
                   max_baseline_days: int) -> list[tuple[GslcScene, GslcScene]]:
    """Select every scene pair whose temporal baseline is in the requested range."""
    if min_baseline_days < 0 or max_baseline_days < min_baseline_days:
        raise ValueError("Require 0 <= min_baseline_days <= max_baseline_days")
    pairs = [
        (first, second)
        for first_index, first in enumerate(scenes[:-1])
        for second in scenes[first_index + 1:]
        if min_baseline_days <= (second.date - first.date).days <= max_baseline_days
    ]
    if not pairs:
        raise ValueError(
            f"No GSLC pairs with temporal baseline {min_baseline_days}–{max_baseline_days} days"
        )
    return pairs


def _load_pair(first: GslcScene, second: GslcScene, config: PixelOffsetConfig):
    if first.format != second.format:
        raise ValueError("Pixel-offset scene formats do not match")
    if first.format == "geotiff":
        with rasterio.open(first.path) as first_file, rasterio.open(second.path) as second_file:
            reference = _read_complex_geotiff(first_file, first.path)
            secondary = _read_complex_geotiff(second_file, second.path)
            if first_file.crs != second_file.crs or first_file.transform != second_file.transform:
                raise ValueError(f"SWEETS SLC grids differ for {first.path.name} and {second.path.name}")
            x_step, y_step = float(first_file.transform.a), float(first_file.transform.e)
        if reference.shape != secondary.shape:
            raise ValueError(f"SWEETS SLC grid shapes differ for {first.path.name} and {second.path.name}")
        return reference, secondary, x_step, y_step, "GeoTIFF band 1 (complex SLC)"
    with h5py.File(first.path, "r") as first_file, h5py.File(second.path, "r") as second_file:
        first_path = locate_complex_dataset(first_file, config.frequency, config.polarization)
        second_path = locate_complex_dataset(second_file, config.frequency, config.polarization)
        reference = first_file[first_path][...]
        secondary = second_file[second_path][...]
        x_step, y_step = coordinate_steps(first_file, first_path)
    if reference.shape != secondary.shape:
        raise ValueError(f"GSLC grids differ for {first.path.name} and {second.path.name}")
    return reference, secondary, x_step, y_step, first_path


def _read_complex_geotiff(dataset: rasterio.io.DatasetReader, path: Path) -> np.ndarray:
    """Read a complex SLC from either one complex band or real/imaginary bands."""
    first_band = dataset.read(1)
    if np.iscomplexobj(first_band):
        return first_band
    if dataset.count >= 2:
        return first_band.astype(np.float32) + 1j * dataset.read(2).astype(np.float32)
    raise ValueError(
        f"SWEETS SLC GeoTIFF must have one complex band or real/imaginary bands: {path}"
    )


def run_pixel_offset(
    gslc_dir: str | Path,
    output_dir: str | Path,
    dataset: str,
    config: PixelOffsetConfig | None = None,
) -> dict:
    config = config or PixelOffsetConfig()
    scenes = discover_scenes(gslc_dir)
    pairs = temporal_pairs(scenes, config.min_baseline_days, config.max_baseline_days)
    output = Path(output_dir) / (dataset if dataset.startswith("Dataset_") else f"Dataset_{dataset}")
    output.mkdir(parents=True, exist_ok=True)
    figures = output / "figures"

    all_x, all_y, all_corr, pair_records = [], [], [], []
    rows = cols = None
    x_step = y_step = None
    dataset_path = None
    for first, second in tqdm(pairs, desc="Pixel-offset tracking", unit="pair", dynamic_ncols=True):
        reference, secondary, pair_x_step, pair_y_step, path = _load_pair(first, second, config)
        if rows is None:
            rows, cols = grid_centers(reference.shape, config.window_pixels, config.search_pixels, config.stride_pixels)
            x_step, y_step, dataset_path = pair_x_step, pair_y_step, path
        elif reference.shape[0] <= int(rows[-1]) or reference.shape[1] <= int(cols[-1]):
            raise ValueError(f"GSLC grid shape changes at {second.path.name}")
        y_offset = np.full((len(rows), len(cols)), np.nan, dtype=np.float32)
        x_offset = np.full_like(y_offset, np.nan)
        correlation = np.full_like(y_offset, np.nan)
        for iy, row in enumerate(rows):
            for ix, col in enumerate(cols):
                dy, dx, peak = normalized_offset(
                    reference, secondary, int(row), int(col),
                    config.window_pixels, config.search_pixels,
                )
                if np.isfinite(peak) and peak >= config.min_correlation:
                    y_offset[iy, ix], x_offset[iy, ix], correlation[iy, ix] = dy, dx, peak
        all_x.append(x_offset); all_y.append(y_offset); all_corr.append(correlation)
        record = {
            "reference_date": first.date.strftime("%Y-%m-%d"),
            "secondary_date": second.date.strftime("%Y-%m-%d"),
            "baseline_days": (second.date - first.date).days,
            "reference_file": str(first.path), "secondary_file": str(second.path),
            "valid_offsets": int(np.count_nonzero(np.isfinite(correlation))),
        }
        pair_records.append(record)
        x_plot = x_offset * pair_x_step if pair_x_step is not None else x_offset
        y_plot = y_offset * pair_y_step if pair_y_step is not None else y_offset
        plot_pair(
            figures / f"{first.date:%Y%m%d}_{second.date:%Y%m%d}_pixel_offset.png",
            f"{first.date:%Y-%m-%d} to {second.date:%Y-%m-%d}",
            x_plot, y_plot, correlation,
            "m" if pair_x_step is not None and pair_y_step is not None else "pixels",
        )

    x_stack, y_stack, correlation_stack = np.stack(all_x), np.stack(all_y), np.stack(all_corr)
    np.save(output / "x_offset_pixels.npy", x_stack)
    np.save(output / "y_offset_pixels.npy", y_stack)
    np.save(output / "correlation.npy", correlation_stack)
    np.save(output / "grid_row_indices.npy", rows)
    np.save(output / "grid_col_indices.npy", cols)
    if x_step is not None and y_step is not None:
        x_stack, y_stack = x_stack * x_step, y_stack * y_step
        np.save(output / "x_offset_m.npy", x_stack)
        np.save(output / "y_offset_m.npy", y_stack)
    plot_summary(figures / "pixel_offset_summary.png", pair_records, x_stack, y_stack,
                 correlation_stack, "m" if x_step is not None and y_step is not None else "pixels")
    (output / "pair_dates.json").write_text(json.dumps(pair_records, indent=2), encoding="utf-8")
    summary = {
        "scene_directory": str(Path(gslc_dir).expanduser().resolve()),
        "input_format": scenes[0].format,
        "output_directory": str(output.resolve()),
        "config": asdict(config), "complex_dataset_path": dataset_path,
        "scene_count": len(scenes), "pair_count": len(pair_records),
        "input_shape": list(reference.shape), "offset_grid_shape": [len(rows), len(cols)],
        "x_coordinate_step_m": x_step, "y_coordinate_step_m": y_step,
        "units": {
            "x_offset_pixels": "geocoded map-grid pixels, positive increasing x",
            "y_offset_pixels": "geocoded map-grid pixels, positive increasing row",
            "correlation": "peak normalized cross-correlation",
        },
    }
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
