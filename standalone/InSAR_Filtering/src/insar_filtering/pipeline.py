"""Batch filtering workflow with unwrapping-compatible output layout."""

from __future__ import annotations

import csv
import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter
from tqdm.auto import tqdm

from .filters import hybrid_filter
from .io import discover_shape, pair_name, read_gamma, write_gamma
from .plots import plot_filter_result, plot_summary


@dataclass
class FilterConfig:
    search_radius: int = 7
    patch_radius: int = 3
    h: float = 0.70
    goldstein_strength: float = 0.40
    make_plots: bool = True


def circular_variance(values: np.ndarray, valid: np.ndarray) -> float:
    unit = np.divide(values, np.abs(values), out=np.zeros_like(values), where=valid)
    local = uniform_filter(unit.real, 3) + 1j * uniform_filter(unit.imag, 3)
    return float(np.nanmedian((1.0 - np.abs(local))[valid])) if np.any(valid) else float("nan")


def _match_files(int_dir: Path) -> list[tuple[Path, Path, tuple[str, str]]]:
    coherence = {pair_name(path): path for path in int_dir.glob("*.coh")}
    result = []
    for phase in sorted(int_dir.glob("*.filt")):
        pair = pair_name(phase)
        if pair not in coherence:
            raise FileNotFoundError(f"No coherence file matches {phase.name}")
        result.append((phase, coherence[pair], pair))
    if not result:
        raise FileNotFoundError(f"No *.filt files in {int_dir}")
    return result


def run_dataset(dataset: str | Path, output: str | Path, config: FilterConfig | None = None,
                selected_pair: str | None = None) -> dict:
    dataset, output = Path(dataset), Path(output)
    config = config or FilterConfig()
    shape = discover_shape(dataset)
    (output / "INT").mkdir(parents=True, exist_ok=True)
    (output / "GEO").mkdir(parents=True, exist_ok=True)
    (output / "diagnostics").mkdir(parents=True, exist_ok=True)
    (output / "figures").mkdir(parents=True, exist_ok=True)
    rows = []
    matched_files = _match_files(dataset / "INT")
    if selected_pair:
        matched_files = [item for item in matched_files if f"{item[2][0]}-{item[2][1]}" == selected_pair]
        if not matched_files:
            raise FileNotFoundError(f"No interferogram found for pair {selected_pair} in {dataset}")
    for phase_file, coherence_file, pair in tqdm(
        matched_files, total=len(matched_files), unit="interferogram",
        desc=f"Filtering {dataset.name}", dynamic_ncols=True,
    ):
        interferogram = read_gamma(phase_file, shape, "c8").astype(np.complex128)
        coherence = read_gamma(coherence_file, shape, "f4").astype(float)
        filtered, concentration, support = hybrid_filter(
            interferogram, coherence, config.search_radius, config.patch_radius,
            config.h, config.goldstein_strength,
        )
        # Preserve filenames so the unwrapping package can consume this directory directly.
        write_gamma(output / "INT" / phase_file.name, filtered, "c8")
        write_gamma(output / "INT" / coherence_file.name, concentration, "f4")
        pair_text = f"{pair[0]}-{pair[1]}"
        np.save(output / "diagnostics" / f"{pair_text}_input_coherence.npy", coherence.astype(np.float32))
        np.save(output / "diagnostics" / f"{pair_text}_phase_concentration.npy", concentration)
        np.save(output / "diagnostics" / f"{pair_text}_nonlocal_support.npy", support)
        valid = np.isfinite(interferogram) & (np.abs(interferogram) > 0)
        row = {
            "pair": pair_text,
            "input_circular_variance": circular_variance(interferogram, valid),
            "filtered_circular_variance": circular_variance(filtered, valid),
            "median_input_coherence": float(np.nanmedian(coherence[valid])),
            "median_phase_concentration": float(np.nanmedian(concentration[valid])),
            "median_absolute_phase_change_rad": float(np.nanmedian(np.abs(np.angle(filtered * np.conj(interferogram)))[valid])),
        }
        rows.append(row)
        if config.make_plots:
            plot_filter_result(output / "figures" / f"{pair_text}_filtering.png", pair_text,
                               interferogram, filtered, coherence, concentration)
    # Copy geometry and crop metadata to keep Dataset_* self-contained and compatible.
    for source in (dataset / "GEO").glob("*"):
        if source.is_file():
            shutil.copy2(source, output / "GEO" / source.name)
    for source in dataset.glob("*.kml"):
        shutil.copy2(source, output / source.name)
    if (dataset / "crop_metadata.json").exists():
        shutil.copy2(dataset / "crop_metadata.json", output / "crop_metadata.json")
    with (output / "filter_quality.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    if config.make_plots:
        plot_summary(output / "figures" / "dataset_filtering_summary.png", rows)
    summary = {"input_directory": str(dataset), "output_directory": str(output),
               "shape_lines_width": list(shape), "interferograms": len(rows),
               "method": "complex_nonlocal_plus_coherence_adaptive_goldstein",
               "config": asdict(config)}
    (output / "filter_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def run_all_datasets(roi_dir: str | Path, output_root: str | Path,
                     config: FilterConfig | None = None,
                     selected_pair: str | None = None) -> list[dict]:
    roi_dir, output_root = Path(roi_dir), Path(output_root)
    datasets = sorted(path for path in roi_dir.glob("Dataset_*") if path.is_dir())
    if not datasets:
        raise FileNotFoundError(f"No Dataset_* directories in {roi_dir}")
    summaries = [run_dataset(path, output_root / path.name, config, selected_pair) for path in tqdm(
        datasets, total=len(datasets), unit="dataset", desc="Filtering datasets",
        dynamic_ncols=True,
    )]
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "batch_summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    return summaries
