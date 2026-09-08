"""Run Cropper -> Filtering -> Unwrapping -> Detrending -> Timeseries."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def dataset_name(value: str | None) -> str | None:
    if not value:
        return None
    return value if value.startswith("Dataset_") else f"Dataset_{value}"


def run(label: str, python: Path, module: str, cwd: Path, arguments: list[str]) -> None:
    if not python.is_file():
        raise FileNotFoundError(f"Python environment not found: {python}")
    command = [str(python), "-m", module, *arguments]
    print(f"\n{'=' * 72}\n{label}\n{'=' * 72}", flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def venv_python(project: Path) -> Path:
    """Return the virtual-environment Python on Windows, macOS, or Linux."""
    candidates = (project / ".venv" / "Scripts" / "python.exe", project / ".venv" / "bin" / "python")
    return next((path for path in candidates if path.is_file()), candidates[0])


def sweets_dem_from_directory(sweets_dir: Path) -> Path | None:
    """Find SWEETS' work-level DEM from an output directory or a descendant."""
    return next(
        (path / "dem.tif" for path in (sweets_dir, *sweets_dir.parents) if (path / "dem.tif").is_file()),
        None,
    )


def saved_sweets_directory(dataset: str | None) -> Path | None:
    """Return the SWEETS source saved by a previous GeoTIFF import, if any."""
    if not dataset:
        return None
    metadata = ROOT / "Data" / "ROI" / dataset / "crop_metadata.json"
    if not metadata.is_file():
        return None
    try:
        value = json.loads(metadata.read_text(encoding="utf-8")).get("sweets_directory")
    except (OSError, json.JSONDecodeError):
        return None
    return Path(value) if value else None


def choose_path(label: str, default: Path | None, *, directory: bool) -> Path:
    """Interactively select an existing directory or file in a terminal."""
    expectation = "directory" if directory else "file"
    while True:
        prompt = f"{label} [{default}]: " if default else f"{label}: "
        value = input(prompt).strip()
        candidate = Path(value).expanduser() if value else default
        if candidate is None:
            print(f"Please enter a {expectation} path.")
            continue
        candidate = candidate.resolve()
        if (candidate.is_dir() if directory else candidate.is_file()):
            return candidate
        print(f"Not an existing {expectation}: {candidate}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the complete InSAR crop, filter, unwrap, detrend, and time-series workflow"
    )
    parser.add_argument("--dataset", help="Optional single Dataset_* directory name for the full workflow")
    parser.add_argument(
        "--source", choices=("gamma", "sweets"), default="gamma",
        help="Input format: existing GAMMA binaries (default) or SWEETS GeoTIFF products.",
    )
    parser.add_argument("--sweets-dir", help="SWEETS work/output directory containing wrapped-phase and coherence GeoTIFFs")
    parser.add_argument(
        "--choose-sweets-paths", action="store_true",
        help="Interactively choose the SWEETS interferogram directory and DEM in the terminal.",
    )
    parser.add_argument(
        "--sweets-phase-pattern", default="*.int.tif",
        help="Recursive glob for SWEETS wrapped-phase GeoTIFFs (default: *.int.tif).",
    )
    parser.add_argument(
        "--sweets-coherence-pattern", default="*.int.cor.tif",
        help="Recursive glob for SWEETS pairwise coherence GeoTIFFs (default: *.int.cor.tif).",
    )
    parser.add_argument(
        "--sweets-skip-unpaired", action="store_true",
        help="Exclude SWEETS interferograms lacking a matching coherence GeoTIFF.",
    )
    crop_choice = parser.add_mutually_exclusive_group()
    crop_choice.add_argument(
        "--crop-kml",
        help="KML used to crop SWEETS input.",
    )
    crop_choice.add_argument(
        "--no-crop", action="store_true",
        help="Process the full SWEETS GeoTIFF extent (the default for --source sweets).",
    )
    parser.add_argument(
        "--sweets-dem",
        help="Optional DEM for SWEETS processing. If omitted, <sweets-dir>/dem.tif is used when present.",
    )
    parser.add_argument(
        "--compute-dataset",
        help=(
            "Optional single Dataset_* directory name for stages after KML cropping "
            "(filtering, unwrapping, detrending, and time-series)."
        ),
    )
    parser.add_argument("--pair", help="Reprocess one YYYYMMDD-YYYYMMDD pair through detrending")
    parser.add_argument("--max-baseline-days", type=int, default=12)
    parser.add_argument(
        "--allow-low-quality-bridges",
        action="store_true",
        dest="allow_low_quality_bridges",
        help=(
            "Allow low-quality interferograms that are structurally required to keep "
            "the time-series network connected. They are included with low weights (default)."
        ),
    )
    parser.add_argument(
        "--no-allow-low-quality-bridges",
        action="store_false",
        dest="allow_low_quality_bridges",
        help="Stop instead of using a low-quality interferogram required for connectivity.",
    )
    parser.set_defaults(allow_low_quality_bridges=True)
    parser.add_argument("--input-source", choices=("filtered", "raw"), default="filtered")
    parser.add_argument(
        "--detrend-degree",
        type=int,
        choices=(0, 1, 2),
        default=1,
        help="Polynomial degree for detrending; default 1 removes a smooth planar ramp.",
    )
    parser.add_argument(
        "--orbit-model",
        choices=(
            "constant", "plane", "range_quadratic", "azimuth_quadratic",
            "range_cubic", "right_edge_quadratic",
            "additive_quadratic", "full_quadratic", "degree",
        ),
        default="plane",
        help=(
            "DEM-independent orbit-ramp model. right_edge_quadratic targets "
            "right-edge residuals; degree preserves the old --detrend-degree behavior."
        ),
    )
    parser.add_argument(
        "--stable-dispersion-percentile",
        type=float,
        default=40.0,
        help="Percentile of temporally quiet pixels used for detrending fit; lower values are more conservative.",
    )
    parser.add_argument(
        "--fit-edge-guard-pixels",
        type=int,
        default=0,
        help="Exclude this many pixels along each image edge from detrending fits while preserving full output coverage.",
    )
    parser.add_argument(
        "--output-edge-guard-pixels",
        type=int,
        default=0,
        help="Mask this many pixels along each image edge in detrended outputs and diagnostic plots.",
    )
    parser.add_argument(
        "--reference-border-fraction",
        type=float,
        default=0.0,
        help=(
            "Optional fractional border width used as the detrending reference; "
            "for example 0.2 fits only the outer 20%% border."
        ),
    )
    parser.add_argument(
        "--no-auto-reference-point",
        action="store_false",
        dest="auto_reference_point",
        help="Disable automatic high-coherence reference-point zeroing.",
    )
    parser.add_argument(
        "--reference-window-pixels",
        type=int,
        default=9,
        help="Odd window size used to validate the representative automatic reference point.",
    )
    parser.add_argument(
        "--reference-pixel-count",
        type=int,
        default=1000,
        help="Maximum number of high-coherence stable pixels used for robust reference zeroing.",
    )
    parser.add_argument(
        "--reference-kml",
        help="Optional KML polygon where zero-deformation reference pixels are allowed to be selected.",
    )
    parser.add_argument(
        "--auto-exclude-from-velocity",
        action="store_true",
        dest="auto_exclude_from_velocity",
        help=(
            "Before detrending, run a preliminary time series directly from unwrapped interferograms "
            "and build an automatic deformation exclude mask from the preliminary velocity map (default)."
        ),
    )
    parser.add_argument(
        "--no-auto-exclude-from-velocity",
        action="store_false",
        dest="auto_exclude_from_velocity",
        help="Disable the default preliminary-velocity deformation exclude mask.",
    )
    parser.set_defaults(auto_exclude_from_velocity=True)
    parser.add_argument("--auto-exclude-velocity-percentile", type=float, default=85.0)
    parser.add_argument("--auto-exclude-min-velocity", type=float, default=0.02)
    parser.add_argument("--auto-exclude-dilation-pixels", type=int, default=6)
    parser.add_argument("--auto-exclude-min-component-pixels", type=int, default=80)
    parser.add_argument("--auto-exclude-max-kml-polygons", type=int, default=20)
    parser.add_argument("--auto-exclude-min-kml-polygon-area-pixels", type=float, default=100.0)
    parser.add_argument(
        "--terrain-correction",
        action="store_true",
        help="Include a low-pass DEM term in detrending; enabled by default.",
    )
    parser.add_argument(
        "--no-terrain-correction",
        action="store_false",
        dest="terrain_correction",
        help="Disable the DEM-related stratified-atmosphere term.",
    )
    parser.set_defaults(terrain_correction=True)
    parser.add_argument("--terrain-degree", type=int, choices=(1, 2), default=2)
    parser.add_argument(
        "--terrain-fit-method",
        choices=("spatial", "hybrid", "local", "binned", "pixel"),
        default="spatial",
        help="Fit DEM-correlated phase from after-orbit residuals using spatially varying DEM terms, a hybrid DEM-bin/local model, local coefficient field, elevation-bin medians, or robust pixels.",
    )
    parser.add_argument("--terrain-bins", type=int, default=30)
    parser.add_argument("--terrain-min-bin-pixels", type=int, default=100)
    parser.add_argument("--terrain-strength", type=float, default=0.3)
    parser.add_argument("--terrain-max-range-fraction", type=float, default=0.6)
    parser.add_argument("--terrain-local-radius-pixels", type=int, default=80)
    parser.add_argument(
        "--terrain-smoothing-pixels",
        type=int,
        default=30,
        help="Low-pass smoothing radius applied to DEM before detrending terrain fitting; use 0 for raw DEM.",
    )
    parser.add_argument(
        "--terrain-smoothing-passes",
        type=int,
        default=1,
        help="Number of repeated DEM smoothing passes.",
    )
    parser.add_argument(
        "--no-turbulent-qa",
        action="store_false",
        dest="turbulent_qa",
        help="Do not write the optional non-DEM low-frequency residual QA layer.",
    )
    parser.set_defaults(turbulent_qa=True)
    parser.add_argument("--turbulent-smoothing-pixels", type=int, default=60)
    parser.add_argument("--margin", type=float, default=0.20)
    parser.add_argument("--skip-cropper", action="store_true")
    parser.add_argument("--skip-filtering", action="store_true")
    parser.add_argument("--skip-unwrapping", action="store_true")
    parser.add_argument("--skip-detrending", action="store_true")
    parser.add_argument("--skip-timeseries", action="store_true")
    args = parser.parse_args()

    crop_dataset = dataset_name(args.dataset)
    saved_dir = saved_sweets_directory(crop_dataset)
    selected_sweets_dir = Path(args.sweets_dir).expanduser() if args.sweets_dir else saved_dir
    if selected_sweets_dir:
        args.sweets_dir = str(selected_sweets_dir)
    if args.choose_sweets_paths:
        selected_sweets_dir = choose_path(
            "SWEETS interferogram directory", selected_sweets_dir, directory=True,
        )
        args.sweets_dir = str(selected_sweets_dir)
        detected_dem = sweets_dem_from_directory(selected_sweets_dir)
        explicit_dem = Path(args.sweets_dem).expanduser() if args.sweets_dem else detected_dem
        args.sweets_dem = str(choose_path("DEM GeoTIFF", explicit_dem, directory=False))
    if args.source == "sweets" and not args.sweets_dir:
        parser.error("--sweets-dir is required when --source sweets (or use --choose-sweets-paths)")
    if args.source == "sweets" and not crop_dataset:
        parser.error("--dataset is required when --source sweets")
    compute_dataset = dataset_name(args.compute_dataset) or crop_dataset
    dataset_args = ["--dataset", compute_dataset] if compute_dataset else []
    pair_args = ["--pair", args.pair] if args.pair else []
    default_dem = ROOT / "Data" / "DEM" / "rasters_USGS10m" / "output_USGS10m.tif"
    sweets_path = Path(args.sweets_dir).expanduser() if args.sweets_dir else saved_dir
    inferred_sweets_dem = sweets_dem_from_directory(sweets_path) if sweets_path else None
    inferred_sweets_dem = inferred_sweets_dem or default_dem
    dem = Path(args.sweets_dem).expanduser() if args.sweets_dem else inferred_sweets_dem

    unwrap_dir = ROOT / "InSAR_Unwrapping"
    cropper_dir = ROOT / "InSAR_KML_Cropper"
    filtering_dir = ROOT / "InSAR_Filtering"
    detrend_dir = ROOT / "InSAR_Detrending"
    timeseries_dir = ROOT / "InSAR_Timeseries"

    if args.source == "sweets":
        run(
            "1/5  SWEETS GeoTIFF Import" + (" + ROI crop" if args.crop_kml else " (no ROI crop)"),
            venv_python(cropper_dir),
            "insar_kml_cropper", cropper_dir,
            [
                "--data-dir", str(ROOT / "Data"),
                "--sweets-dir", str(args.sweets_dir),
                "--dataset", crop_dataset,
                "--phase-pattern", args.sweets_phase_pattern,
                "--coherence-pattern", args.sweets_coherence_pattern,
                "--margin", str(args.margin),
                "--overwrite",
                *(["--skip-unpaired"] if args.sweets_skip_unpaired else []),
                *( ["--no-crop"] if args.no_crop else [] ),
                *(["--crop-kml", str(args.crop_kml)] if args.crop_kml else []),
            ],
        )
    elif not args.skip_cropper:
        crop_arguments = [
            "--data-dir", str(ROOT / "Data"),
            "--margin", str(args.margin),
            "--overwrite",
            *pair_args,
        ]
        if crop_dataset:
            roi_name = crop_dataset.removeprefix("Dataset_")
            crop_arguments.extend(["--kml", str(ROOT / "Data" / "ROI" / f"{roi_name}.kml")])
        run(
            "1/5  InSAR KML Cropper",
            venv_python(cropper_dir),
            "insar_kml_cropper", cropper_dir, crop_arguments,
        )

    if not args.skip_filtering:
        run(
            "2/5  InSAR Filtering",
            venv_python(filtering_dir),
            "insar_filtering", filtering_dir,
            [
                "--roi-dir", str(ROOT / "Data" / "ROI"),
                "--output-dir", str(filtering_dir / "OUTPUT"),
                *dataset_args,
                *pair_args,
            ],
        )

    if not args.skip_unwrapping:
        unwrap_input = (
            ROOT / "InSAR_Filtering" / "OUTPUT"
            if args.input_source == "filtered"
            else ROOT / "Data" / "ROI"
        )
        input_option = "--filtered-root" if args.input_source == "filtered" else "--roi-dir"
        run(
            "3/5  InSAR Unwrapping",
            venv_python(unwrap_dir),
            "insar_unwrapping", unwrap_dir,
            [
                "--input-source", args.input_source,
                input_option, str(unwrap_input),
                "--output-dir", str(unwrap_dir / "OUTPUT"),
                "--dem", str(dem),
                *dataset_args,
                *pair_args,
            ],
        )

    if not args.skip_detrending:
        exclude_masks: dict[str, Path] = {}
        if args.auto_exclude_from_velocity:
            preliminary_output = timeseries_dir / "OUTPUT_preliminary_unwrapped"
            auto_mask_output = detrend_dir / "OUTPUT_auto_masks"
            run(
                "3.5/5  Preliminary Time Series for Auto Exclude Mask",
                venv_python(timeseries_dir),
                "insar_timeseries", timeseries_dir,
                [
                    "--input-source", "unwrapping",
                    "--unwrapping-root", str(unwrap_dir / "OUTPUT"),
                    "--output-dir", str(preliminary_output),
                    "--max-baseline-days", str(args.max_baseline_days),
                    *(["--allow-low-quality-bridges"] if args.allow_low_quality_bridges else []),
                    *dataset_args,
                ],
            )
            mask_datasets = ([compute_dataset] if compute_dataset else sorted(
                path.name for path in preliminary_output.glob("Dataset_*") if path.is_dir()
            ))
            if not mask_datasets:
                raise FileNotFoundError(f"No preliminary Dataset_* outputs in {preliminary_output}")
            for mask_dataset in mask_datasets:
                run(
                    f"3.6/5  Auto Deformation Exclude Mask ({mask_dataset})",
                    venv_python(detrend_dir),
                    "insar_detrending.auto_mask_cli", detrend_dir,
                    [
                        "--timeseries-root", str(preliminary_output),
                        "--filtering-root", str(filtering_dir / "OUTPUT"),
                        "--output-dir", str(auto_mask_output),
                        "--dataset", mask_dataset,
                        "--velocity-percentile", str(args.auto_exclude_velocity_percentile),
                        "--min-velocity", str(args.auto_exclude_min_velocity),
                        "--dilation-pixels", str(args.auto_exclude_dilation_pixels),
                        "--min-component-pixels", str(args.auto_exclude_min_component_pixels),
                        "--max-kml-polygons", str(args.auto_exclude_max_kml_polygons),
                        "--min-kml-polygon-area-pixels", str(args.auto_exclude_min_kml_polygon_area_pixels),
                    ],
                )
                exclude_masks[mask_dataset] = auto_mask_output / mask_dataset / "auto_deformation_exclude_mask.npy"

        detrend_datasets = [compute_dataset] if compute_dataset else (list(exclude_masks) or [None])
        for detrend_dataset in detrend_datasets:
            detrend_dataset_args = ["--dataset", detrend_dataset] if detrend_dataset else []
            exclude_mask_args = (
                ["--exclude-mask", str(exclude_masks[detrend_dataset])]
                if detrend_dataset in exclude_masks else []
            )
            run(
                f"4/5  InSAR Detrending{f' ({detrend_dataset})' if detrend_dataset else ''}",
                venv_python(detrend_dir),
                "insar_detrending", detrend_dir,
                [
                "--unwrapping-root", str(unwrap_dir / "OUTPUT"),
                "--geo-root", str(ROOT / "InSAR_Filtering" / "OUTPUT"),
                "--dem", str(dem),
                "--output-dir", str(detrend_dir / "OUTPUT"),
                "--degree", str(args.detrend_degree),
                "--orbit-model", args.orbit_model,
                "--stable-dispersion-percentile", str(args.stable_dispersion_percentile),
                "--fit-edge-guard-pixels", str(args.fit_edge_guard_pixels),
                "--output-edge-guard-pixels", str(args.output_edge_guard_pixels),
                "--reference-border-fraction", str(args.reference_border_fraction),
                "--reference-window-pixels", str(args.reference_window_pixels),
                "--reference-pixel-count", str(args.reference_pixel_count),
                *(["--reference-kml", str(args.reference_kml)] if args.reference_kml else []),
                *exclude_mask_args,
                *(["--no-auto-reference-point"] if not args.auto_reference_point else []),
                "--terrain-degree", str(args.terrain_degree),
                "--terrain-fit-method", args.terrain_fit_method,
                "--terrain-bins", str(args.terrain_bins),
                "--terrain-min-bin-pixels", str(args.terrain_min_bin_pixels),
                "--terrain-local-radius-pixels", str(args.terrain_local_radius_pixels),
                "--terrain-strength", str(args.terrain_strength),
                "--terrain-max-range-fraction", str(args.terrain_max_range_fraction),
                "--terrain-smoothing-pixels", str(args.terrain_smoothing_pixels),
                "--terrain-smoothing-passes", str(args.terrain_smoothing_passes),
                "--turbulent-smoothing-pixels", str(args.turbulent_smoothing_pixels),
                *(["--no-turbulent-qa"] if not args.turbulent_qa else []),
                *(["--terrain-correction"] if args.terrain_correction else []),
                *detrend_dataset_args,
                *pair_args,
                ],
            )

    if not args.skip_timeseries:
        run(
            "5/5  InSAR Timeseries",
            venv_python(timeseries_dir),
            "insar_timeseries", timeseries_dir,
            [
                "--input-source", "detrending",
                "--detrending-root", str(detrend_dir / "OUTPUT"),
                "--output-dir", str(timeseries_dir / "OUTPUT"),
                "--max-baseline-days", str(args.max_baseline_days),
                *(["--allow-low-quality-bridges"] if args.allow_low_quality_bridges else []),
                *dataset_args,
            ],
        )

    print("\nComplete InSAR workflow finished successfully.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(f"\nPipeline stopped because a stage failed (exit code {error.returncode}).", file=sys.stderr)
        raise SystemExit(error.returncode)
