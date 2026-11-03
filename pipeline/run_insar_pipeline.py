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


def saved_sweets_directory(dataset: str | None, data_root: Path) -> Path | None:
    """Return the SWEETS source saved by a previous GeoTIFF import, if any."""
    if not dataset:
        return None
    metadata = data_root / "ROI" / dataset / "crop_metadata.json"
    if not metadata.is_file():
        return None
    try:
        value = json.loads(metadata.read_text(encoding="utf-8")).get("sweets_directory")
    except (OSError, json.JSONDecodeError):
        return None
    return Path(value) if value else None


def project_label(dataset: str) -> str:
    """Return the project-folder name from either NAME or Dataset_NAME."""
    return dataset.removeprefix("Dataset_")


def default_landfill_project_root(dataset: str) -> Path:
    """Infer the shared CA_Landfill project directory beside the ZoomInSAR clone."""
    return ROOT.parent.parent / "CA_Landfill" / project_label(dataset)


def discover_sweets_slc_directory(project_root: Path) -> Path:
    """Find the one usable SWEETS complex-SLC stack for automatic POT."""
    scratch = project_root / "work" / "gslcs" / "scratch"
    if not scratch.is_dir():
        raise FileNotFoundError(f"SWEETS SLC scratch directory not found: {scratch}")
    candidates = [
        directory for directory in sorted(path for path in scratch.iterdir() if path.is_dir())
        if sum(1 for _ in directory.rglob("*.slc.tif")) >= 2
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"No stack with at least two SWEETS *.slc.tif files found below {scratch}")
    found = ", ".join(str(path) for path in candidates)
    raise ValueError(
        "More than one SWEETS SLC stack was found. Select one with --gslc-dir. "
        f"Candidates: {found}"
    )


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
        "--output-root",
        help=(
            "Directory for imported data and all generated outputs. Defaults to the "
            "pipeline directory; code and virtual environments remain in place."
        ),
    )
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
    parser.add_argument("--pixel-offset", action="store_true",
                        help="Run optional pixel-offset tracking after the phase-InSAR workflow.")
    parser.add_argument("--pixel-offset-only", action="store_true",
                        help="Run only pixel-offset tracking; do not run phase-InSAR stages.")
    parser.add_argument("--gslc-dir",
                        help="Folder containing NISAR GSLC .h5 or SWEETS Sentinel complex *.slc.tif files")
    parser.add_argument(
        "--project-root",
        help=(
            "Optional project folder containing work/. With --pixel-offset-only, "
            "omit this and --gslc-dir to infer CA_Landfill/<dataset> automatically."
        ),
    )
    parser.add_argument("--choose-gslc-dir", action="store_true",
                        help="Prompt for the NISAR GSLC or SWEETS Sentinel SLC directory in the terminal.")
    parser.add_argument("--pot-frequency", choices=("A", "B"), default="A")
    parser.add_argument("--pot-polarization", default="HH")
    parser.add_argument("--pot-window-pixels", type=int, default=64)
    parser.add_argument("--pot-search-pixels", type=int, default=12)
    parser.add_argument("--pot-stride-pixels", type=int, default=32)
    parser.add_argument("--pot-min-correlation", type=float, default=0.25)
    parser.add_argument("--pot-min-baseline-days", type=int, default=0)
    parser.add_argument("--pot-max-baseline-days", type=int, default=24)
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
        "--timeseries-wavelength",
        type=float,
        help=(
            "Radar wavelength in metres for phase-to-LOS conversion. Defaults to "
            "0.238 m for --source sweets (the configured NISAR L-band workflow) "
            "and 0.056 m for GAMMA/Sentinel-1 input."
        ),
    )
    parser.add_argument(
        "--timeseries-phase-sign",
        type=float,
        choices=(-1.0, 1.0),
        help=(
            "LOS sign applied when converting unwrapped phase to displacement. Defaults to "
            "-1 for the configured SWEETS/NISAR workflow and +1 for GAMMA/Sentinel-1 input."
        ),
    )
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

    if args.pixel_offset_only:
        args.pixel_offset = True
        args.skip_cropper = args.skip_filtering = args.skip_unwrapping = True
        args.skip_detrending = args.skip_timeseries = True

    # SWEETS is used here for NISAR GSLC-derived interferograms. Keep legacy
    # GAMMA/Sentinel-1 values unchanged, while allowing explicit overrides for
    # a different SWEETS sensor or independently validated phase convention.
    if args.timeseries_wavelength is None:
        args.timeseries_wavelength = 0.238 if args.source == "sweets" else 0.056
    if args.timeseries_phase_sign is None:
        args.timeseries_phase_sign = -1.0 if args.source == "sweets" else 1.0
    if args.timeseries_wavelength <= 0:
        parser.error("--timeseries-wavelength must be positive")

    crop_dataset = dataset_name(args.dataset)
    automatic_project_root = (
        Path(args.project_root).expanduser().resolve()
        if args.project_root else default_landfill_project_root(crop_dataset)
        if args.pixel_offset_only and crop_dataset else None
    )
    if args.output_root:
        output_root = Path(args.output_root).expanduser().resolve()
    elif automatic_project_root is not None and automatic_project_root.is_dir():
        output_root = automatic_project_root / "ZoomInSAR_Results"
        print(f"Using project output directory: {output_root}", flush=True)
    else:
        output_root = ROOT
    data_root = output_root / "Data"
    saved_dir = saved_sweets_directory(crop_dataset, data_root)
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
    if args.pixel_offset and not compute_dataset:
        parser.error("--pixel-offset requires --dataset or --compute-dataset")
    if args.pixel_offset and args.choose_gslc_dir:
        args.gslc_dir = str(choose_path(
            "SLC/GSLC directory", Path(args.gslc_dir).expanduser() if args.gslc_dir else None,
            directory=True,
        ))
    if args.pixel_offset and not args.gslc_dir and automatic_project_root is not None:
        args.gslc_dir = str(discover_sweets_slc_directory(automatic_project_root))
        print(f"Using discovered SLC directory: {args.gslc_dir}", flush=True)
    if args.pixel_offset and not args.gslc_dir:
        parser.error("--gslc-dir is required with --pixel-offset (or use --choose-gslc-dir)")
    dataset_args = ["--dataset", compute_dataset] if compute_dataset else []
    pair_args = ["--pair", args.pair] if args.pair else []
    output_dem = data_root / "DEM" / "rasters_USGS10m" / "output_USGS10m.tif"
    bundled_dem = ROOT / "Data" / "DEM" / "rasters_USGS10m" / "output_USGS10m.tif"
    default_dem = output_dem if output_dem.is_file() else bundled_dem
    sweets_path = Path(args.sweets_dir).expanduser() if args.sweets_dir else saved_dir
    inferred_sweets_dem = sweets_dem_from_directory(sweets_path) if sweets_path else None
    inferred_sweets_dem = inferred_sweets_dem or default_dem
    dem = Path(args.sweets_dem).expanduser() if args.sweets_dem else inferred_sweets_dem

    unwrap_dir = ROOT / "InSAR_Unwrapping"
    cropper_dir = ROOT / "InSAR_KML_Cropper"
    filtering_dir = ROOT / "InSAR_Filtering"
    detrend_dir = ROOT / "InSAR_Detrending"
    timeseries_dir = ROOT / "InSAR_Timeseries"
    pixel_offset_dir = ROOT / "InSAR_PixelOffset"
    filtering_output = output_root / "InSAR_Filtering" / "OUTPUT"
    unwrapping_output = output_root / "InSAR_Unwrapping" / "OUTPUT"
    detrending_output = output_root / "InSAR_Detrending" / "OUTPUT"
    timeseries_output = output_root / "InSAR_Timeseries" / "OUTPUT"
    pixel_offset_output = output_root / "InSAR_PixelOffset" / "OUTPUT"

    if args.source == "sweets":
        run(
            "1/5  SWEETS GeoTIFF Import" + (" + ROI crop" if args.crop_kml else " (no ROI crop)"),
            venv_python(cropper_dir),
            "insar_kml_cropper", cropper_dir,
            [
                "--data-dir", str(data_root),
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
            "--data-dir", str(data_root),
            "--margin", str(args.margin),
            "--overwrite",
            *pair_args,
        ]
        if crop_dataset:
            roi_name = crop_dataset.removeprefix("Dataset_")
            crop_arguments.extend(["--kml", str(data_root / "ROI" / f"{roi_name}.kml")])
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
                "--roi-dir", str(data_root / "ROI"),
                "--output-dir", str(filtering_output),
                *dataset_args,
                *pair_args,
            ],
        )

    if not args.skip_unwrapping:
        unwrap_input = (
            filtering_output
            if args.input_source == "filtered"
            else data_root / "ROI"
        )
        input_option = "--filtered-root" if args.input_source == "filtered" else "--roi-dir"
        run(
            "3/5  InSAR Unwrapping",
            venv_python(unwrap_dir),
            "insar_unwrapping", unwrap_dir,
            [
                "--input-source", args.input_source,
                input_option, str(unwrap_input),
                "--output-dir", str(unwrapping_output),
                "--dem", str(dem),
                *dataset_args,
                *pair_args,
            ],
        )

    if not args.skip_detrending:
        exclude_masks: dict[str, Path] = {}
        if args.auto_exclude_from_velocity:
            preliminary_output = output_root / "InSAR_Timeseries" / "OUTPUT_preliminary_unwrapped"
            auto_mask_output = output_root / "InSAR_Detrending" / "OUTPUT_auto_masks"
            run(
                "3.5/5  Preliminary Time Series for Auto Exclude Mask",
                venv_python(timeseries_dir),
                "insar_timeseries", timeseries_dir,
                [
                    "--input-source", "unwrapping",
                    "--unwrapping-root", str(unwrapping_output),
                    "--output-dir", str(preliminary_output),
                    "--max-baseline-days", str(args.max_baseline_days),
                    "--wavelength", str(args.timeseries_wavelength),
                    "--phase-sign", str(args.timeseries_phase_sign),
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
                        "--filtering-root", str(filtering_output),
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
                "--unwrapping-root", str(unwrapping_output),
                "--geo-root", str(filtering_output),
                "--dem", str(dem),
                "--output-dir", str(detrending_output),
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
                "--detrending-root", str(detrending_output),
                "--output-dir", str(timeseries_output),
                "--max-baseline-days", str(args.max_baseline_days),
                "--wavelength", str(args.timeseries_wavelength),
                "--phase-sign", str(args.timeseries_phase_sign),
                *(["--allow-low-quality-bridges"] if args.allow_low_quality_bridges else []),
                *dataset_args,
            ],
        )

    if args.pixel_offset:
        run(
            "6/6  SLC/GSLC Pixel-Offset Tracking",
            venv_python(pixel_offset_dir),
            "insar_pixel_offset", pixel_offset_dir,
            [
                "--gslc-dir", str(Path(args.gslc_dir).expanduser().resolve()),
                "--output-dir", str(pixel_offset_output),
                "--dataset", compute_dataset,
                "--frequency", args.pot_frequency,
                "--polarization", args.pot_polarization,
                "--window-pixels", str(args.pot_window_pixels),
                "--search-pixels", str(args.pot_search_pixels),
                "--stride-pixels", str(args.pot_stride_pixels),
                "--min-correlation", str(args.pot_min_correlation),
                "--min-baseline-days", str(args.pot_min_baseline_days),
                "--max-baseline-days", str(args.pot_max_baseline_days),
            ],
        )

    print("\nComplete InSAR workflow finished successfully.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(f"\nPipeline stopped because a stage failed (exit code {error.returncode}).", file=sys.stderr)
        raise SystemExit(error.returncode)
