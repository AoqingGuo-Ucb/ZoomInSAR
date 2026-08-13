"""Command-line interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .pipeline import DetrendConfig, run_all_datasets, run_dataset


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Robust polynomial detrending of unwrapped InSAR")
    command.add_argument("--unwrapping-root", default="../InSAR_Unwrapping/OUTPUT")
    command.add_argument("--output-dir", default="OUTPUT")
    command.add_argument("--dataset")
    command.add_argument("--pair", help="Process only one YYYYMMDD-YYYYMMDD interferogram")
    command.add_argument("--degree", type=int, choices=(0, 1, 2), default=1)
    command.add_argument(
        "--orbit-model",
        choices=(
            "constant", "plane", "range_quadratic", "azimuth_quadratic",
            "range_cubic", "right_edge_quadratic",
            "additive_quadratic", "full_quadratic", "degree",
        ),
        default="plane",
        help=(
            "DEM-independent orbit-ramp model. Use right_edge_quadratic for "
            "right-edge residuals, or range_quadratic for conservative "
            "one-direction curvature; degree preserves the old --degree behavior."
        ),
    )
    command.add_argument("--mad-scale", type=float, default=3.5)
    command.add_argument("--max-iterations", type=int, default=5)
    command.add_argument("--stable-valid-fraction", type=float, default=0.8)
    command.add_argument("--stable-dispersion-percentile", type=float, default=40.0)
    command.add_argument(
        "--fit-edge-guard-pixels",
        type=int,
        default=0,
        help="Exclude this many pixels along each image edge from detrending fits while preserving full output coverage.",
    )
    command.add_argument(
        "--output-edge-guard-pixels",
        type=int,
        default=0,
        help="Mask this many pixels along each image edge in the final output and diagnostic plots.",
    )
    command.add_argument(
        "--reference-border-fraction",
        type=float,
        default=0.0,
        help="Fit detrending only from this fractional image border width, e.g. 0.2; use 0 for all stable pixels.",
    )
    command.add_argument(
        "--no-auto-reference-point",
        action="store_false",
        dest="auto_reference_point",
        help="Disable automatic high-coherence reference-point zeroing.",
    )
    command.add_argument(
        "--reference-window-pixels",
        type=int,
        default=9,
        help="Odd window size used to validate the representative automatic reference point.",
    )
    command.add_argument(
        "--reference-pixel-count",
        type=int,
        default=1000,
        help="Maximum number of high-coherence stable pixels used for robust reference zeroing.",
    )
    command.add_argument(
        "--reference-kml",
        help="Optional KML polygon where zero-deformation reference pixels are allowed to be selected.",
    )
    command.add_argument("--preserve-median", action="store_true")
    command.add_argument("--exclude-mask", help="Optional .npy bool mask excluded from fitting")
    command.add_argument(
        "--terrain-correction",
        action="store_true",
        help="Include a low-pass DEM term in the fitted correction; enabled by default.",
    )
    command.add_argument(
        "--no-terrain-correction",
        action="store_false",
        dest="terrain_correction",
        help="Disable the DEM-related stratified-atmosphere term.",
    )
    command.set_defaults(terrain_correction=True)
    command.add_argument("--terrain-degree", type=int, choices=(1, 2), default=2)
    command.add_argument(
        "--terrain-fit-method",
        choices=("spatial", "hybrid", "local", "binned", "pixel"),
        default="spatial",
        help="Fit DEM-correlated phase from after-orbit residuals using spatially varying DEM terms, a hybrid DEM-bin/local model, local coefficient field, elevation-bin medians, or robust pixels.",
    )
    command.add_argument("--terrain-bins", type=int, default=30)
    command.add_argument("--terrain-min-bin-pixels", type=int, default=100)
    command.add_argument(
        "--terrain-strength",
        type=float,
        default=0.3,
        help="Fraction of the accepted DEM-correlated term subtracted from the final interferogram.",
    )
    command.add_argument(
        "--terrain-max-range-fraction",
        type=float,
        default=0.6,
        help="Maximum DEM-term robust range as a fraction of the after-orbit residual robust range.",
    )
    command.add_argument(
        "--terrain-local-guard-pixels",
        type=int,
        default=96,
        help="Block size for rejecting DEM corrections that locally increase residual scatter.",
    )
    command.add_argument(
        "--terrain-local-guard-tolerance",
        type=float,
        default=0.15,
        help="Reject DEM correction if any valid local block increases robust scatter by more than this fraction.",
    )
    command.add_argument(
        "--terrain-local-guard-min-pixels",
        type=int,
        default=200,
        help="Minimum valid pixels required for a local DEM-correction guard block.",
    )
    command.add_argument(
        "--terrain-local-radius-pixels",
        type=int,
        default=80,
        help="Smoothing radius used to estimate the spatially varying local DEM coefficient.",
    )
    command.add_argument(
        "--terrain-smoothing-pixels",
        type=int,
        default=30,
        help="Low-pass smoothing radius applied to DEM before terrain fitting; use 0 to keep raw DEM.",
    )
    command.add_argument(
        "--terrain-smoothing-passes",
        type=int,
        default=1,
        help="Number of repeated DEM smoothing passes; default 1 keeps more terrain structure than the previous strong smoothing.",
    )
    command.add_argument(
        "--no-turbulent-qa",
        action="store_false",
        dest="turbulent_qa",
        help="Do not write the optional non-DEM low-frequency residual QA layer.",
    )
    command.set_defaults(turbulent_qa=True)
    command.add_argument(
        "--turbulent-smoothing-pixels",
        type=int,
        default=60,
        help="Low-pass smoothing radius for the diagnostic turbulent-atmosphere QA layer.",
    )
    command.add_argument("--geo-root", default="../InSAR_Filtering/OUTPUT")
    command.add_argument("--dem", default="../Data/DEM/rasters_USGS10m/output_USGS10m.tif")
    command.add_argument("--no-plots", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    config = DetrendConfig(
        degree=args.degree, orbit_model=args.orbit_model,
        mad_scale=args.mad_scale, max_iterations=args.max_iterations,
        stable_valid_fraction=args.stable_valid_fraction,
        stable_dispersion_percentile=args.stable_dispersion_percentile,
        fit_edge_guard_pixels=args.fit_edge_guard_pixels,
        output_edge_guard_pixels=args.output_edge_guard_pixels,
        reference_border_fraction=args.reference_border_fraction,
        auto_reference_point=args.auto_reference_point,
        reference_window_pixels=args.reference_window_pixels,
        reference_pixel_count=args.reference_pixel_count,
        reference_kml=args.reference_kml,
        preserve_median=args.preserve_median, exclude_mask=args.exclude_mask,
        make_plots=not args.no_plots,
        terrain_correction=args.terrain_correction,
        terrain_degree=args.terrain_degree,
        terrain_fit_method=args.terrain_fit_method,
        terrain_bins=args.terrain_bins,
        terrain_min_bin_pixels=args.terrain_min_bin_pixels,
        terrain_local_radius_pixels=args.terrain_local_radius_pixels,
        terrain_strength=args.terrain_strength,
        terrain_max_range_fraction=args.terrain_max_range_fraction,
        terrain_local_guard_pixels=args.terrain_local_guard_pixels,
        terrain_local_guard_tolerance=args.terrain_local_guard_tolerance,
        terrain_local_guard_min_pixels=args.terrain_local_guard_min_pixels,
        terrain_smoothing_pixels=args.terrain_smoothing_pixels,
        terrain_smoothing_passes=args.terrain_smoothing_passes,
        turbulent_qa=args.turbulent_qa,
        turbulent_smoothing_pixels=args.turbulent_smoothing_pixels,
        geo_root=args.geo_root, dem_path=args.dem,
        selected_pair=args.pair,
    )
    if args.dataset:
        result = run_dataset(Path(args.unwrapping_root) / args.dataset,
                             Path(args.output_dir) / args.dataset, config)
    else:
        result = run_all_datasets(args.unwrapping_root, args.output_dir, config)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
