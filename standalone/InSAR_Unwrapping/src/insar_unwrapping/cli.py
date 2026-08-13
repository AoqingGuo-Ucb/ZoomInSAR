"""Command-line interface."""

from __future__ import annotations

import argparse
import json

from .pipeline import PipelineConfig, run_all_datasets, run_pipeline


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="MCF unwrap and closure-correct GAMMA interferograms")
    command.add_argument("--input-source", choices=("raw", "filtered"), default="raw",
                         help="Use raw cropped data or InSAR_Filtering output")
    command.add_argument("--roi-dir", default="../Data/ROI", help="Raw parent directory containing Dataset_* inputs")
    command.add_argument("--filtered-root", default="../InSAR_Filtering/OUTPUT",
                         help="Parent directory containing filtered Dataset_* inputs")
    command.add_argument("--output-dir", default="OUTPUT", help="Root directory for all generated products")
    command.add_argument("--dataset", help="Process one Dataset_* name instead of every dataset")
    command.add_argument("--pair", help="Process only one YYYYMMDD-YYYYMMDD interferogram")
    command.add_argument("--lines", type=int)
    command.add_argument("--width", type=int)
    command.add_argument("--coherence-threshold", type=float, default=0.05)
    command.add_argument("--max-branch-cut-length", type=float, default=15.0,
                         help="Maximum residue-pair or residue-boundary cut length in pixels")
    command.add_argument("--no-spatial-repair", action="store_true")
    command.add_argument("--spatial-sigma", type=float, default=12.0)
    command.add_argument("--spatial-min-pixels", type=int, default=40)
    command.add_argument("--no-region-graph-repair", action="store_true",
                         help="Disable image-based spatial integer-cycle region repair")
    command.add_argument("--region-edge-sigma", type=float, default=1.5)
    command.add_argument("--region-edge-threshold-cycles", type=float, default=0.25)
    command.add_argument("--region-boundary-width", type=int, default=5,
                         help="Inner/outer multi-pixel bands for spatial region comparison")
    command.add_argument("--region-max-iterations", type=int, default=3)
    command.add_argument("--region-open-boundary-max-gap", type=int, default=12)
    command.add_argument("--region-open-boundary-snap", type=int, default=10)
    command.add_argument("--no-edge-strip-repair",action="store_true")
    command.add_argument("--edge-strip-min-pixels",type=int,default=40)
    command.add_argument("--edge-strip-boundary-width",type=int,default=5)
    command.add_argument("--no-closure-repair", action="store_true")
    command.add_argument("--closure-iterations", type=int, default=8)
    command.add_argument("--closure-min-pixels", type=int, default=20)
    command.add_argument("--closure-boundary-width", type=int, default=3,
                         help="Inner/outer boundary-band width used to validate each closure repair")
    command.add_argument("--no-auto-exclude-bad", action="store_true",
                         help="Keep severe spatial outliers even when safely removable from the date network")
    command.add_argument("--qc-max-residual-fraction", type=float, default=0.04)
    command.add_argument("--qc-max-phase-range-cycles", type=float, default=2.0)
    command.add_argument("--no-plots", action="store_true", help="Do not create PNG diagnostic figures")
    command.add_argument("--no-water-mask", action="store_true",
                         help="Disable the default longitude/latitude water mask")
    command.add_argument("--dem", default="../Data/DEM/rasters_USGS10m/output_USGS10m.tif",
                         help="DEM sampled at lon/lat pixels for water masking")
    command.add_argument("--water-max-elevation", type=float, default=0.0,
                         help="DEM elevation at or below this value is water (default: 0 m)")
    return command


def main() -> None:
    args = parser().parse_args()
    if (args.lines is None) != (args.width is None):
        raise SystemExit("--lines and --width must be supplied together")
    shape = (args.lines, args.width) if args.lines else None
    config = PipelineConfig(
        coherence_threshold=args.coherence_threshold,
        max_branch_cut_length=args.max_branch_cut_length,
        spatial_repair=not args.no_spatial_repair,
        spatial_sigma=args.spatial_sigma,
        spatial_min_pixels=args.spatial_min_pixels,
        region_graph_repair=not args.no_region_graph_repair,
        region_edge_sigma=args.region_edge_sigma,
        region_edge_threshold_cycles=args.region_edge_threshold_cycles,
        region_boundary_width=args.region_boundary_width,
        region_max_iterations=args.region_max_iterations,
        region_open_boundary_max_gap=args.region_open_boundary_max_gap,
        region_open_boundary_snap=args.region_open_boundary_snap,
        edge_strip_repair=not args.no_edge_strip_repair,
        edge_strip_min_pixels=args.edge_strip_min_pixels,
        edge_strip_boundary_width=args.edge_strip_boundary_width,
        closure_repair=not args.no_closure_repair,
        closure_iterations=args.closure_iterations,
        closure_min_pixels=args.closure_min_pixels,
        closure_boundary_width=args.closure_boundary_width,
        auto_exclude_bad=not args.no_auto_exclude_bad,
        qc_max_residual_fraction=args.qc_max_residual_fraction,
        qc_max_phase_range_cycles=args.qc_max_phase_range_cycles,
        make_plots=not args.no_plots,
        mask_water=not args.no_water_mask,
        dem_path=args.dem,
        water_max_elevation=args.water_max_elevation,
        selected_pair=args.pair,
    )
    input_root = args.filtered_root if args.input_source == "filtered" else args.roi_dir
    if args.dataset:
        data_dir = __import__("pathlib").Path(input_root) / args.dataset
        summary = run_pipeline(data_dir, __import__("pathlib").Path(args.output_dir) / args.dataset, shape, config)
    else:
        if shape:
            raise SystemExit("--lines/--width are supported only together with --dataset")
        summary = run_all_datasets(input_root, args.output_dir, config)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
