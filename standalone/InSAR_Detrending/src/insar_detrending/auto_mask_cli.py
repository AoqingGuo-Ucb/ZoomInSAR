"""Command-line interface for automatic deformation-mask generation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .auto_mask import AutoMaskConfig, build_mask_from_velocity


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description="Build an automatic deformation exclusion mask from a preliminary velocity map"
    )
    command.add_argument("--timeseries-root", default="../InSAR_Timeseries/OUTPUT_preliminary_unwrapped")
    command.add_argument("--filtering-root", default="../InSAR_Filtering/OUTPUT")
    command.add_argument("--output-dir", default="OUTPUT_auto_masks")
    command.add_argument("--dataset", required=True, help="Dataset_* name to process")
    command.add_argument("--velocity-percentile", type=float, default=85.0)
    command.add_argument("--min-velocity", type=float, default=0.02, help="Minimum |velocity| threshold in m/year")
    command.add_argument("--dilation-pixels", type=int, default=6)
    command.add_argument("--min-component-pixels", type=int, default=80)
    command.add_argument("--max-kml-polygons", type=int, default=20)
    command.add_argument("--min-kml-polygon-area-pixels", type=float, default=100.0)
    return command


def main() -> None:
    args = parser().parse_args()
    config = AutoMaskConfig(
        velocity_percentile=args.velocity_percentile,
        min_velocity_m_per_year=args.min_velocity,
        dilation_pixels=args.dilation_pixels,
        min_component_pixels=args.min_component_pixels,
        max_kml_polygons=args.max_kml_polygons,
        min_kml_polygon_area_pixels=args.min_kml_polygon_area_pixels,
    )
    result = build_mask_from_velocity(
        Path(args.timeseries_root) / args.dataset,
        Path(args.filtering_root) / args.dataset,
        Path(args.output_dir) / args.dataset,
        config,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
