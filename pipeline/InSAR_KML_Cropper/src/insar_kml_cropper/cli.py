"""Command line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from .crop import NoOverlapError, crop_dataset, crop_one_roi
from .sweets import import_sweets_dataset


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Crop GAMMA rasters or import SWEETS GeoTIFF interferograms")
    parser.add_argument(
        "--data-dir",
        default="../Data",
        help="Directory containing INT, GEO and ROI (default: ../Data from the package root)",
    )
    parser.add_argument("--kml", help="Process one KML; default: all ROI/*.kml")
    parser.add_argument("--margin", type=float, default=0.20, help="Expansion per side (default: 0.20)")
    parser.add_argument("--lines", type=int, help="Source rows if no GEO .par is available")
    parser.add_argument("--width", type=int, help="Source columns if no GEO .par is available")
    parser.add_argument("--overwrite", action="store_true", help="Allow replacing existing output files")
    parser.add_argument("--pair", help="Process only one YYYYMMDD-YYYYMMDD interferogram")
    parser.add_argument("--sweets-dir", help="Import SWEETS GeoTIFFs from this directory instead of GAMMA files")
    parser.add_argument("--dataset", help="Dataset name used with --sweets-dir")
    parser.add_argument("--phase-pattern", default="*.tif*", help="Recursive GeoTIFF glob for SWEETS wrapped phase")
    parser.add_argument("--coherence-pattern", default="*.tif*", help="Recursive GeoTIFF glob for SWEETS coherence")
    crop_choice = parser.add_mutually_exclusive_group()
    crop_choice.add_argument("--crop-kml", help="KML used to crop SWEETS input")
    crop_choice.add_argument("--no-crop", action="store_true", help="Import the complete SWEETS grid (the default)")
    parser.add_argument(
        "--skip-unpaired", action="store_true",
        help="Exclude SWEETS interferograms that do not have both phase and coherence GeoTIFFs.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.sweets_dir:
        if not args.dataset:
            raise SystemExit("--dataset is required with --sweets-dir")
        output = import_sweets_dataset(
            args.sweets_dir, Path(args.data_dir) / "ROI", args.dataset,
            phase_pattern=args.phase_pattern, coherence_pattern=args.coherence_pattern,
            crop_kml=args.crop_kml, margin=args.margin, overwrite=args.overwrite,
            skip_unpaired=args.skip_unpaired,
        )
        print(f"Created: {output}")
        return
    if (args.lines is None) != (args.width is None):
        raise SystemExit("--lines and --width must be provided together")
    shape = (args.lines, args.width) if args.lines else None
    try:
        outputs = (
            [crop_one_roi(args.data_dir, args.kml, args.margin, shape, args.overwrite, args.pair)]
            if args.kml
            else crop_dataset(args.data_dir, args.margin, shape, args.overwrite, args.pair)
        )
    except NoOverlapError as exc:
        raise SystemExit(str(exc)) from exc
    for output in outputs:
        print(f"Created: {output}")


if __name__ == "__main__":
    main()
