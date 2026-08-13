"""Command line interface."""

from __future__ import annotations

import argparse

from .crop import NoOverlapError, crop_dataset, crop_one_roi


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Crop GAMMA InSAR files using KML bounds")
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
    return parser


def main() -> None:
    args = build_parser().parse_args()
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
