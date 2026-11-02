"""Command-line interface for GSLC pixel-offset tracking."""

from __future__ import annotations

import argparse
import json

from .pipeline import PixelOffsetConfig, run_pixel_offset


def main() -> None:
    parser = argparse.ArgumentParser(description="Track map-grid offsets in consecutive NISAR GSLC products")
    parser.add_argument("--gslc-dir", required=True, help="Folder containing original NISAR GSLC .h5 files")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--frequency", default="A", choices=("A", "B", "a", "b"))
    parser.add_argument("--polarization", default="HH", help="Complex GSLC polarization layer, e.g. HH or VV")
    parser.add_argument("--window-pixels", type=int, default=64)
    parser.add_argument("--search-pixels", type=int, default=12)
    parser.add_argument("--stride-pixels", type=int, default=32)
    parser.add_argument("--min-correlation", type=float, default=0.25)
    parser.add_argument("--max-baseline-days", type=int, default=24)
    args = parser.parse_args()
    if not 0 <= args.min_correlation <= 1:
        parser.error("--min-correlation must be between 0 and 1")
    result = run_pixel_offset(
        args.gslc_dir, args.output_dir, args.dataset,
        PixelOffsetConfig(
            frequency=args.frequency.upper(), polarization=args.polarization.upper(),
            window_pixels=args.window_pixels, search_pixels=args.search_pixels,
            stride_pixels=args.stride_pixels, min_correlation=args.min_correlation,
            max_baseline_days=args.max_baseline_days,
        ),
    )
    print(json.dumps(result, indent=2))
