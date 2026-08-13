"""Command-line interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .pipeline import TimeseriesConfig, run_all_datasets, run_dataset


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Connectivity-preserving SBAS time-series inversion")
    command.add_argument("--input-source", choices=("unwrapping", "detrending"), default="detrending",
                         help="Use detrended output by default, or direct unwrapping output")
    command.add_argument("--unwrapping-root", default="../InSAR_Unwrapping/OUTPUT",
                         help="InSAR_Unwrapping OUTPUT containing Dataset_* directories")
    command.add_argument("--detrending-root", default="../InSAR_Detrending/OUTPUT",
                         help="InSAR_Detrending OUTPUT containing Dataset_* directories")
    command.add_argument("--output-dir", default="OUTPUT")
    command.add_argument("--dataset", help="Process one Dataset_* only")
    command.add_argument("--max-baseline-days", type=int, default=12)
    command.add_argument("--wavelength", type=float, default=0.056, help="Radar wavelength in meters")
    command.add_argument("--phase-sign", type=float, choices=(-1.0, 1.0), default=1.0)
    command.add_argument("--no-plots", action="store_true")
    command.add_argument("--allow-low-quality-bridges", action="store_true",
                         dest="allow_low_quality_bridges",
                         help="Explicitly allow QC-bad necessary bridges with very low weights")
    command.add_argument("--no-allow-low-quality-bridges", action="store_false",
                         dest="allow_low_quality_bridges",
                         help="Stop if a QC-bad bridge is required for network connectivity")
    command.set_defaults(allow_low_quality_bridges=True)
    command.add_argument("--pick-point", action="store_true",
                         help="Open an existing dataset mean-velocity map and click a time-series point")
    return command


def main() -> None:
    args = parser().parse_args()
    if args.pick_point:
        if not args.dataset:
            raise SystemExit("--pick-point requires --dataset Dataset_NAME")
        from .interactive import interactive_point_timeseries
        interactive_point_timeseries(args.dataset, args.output_dir)
        return
    config = TimeseriesConfig(
        args.max_baseline_days, args.wavelength, args.phase_sign,
        not args.no_plots, not args.allow_low_quality_bridges,
    )
    input_root = args.detrending_root if args.input_source == "detrending" else args.unwrapping_root
    if args.dataset:
        result = run_dataset(Path(input_root) / args.dataset,
                             Path(args.output_dir) / args.dataset, config)
    else:
        result = run_all_datasets(input_root, args.output_dir, config)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
