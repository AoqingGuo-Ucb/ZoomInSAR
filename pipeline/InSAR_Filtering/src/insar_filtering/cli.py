"""Command-line interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .pipeline import FilterConfig, run_all_datasets, run_dataset


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Filter cropped complex InSAR interferograms")
    command.add_argument("--roi-dir", default="../Data/ROI")
    command.add_argument("--output-dir", default="OUTPUT")
    command.add_argument("--dataset", help="Process only one Dataset_* directory")
    command.add_argument("--pair", help="Process only one YYYYMMDD-YYYYMMDD interferogram")
    command.add_argument("--search-radius", type=int, default=7)
    command.add_argument("--patch-radius", type=int, default=3)
    command.add_argument("--bandwidth", type=float, default=0.70)
    command.add_argument("--goldstein-strength", type=float, default=0.40)
    command.add_argument("--no-plots", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    config = FilterConfig(args.search_radius, args.patch_radius, args.bandwidth,
                          args.goldstein_strength, not args.no_plots)
    if args.dataset:
        result = run_dataset(Path(args.roi_dir) / args.dataset,
                             Path(args.output_dir) / args.dataset, config, args.pair)
    else:
        result = run_all_datasets(args.roi_dir, args.output_dir, config, args.pair)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
