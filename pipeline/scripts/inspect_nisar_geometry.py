#!/usr/bin/env python3
"""List NISAR/CSLC HDF5 geometry datasets below a project directory."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py


KEYWORDS = ("los", "look", "incidence", "heading", "azimuth", "geolocation", "radargrid")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.expanduser().resolve()
    files = sorted({*root.rglob("*.h5"), *root.rglob("*.H5"), *root.rglob("*.hdf5"), *root.rglob("*.HDF5")})
    if not files:
        raise SystemExit(f"No HDF5 files found below {root}")
    print(f"Found {len(files)} HDF5 file(s) below {root}")
    matched_files = 0
    for path in files:
        matches: list[tuple[str, tuple[int, ...]]] = []
        try:
            with h5py.File(path, "r") as handle:
                def visit(name: str, item: object) -> None:
                    if isinstance(item, h5py.Dataset) and any(word in name.lower() for word in KEYWORDS):
                        matches.append((name, tuple(item.shape)))
                handle.visititems(visit)
        except OSError as error:
            print(f"Unreadable: {path} ({error})")
            continue
        if matches:
            matched_files += 1
            print(f"\n{path}")
            for name, shape in matches:
                print(f"  {name}  shape={shape}")
    if not matched_files:
        print("No LOS/look/incidence/heading/azimuth datasets were found in these files.")


if __name__ == "__main__":
    main()
