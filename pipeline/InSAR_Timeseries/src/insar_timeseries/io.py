"""Unwrapped interferogram and metadata input/output."""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np


DATE_PAIR = re.compile(r"(?<!\d)(\d{8})[-_](\d{8})(?!\d)")


def pair_from_name(path: str | Path) -> tuple[datetime, datetime]:
    match = DATE_PAIR.search(Path(path).name)
    if not match:
        raise ValueError(f"No YYYYMMDD-YYYYMMDD pair in {path}")
    first, second = (datetime.strptime(value, "%Y%m%d") for value in match.groups())
    if second <= first:
        raise ValueError(f"Date pair is not chronological: {path}")
    return first, second


def read_shape(dataset_output: Path) -> tuple[int, int]:
    summary = json.loads((dataset_output / "run_summary.json").read_text(encoding="utf-8"))
    values = summary["shape_lines_width"]
    return int(values[0]), int(values[1])


def discover_unwrapped(dataset_output: Path) -> list[tuple[Path, datetime, datetime]]:
    result = [(path, *pair_from_name(path)) for path in sorted((dataset_output / "unwrapped").glob("*.unw"))]
    if not result:
        raise FileNotFoundError(f"No unwrapped/*.unw files in {dataset_output}")
    pairs = [(a, b) for _, a, b in result]
    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate date pairs in unwrapped inputs")
    return result


def read_unwrapped(path: Path, shape: tuple[int, int]) -> np.ndarray:
    expected = int(np.prod(shape)) * 4
    if path.stat().st_size != expected:
        raise ValueError(f"{path}: {path.stat().st_size} bytes, expected {expected}")
    return np.asarray(np.memmap(path, dtype=">f4", mode="r", shape=shape)).astype(float)
