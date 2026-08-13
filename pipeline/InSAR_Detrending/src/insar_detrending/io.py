"""Read and write InSAR_Unwrapping-compatible products."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def read_shape(dataset: Path) -> tuple[int, int]:
    values = json.loads((dataset / "run_summary.json").read_text(encoding="utf-8"))["shape_lines_width"]
    return int(values[0]), int(values[1])


def discover_unwrapped(dataset: Path) -> list[Path]:
    files = sorted((dataset / "unwrapped").glob("*.unw"))
    if not files:
        raise FileNotFoundError(f"No unwrapped/*.unw in {dataset}")
    return files


def read_unwrapped(path: Path, shape: tuple[int, int]) -> np.ndarray:
    expected = int(np.prod(shape)) * 4
    if path.stat().st_size != expected:
        raise ValueError(f"{path}: {path.stat().st_size} bytes, expected {expected}")
    return np.asarray(np.memmap(path, dtype=">f4", mode="r", shape=shape)).astype(float)


def write_unwrapped(path: Path, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.asarray(values, dtype=">f4").tofile(path)
