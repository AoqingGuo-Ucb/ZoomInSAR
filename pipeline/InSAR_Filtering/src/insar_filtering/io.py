"""GAMMA binary input/output and cropped-dataset discovery."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np


def discover_shape(dataset: Path) -> tuple[int, int]:
    metadata = dataset / "crop_metadata.json"
    if metadata.exists():
        values = json.loads(metadata.read_text(encoding="utf-8"))["output_shape_lines_width"]
        return int(values[0]), int(values[1])
    par_files = sorted((dataset / "GEO").glob("*.par"))
    if not par_files:
        raise FileNotFoundError(f"No crop_metadata.json or GEO/*.par in {dataset}")
    text = par_files[0].read_text(encoding="utf-8", errors="replace")
    lines = re.search(r"^azimuth_lines:\s+(\d+)", text, re.MULTILINE)
    width = re.search(r"^range_samples:\s+(\d+)", text, re.MULTILINE)
    if not lines or not width:
        raise ValueError(f"Cannot parse shape from {par_files[0]}")
    return int(lines.group(1)), int(width.group(1))


def read_gamma(path: Path, shape: tuple[int, int], dtype: str) -> np.ndarray:
    disk_dtype = np.dtype(dtype).newbyteorder(">")
    expected = int(np.prod(shape)) * disk_dtype.itemsize
    if path.stat().st_size != expected:
        raise ValueError(f"{path}: {path.stat().st_size} bytes, expected {expected}")
    return np.asarray(np.memmap(path, dtype=disk_dtype, mode="r", shape=shape))


def write_gamma(path: Path, values: np.ndarray, dtype: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.asarray(values, dtype=np.dtype(dtype).newbyteorder(">")).tofile(path)


def pair_name(path: Path) -> tuple[str, str]:
    match = re.search(r"(?<!\d)(\d{8})[-_](\d{8})(?!\d)", path.name)
    if not match:
        raise ValueError(f"No date pair in {path.name}")
    return match.group(1), match.group(2)
