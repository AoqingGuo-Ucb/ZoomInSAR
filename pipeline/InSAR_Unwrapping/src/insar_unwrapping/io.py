"""GAMMA binary raster input/output."""

from __future__ import annotations

import re
import json
from pathlib import Path

import numpy as np


def parse_shape(par_file: str | Path) -> tuple[int, int]:
    text = Path(par_file).read_text(encoding="utf-8", errors="replace")
    width = re.search(r"^range_samples:\s+(\d+)", text, re.MULTILINE)
    lines = re.search(r"^azimuth_lines:\s+(\d+)", text, re.MULTILINE)
    if not width or not lines:
        raise ValueError(f"Cannot read range_samples/azimuth_lines from {par_file}")
    return int(lines.group(1)), int(width.group(1))


def discover_shape(data_dir: Path, explicit: tuple[int, int] | None = None) -> tuple[int, int]:
    if explicit:
        return explicit
    candidates = sorted((data_dir / "GEO").glob("*.par")) + sorted(data_dir.glob("*.par"))
    metadata = data_dir / "crop_metadata.json"
    if not candidates and metadata.exists():
        values = json.loads(metadata.read_text(encoding="utf-8"))["output_shape_lines_width"]
        return int(values[0]), int(values[1])
    if not candidates:
        raise FileNotFoundError("No GAMMA .par file found; provide --lines and --width")
    shapes = {parse_shape(path) for path in candidates}
    if len(shapes) != 1:
        raise ValueError(f"Inconsistent raster shapes in parameter files: {sorted(shapes)}")
    return shapes.pop()


def read_gamma(path: str | Path, shape: tuple[int, int], dtype: str) -> np.ndarray:
    path = Path(path)
    disk_dtype = np.dtype(dtype).newbyteorder(">")
    expected = int(np.prod(shape)) * disk_dtype.itemsize
    if path.stat().st_size != expected:
        raise ValueError(f"{path} has {path.stat().st_size} bytes; expected {expected}")
    return np.asarray(np.memmap(path, dtype=disk_dtype, mode="r", shape=shape))


def read_wrapped(path: str | Path, shape: tuple[int, int]) -> np.ndarray:
    """Read a big-endian complex64 interferogram and return phase in radians."""
    return np.angle(read_gamma(path, shape, "c8")).astype(np.float64)


def read_coherence(path: str | Path, shape: tuple[int, int]) -> np.ndarray:
    return read_gamma(path, shape, "f4").astype(np.float64)


def write_unwrapped(path: str | Path, phase: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.asarray(phase, dtype=">f4").tofile(path)


def pair_name(path: str | Path) -> tuple[str, str]:
    match = re.search(r"(?<!\d)(\d{8})[-_](\d{8})(?!\d)", Path(path).name)
    if not match:
        raise ValueError(f"Filename does not contain YYYYMMDD-YYYYMMDD: {path}")
    return match.group(1), match.group(2)
