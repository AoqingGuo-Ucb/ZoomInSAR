"""Binary raster and GAMMA parameter-file helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import BinaryIO, Iterator

import numpy as np


def parse_gamma_shape(par_file: Path) -> tuple[int, int]:
    """Return (azimuth_lines, range_samples) from a GAMMA .par file."""
    text = par_file.read_text(encoding="utf-8", errors="replace")
    width = re.search(r"^range_samples:\s+(\d+)", text, re.MULTILINE)
    lines = re.search(r"^azimuth_lines:\s+(\d+)", text, re.MULTILINE)
    if not width or not lines:
        raise ValueError(f"Cannot find range_samples/azimuth_lines in {par_file}")
    return int(lines.group(1)), int(width.group(1))


def discover_shape(geo_dir: Path, shape: tuple[int, int] | None = None) -> tuple[int, int]:
    if shape is not None:
        return shape
    par_files = sorted(geo_dir.glob("*.par"))
    if not par_files:
        raise FileNotFoundError(
            f"No .par file found in {geo_dir}; provide --lines and --width"
        )
    shapes = {parse_gamma_shape(path) for path in par_files}
    if len(shapes) != 1:
        raise ValueError(f"GEO parameter files contain inconsistent shapes: {shapes}")
    return shapes.pop()


def raster_memmap(path: Path, shape: tuple[int, int], dtype: str) -> np.memmap:
    expected = int(np.prod(shape)) * np.dtype(dtype).itemsize
    actual = path.stat().st_size
    if actual != expected:
        raise ValueError(f"Unexpected size for {path}: {actual} bytes; expected {expected}")
    return np.memmap(path, dtype=dtype, mode="r", shape=shape, order="C")


def freadbkB(
    infile: str | Path,
    lines: int,
    bkformat: str = "float32",
    r0: int = 0,
    rN: int = 0,
    c0: int = 0,
    cN: int = 0,
) -> tuple[np.ndarray, int]:
    """Read a big-endian GAMMA raster (compatible with the original freadbkB API).

    Indices are one-based and inclusive, as in the original MATLAB/Python helper;
    zero selects the first/last available index.
    """
    path = Path(infile)
    fmt = bkformat.lower()
    if fmt == "mph":
        fmt = "cpxfloat32"
    if fmt == "hgt":
        raise ValueError("Use a dedicated height reader for 'hgt' format")
    is_complex = fmt.startswith("cpx")
    # NumPy accepts aliases such as ``float32`` but not the combined string
    # ``>float32``; construct first, then apply GAMMA's big-endian byte order.
    scalar_dtype = np.dtype(fmt[3:] if is_complex else fmt).newbyteorder(">")
    bytes_per_pixel = scalar_dtype.itemsize * (2 if is_complex else 1)
    pixels = path.stat().st_size // bytes_per_pixel
    if pixels * bytes_per_pixel != path.stat().st_size or pixels % lines:
        raise ValueError("File size is incompatible with lines and format")
    width = pixels // lines
    r0, rN = (r0 or 1), (rN or lines)
    c0, cN = (c0 or 1), (cN or width)
    if not (1 <= r0 <= rN <= lines and 1 <= c0 <= cN <= width):
        raise ValueError("Invalid one-based read range")
    raw_shape = (lines, width * (2 if is_complex else 1))
    raw = np.memmap(path, dtype=scalar_dtype, mode="r", shape=raw_shape)
    if is_complex:
        part = np.asarray(raw[r0 - 1 : rN, 2 * (c0 - 1) : 2 * cN])
        data = part[:, 0::2] + 1j * part[:, 1::2]
    else:
        data = np.asarray(raw[r0 - 1 : rN, c0 - 1 : cN])
    return data, data.size


def write_window(
    source: Path,
    destination: Path,
    shape: tuple[int, int],
    dtype: str,
    window: tuple[int, int, int, int],
) -> None:
    """Write a rectangular row/column window while preserving binary dtype/order."""
    r0, r1, c0, c1 = window
    source_array = raster_memmap(source, shape, dtype)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as stream:
        for row in range(r0, r1):
            np.asarray(source_array[row, c0:c1]).tofile(stream)


def iter_chunks(stream: BinaryIO, size: int = 1024 * 1024) -> Iterator[bytes]:
    while chunk := stream.read(size):
        yield chunk
