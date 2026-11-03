"""Discovery and safe reading of geocoded complex SLC/GSLC scenes."""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np


_DATE = re.compile(r"(20\d{6})(?:T\d{6})?")


@dataclass(frozen=True)
class GslcScene:
    path: Path
    date: datetime
    format: str


def acquisition_date(path: Path) -> datetime:
    match = _DATE.search(path.name)
    if not match:
        raise ValueError(f"Could not find a YYYYMMDD acquisition date in {path.name}")
    return datetime.strptime(match.group(1), "%Y%m%d")


def discover_scenes(directory: str | Path) -> list[GslcScene]:
    directory = Path(directory).expanduser()
    scenes = []
    hdf_files = {*directory.rglob("*.h5"), *directory.rglob("*.hdf5")}
    # SWEETS writes the geocoded, complex Sentinel-1 SLC to these GeoTIFFs.
    slc_tifs = {*directory.rglob("*.slc.tif"), *directory.rglob("*.SLC.tif")}
    for path in sorted(hdf_files | slc_tifs):
        try:
            scenes.append(GslcScene(
                path, acquisition_date(path),
                "geotiff" if path in slc_tifs else "hdf5",
            ))
        except ValueError:
            continue
    scenes.sort(key=lambda scene: (scene.date, scene.path.name))
    if len(scenes) < 2:
        raise FileNotFoundError(
            f"Need at least two dated NISAR GSLC .h5 or SWEETS *.slc.tif files below {directory}"
        )
    formats = {scene.format for scene in scenes}
    if len(formats) != 1:
        raise ValueError(
            "Do not mix NISAR GSLC HDF5 files and SWEETS Sentinel SLC GeoTIFFs in one run; "
            "choose one scene directory."
        )
    return scenes


def _complex_paths(handle: h5py.File, frequency: str, polarization: str) -> list[str]:
    paths: list[str] = []
    freq = f"frequency{frequency.upper()}"
    pol = polarization.upper()

    def visit(name: str, item: h5py.Dataset | h5py.Group) -> None:
        if not isinstance(item, h5py.Dataset) or item.ndim != 2:
            return
        if not np.issubdtype(item.dtype, np.complexfloating):
            return
        parts = name.split("/")
        if freq in parts and parts[-1].upper() == pol:
            paths.append(name)

    handle.visititems(visit)
    return sorted(paths)


def locate_complex_dataset(handle: h5py.File, frequency: str, polarization: str) -> str:
    candidates = _complex_paths(handle, frequency, polarization)
    if not candidates:
        raise ValueError(
            f"No complex frequency{frequency.upper()}/{polarization.upper()} GSLC layer in {handle.filename}"
        )
    return candidates[0]


def coordinate_steps(handle: h5py.File, dataset_path: str) -> tuple[float | None, float | None]:
    """Return x/y map-coordinate spacing for a GSLC image, when provided."""
    group = handle[dataset_path].parent
    x = group.get("xCoordinates")
    y = group.get("yCoordinates")
    if x is None or y is None or len(x) < 2 or len(y) < 2:
        return None, None
    return float(np.nanmedian(np.diff(x[...]))), float(np.nanmedian(np.diff(y[...])))
