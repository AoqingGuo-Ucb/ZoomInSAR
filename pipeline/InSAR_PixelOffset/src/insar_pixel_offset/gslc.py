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


def _all_complex_paths(handle: h5py.File) -> list[str]:
    """Return every 2-D complex dataset in an HDF5 SLC/CSLC product."""
    paths: list[str] = []

    def visit(name: str, item: h5py.Dataset | h5py.Group) -> None:
        if not isinstance(item, h5py.Dataset) or item.ndim != 2:
            return
        if not np.issubdtype(item.dtype, np.complexfloating):
            return
        paths.append(name)

    handle.visititems(visit)
    return sorted(paths)


def _complex_paths(handle: h5py.File, frequency: str, polarization: str) -> list[str]:
    freq = f"frequency{frequency.upper()}"
    pol = polarization.upper()
    return [
        name for name in _all_complex_paths(handle)
        if freq in name.split("/") and name.split("/")[-1].upper() == pol
    ]


def locate_complex_dataset(handle: h5py.File, frequency: str, polarization: str) -> str:
    candidates = _complex_paths(handle, frequency, polarization)
    if candidates:
        return candidates[0]

    # Sentinel-1 SWEETS/COMPASS CSLC HDF5 products need not use the NISAR
    # frequencyA/HH group convention. Prefer the requested polarization, then
    # fall back only when the product has one unambiguous complex image layer.
    all_candidates = _all_complex_paths(handle)
    pol = polarization.upper()
    pol_candidates = [name for name in all_candidates if name.split("/")[-1].upper() == pol]
    if len(pol_candidates) == 1:
        return pol_candidates[0]
    if len(all_candidates) == 1:
        return all_candidates[0]
    available = ", ".join(all_candidates[:8]) or "none"
    raise ValueError(
        f"No unambiguous complex layer for frequency{frequency.upper()}/{pol} in {handle.filename}. "
        f"Available complex datasets: {available}"
    )


def coordinate_steps(handle: h5py.File, dataset_path: str) -> tuple[float | None, float | None]:
    """Return x/y map-coordinate spacing for a GSLC image, when provided."""
    candidate = handle[dataset_path].parent
    while True:
        x = candidate.get("xCoordinates")
        y = candidate.get("yCoordinates")
        if x is None:
            x = candidate.get("x_coordinates")
        if y is None:
            y = candidate.get("y_coordinates")
        if x is not None and y is not None and len(x) >= 2 and len(y) >= 2:
            return float(np.nanmedian(np.diff(x[...]))), float(np.nanmedian(np.diff(y[...])))
        if candidate.name == "/":
            break
        candidate = candidate.parent
    return None, None
