"""Import georeferenced SWEETS interferogram GeoTIFFs into ZoomInSAR datasets."""

from __future__ import annotations

import json
import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import transform

from .kml import expand_bounds, kml_bounds


_PAIR = re.compile(r"(?<!\d)(\d{8})[-_](\d{8})(?!\d)")
_COHERENCE = re.compile(
    r"(?:^|[_\-.])(?:coh(?:erence)?|corr(?:elation)?|cor)(?=[_\-.]|$)",
    re.IGNORECASE,
)
_PHASE = re.compile(r"(?:^|[_\-.])(int|ifg|interferogram|wrapped|phase)(?=[_\-.]|$)", re.IGNORECASE)
_AUXILIARY = re.compile(
    r"(?:^|[_\-.])(?:temporal[_\-.]?coherence|similarity|shp[_\-.]?counts|mask)(?=[_\-.]|$)",
    re.IGNORECASE,
)


def _pair(path: Path) -> tuple[str, str]:
    match = _PAIR.search(path.name)
    if not match:
        raise ValueError(f"No YYYYMMDD-YYYYMMDD date pair in {path.name}")
    return match.group(1), match.group(2)


def _discover(root: Path, pattern: str, kind: str) -> dict[tuple[str, str], Path]:
    candidates = sorted(root.rglob(pattern))
    selected: dict[tuple[str, str], Path] = {}
    for path in candidates:
        name = path.name
        # Dolphin/SWEETS also writes masks and quality layers with date pairs.
        # They are not pairwise coherence rasters and must never enter this workflow.
        if _AUXILIARY.search(name):
            continue
        # SWEETS commonly writes ``<pair>.int.tif`` and
        # ``<pair>.int.cor.tif``.  Test the compound suffix first so the
        # ``int`` token in the coherence filename can never be mistaken for
        # a second phase product.
        suffixes = "".join(path.suffixes).lower()
        is_coherence = (
            suffixes.endswith(".int.cor.tif")
            or suffixes.endswith(".int.cor.tiff")
            or bool(_COHERENCE.search(name))
        )
        is_phase = bool(_PHASE.search(name)) and not is_coherence
        if (kind == "coherence" and not is_coherence) or (kind == "phase" and not is_phase):
            continue
        pair = _pair(path)
        if pair in selected:
            raise ValueError(f"More than one {kind} GeoTIFF found for {pair}: {selected[pair]} and {path}")
        selected[pair] = path
    if not selected:
        hint = "--phase-pattern" if kind == "phase" else "--coherence-pattern"
        raise FileNotFoundError(f"No SWEETS {kind} GeoTIFFs found below {root}; adjust {hint} if needed")
    return selected


def _read_phase(path: Path) -> tuple[np.ndarray, dict]:
    with rasterio.open(path) as source:
        if source.count != 1:
            raise ValueError(f"{path} must have one phase band, found {source.count}")
        values = source.read(1)
        profile = {"crs": source.crs, "transform": source.transform, "shape": source.shape, "nodata": source.nodata}
    if np.iscomplexobj(values):
        phase = np.angle(values)
        valid = np.isfinite(values.real) & np.isfinite(values.imag) & (np.abs(values) > 0)
    else:
        phase = values.astype(float)
        valid = np.isfinite(phase)
    if profile["nodata"] is not None:
        valid &= ~np.isclose(values, profile["nodata"])
    return np.where(valid, np.exp(1j * phase), 0).astype(np.complex64), profile


def _read_coherence(path: Path, profile: dict) -> np.ndarray:
    with rasterio.open(path) as source:
        if source.count != 1:
            raise ValueError(f"{path} must have one coherence band, found {source.count}")
        if source.shape != profile["shape"] or source.crs != profile["crs"] or source.transform != profile["transform"]:
            raise ValueError(f"{path} is not on the same grid as its phase GeoTIFF")
        values = source.read(1).astype(np.float32)
        nodata = source.nodata
    valid = np.isfinite(values)
    if nodata is not None:
        valid &= ~np.isclose(values, nodata)
    # Numerical noise outside the physical coherence interval should not break filtering.
    return np.where(valid, np.clip(values, 0.0, 1.0), 0.0).astype(np.float32)


def _lon_lat(profile: dict) -> tuple[np.ndarray, np.ndarray]:
    crs = profile["crs"]
    if crs is None:
        raise ValueError("SWEETS GeoTIFF has no CRS; cannot create required longitude/latitude geometry")
    rows, cols = np.indices(profile["shape"])
    xs, ys = rasterio.transform.xy(profile["transform"], rows, cols, offset="center")
    lon, lat = transform(crs, "EPSG:4326", xs, ys)
    return np.asarray(lon, dtype=np.float32).reshape(profile["shape"]), np.asarray(lat, dtype=np.float32).reshape(profile["shape"])


def _crop_mask(longitude: np.ndarray, latitude: np.ndarray, kml: Path, margin: float) -> tuple[slice, slice]:
    west, east, south, north = expand_bounds(kml_bounds(kml), margin)
    inside = (longitude >= west) & (longitude <= east) & (latitude >= south) & (latitude <= north)
    rows, cols = np.nonzero(inside)
    if not rows.size:
        raise ValueError(f"KML does not overlap SWEETS GeoTIFF: {kml}")
    return slice(rows.min(), rows.max() + 1), slice(cols.min(), cols.max() + 1)


def import_sweets_dataset(
    sweets_dir: str | Path,
    output_root: str | Path,
    dataset_name: str,
    *,
    phase_pattern: str = "*.int.tif",
    coherence_pattern: str = "*.int.cor.tif",
    crop_kml: str | Path | None = None,
    margin: float = 0.20,
    overwrite: bool = False,
    skip_unpaired: bool = False,
) -> Path:
    """Create a ``Dataset_*`` directory consumable by existing ZoomInSAR stages.

    Scalar phase GeoTIFFs are converted to unit-magnitude complex interferograms;
    complex GeoTIFFs retain their wrapped phase.  Cropping is optional: without a
    KML the full SWEETS grid is imported unchanged.
    """
    sweets_dir, output_root = Path(sweets_dir), Path(output_root)
    if not sweets_dir.is_dir():
        raise FileNotFoundError(f"SWEETS directory not found: {sweets_dir}")
    name = dataset_name if dataset_name.startswith("Dataset_") else f"Dataset_{dataset_name}"
    output = output_root / name
    if output.exists() and any(output.iterdir()):
        if not overwrite:
            raise FileExistsError(f"Output exists: {output}; use --overwrite to replace it")
        shutil.rmtree(output)
    phase_files = _discover(sweets_dir, phase_pattern, "phase")
    coherence_files = _discover(sweets_dir, coherence_pattern, "coherence")
    missing = sorted(set(phase_files) ^ set(coherence_files))
    if missing and not skip_unpaired:
        raise FileNotFoundError(
            "SWEETS phase/coherence pairs do not match: "
            f"{missing}. Use --skip-unpaired only when intentionally excluding "
            "interferograms without coherence."
        )
    if missing:
        warnings.warn(
            f"Excluding {len(missing)} SWEETS interferogram(s) without matching "
            "phase/coherence products; details are recorded in crop_metadata.json.",
            stacklevel=2,
        )
        common_pairs = set(phase_files) & set(coherence_files)
        phase_files = {pair: path for pair, path in phase_files.items() if pair in common_pairs}
        coherence_files = {pair: path for pair, path in coherence_files.items() if pair in common_pairs}
    if not phase_files:
        raise FileNotFoundError("No complete SWEETS phase/coherence pairs were found")

    first_pair = sorted(phase_files)[0]
    first_phase, profile = _read_phase(phase_files[first_pair])
    longitude, latitude = _lon_lat(profile)
    row_slice, col_slice = (slice(None), slice(None))
    crop_path = Path(crop_kml) if crop_kml else None
    if crop_path:
        row_slice, col_slice = _crop_mask(longitude, latitude, crop_path, margin)
    output_shape = first_phase[row_slice, col_slice].shape
    (output / "INT").mkdir(parents=True)
    (output / "GEO").mkdir(parents=True)

    for pair in sorted(phase_files):
        phase, pair_profile = _read_phase(phase_files[pair])
        if pair_profile["shape"] != profile["shape"] or pair_profile["crs"] != profile["crs"] or pair_profile["transform"] != profile["transform"]:
            raise ValueError(f"{phase_files[pair]} is not on the common SWEETS grid")
        coherence = _read_coherence(coherence_files[pair], profile)
        label = f"{pair[0]}-{pair[1]}"
        np.asarray(phase[row_slice, col_slice], dtype=">c8").tofile(output / "INT" / f"{label}.tflt.filt")
        np.asarray(coherence[row_slice, col_slice], dtype=">f4").tofile(output / "INT" / f"{label}.tflt.coh")

    np.asarray(longitude[row_slice, col_slice], dtype=">f4").tofile(output / "GEO" / "sweets.lon")
    np.asarray(latitude[row_slice, col_slice], dtype=">f4").tofile(output / "GEO" / "sweets.lat")
    if crop_path:
        shutil.copy2(crop_path, output / crop_path.name)
    metadata = {
        "source": "SWEETS GeoTIFF",
        "sweets_directory": str(sweets_dir.resolve()),
        "phase_pattern": phase_pattern,
        "coherence_pattern": coherence_pattern,
        "skipped_unpaired_pairs": [f"{first}-{second}" for first, second in missing],
        "source_shape_lines_width": list(profile["shape"]),
        "output_shape_lines_width": list(output_shape),
        "crop_kml": str(crop_path.resolve()) if crop_path else None,
        "margin_fraction_per_side": margin if crop_path else None,
        "formats": {"filt": "big-endian complex64", "coh": "big-endian float32", "lon_lat": "big-endian float32"},
    }
    (output / "crop_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return output
