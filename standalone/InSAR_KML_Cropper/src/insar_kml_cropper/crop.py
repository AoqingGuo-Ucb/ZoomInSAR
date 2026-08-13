"""End-to-end dataset cropping."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
from tqdm.auto import tqdm

from .io import discover_shape, raster_memmap, write_window
from .kml import expand_bounds, kml_bounds
from .plots import plot_kml_crop_overlay


class NoOverlapError(ValueError):
    """Raised when an ROI does not intersect valid geocoded pixels."""


def _find_lon_lat(geo_dir: Path) -> tuple[Path, Path]:
    lon_files, lat_files = sorted(geo_dir.glob("*.lon")), sorted(geo_dir.glob("*.lat"))
    if len(lon_files) != 1 or len(lat_files) != 1:
        raise ValueError("GEO must contain exactly one .lon file and one .lat file")
    return lon_files[0], lat_files[0]


def _window_from_coordinates(
    lon: np.ndarray, lat: np.ndarray, bounds: tuple[float, float, float, float]
) -> tuple[int, int, int, int]:
    west, east, south, north = bounds
    valid = np.isfinite(lon) & np.isfinite(lat) & (lon != 0) & (lat != 0)
    inside = valid & (lon >= west) & (lon <= east) & (lat >= south) & (lat <= north)
    rows, cols = np.nonzero(inside)
    if not rows.size:
        valid_rows, valid_cols = np.nonzero(valid)
        if valid_rows.size:
            lon_min, lon_max = np.nanmin(lon[valid]), np.nanmax(lon[valid])
            lat_min, lat_max = np.nanmin(lat[valid]), np.nanmax(lat[valid])
            raise NoOverlapError(
                "Expanded KML bounds do not overlap valid longitude/latitude pixels. "
                f"KML bounds west/east/south/north={bounds}; "
                f"GEO coverage west/east/south/north=({lon_min}, {lon_max}, {lat_min}, {lat_max})"
            )
        raise NoOverlapError("No valid longitude/latitude pixels found in GEO rasters")
    return int(rows.min()), int(rows.max() + 1), int(cols.min()), int(cols.max() + 1)


def crop_one_roi(
    data_dir: str | Path,
    kml_file: str | Path,
    margin: float = 0.20,
    shape: tuple[int, int] | None = None,
    overwrite: bool = False,
    pair: str | None = None,
) -> Path:
    data_dir, kml_file = Path(data_dir), Path(kml_file)
    int_dir, geo_dir = data_dir / "INT", data_dir / "GEO"
    rows, cols = discover_shape(geo_dir, shape)
    lon_file, lat_file = _find_lon_lat(geo_dir)
    lon = raster_memmap(lon_file, (rows, cols), ">f4")
    lat = raster_memmap(lat_file, (rows, cols), ">f4")
    original_bounds = kml_bounds(kml_file)
    expanded = expand_bounds(original_bounds, margin)
    window = _window_from_coordinates(lon, lat, expanded)

    output = data_dir / "ROI" / f"Dataset_{kml_file.stem}"
    if output.exists() and any(output.iterdir()) and not overwrite:
        raise FileExistsError(f"Output exists: {output}; use --overwrite to replace files")
    (output / "INT").mkdir(parents=True, exist_ok=True)
    (output / "GEO").mkdir(parents=True, exist_ok=True)

    phase_files = sorted(int_dir.glob("*.filt"))
    coherence_files = sorted(int_dir.glob("*.coh"))
    if pair:
        phase_files = [path for path in phase_files if pair in path.name]
        coherence_files = [path for path in coherence_files if pair in path.name]
        if not phase_files or not coherence_files:
            raise FileNotFoundError(f"No .filt/.coh inputs found for pair {pair}")
    for source in tqdm(
        phase_files, total=len(phase_files), unit="interferogram",
        desc=f"Cropping {kml_file.stem}", dynamic_ncols=True,
    ):
        write_window(source, output / "INT" / source.name, (rows, cols), ">c8", window)
    for source in coherence_files:
        write_window(source, output / "INT" / source.name, (rows, cols), ">f4", window)
    for source in (lon_file, lat_file):
        write_window(source, output / "GEO" / source.name, (rows, cols), ">f4", window)
    shutil.copy2(kml_file, output / kml_file.name)

    r0, r1, c0, c1 = window
    cropped_lon = np.asarray(lon[r0:r1, c0:c1])
    cropped_lat = np.asarray(lat[r0:r1, c0:c1])
    first_coherence = next(iter(sorted((output / "INT").glob("*.coh"))), None)
    cropped_coherence = (
        np.asarray(np.memmap(first_coherence, dtype=">f4", mode="r", shape=cropped_lon.shape))
        if first_coherence is not None else None
    )
    plot_kml_crop_overlay(
        output / "kml_crop_overlay.png", kml_file, cropped_lon, cropped_lat,
        cropped_coherence, expanded,
    )

    metadata = {
        "kml": kml_file.name,
        "margin_fraction_per_side": margin,
        "kml_bounds_west_east_south_north": original_bounds,
        "expanded_bounds_west_east_south_north": expanded,
        "source_shape_lines_width": [rows, cols],
        "crop_window_zero_based_end_exclusive": {"row_start": r0, "row_end": r1, "col_start": c0, "col_end": c1},
        "output_shape_lines_width": [r1 - r0, c1 - c0],
        "dtype_byte_order": "big-endian",
        "formats": {"filt": "complex64", "coh": "float32", "lon_lat": "float32"},
    }
    (output / "crop_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return output


def crop_dataset(
    data_dir: str | Path,
    margin: float = 0.20,
    shape: tuple[int, int] | None = None,
    overwrite: bool = False,
    pair: str | None = None,
) -> list[Path]:
    data_dir = Path(data_dir)
    kml_files = sorted((data_dir / "ROI").glob("*.kml"))
    if not kml_files:
        raise FileNotFoundError(f"No KML files found in {data_dir / 'ROI'}")
    outputs: list[Path] = []
    skipped: list[dict[str, str]] = []
    for path in tqdm(
        kml_files, total=len(kml_files), unit="ROI", desc="Cropping ROIs",
        dynamic_ncols=True,
    ):
        try:
            outputs.append(crop_one_roi(data_dir, path, margin, shape, overwrite, pair))
        except NoOverlapError as exc:
            skipped.append({"kml": path.name, "reason": str(exc)})
            tqdm.write(f"Skipped ROI without GEO overlap: {path.name}")
    if skipped:
        skipped_path = data_dir / "ROI" / "skipped_no_overlap.json"
        skipped_path.write_text(
            json.dumps(skipped, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        tqdm.write(f"Skipped {len(skipped)} ROI(s); details written to {skipped_path}")
    if not outputs:
        raise NoOverlapError("No ROI KML overlaps the valid longitude/latitude coverage")
    return outputs
