"""Create the compact PBLC example shipped with the pipeline repository."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import from_bounds


PAIRS = (
    "20220106-20220118",
    "20220118-20220130",
    "20220106-20220130",
)
SOURCE_SHAPE = (410, 462)
WINDOW = (221, 317, 106, 202)  # row start/end, column start/end


def subset_binary(source: Path, destination: Path, dtype: str) -> np.ndarray:
    r0, r1, c0, c1 = WINDOW
    array = np.memmap(source, dtype=np.dtype(dtype).newbyteorder(">"), mode="r", shape=SOURCE_SHAPE)
    subset = np.asarray(array[r0:r1, c0:c1])
    destination.parent.mkdir(parents=True, exist_ok=True)
    subset.astype(np.dtype(dtype).newbyteorder(">"), copy=False).tofile(destination)
    return subset


def crop_dem(source: Path, destination: Path, bounds: tuple[float, float, float, float]) -> None:
    west, east, south, north = bounds
    with rasterio.open(source) as src:
        window = from_bounds(west, south, east, north, src.transform).round_offsets().round_lengths()
        data = src.read(window=window, boundless=True, fill_value=src.nodata)
        profile = src.profile.copy()
        profile.update(width=data.shape[2], height=data.shape[1], transform=src.window_transform(window))
        destination.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(destination, "w", **profile) as dst:
            dst.write(data)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Full Dataset_PBLC directory")
    parser.add_argument("destination", type=Path, help="Pipeline Data directory")
    parser.add_argument("--dem", type=Path, required=True, help="Source geographic DEM")
    args = parser.parse_args()

    dataset = args.destination / "ROI" / "Dataset_PBLC_Demo"
    if dataset.exists():
        shutil.rmtree(dataset)
    geo = dataset / "GEO"
    ints = dataset / "INT"

    lon = subset_binary(args.source / "GEO" / "sla.lon", geo / "pblc_demo.lon", "f4")
    lat = subset_binary(args.source / "GEO" / "sla.lat", geo / "pblc_demo.lat", "f4")
    for pair in PAIRS:
        subset_binary(args.source / "INT" / f"{pair}.tflt.coh", ints / f"{pair}.tflt.coh", "f4")
        subset_binary(args.source / "INT" / f"{pair}.tflt.filt", ints / f"{pair}.tflt.filt", "c8")

    par = (args.source / "GEO" / "sla.mli.par").read_text(encoding="utf-8", errors="replace")
    par = re.sub(r"(?m)^range_samples:\s+\d+", "range_samples:                    96", par)
    par = re.sub(r"(?m)^azimuth_lines:\s+\d+", "azimuth_lines:                    96", par)
    (geo / "pblc_demo.mli.par").write_text(par, encoding="utf-8")

    valid = np.isfinite(lon) & np.isfinite(lat)
    bounds = (float(lon[valid].min()), float(lon[valid].max()), float(lat[valid].min()), float(lat[valid].max()))
    metadata = {
        "demo": "PBLC compact three-interferogram example",
        "source_shape_lines_width": list(SOURCE_SHAPE),
        "crop_window_zero_based_end_exclusive": dict(zip(("row_start", "row_end", "col_start", "col_end"), WINDOW)),
        "output_shape_lines_width": [96, 96],
        "geographic_bounds_west_east_south_north": list(bounds),
        "dtype_byte_order": "big-endian",
        "formats": {"filt": "complex64", "coh": "float32", "lon_lat": "float32"},
        "pairs": list(PAIRS),
    }
    (dataset / "crop_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    west, east, south, north = bounds
    kml = f'''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>PBLC Demo</name>
<Placemark><name>PBLC Demo extent</name><Polygon><outerBoundaryIs><LinearRing><coordinates>
{west},{south},0 {east},{south},0 {east},{north},0 {west},{north},0 {west},{south},0
</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark></Document></kml>'''
    (dataset / "PBLC_Demo.kml").write_text(kml, encoding="utf-8")

    readme = """# PBLC compact demo

This is a 96 x 96 pixel real-data excerpt with three interferograms connecting
2022-01-06, 2022-01-18, and 2022-01-30 as one closed triangle. Files use GAMMA
big-endian conventions and include complex filtered phase, float32 coherence,
longitude, latitude, metadata, and a matching DEM.

Run from the pipeline directory:

```powershell
.\\run_insar_pipeline.cmd --dataset PBLC_Demo --skip-cropper
```

The example starts after cropping to remain compact. Three interferograms can
verify the software workflow but cannot provide a scientifically meaningful
long-term velocity estimate. Before public redistribution, confirm that the
underlying PBLC data license permits it; the software MIT License does not
automatically relicense research data.
"""
    (dataset / "README.md").write_text(readme, encoding="utf-8")

    crop_dem(args.dem, args.destination / "DEM" / "rasters_USGS10m" / "output_USGS10m.tif", bounds)


if __name__ == "__main__":
    main()
