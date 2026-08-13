import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from insar_unwrapping.water import water_mask_from_dem


def test_dem_water_mask_at_lon_lat_pixels(tmp_path: Path):
    dataset = tmp_path / "Dataset_Test"
    (dataset / "GEO").mkdir(parents=True)
    longitude = np.array([[0.5, 1.5], [0.5, 1.5]], dtype=">f4")
    latitude = np.array([[1.5, 1.5], [0.5, 0.5]], dtype=">f4")
    longitude.tofile(dataset / "GEO" / "x.lon")
    latitude.tofile(dataset / "GEO" / "x.lat")
    (dataset / "crop_metadata.json").write_text(
        json.dumps({"output_shape_lines_width": [2, 2]}), encoding="utf-8"
    )
    dem_path = tmp_path / "dem.tif"
    with rasterio.open(
        dem_path, "w", driver="GTiff", height=2, width=2, count=1,
        dtype="float32", crs="EPSG:4326", transform=from_origin(0, 2, 1, 1),
        nodata=-9999,
    ) as dem:
        dem.write(np.array([[10, -1], [5, 0]], dtype=np.float32), 1)
    water, valid, _, _, elevation = water_mask_from_dem(dataset, (2, 2), dem_path)
    np.testing.assert_array_equal(water, [[False, True], [False, True]])
    assert np.all(valid) and np.isfinite(elevation).all()
