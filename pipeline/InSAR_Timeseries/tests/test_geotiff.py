import numpy as np
import rasterio

from insar_timeseries.geotiff import export_velocity_geotiff
from insar_timeseries.plots import select_representative_points


def test_exports_qgis_ready_velocity_and_point(tmp_path):
    yy, xx = np.mgrid[:20, :30]
    longitude = -119.0 + xx * 0.001 + yy * 0.00005
    latitude = 34.0 - yy * 0.001
    velocity = (xx - 15) * 0.001
    path = tmp_path / "velocity.tif"
    export_velocity_geotiff(path, velocity, longitude, latitude)
    with rasterio.open(path) as source:
        assert source.crs.to_epsg() == 4326
        assert source.nodata == -9999.0
        assert source.tags(1)["units"] == "mm/year"
        assert source.read(1).shape == velocity.shape
    points = select_representative_points(velocity, longitude, latitude)
    assert points[0]["id"] == "P1"
    assert np.isfinite(points[0]["longitude"])
