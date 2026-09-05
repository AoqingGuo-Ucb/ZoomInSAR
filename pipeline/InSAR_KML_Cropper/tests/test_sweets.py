from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from insar_kml_cropper.sweets import import_sweets_dataset


def _write(path: Path, values: np.ndarray) -> None:
    with rasterio.open(
        path, "w", driver="GTiff", height=values.shape[0], width=values.shape[1],
        count=1, dtype=values.dtype, crs="EPSG:4326", transform=from_origin(-118, 35, 0.01, 0.01),
    ) as destination:
        destination.write(values, 1)


def test_import_sweets_phase_and_coherence_without_crop(tmp_path: Path):
    source = tmp_path / "sweets"
    source.mkdir()
    phase = np.array([[0.0, np.pi / 2], [-np.pi / 2, np.pi]], dtype=np.float32)
    coherence = np.array([[0.2, 0.5], [0.8, 1.0]], dtype=np.float32)
    _write(source / "20200101_20200113_int.tif", phase)
    _write(source / "20200101_20200113_int.cor.tif", coherence)

    dataset = import_sweets_dataset(source, tmp_path / "ROI", "example")

    result_phase = np.fromfile(dataset / "INT" / "20200101-20200113.tflt.filt", dtype=">c8").reshape(2, 2)
    result_coherence = np.fromfile(dataset / "INT" / "20200101-20200113.tflt.coh", dtype=">f4").reshape(2, 2)
    longitude = np.fromfile(dataset / "GEO" / "sweets.lon", dtype=">f4").reshape(2, 2)
    np.testing.assert_allclose(result_phase, np.exp(1j * phase))
    np.testing.assert_allclose(result_coherence, coherence)
    np.testing.assert_allclose(longitude[0], [-117.995, -117.985])
