from pathlib import Path

import numpy as np

from insar_filtering.io import read_gamma, write_gamma


def test_big_endian_complex_roundtrip(tmp_path: Path):
    values = (np.arange(12).reshape(3, 4) + 1j).astype(np.complex64)
    path = tmp_path / "test.filt"
    write_gamma(path, values, "c8")
    np.testing.assert_array_equal(read_gamma(path, values.shape, "c8"), values)
