from pathlib import Path

import numpy as np

from insar_kml_cropper.io import freadbkB


def test_freadbk_complex_window(tmp_path: Path):
    expected = (np.arange(12).reshape(3, 4) + 1j * np.arange(20, 32).reshape(3, 4)).astype(">c8")
    path = tmp_path / "test.filt"
    expected.tofile(path)
    actual, count = freadbkB(path, lines=3, bkformat="cpxfloat32", r0=2, rN=3, c0=2, cN=4)
    np.testing.assert_array_equal(actual, expected[1:3, 1:4])
    assert count == 6
