import numpy as np

from insar_pixel_offset.tracking import grid_centers, normalized_offset


def test_recovers_known_integer_offset():
    rng = np.random.default_rng(4)
    reference = rng.normal(size=(128, 128)).astype(np.complex64)
    secondary = np.zeros_like(reference)
    # The matching feature moves 3 rows down and 2 columns right in secondary.
    secondary[3:, 2:] = reference[:-3, :-2]
    dy, dx, correlation = normalized_offset(reference, secondary, 64, 64, 32, 6)
    assert (dy, dx) == (3.0, 2.0)
    assert correlation > 0.99


def test_grid_keeps_windows_inside_image():
    rows, cols = grid_centers((200, 220), window=64, search=12, stride=32)
    assert rows[0] >= 44 and cols[0] >= 44
    assert rows[-1] < 156 and cols[-1] < 176
