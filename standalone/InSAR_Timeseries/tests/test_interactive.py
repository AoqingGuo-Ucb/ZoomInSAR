import numpy as np

from insar_timeseries.interactive import _nearest_valid_pixel


def test_nearest_valid_pixel_uses_clicked_pixel():
    velocity = np.full((4, 5), np.nan)
    velocity[2, 3] = 0.1
    assert _nearest_valid_pixel(velocity, 3.1, 1.9) == (2, 3)


def test_nearest_valid_pixel_moves_off_nan():
    velocity = np.full((5, 5), np.nan)
    velocity[1, 1] = 0.1
    velocity[4, 4] = 0.2
    assert _nearest_valid_pixel(velocity, 2.0, 2.0) == (1, 1)
