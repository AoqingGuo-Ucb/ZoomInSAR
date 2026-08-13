from datetime import datetime, timedelta

import numpy as np

from insar_timeseries.solver import invert_timeseries


def test_recovers_known_cumulative_phase():
    start = datetime(2020, 1, 1)
    dates = [start + timedelta(days=value) for value in (0, 12, 24)]
    pairs = [(dates[0], dates[1]), (dates[1], dates[2]), (dates[0], dates[2])]
    truth = np.array([0.0, 1.5, 4.0])
    observations = np.array([1.5, 2.5, 4.0])[:, None, None]
    displacement, _, residual, stats = invert_timeseries(
        observations, pairs, dates, wavelength_m=4 * np.pi
    )
    np.testing.assert_allclose(displacement[:, 0, 0], truth)
    assert residual[0, 0] < 1e-12 and stats.pixels_solved == 1
