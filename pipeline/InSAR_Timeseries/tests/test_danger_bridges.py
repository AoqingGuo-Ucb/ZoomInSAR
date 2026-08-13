from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import numpy as np

from insar_timeseries.plots import _plot_danger_bridges


def test_low_quality_bridge_is_drawn_red():
    start = datetime(2024, 1, 1)
    dates = [start + timedelta(days=12 * index) for index in range(3)]
    fig, axis = plt.subplots()
    _plot_danger_bridges(axis, dates, np.array([0.0, 1.0, 2.0]), [(dates[1], dates[2])])
    red = [line for line in axis.lines if line.get_color() == "red"]
    assert len(red) == 1
    assert red[0].get_label() == "Low-quality connectivity bridge (danger)"
    plt.close(fig)
