from datetime import datetime, timedelta

from insar_timeseries.network import design_matrix, select_connected_network
from insar_timeseries.quality import EdgeQuality


def test_adds_shortest_bridge_when_threshold_network_is_disconnected():
    start = datetime(2020, 1, 1)
    dates = [start + timedelta(days=value) for value in (0, 12, 40, 52)]
    pairs = [(dates[0], dates[1]), (dates[1], dates[2]), (dates[2], dates[3]), (dates[0], dates[2])]
    selected, details, output_dates = select_connected_network(pairs, 12)
    assert selected == [0, 1, 2]
    assert "connectivity_bridge" in details[1].reason
    assert design_matrix(output_dates, [pairs[i] for i in selected]).shape == (3, 3)


def test_threshold_24_keeps_redundant_cycle_edges():
    start = datetime(2020, 1, 1)
    dates = [start + timedelta(days=value) for value in (0, 12, 24)]
    pairs = [(dates[0], dates[1]), (dates[1], dates[2]), (dates[0], dates[2])]
    selected, _, _ = select_connected_network(pairs, 24)
    assert selected == [0, 1, 2]


def test_prefers_longer_good_bridge_over_short_bad_bridge():
    start = datetime(2020, 1, 1)
    dates = [start + timedelta(days=value) for value in (0, 12, 24)]
    pairs = [(dates[0], dates[1]), (dates[1], dates[2]), (dates[0], dates[2])]
    quality = [EdgeQuality(), EdgeQuality(3, 0.03, "network_protected", False), EdgeQuality()]
    selected, details, _ = select_connected_network(pairs, 12, quality)
    assert selected == [0, 2]
    assert not details[1].selected
    assert details[2].reason == "quality_connectivity_bridge"


def test_keeps_only_available_bad_bridge_with_low_weight():
    start = datetime(2020, 1, 1)
    dates = [start + timedelta(days=value) for value in (0, 12, 24)]
    pairs = [(dates[0], dates[1]), (dates[1], dates[2])]
    quality = [EdgeQuality(), EdgeQuality(3, 0.03, "network_protected", False)]
    selected, details, _ = select_connected_network(pairs, 12, quality)
    assert selected == [0, 1]
    assert details[1].reason == "low_quality_connectivity_bridge"
    assert details[1].inversion_weight == 0.03
