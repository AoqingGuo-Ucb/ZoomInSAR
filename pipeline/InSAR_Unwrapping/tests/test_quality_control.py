import numpy as np

from insar_unwrapping.quality_control import select_network_safe_interferograms


def test_excludes_bad_redundant_edge_but_keeps_network_connected():
    pairs = [("1", "2"), ("2", "3"), ("1", "3")]
    stack = np.zeros((3, 40, 40), dtype=float)
    stack[2, :20, :20] = 6 * np.pi
    keep, rows = select_network_safe_interferograms(
        stack, pairs, sigma=3, max_residual_fraction=0.01,
        max_phase_range_cycles=1.0,
    )
    assert keep == [0, 1]
    assert rows[2]["status"] == "excluded_spatial_outlier"


def test_protects_bad_bridge_edge():
    pairs = [("1", "2"), ("2", "3")]
    stack = np.zeros((2, 40, 40), dtype=float)
    stack[0, :20, :20] = 6 * np.pi
    keep, rows = select_network_safe_interferograms(
        stack, pairs, sigma=3, max_residual_fraction=0.01,
        max_phase_range_cycles=1.0,
    )
    assert keep == [0, 1]
    assert rows[0]["status"] == "network_protected"
