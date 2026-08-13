import numpy as np

from insar_unwrapping.corrections import TWOPI, repair_edge_connected_cycles, repair_region_graph_cycles


def test_repairs_closed_integer_cycle_island():
    phase = np.zeros((80, 90), dtype=float)
    phase[20:60, 25:70] = 2 * np.pi
    coherence = np.full_like(phase, 0.8)
    corrected, cycles, stats = repair_region_graph_cycles(
        phase, coherence, min_pixels=30
    )
    assert stats.regions_repaired >= 1
    assert np.count_nonzero(np.abs(corrected) > 1e-8) <= 0.01 * (40 * 45)
    assert np.count_nonzero(cycles) >= 0.99 * (40 * 45)


def test_repairs_two_regions_with_different_cycle_offsets():
    phase = np.zeros((100, 120), dtype=float)
    phase[10:42, 12:48] = 2 * np.pi
    phase[55:90, 68:110] = -4 * np.pi
    corrected, cycles, stats = repair_region_graph_cycles(phase, min_pixels=30)
    assert stats.regions_repaired >= 2
    changed_pixels = 32 * 36 + 35 * 42
    assert np.count_nonzero(np.abs(corrected) > 1e-8) <= 0.01 * changed_pixels
    assert np.min(cycles) <= -1
    assert np.max(cycles) >= 2


def test_does_not_remove_noninteger_real_step():
    phase = np.zeros((80, 90), dtype=float)
    phase[20:60, 25:70] = 0.65 * 2 * np.pi
    corrected, cycles, stats = repair_region_graph_cycles(phase, min_pixels=30)
    assert stats.regions_repaired == 0
    assert np.array_equal(corrected, phase)
    assert not np.any(cycles)


def test_local_pairing_cancels_a_long_wavelength_ramp():
    yy, xx = np.mgrid[:100, :120]
    ramp = 0.035 * xx - 0.018 * yy
    phase = ramp.copy()
    phase[25:75, 35:90] += 2 * np.pi
    corrected, cycles, stats = repair_region_graph_cycles(phase, min_pixels=30)
    assert stats.regions_repaired >= 1
    residual = corrected - ramp
    assert np.count_nonzero(np.abs(residual) > 1e-6) <= 0.01 * (50 * 55)
    assert np.min(cycles) <= -1


def test_repairs_integer_region_closed_by_image_boundary():
    phase = np.zeros((90, 110), dtype=float)
    # The top image edge is one side of the region boundary; only the two side
    # walls and lower edge are present inside the valid raster.
    phase[:48, 25:82] = 2 * np.pi
    corrected, cycles, stats = repair_region_graph_cycles(phase, min_pixels=30)
    assert stats.regions_repaired >= 1
    assert np.count_nonzero(np.abs(corrected) > 1e-8) <= 0.01 * (48 * 57)
    assert np.min(cycles) <= -1


def test_edge_strip_repairs_top_plateau_over_quadratic_trend():
    yy,xx=np.mgrid[:100,:120]
    trend=0.012*xx+0.00015*xx*xx-0.008*yy
    phase=trend.copy();phase[:38,72:112]+=2*np.pi
    corrected,cycles,stats=repair_edge_connected_cycles(phase,min_pixels=30)
    assert stats.components_repaired>=1
    assert np.count_nonzero(np.abs(corrected-trend)>1e-6)<=0.06*(38*40)
    assert np.min(cycles)<=-1


def test_edge_strip_rejects_noninteger_boundary_signal():
    phase=np.zeros((90,110));phase[:35,60:100]+=0.60*TWOPI
    corrected,cycles,stats=repair_edge_connected_cycles(phase,min_pixels=30)
    assert stats.components_repaired==0
    assert np.array_equal(corrected,phase)
    assert not np.any(cycles)
