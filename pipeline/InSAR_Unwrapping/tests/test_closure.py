import numpy as np

from insar_unwrapping.closure import repair_phase_closure


def test_repairs_integer_closure_region_on_low_quality_edge():
    shape = (30, 30)
    stack = np.zeros((3, *shape), dtype=float)
    stack[2, 5:25, 5:25] = 2 * np.pi
    coherence = np.ones_like(stack)
    coherence[2] = 0.2
    result, correction, stats = repair_phase_closure(
        stack, [("1", "2"), ("2", "3"), ("1", "3")], coherence,
        min_component_pixels=10,
    )
    assert stats.triangles == 1
    assert stats.components_repaired >= 1
    assert np.nanmax(np.abs(result[0] + result[1] - result[2])) < 1e-10
    assert np.count_nonzero(correction) >= 400


def test_selects_spatially_correct_edge_instead_of_lowest_quality_edge():
    shape = (30, 30)
    stack = np.zeros((3, *shape), dtype=float)
    # AB contains the actual spatial step, but BC is deliberately assigned the
    # lowest quality. Correcting BC would close the triangle while creating a
    # new 2-pi spatial discontinuity on an otherwise continuous interferogram.
    stack[0, 5:25, 5:25] = 2 * np.pi
    coherence = np.ones_like(stack)
    coherence[1] = 0.2
    result, correction, stats = repair_phase_closure(
        stack, [("1", "2"), ("2", "3"), ("1", "3")], coherence,
        min_component_pixels=10,
    )
    assert stats.components_repaired >= 1
    assert np.count_nonzero(correction[0]) == 400
    assert not np.any(correction[1:])
    assert np.nanmax(np.abs(result[0] + result[1] - result[2])) < 1e-10


def test_rejects_region_without_a_measurable_outer_boundary():
    shape = (30, 30)
    stack = np.zeros((3, *shape), dtype=float)
    stack[0] = 2 * np.pi
    result, correction, stats = repair_phase_closure(
        stack, [("1", "2"), ("2", "3"), ("1", "3")],
        min_component_pixels=10,
    )
    assert stats.components_repaired == 0
    assert stats.spatially_rejected_components >= 1
    assert np.array_equal(result, stack)
    assert not np.any(correction)


def test_peels_nested_two_cycle_error_across_different_edges():
    shape = (50, 50)
    stack = np.zeros((3, *shape), dtype=float)
    # AB has a compact +1-cycle island. AC has a broad -1-cycle region.
    # Closure is therefore +1 in the broad region and +2 in its central island.
    stack[0, 15:35, 15:35] = 2 * np.pi
    stack[2, 5:45, 5:45] = -2 * np.pi
    result, correction, stats = repair_phase_closure(
        stack, [("1", "2"), ("2", "3"), ("1", "3")],
        max_iterations=3, min_component_pixels=10,
    )
    assert stats.components_repaired >= 2
    assert np.count_nonzero(correction[0]) == 400
    assert np.count_nonzero(correction[2]) == 1600
    assert np.nanmax(np.abs(result[0] + result[1] - result[2])) < 1e-10
