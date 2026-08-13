import numpy as np

from insar_unwrapping.residues import compute_residues, minimum_cost_branch_cuts, wrap


def test_wrap_range():
    result = wrap(np.array([-4 * np.pi, -1.5 * np.pi, 0, 1.5 * np.pi, 4 * np.pi]))
    assert np.all(result >= -np.pi) and np.all(result <= np.pi)


def test_smooth_plane_has_no_residues():
    y, x = np.mgrid[:20, :30]
    wrapped = wrap(0.2 * x + 0.1 * y)
    assert np.count_nonzero(compute_residues(wrapped)) == 0


def test_rejects_overlong_branch_cut():
    residues = np.zeros((30, 30), dtype=np.int8)
    residues[2, 2] = 1
    residues[25, 25] = -1
    cuts, blocked_lr, blocked_ud, edges, rejected = minimum_cost_branch_cuts(
        residues, max_length=10
    )
    assert rejected >= 1
    assert edges
    assert all(np.hypot(a[0] - b[0], a[1] - b[1]) <= 10 for a, b in edges)
    assert np.any(cuts) and (np.any(blocked_lr) or np.any(blocked_ud))
