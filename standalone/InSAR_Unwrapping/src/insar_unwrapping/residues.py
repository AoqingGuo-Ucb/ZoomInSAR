"""Phase residues and minimum-cost branch cuts."""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

TWOPI = 2.0 * np.pi
MAX_DENSE_PAIR_COSTS = 10_000_000


def wrap(phase: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * phase))


def compute_residues(wrapped: np.ndarray, valid: np.ndarray | None = None) -> np.ndarray:
    """Calculate integer residues around each 2x2 plaquette."""
    phase = np.asarray(wrapped, dtype=float)
    if valid is None:
        valid = np.isfinite(phase)
    cell_valid = valid[:-1, :-1] & valid[:-1, 1:] & valid[1:, 1:] & valid[1:, :-1]
    circulation = (
        wrap(phase[:-1, 1:] - phase[:-1, :-1])
        + wrap(phase[1:, 1:] - phase[:-1, 1:])
        + wrap(phase[1:, :-1] - phase[1:, 1:])
        + wrap(phase[:-1, :-1] - phase[1:, :-1])
    )
    result = np.zeros(phase.shape, dtype=np.int8)
    result[:-1, :-1][cell_valid] = np.rint(circulation[cell_valid] / TWOPI).astype(np.int8)
    return result


def _line(y0: int, x0: int, y1: int, x1: int) -> tuple[np.ndarray, np.ndarray]:
    n = max(abs(y1 - y0), abs(x1 - x0)) + 1
    return np.rint(np.linspace(y0, y1, n)).astype(int), np.rint(np.linspace(x0, x1, n)).astype(int)


def _coherence_pair_factor(
    coherence: np.ndarray | None,
    positive: np.ndarray,
    negative: np.ndarray,
    p_index: int,
    n_index: int,
) -> float:
    if coherence is None:
        return 1.0
    quality = (
        coherence[positive[p_index, 0], positive[p_index, 1]]
        + coherence[negative[n_index, 0], negative[n_index, 1]]
    ) / 2.0
    return float(0.75 + 0.5 * np.nan_to_num(quality, nan=0.0))


def _match_residues(
    positive: np.ndarray,
    negative: np.ndarray,
    coherence: np.ndarray | None,
    max_length: float,
) -> tuple[list[tuple[int, int, float]], int]:
    """Match opposite residues without constructing an unbounded distance matrix."""
    pair_count = int(len(positive) * len(negative))
    if pair_count <= MAX_DENSE_PAIR_COSTS:
        dy = positive[:, 0][:, None] - negative[:, 0][None, :]
        dx = positive[:, 1][:, None] - negative[:, 1][None, :]
        distance = np.hypot(dy, dx)
        cost = distance.copy()
        if coherence is not None:
            endpoint_quality = (
                coherence[positive[:, 0], positive[:, 1]][:, None]
                + coherence[negative[:, 0], negative[:, 1]][None, :]
            ) / 2.0
            cost *= 0.75 + 0.5 * np.nan_to_num(endpoint_quality, nan=0.0)
        pi, ni = linear_sum_assignment(cost)
        matches = []
        rejected = 0
        for p_index, n_index in zip(pi, ni):
            length = float(distance[p_index, n_index])
            if length > max_length:
                rejected += 1
                continue
            matches.append((int(p_index), int(n_index), length))
        return matches, rejected

    tree = cKDTree(negative.astype(float))
    candidate_edges: list[tuple[float, int, int, float]] = []
    for p_index, point in enumerate(positive.astype(float)):
        nearby = tree.query_ball_point(point, r=float(max_length))
        if not nearby:
            continue
        nearby_array = np.asarray(nearby, dtype=int)
        delta = negative[nearby_array].astype(float) - point
        lengths = np.hypot(delta[:, 0], delta[:, 1])
        for n_index, length in zip(nearby_array, lengths):
            factor = _coherence_pair_factor(coherence, positive, negative, p_index, int(n_index))
            candidate_edges.append((float(length) * factor, int(p_index), int(n_index), float(length)))
    candidate_edges.sort(key=lambda item: item[0])
    used_pos: set[int] = set()
    used_neg: set[int] = set()
    matches: list[tuple[int, int, float]] = []
    for _, p_index, n_index, length in candidate_edges:
        if p_index in used_pos or n_index in used_neg:
            continue
        used_pos.add(p_index)
        used_neg.add(n_index)
        matches.append((p_index, n_index, length))
    return matches, 0


def minimum_cost_branch_cuts(
    residues: np.ndarray,
    coherence: np.ndarray | None = None,
    max_length: float = 15.0,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    list[tuple[tuple[int, int], tuple[int, int]]],
    int,
]:
    """Pair opposite residues by minimum total cost and rasterize branch cuts.

    Euclidean distance is the base cost. Low mean coherence slightly lowers a
    path's cost, encouraging cuts to pass through less reliable pixels.
    Unbalanced residues are connected to their nearest image boundary.
    """
    positive = np.column_stack(np.where(residues > 0))
    negative = np.column_stack(np.where(residues < 0))
    cuts = np.zeros(residues.shape, dtype=bool)
    # blocked_lr[y, x] blocks movement between pixels (y,x) and (y,x+1).
    # blocked_ud[y, x] blocks movement between pixels (y,x) and (y+1,x).
    blocked_lr = np.zeros((residues.shape[0], residues.shape[1] - 1), dtype=bool)
    blocked_ud = np.zeros((residues.shape[0] - 1, residues.shape[1]), dtype=bool)
    edges: list[tuple[tuple[int, int], tuple[int, int]]] = []
    used_pos: set[int] = set()
    used_neg: set[int] = set()
    rejected = 0

    def add_dual_path(start: np.ndarray | tuple[int, int], end: np.ndarray | tuple[int, int]) -> None:
        """Rasterize a cut on the dual grid and mark only crossed primal edges."""
        nonlocal cuts, blocked_lr, blocked_ud
        yy, xx = _line(int(start[0]), int(start[1]), int(end[0]), int(end[1]))
        cuts[yy, xx] = True
        for y0, x0, y1, x1 in zip(yy[:-1], xx[:-1], yy[1:], xx[1:]):
            # Split a diagonal dual-grid step deterministically into two crossings.
            if x1 != x0:
                row = int(np.clip(min(y0, residues.shape[0] - 2), 0, blocked_ud.shape[0] - 1))
                col = int(np.clip(max(x0, x1), 0, blocked_ud.shape[1] - 1))
                blocked_ud[row, col] = True
            if y1 != y0:
                row = int(np.clip(max(y0, y1), 0, blocked_lr.shape[0] - 1))
                col = int(np.clip(min(x0, residues.shape[1] - 2), 0, blocked_lr.shape[1] - 1))
                blocked_lr[row, col] = True

    if len(positive) and len(negative):
        matches, rejected_pairs = _match_residues(positive, negative, coherence, max_length)
        rejected += rejected_pairs
        for p_index, n_index, length in matches:
            p, n = positive[p_index], negative[n_index]
            used_pos.add(int(p_index)); used_neg.add(int(n_index))
            add_dual_path(p, n)
            edges.append(((int(p[0]), int(p[1])), (int(n[0]), int(n[1]))))
    height, width = residues.shape
    for points, used in ((positive, used_pos), (negative, used_neg)):
        for index, point in enumerate(points):
            if index in used:
                continue
            y, x = map(int, point)
            boundaries = [(y, 0), (y, width - 1), (0, x), (height - 1, x)]
            target = min(boundaries, key=lambda q: (q[0] - y) ** 2 + (q[1] - x) ** 2)
            if float(np.hypot(target[0] - y, target[1] - x)) > max_length:
                rejected += 1
                continue
            add_dual_path((y, x), target)
            edges.append(((y, x), target))
    return cuts, blocked_lr, blocked_ud, edges, rejected
