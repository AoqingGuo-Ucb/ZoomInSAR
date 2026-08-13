"""Minimum-cost branch-cut phase unwrapping."""

from __future__ import annotations

import heapq
from dataclasses import dataclass

import numpy as np

from .residues import compute_residues, minimum_cost_branch_cuts, wrap


@dataclass
class UnwrapDiagnostics:
    positive_residues: int
    negative_residues: int
    branch_cut_pixels: int
    branch_cut_edges: int
    valid_pixels: int
    unwrapped_pixels: int
    components: int
    rejected_long_cuts: int
    maximum_cut_length: float


def _can_move(y: int, x: int, ny: int, nx: int, blocked_lr: np.ndarray, blocked_ud: np.ndarray) -> bool:
    if y == ny:
        return not blocked_lr[y, min(x, nx)]
    return not blocked_ud[min(y, ny), x]


def quality_guided_unwrap(
    wrapped: np.ndarray,
    coherence: np.ndarray,
    valid: np.ndarray,
    blocked_lr: np.ndarray,
    blocked_ud: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Grow from high-coherence pixels without crossing branch-cut edges."""
    height, width = wrapped.shape
    result = np.full(wrapped.shape, np.nan, dtype=float)
    visited = np.zeros(wrapped.shape, dtype=bool)
    labels = np.full(wrapped.shape, -1, dtype=np.int32)
    components = 0
    score = np.where(valid, np.nan_to_num(coherence, nan=0.0), -np.inf)
    while np.any(valid & ~visited):
        remaining = valid & ~visited
        seed_flat = int(np.argmax(np.where(remaining, score, -np.inf)))
        sy, sx = np.unravel_index(seed_flat, wrapped.shape)
        result[sy, sx] = wrapped[sy, sx]
        visited[sy, sx] = True
        labels[sy, sx] = components
        components += 1
        queue: list[tuple[float, int, int]] = [(-float(score[sy, sx]), sy, sx)]
        while queue:
            _, y, x = heapq.heappop(queue)
            for ny, nx in ((y - 1, x), (y + 1, x), (y, x - 1), (y, x + 1)):
                if not (0 <= ny < height and 0 <= nx < width):
                    continue
                if visited[ny, nx] or not valid[ny, nx] or not _can_move(y, x, ny, nx, blocked_lr, blocked_ud):
                    continue
                result[ny, nx] = result[y, x] + wrap(wrapped[ny, nx] - wrapped[y, x])
                visited[ny, nx] = True
                labels[ny, nx] = components - 1
                heapq.heappush(queue, (-float(score[ny, nx]), ny, nx))
    return result, labels, components


def align_component_cycles(
    wrapped: np.ndarray,
    unwrapped: np.ndarray,
    labels: np.ndarray,
    blocked_lr: np.ndarray,
    blocked_ud: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Align independently grown regions using phase differences across cut edges."""
    n_components = int(labels.max() + 1)
    if n_components <= 1:
        return unwrapped, 0
    observations: dict[tuple[int, int], list[float]] = {}

    def observe(y0: int, x0: int, y1: int, x1: int) -> None:
        a, b = int(labels[y0, x0]), int(labels[y1, x1])
        if a < 0 or b < 0 or a == b:
            return
        # Integer cycles to add to component b relative to component a.
        cycles = (unwrapped[y0, x0] + wrap(wrapped[y1, x1] - wrapped[y0, x0]) - unwrapped[y1, x1]) / (2 * np.pi)
        if a < b:
            observations.setdefault((a, b), []).append(float(cycles))
        else:
            observations.setdefault((b, a), []).append(float(-cycles))

    for y, x in np.column_stack(np.where(blocked_lr)):
        observe(int(y), int(x), int(y), int(x + 1))
    for y, x in np.column_stack(np.where(blocked_ud)):
        observe(int(y), int(x), int(y + 1), int(x))

    graph: dict[int, list[tuple[int, int, int]]] = {i: [] for i in range(n_components)}
    for (a, b), values in observations.items():
        cycle = int(np.rint(np.median(values)))
        weight = len(values)
        graph[a].append((b, cycle, weight))
        graph[b].append((a, -cycle, weight))
    offsets: dict[int, int] = {}
    # Start with the largest component so its reference remains unchanged.
    sizes = np.bincount(labels[labels >= 0], minlength=n_components)
    for root in np.argsort(sizes)[::-1]:
        root = int(root)
        if root in offsets:
            continue
        offsets[root] = 0
        queue = [root]
        while queue:
            current = queue.pop(0)
            for neighbor, relative, _ in sorted(graph[current], key=lambda item: -item[2]):
                if neighbor not in offsets:
                    offsets[neighbor] = offsets[current] + relative
                    queue.append(neighbor)
    aligned = unwrapped.copy()
    changed = 0
    for component, cycles in offsets.items():
        if cycles:
            mask = labels == component
            aligned[mask] += cycles * 2 * np.pi
            changed += int(np.count_nonzero(mask))
    return aligned, changed


def unwrap_mcf(
    wrapped: np.ndarray,
    coherence: np.ndarray | None = None,
    coherence_threshold: float = 0.05,
    max_branch_cut_length: float = 15.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, UnwrapDiagnostics]:
    """Unwrap one interferogram using residues, minimum-cost cuts, and quality growth."""
    wrapped = np.asarray(wrapped, dtype=float)
    if coherence is None:
        coherence = np.ones(wrapped.shape, dtype=float)
    coherence = np.asarray(coherence, dtype=float)
    valid = np.isfinite(wrapped) & np.isfinite(coherence) & (coherence > coherence_threshold)
    residues = compute_residues(wrapped, valid)
    cuts, blocked_lr, blocked_ud, edges, rejected = minimum_cost_branch_cuts(
        residues, coherence, max_length=max_branch_cut_length
    )
    # A geometric cut may cross an invalid/water cell while connecting two land
    # residues. Remove those portions so diagnostics and graph alignment never
    # treat masked water as part of the branch-cut network.
    cuts[~valid] = False
    blocked_lr &= valid[:, :-1] & valid[:, 1:]
    blocked_ud &= valid[:-1, :] & valid[1:, :]
    # Invalid pixels are excluded by valid; they are not converted into cut pixels.
    unwrapped, labels, components = quality_guided_unwrap(
        wrapped, coherence, valid, blocked_lr, blocked_ud
    )
    unwrapped, _ = align_component_cycles(
        wrapped, unwrapped, labels, blocked_lr, blocked_ud
    )
    offsets = np.full(wrapped.shape, 0, dtype=np.int32)
    finite = np.isfinite(unwrapped)
    offsets[finite] = np.rint((unwrapped[finite] - wrapped[finite]) / (2.0 * np.pi)).astype(np.int32)
    diagnostics = UnwrapDiagnostics(
        positive_residues=int(np.count_nonzero(residues > 0)),
        negative_residues=int(np.count_nonzero(residues < 0)),
        branch_cut_pixels=int(np.count_nonzero(cuts & valid)),
        branch_cut_edges=len(edges),
        valid_pixels=int(np.count_nonzero(valid)),
        unwrapped_pixels=int(np.count_nonzero(finite)),
        components=components,
        rejected_long_cuts=rejected,
        maximum_cut_length=float(max_branch_cut_length),
    )
    return unwrapped, offsets, cuts, diagnostics
