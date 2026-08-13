"""Temporal-baseline filtering with mandatory network connectivity."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

import numpy as np


@dataclass(frozen=True)
class NetworkEdge:
    index: int
    first: datetime
    second: datetime
    baseline_days: int
    selected: bool
    reason: str
    quality_tier: int = 0
    inversion_weight: float = 1.0
    qc_status: str = "quality_metadata_unavailable"
    detrending_accepted: bool | None = None


class _DisjointSet:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, first, second) -> bool:
        a, b = self.find(first), self.find(second)
        if a == b:
            return False
        self.parent[b] = a
        return True


def design_matrix(
    dates: list[datetime], pairs: list[tuple[datetime, datetime]]
) -> np.ndarray:
    """Build x(t2)-x(t1)=observation with the first date fixed to zero."""
    index = {date: i for i, date in enumerate(dates)}
    matrix = np.zeros((len(pairs), len(dates) - 1), dtype=float)
    for row, (first, second) in enumerate(pairs):
        i, j = index[first], index[second]
        if i:
            matrix[row, i - 1] = -1.0
        if j:
            matrix[row, j - 1] = 1.0
    return matrix


def select_connected_network(
    pairs: list[tuple[datetime, datetime]], max_baseline_days: int,
    quality: list[object] | None = None,
) -> tuple[list[int], list[NetworkEdge], list[datetime]]:
    """Select a connected network, preferring quality over temporal length."""
    if max_baseline_days <= 0:
        raise ValueError("max_baseline_days must be positive")
    dates = sorted({date for pair in pairs for date in pair})
    dsu = _DisjointSet(dates)
    selected: set[int] = set()
    reasons: dict[int, str] = {}
    baselines = [(second - first).days for first, second in pairs]
    if quality is None:
        quality = [None] * len(pairs)
    if len(quality) != len(pairs):
        raise ValueError("quality must contain one item per interferogram")
    tiers = [int(getattr(item, "tier", 0)) if item is not None else 0 for item in quality]
    # Keep normal-quality short edges, including useful redundant cycles. Bad
    # short edges are deferred so a longer but reliable bridge can replace them.
    for index, ((first, second), baseline) in enumerate(zip(pairs, baselines)):
        if baseline <= max_baseline_days and tiers[index] <= 1:
            selected.add(index); reasons[index] = "within_threshold_quality_ok"; dsu.union(first, second)
    # Quality-aware Kruskal: choose a longer reliable bridge before a shorter
    # network-protected/detrending-failed edge.
    for index in sorted(range(len(pairs)), key=lambda i: (tiers[i], baselines[i], pairs[i][0], pairs[i][1])):
        if index in selected:
            continue
        first, second = pairs[index]
        if dsu.union(first, second):
            selected.add(index)
            reasons[index] = "low_quality_connectivity_bridge" if tiers[index] >= 2 else "quality_connectivity_bridge"
    if len({dsu.find(date) for date in dates}) != 1:
        raise ValueError("Available interferograms cannot form a connected acquisition network")
    chosen = sorted(selected)
    chosen_pairs = [pairs[index] for index in chosen]
    rank = int(np.linalg.matrix_rank(design_matrix(dates, chosen_pairs)))
    if rank != len(dates) - 1:
        raise ValueError(f"Selected network rank is {rank}; expected {len(dates) - 1}")
    details = [
        NetworkEdge(
            index, first, second, baselines[index], index in selected,
            reasons.get(index, "not_required"), tiers[index],
            float(getattr(quality[index], "weight", 1.0)) if quality[index] is not None else 1.0,
            str(getattr(quality[index], "qc_status", "quality_metadata_unavailable")) if quality[index] is not None else "quality_metadata_unavailable",
            getattr(quality[index], "detrending_accepted", None) if quality[index] is not None else None,
        )
        for index, (first, second) in enumerate(pairs)
    ]
    return chosen, details, dates
