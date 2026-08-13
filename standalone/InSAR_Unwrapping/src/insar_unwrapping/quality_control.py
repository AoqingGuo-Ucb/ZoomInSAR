"""Automatic interferogram rejection without breaking the date network."""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .corrections import normalized_gaussian

TWOPI = 2.0 * np.pi


@dataclass
class InterferogramQC:
    pair: str
    spatial_residual_fraction: float
    robust_phase_range_cycles: float
    bad_candidate: bool
    included: bool
    status: str


def _connected(pairs: list[tuple[str, str]], keep: set[int]) -> bool:
    dates = {date for pair in pairs for date in pair}
    if not dates:
        return True
    adjacency = {date: set() for date in dates}
    for index in keep:
        a, b = pairs[index]
        adjacency[a].add(b); adjacency[b].add(a)
    seen = set()
    pending = [next(iter(dates))]
    while pending:
        date = pending.pop()
        if date in seen:
            continue
        seen.add(date); pending.extend(adjacency[date] - seen)
    return seen == dates


def select_network_safe_interferograms(
    stack: np.ndarray,
    pairs: list[tuple[str, str]],
    sigma: float = 12.0,
    residual_threshold_cycles: float = 0.45,
    max_residual_fraction: float = 0.04,
    max_phase_range_cycles: float = 2.0,
) -> tuple[list[int], list[dict]]:
    """Reject severe spatial outliers only when their removal preserves connectivity."""
    metrics = []
    for index, phase in enumerate(stack):
        valid = np.isfinite(phase)
        if not np.any(valid):
            fraction = 1.0
            phase_range = float("inf")
        else:
            smooth = normalized_gaussian(phase, valid, sigma)
            residual = np.abs((phase - smooth) / TWOPI)
            fraction = float(np.mean(residual[valid] >= residual_threshold_cycles))
            low, high = np.nanpercentile(phase[valid], [1, 99])
            phase_range = float((high - low) / TWOPI)
        bad = fraction >= max_residual_fraction and phase_range >= max_phase_range_cycles
        metrics.append((index, fraction, phase_range, bad))

    keep = set(range(len(pairs)))
    status = {index: "included" for index in keep}
    candidates = sorted(
        (item for item in metrics if item[3]),
        key=lambda item: item[1] * item[2], reverse=True,
    )
    for index, _, _, _ in candidates:
        trial = keep - {index}
        if _connected(pairs, trial):
            keep = trial
            status[index] = "excluded_spatial_outlier"
        else:
            status[index] = "network_protected"

    rows = [
        asdict(InterferogramQC(
            pair=f"{pairs[index][0]}-{pairs[index][1]}",
            spatial_residual_fraction=fraction,
            robust_phase_range_cycles=phase_range,
            bad_candidate=bad,
            included=index in keep,
            status=status[index],
        ))
        for index, fraction, phase_range, bad in metrics
    ]
    return sorted(keep), rows
