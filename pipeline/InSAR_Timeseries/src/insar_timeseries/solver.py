"""Cached least-squares SBAS inversion for raster stacks."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

import numpy as np

from .network import design_matrix


@dataclass
class InversionStats:
    pixels_total: int
    pixels_solved: int
    validity_patterns: int
    dates: int
    interferograms: int


def invert_timeseries(
    phase_stack: np.ndarray,
    pairs: list[tuple[datetime, datetime]],
    dates: list[datetime],
    wavelength_m: float = 0.056,
    phase_sign: float = 1.0,
    edge_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, InversionStats]:
    """Invert phase differences to cumulative LOS displacement and residual RMS.

    Pixels sharing the same valid-edge pattern reuse one weighted pseudoinverse.
    A pixel is solved only when its valid subnetwork retains full temporal rank.
    """
    phase_stack = np.asarray(phase_stack, dtype=float)
    n_edges, height, width = phase_stack.shape
    matrix = design_matrix(dates, pairs)
    if edge_weights is None:
        edge_weights = np.ones(n_edges, dtype=float)
    edge_weights = np.asarray(edge_weights, dtype=float)
    if edge_weights.shape != (n_edges,) or np.any(edge_weights <= 0):
        raise ValueError("edge_weights must be one positive value per interferogram")
    flat = phase_stack.reshape(n_edges, -1)
    valid = np.isfinite(flat)
    cumulative_phase = np.full((len(dates), flat.shape[1]), np.nan, dtype=float)
    residual_rms = np.full(flat.shape[1], np.nan, dtype=float)
    # Pack validity masks into bytes so repeated coastal/water patterns are grouped.
    packed = np.packbits(valid.T, axis=1)
    unique, inverse = np.unique(packed, axis=0, return_inverse=True)
    solved = 0
    for pattern_index in range(len(unique)):
        pixels = np.where(inverse == pattern_index)[0]
        edge_mask = valid[:, pixels[0]]
        if np.count_nonzero(edge_mask) < len(dates) - 1:
            continue
        local_matrix = matrix[edge_mask]
        if np.linalg.matrix_rank(local_matrix) != len(dates) - 1:
            continue
        sqrt_weight = np.sqrt(edge_weights[edge_mask])
        weighted_matrix = local_matrix * sqrt_weight[:, None]
        observations = flat[edge_mask][:, pixels] * sqrt_weight[:, None]
        solution = np.linalg.pinv(weighted_matrix) @ observations
        cumulative_phase[0, pixels] = 0.0
        cumulative_phase[1:, pixels] = solution
        prediction = local_matrix @ solution
        residual = flat[edge_mask][:, pixels] - prediction
        residual_rms[pixels] = np.sqrt(np.mean(residual * residual, axis=0))
        solved += len(pixels)
    displacement = cumulative_phase * (phase_sign * wavelength_m / (4.0 * np.pi))
    elapsed_years = np.array([(date - dates[0]).days / 365.25 for date in dates])
    centered_time = elapsed_years - elapsed_years.mean()
    denom = float(np.sum(centered_time ** 2))
    mean_velocity = np.full(flat.shape[1], np.nan)
    finite_all = np.all(np.isfinite(displacement), axis=0)
    if denom > 0 and np.any(finite_all):
        centered_disp = displacement[:, finite_all] - displacement[:, finite_all].mean(axis=0)
        mean_velocity[finite_all] = centered_time @ centered_disp / denom
    return (
        displacement.reshape(len(dates), height, width),
        mean_velocity.reshape(height, width),
        residual_rms.reshape(height, width),
        InversionStats(height * width, solved, len(unique), len(dates), n_edges),
    )
