"""Edge-preserving complex-domain phase filters."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter


def _shift(values: np.ndarray, dy: int, dx: int) -> tuple[np.ndarray, np.ndarray]:
    shifted = np.roll(values, (dy, dx), axis=(0, 1))
    valid = np.ones(values.shape, dtype=bool)
    if dy > 0:
        valid[:dy] = False
    elif dy < 0:
        valid[dy:] = False
    if dx > 0:
        valid[:, :dx] = False
    elif dx < 0:
        valid[:, dx:] = False
    return shifted, valid


def nonlocal_phase_filter(
    interferogram: np.ndarray,
    coherence: np.ndarray,
    search_radius: int = 7,
    patch_radius: int = 3,
    h: float = 0.70,
    spatial_sigma: float = 3.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Filter wrapped phase by coherence-adaptive nonlocal complex averaging.

    Patch distance uses circular complex similarity, so phase values around
    -pi/+pi remain neighbors. Low-coherence pixels receive a broader bandwidth;
    high-coherence pixels retain stronger edge and fringe selectivity.
    """
    values = np.asarray(interferogram, dtype=np.complex128)
    coherence = np.clip(np.nan_to_num(coherence, nan=0.0), 0.0, 1.0)
    amplitude = np.abs(values)
    valid = np.isfinite(values.real) & np.isfinite(values.imag) & (amplitude > 0)
    unit = np.divide(values, amplitude, out=np.zeros_like(values), where=valid)
    numerator = np.zeros_like(unit)
    denominator = np.zeros(unit.shape, dtype=float)
    patch_size = 2 * patch_radius + 1
    bandwidth = np.maximum(0.08, h * (0.35 + 0.90 * (1.0 - coherence)))
    for dy in range(-search_radius, search_radius + 1):
        for dx in range(-search_radius, search_radius + 1):
            shifted, in_bounds = _shift(unit, dy, dx)
            shifted_coh, _ = _shift(coherence, dy, dx)
            shifted_valid, _ = _shift(valid, dy, dx)
            pair_valid = valid & shifted_valid & in_bounds
            circular_distance = np.where(pair_valid, 1.0 - np.real(unit * np.conj(shifted)), 1.0)
            patch_distance = uniform_filter(circular_distance, size=patch_size, mode="nearest")
            spatial_weight = np.exp(-(dy * dy + dx * dx) / (2.0 * spatial_sigma * spatial_sigma))
            reliability = np.sqrt(np.maximum(coherence * shifted_coh, 1e-3))
            weight = np.exp(-patch_distance / (bandwidth * bandwidth)) * spatial_weight
            weight *= (0.20 + 0.80 * reliability) * pair_valid
            if dy == 0 and dx == 0:
                weight += valid.astype(float)
            numerator += weight * shifted
            denominator += weight
    estimate = np.divide(numerator, denominator, out=unit.copy(), where=denominator > 1e-12)
    concentration = np.clip(np.abs(estimate), 0.0, 1.0)
    filtered_unit = np.divide(estimate, np.abs(estimate), out=unit.copy(), where=np.abs(estimate) > 1e-12)
    filtered = amplitude * filtered_unit
    filtered[~valid] = 0.0
    return filtered.astype(np.complex64), concentration.astype(np.float32), denominator.astype(np.float32)


def adaptive_goldstein(
    interferogram: np.ndarray,
    coherence: np.ndarray,
    strength: float = 0.40,
    sigma: float = 1.0,
) -> np.ndarray:
    """Apply a mild coherence-adaptive spectral enhancement to a full raster."""
    values = np.asarray(interferogram, dtype=np.complex128)
    spectrum = np.fft.fft2(values)
    magnitude = gaussian_filter(np.abs(spectrum), sigma=sigma, mode="wrap")
    normalized = magnitude / max(float(np.nanmax(magnitude)), 1e-12)
    alpha = float(np.clip(strength * np.nanmedian(1.0 - coherence), 0.0, 0.8))
    enhanced = spectrum * np.power(np.maximum(normalized, 1e-8), alpha)
    result = np.fft.ifft2(enhanced)
    original_amplitude = np.abs(values)
    result_unit = np.divide(result, np.abs(result), out=np.zeros_like(result), where=np.abs(result) > 0)
    return (original_amplitude * result_unit).astype(np.complex64)


def hybrid_filter(
    interferogram: np.ndarray,
    coherence: np.ndarray,
    search_radius: int = 7,
    patch_radius: int = 3,
    h: float = 0.70,
    goldstein_strength: float = 0.40,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nonlocal_result, concentration, support = nonlocal_phase_filter(
        interferogram, coherence, search_radius, patch_radius, h
    )
    filtered = adaptive_goldstein(nonlocal_result, coherence, goldstein_strength)
    return filtered, concentration, support
