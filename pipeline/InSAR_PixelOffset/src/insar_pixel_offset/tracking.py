"""Amplitude normalized-cross-correlation pixel-offset tracker."""

from __future__ import annotations

import numpy as np
from scipy.signal import correlate2d


def grid_centers(shape: tuple[int, int], window: int, search: int, stride: int) -> tuple[np.ndarray, np.ndarray]:
    if window < 8 or window % 2:
        raise ValueError("window must be an even integer of at least 8 pixels")
    if search < 1 or stride < 1:
        raise ValueError("search and stride must be positive")
    margin = window // 2 + search
    rows = np.arange(margin, shape[0] - margin, stride, dtype=int)
    cols = np.arange(margin, shape[1] - margin, stride, dtype=int)
    if not rows.size or not cols.size:
        raise ValueError("GSLC image is too small for the selected window and search radius")
    return rows, cols


def normalized_offset(reference: np.ndarray, secondary: np.ndarray, row: int, col: int,
                      window: int, search: int) -> tuple[float, float, float]:
    """Return secondary-minus-reference (row, col) offset and NCC peak."""
    half = window // 2
    template = np.abs(reference[row - half:row + half, col - half:col + half]).astype(float)
    chip = np.abs(secondary[
        row - half - search:row + half + search,
        col - half - search:col + half + search,
    ]).astype(float)
    if template.shape != (window, window) or chip.shape != (window + 2 * search, window + 2 * search):
        return np.nan, np.nan, np.nan
    if not (np.all(np.isfinite(template)) and np.all(np.isfinite(chip))):
        return np.nan, np.nan, np.nan
    template -= template.mean()
    template_energy = float(np.sum(template * template))
    if template_energy <= 1e-12:
        return np.nan, np.nan, np.nan
    count = float(window * window)
    sums = correlate2d(chip, np.ones((window, window)), mode="valid")
    sums_sq = correlate2d(chip * chip, np.ones((window, window)), mode="valid")
    energy = np.maximum(sums_sq - sums * sums / count, 0.0)
    numerator = correlate2d(chip, template, mode="valid")
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = numerator / np.sqrt(energy * template_energy)
    if not np.any(np.isfinite(correlation)):
        return np.nan, np.nan, np.nan
    peak_row, peak_col = np.unravel_index(np.nanargmax(correlation), correlation.shape)
    return float(peak_row - search), float(peak_col - search), float(correlation[peak_row, peak_col])
