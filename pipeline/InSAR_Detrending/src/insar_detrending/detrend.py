"""Robust polynomial surface estimation adapted from the existing SBAS workflow."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np


@dataclass
class DetrendStats:
    accepted: bool
    orbit_accepted: bool
    terrain_accepted: bool
    valid_pixels: int
    initial_fit_pixels: int
    final_fit_pixels: int
    orbit_iterations: int
    terrain_iterations: int
    terrain_fit_method: str
    terrain_bins_used: int
    terrain_local_worst_ratio: float
    terrain_local_bad_blocks: int
    terrain_local_blocks: int
    before_robust_std_rad: float
    after_robust_std_rad: float
    correction_robust_range_rad: float
    orbit_robust_range_rad: float
    terrain_robust_range_rad: float
    turbulent_qa_robust_range_rad: float
    terrain_enabled: bool
    terrain_degree: int
    elevation_phase_correlation_before: float
    elevation_phase_correlation_after: float


def polynomial_terms(x: np.ndarray, y: np.ndarray, degree: int) -> np.ndarray:
    if degree not in (0, 1, 2):
        raise ValueError("Polynomial degree must be 0, 1, or 2")
    terms = [np.ones_like(x)]
    if degree >= 1:
        terms.extend([x, y])
    if degree >= 2:
        terms.extend([x * x, x * y, y * y])
    return np.column_stack(terms)


def orbit_terms(x: np.ndarray, y: np.ndarray, model: str, degree: int) -> np.ndarray:
    """Build DEM-independent orbit-ramp terms.

    The named models keep the ramp geometry more constrained than a full
    quadratic surface, which helps avoid removing broad landslide deformation.
    """
    if model == "degree":
        return polynomial_terms(x, y, degree)
    terms = [np.ones_like(x)]
    if model == "constant":
        return np.column_stack(terms)
    terms.extend([x, y])
    if model == "plane":
        return np.column_stack(terms)
    if model == "range_quadratic":
        terms.append(x * x)
    elif model == "range_cubic":
        terms.extend([x * x, x * x * x])
    elif model == "right_edge_quadratic":
        edge = np.clip((x - 0.35) / 0.65, 0.0, None)
        terms.extend([edge * edge, edge * edge * edge])
    elif model == "azimuth_quadratic":
        terms.append(y * y)
    elif model == "additive_quadratic":
        terms.extend([x * x, y * y])
    elif model == "full_quadratic":
        terms.extend([x * x, x * y, y * y])
    else:
        raise ValueError(f"Unsupported orbit model: {model}")
    return np.column_stack(terms)


def robust_std(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if not values.size:
        return float("nan")
    median = np.median(values)
    return float(1.4826 * np.median(np.abs(values - median)))


def _robust_surface_fit(
    values: np.ndarray,
    design: np.ndarray,
    candidate_mask: np.ndarray,
    mad_scale: float = 3.5,
    max_iterations: int = 5,
) -> tuple[np.ndarray, np.ndarray, int]:
    keep = candidate_mask.ravel().copy()
    initial_count = int(np.count_nonzero(keep))
    if initial_count < design.shape[1] * 3:
        raise ValueError("Too few stable pixels for polynomial detrending")
    iterations = 0
    coefficients = None
    for iteration in range(max_iterations):
        coefficients, *_ = np.linalg.lstsq(design[keep], values.ravel()[keep], rcond=None)
        surface_flat = design @ coefficients
        residual = values.ravel() - surface_flat
        center = np.nanmedian(residual[keep])
        scale = robust_std(residual[keep])
        new_keep = candidate_mask.ravel() & np.isfinite(residual)
        if np.isfinite(scale) and scale > 1e-10:
            new_keep &= np.abs(residual - center) <= mad_scale * scale
        iterations = iteration + 1
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    coefficients, *_ = np.linalg.lstsq(design[keep], values.ravel()[keep], rcond=None)
    surface = (design @ coefficients).reshape(values.shape)
    return surface, keep.reshape(values.shape), iterations


def _robust_range(values: np.ndarray, valid: np.ndarray) -> float:
    if not np.any(valid):
        return float("nan")
    q02, q98 = np.nanpercentile(values[valid], [2, 98])
    return float(q98 - q02)


def _local_worsening_stats(
    before: np.ndarray,
    after: np.ndarray,
    valid_mask: np.ndarray,
    block_pixels: int,
    tolerance: float,
    min_pixels: int,
) -> tuple[float, int, int]:
    """Measure whether a candidate correction makes local residuals worse."""
    block_pixels = max(1, int(block_pixels))
    min_pixels = max(1, int(min_pixels))
    tolerance = max(0.0, float(tolerance))
    height, width = before.shape
    worst_ratio = 1.0
    bad_blocks = 0
    total_blocks = 0
    for row0 in range(0, height, block_pixels):
        row1 = min(row0 + block_pixels, height)
        for col0 in range(0, width, block_pixels):
            col1 = min(col0 + block_pixels, width)
            block_mask = (
                valid_mask[row0:row1, col0:col1]
                & np.isfinite(before[row0:row1, col0:col1])
                & np.isfinite(after[row0:row1, col0:col1])
            )
            if int(np.count_nonzero(block_mask)) < min_pixels:
                continue
            before_std = robust_std(before[row0:row1, col0:col1][block_mask])
            after_std = robust_std(after[row0:row1, col0:col1][block_mask])
            if not (np.isfinite(before_std) and np.isfinite(after_std)) or before_std <= 1e-10:
                continue
            ratio = float(after_std / before_std)
            worst_ratio = max(worst_ratio, ratio)
            total_blocks += 1
            if ratio > 1.0 + tolerance:
                bad_blocks += 1
    return worst_ratio, bad_blocks, total_blocks


def _terrain_design(elevation_normalized: np.ndarray, terrain_degree: int) -> np.ndarray:
    columns = [elevation_normalized.ravel()]
    if terrain_degree >= 2:
        columns.append(elevation_normalized.ravel() ** 2)
    return np.column_stack(columns)


def _spatial_terrain_design(
    elevation_normalized: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    terrain_degree: int,
) -> np.ndarray:
    """DEM-correlated terms whose coefficient varies smoothly in space."""
    h = elevation_normalized.ravel()
    xr = x.ravel()
    yr = y.ravel()
    columns = [h, h * xr, h * yr]
    if terrain_degree >= 2:
        columns.extend([h * xr * xr, h * xr * yr, h * yr * yr, h * h])
    return np.column_stack(columns)


def _binned_terrain_fit(
    values: np.ndarray,
    elevation_normalized: np.ndarray,
    candidate_mask: np.ndarray,
    terrain_degree: int,
    bins: int,
    min_bin_pixels: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Fit residual phase as a function of DEM using elevation-bin medians."""
    candidate = (
        np.asarray(candidate_mask, dtype=bool)
        & np.isfinite(values) & np.isfinite(elevation_normalized)
    )
    h = elevation_normalized[candidate]
    phase = values[candidate]
    if h.size < max(10, (terrain_degree + 1) * min_bin_pixels):
        raise ValueError("Too few stable pixels for binned terrain fitting")
    edges = np.nanpercentile(h, np.linspace(0, 100, max(3, int(bins)) + 1))
    edges = np.unique(edges[np.isfinite(edges)])
    if edges.size < terrain_degree + 2:
        raise ValueError("DEM elevation has insufficient bin variation for terrain fitting")
    centers, medians, weights = [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        in_bin = (h >= left) & (h <= right if right == edges[-1] else h < right)
        count = int(np.count_nonzero(in_bin))
        if count < int(min_bin_pixels):
            continue
        centers.append(float(np.nanmedian(h[in_bin])))
        medians.append(float(np.nanmedian(phase[in_bin])))
        weights.append(float(np.sqrt(count)))
    if len(centers) < terrain_degree + 1:
        raise ValueError("Too few populated DEM bins for terrain fitting")
    centers_array = np.asarray(centers)
    medians_array = np.asarray(medians)
    weights_array = np.asarray(weights)
    design_columns = [centers_array]
    if terrain_degree >= 2:
        design_columns.append(centers_array ** 2)
    design = np.column_stack(design_columns)
    coefficients, *_ = np.linalg.lstsq(
        design * weights_array[:, None], medians_array * weights_array, rcond=None
    )
    full_design = _terrain_design(elevation_normalized, terrain_degree)
    surface = (full_design @ coefficients).reshape(values.shape)
    return surface, candidate, len(centers)


def _binned_terrain_lookup(
    values: np.ndarray,
    elevation_normalized: np.ndarray,
    candidate_mask: np.ndarray,
    bins: int,
    min_bin_pixels: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Estimate a non-parametric residual-vs-DEM relation from after-orbit data."""
    candidate = (
        np.asarray(candidate_mask, dtype=bool)
        & np.isfinite(values) & np.isfinite(elevation_normalized)
    )
    h = elevation_normalized[candidate]
    phase = values[candidate]
    if h.size < max(10, int(min_bin_pixels) * 3):
        raise ValueError("Too few stable pixels for DEM-bin terrain fitting")
    edges = np.nanpercentile(h, np.linspace(0, 100, max(3, int(bins)) + 1))
    edges = np.unique(edges[np.isfinite(edges)])
    if edges.size < 4:
        raise ValueError("DEM elevation has insufficient bin variation for terrain lookup")
    centers, medians = [], []
    for left, right in zip(edges[:-1], edges[1:]):
        in_bin = (h >= left) & (h <= right if right == edges[-1] else h < right)
        if int(np.count_nonzero(in_bin)) < int(min_bin_pixels):
            continue
        centers.append(float(np.nanmedian(h[in_bin])))
        medians.append(float(np.nanmedian(phase[in_bin])))
    if len(centers) < 3:
        raise ValueError("Too few populated DEM bins for terrain lookup")
    centers_array = np.asarray(centers)
    medians_array = np.asarray(medians)
    order = np.argsort(centers_array)
    centers_array = centers_array[order]
    medians_array = medians_array[order]
    surface = np.interp(
        elevation_normalized.ravel(),
        centers_array,
        medians_array,
        left=medians_array[0],
        right=medians_array[-1],
    ).reshape(values.shape)
    surface -= float(np.nanmedian(surface[candidate]))
    surface[~np.isfinite(values)] = np.nan
    return surface, candidate, len(centers)


def _nanmean_filter_axis(values: np.ndarray, radius: int, axis: int) -> np.ndarray:
    if radius <= 0:
        return values.copy()
    moved = np.moveaxis(values, axis, 0)
    valid = np.isfinite(moved)
    filled = np.where(valid, moved, 0.0)
    weights = valid.astype(float)
    pad_width = [(0, 0)] * filled.ndim
    pad_width[0] = (radius, radius)
    filled_pad = np.pad(filled, pad_width, mode="constant", constant_values=0.0)
    weights_pad = np.pad(weights, pad_width, mode="constant", constant_values=0.0)
    leading = np.zeros((1, *filled_pad.shape[1:]), dtype=float)
    sums = np.concatenate([leading, np.cumsum(filled_pad, axis=0)], axis=0)
    counts = np.concatenate([leading, np.cumsum(weights_pad, axis=0)], axis=0)
    window = 2 * radius + 1
    sums = sums[window:] - sums[:-window]
    counts = counts[window:] - counts[:-window]
    with np.errstate(invalid="ignore", divide="ignore"):
        filtered = sums / counts
    filtered[counts == 0] = np.nan
    return np.moveaxis(filtered, 0, axis)


def _local_terrain_fit(
    values: np.ndarray,
    elevation_normalized: np.ndarray,
    candidate_mask: np.ndarray,
    radius_pixels: int,
    coefficient_clip_percentile: float = 98.0,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Estimate a spatially varying DEM coefficient from after-orbit residuals.

    The coefficient field is the local least-squares slope between residual
    phase and normalized elevation, estimated only from stable fit pixels:
    beta(x, y) = smooth(H * phase) / smooth(H * H).
    """
    radius_pixels = max(1, int(radius_pixels))
    candidate = (
        np.asarray(candidate_mask, dtype=bool)
        & np.isfinite(values) & np.isfinite(elevation_normalized)
    )
    if np.count_nonzero(candidate) < 100:
        raise ValueError("Too few stable pixels for local terrain fitting")
    numerator = np.where(candidate, elevation_normalized * values, np.nan)
    denominator = np.where(candidate, elevation_normalized * elevation_normalized, np.nan)
    for _ in range(2):
        numerator = _nanmean_filter_axis(numerator, radius_pixels, axis=0)
        numerator = _nanmean_filter_axis(numerator, radius_pixels, axis=1)
        denominator = _nanmean_filter_axis(denominator, radius_pixels, axis=0)
        denominator = _nanmean_filter_axis(denominator, radius_pixels, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = numerator / denominator
    valid_beta = np.isfinite(beta) & (np.abs(denominator) > 1e-10)
    if np.count_nonzero(valid_beta) < 100:
        raise ValueError("Local terrain fitting produced too few finite coefficients")
    clip = float(np.nanpercentile(np.abs(beta[valid_beta]), coefficient_clip_percentile))
    if np.isfinite(clip) and clip > 0:
        beta = np.clip(beta, -clip, clip)
    beta = _nanmean_filter_axis(beta, radius_pixels, axis=0)
    beta = _nanmean_filter_axis(beta, radius_pixels, axis=1)
    surface = beta * elevation_normalized
    surface[~np.isfinite(values)] = np.nan
    return surface, candidate, int(np.count_nonzero(valid_beta))


def _hybrid_terrain_fit(
    values: np.ndarray,
    elevation_normalized: np.ndarray,
    candidate_mask: np.ndarray,
    bins: int,
    min_bin_pixels: int,
    radius_pixels: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Combine a DEM-bin residual curve with local after-orbit modulation.

    The binned curve captures the stable-pixel phase that is systematically
    organized by elevation. The local term is estimated only from the residual
    that remains after that DEM-bin correction, so it can absorb broad spatial
    changes in the DEM coefficient without becoming a generic low-pass filter.
    """
    binned_surface, candidate, bins_used = _binned_terrain_lookup(
        values, elevation_normalized, candidate_mask, bins, min_bin_pixels
    )
    residual = values - binned_surface
    local_surface, local_candidate, local_support = _local_terrain_fit(
        residual, elevation_normalized, candidate, radius_pixels
    )
    surface = binned_surface + local_surface
    surface -= float(np.nanmedian(surface[candidate]))
    surface[~np.isfinite(values)] = np.nan
    return surface, candidate & local_candidate, max(bins_used, local_support)


def lowpass_residual(values: np.ndarray, fit_mask: np.ndarray,
                     radius_pixels: int, passes: int = 3) -> np.ndarray:
    """Estimate a non-DEM long-wavelength residual from stable pixels only.

    This field is diagnostic by default. It is not subtracted from the final
    interferogram because low-frequency atmospheric phase and broad landslide
    deformation cannot be separated reliably from one interferogram alone.
    """
    residual = np.where(fit_mask, values, np.nan).astype(float)
    for _ in range(max(1, int(passes))):
        residual = _nanmean_filter_axis(residual, int(radius_pixels), axis=0)
        residual = _nanmean_filter_axis(residual, int(radius_pixels), axis=1)
    residual[~np.isfinite(values)] = np.nan
    return residual


def robust_polynomial_detrend(
    image: np.ndarray,
    degree: int = 2,
    fit_mask: np.ndarray | None = None,
    exclude_mask: np.ndarray | None = None,
    mad_scale: float = 3.5,
    max_iterations: int = 5,
    preserve_median: bool = False,
    elevation: np.ndarray | None = None,
    terrain_degree: int = 2,
    terrain_fit_method: str = "spatial",
    terrain_bins: int = 30,
    terrain_min_bin_pixels: int = 100,
    terrain_local_radius_pixels: int = 80,
    terrain_strength: float = 0.3,
    terrain_max_range_fraction: float = 0.6,
    terrain_local_guard_pixels: int = 96,
    terrain_local_guard_tolerance: float = 0.15,
    terrain_local_guard_min_pixels: int = 200,
    orbit_model: str = "degree",
    turbulent_qa: bool = True,
    turbulent_smoothing_pixels: int = 60,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, DetrendStats, dict[str, np.ndarray]]:
    """Sequentially estimate orbit and DEM-correlated phase corrections.

    The removed correction is split into (1) a DEM-independent polynomial orbit
    ramp and (2) a low-pass DEM-correlated term. A non-DEM low-frequency
    residual is optionally saved only as a QA layer, not removed by default.
    """
    image = np.asarray(image, dtype=float)
    if terrain_degree not in (0, 1, 2):
        raise ValueError("Terrain degree must be 0, 1, or 2")
    valid = np.isfinite(image)
    if fit_mask is None:
        fit_mask = valid.copy()
    else:
        fit_mask = np.asarray(fit_mask, dtype=bool) & valid
    if exclude_mask is not None:
        fit_mask &= ~np.asarray(exclude_mask, dtype=bool)
    terrain_enabled = elevation is not None and terrain_degree > 0
    if terrain_enabled:
        elevation = np.asarray(elevation, dtype=float)
        if elevation.shape != image.shape:
            raise ValueError("Elevation and phase shapes do not match")
        fit_mask &= np.isfinite(elevation)
    height, width = image.shape
    yy, xx = np.mgrid[:height, :width]
    x = 2.0 * xx / max(width - 1, 1) - 1.0
    y = 2.0 * yy / max(height - 1, 1) - 1.0
    orbit_design = orbit_terms(x.ravel(), y.ravel(), orbit_model, degree)
    initial_count = int(np.count_nonzero(fit_mask))

    orbit_surface, orbit_mask, orbit_iterations = _robust_surface_fit(
        image, orbit_design, fit_mask, mad_scale, max_iterations
    )
    if preserve_median:
        orbit_surface -= np.nanmedian(orbit_surface[orbit_mask])
    orbit_corrected = image - orbit_surface
    before_scale = robust_std(image[orbit_mask])
    orbit_after_scale = robust_std(orbit_corrected[orbit_mask])
    orbit_accepted = bool(
        np.isfinite(before_scale) and np.isfinite(orbit_after_scale)
        and orbit_after_scale < before_scale
    )
    if not orbit_accepted:
        orbit_surface = np.zeros_like(image)
        orbit_corrected = image.copy()
        orbit_mask = fit_mask.copy()
        orbit_after_scale = before_scale
    after_orbit = orbit_corrected.copy()

    terrain_surface = np.zeros_like(image)
    terrain_iterations = 0
    terrain_bins_used = 0
    terrain_local_worst_ratio = float("nan")
    terrain_local_bad_blocks = 0
    terrain_local_blocks = 0
    terrain_accepted = False
    final_mask = orbit_mask.copy()
    elevation_normalized = None
    corrected = orbit_corrected.copy()
    if terrain_enabled:
        center_h = float(np.nanmedian(elevation[fit_mask]))
        scale_h = float(np.nanstd(elevation[fit_mask]))
        if not np.isfinite(scale_h) or scale_h < 1e-10:
            raise ValueError("DEM elevation has insufficient variation for terrain fitting")
        elevation_normalized = (elevation - center_h) / scale_h
        if terrain_fit_method == "binned":
            terrain_surface_candidate, terrain_mask, terrain_bins_used = _binned_terrain_fit(
                orbit_corrected, elevation_normalized, orbit_mask, terrain_degree,
                terrain_bins, terrain_min_bin_pixels,
            )
            terrain_iterations = 1
        elif terrain_fit_method == "spatial":
            terrain_design = _spatial_terrain_design(
                elevation_normalized, x, y, terrain_degree
            )
            terrain_surface_candidate, terrain_mask, terrain_iterations = _robust_surface_fit(
                orbit_corrected, terrain_design, orbit_mask, mad_scale, max_iterations
            )
        elif terrain_fit_method == "local":
            terrain_surface_candidate, terrain_mask, terrain_bins_used = _local_terrain_fit(
                orbit_corrected, elevation_normalized, orbit_mask,
                terrain_local_radius_pixels,
            )
            terrain_iterations = 1
        elif terrain_fit_method == "hybrid":
            terrain_surface_candidate, terrain_mask, terrain_bins_used = _hybrid_terrain_fit(
                orbit_corrected, elevation_normalized, orbit_mask,
                terrain_bins, terrain_min_bin_pixels, terrain_local_radius_pixels,
            )
            terrain_iterations = 2
        elif terrain_fit_method == "pixel":
            terrain_design = _terrain_design(elevation_normalized, terrain_degree)
            terrain_surface_candidate, terrain_mask, terrain_iterations = _robust_surface_fit(
                orbit_corrected, terrain_design, orbit_mask, mad_scale, max_iterations
            )
        else:
            raise ValueError(f"Unsupported terrain fit method: {terrain_fit_method}")
        terrain_scale_factor = float(terrain_strength)
        if not np.isfinite(terrain_scale_factor) or terrain_scale_factor < 0:
            raise ValueError("terrain_strength must be a finite non-negative value")
        terrain_scale_factor = min(terrain_scale_factor, 1.0)
        candidate_range = _robust_range(terrain_surface_candidate, terrain_mask)
        residual_range = _robust_range(orbit_corrected, terrain_mask)
        if (
            np.isfinite(candidate_range) and candidate_range > 1e-10
            and np.isfinite(residual_range) and residual_range > 0
            and np.isfinite(terrain_max_range_fraction)
            and terrain_max_range_fraction > 0
        ):
            terrain_scale_factor = min(
                terrain_scale_factor,
                float(terrain_max_range_fraction) * residual_range / candidate_range,
            )
        terrain_surface_candidate = terrain_surface_candidate * terrain_scale_factor
        terrain_corrected = orbit_corrected - terrain_surface_candidate
        terrain_before_scale = robust_std(orbit_corrected[terrain_mask])
        terrain_after_scale = robust_std(terrain_corrected[terrain_mask])
        terrain_guard_mask = valid.copy()
        if exclude_mask is not None:
            terrain_guard_mask &= ~np.asarray(exclude_mask, dtype=bool)
        terrain_local_worst_ratio, terrain_local_bad_blocks, terrain_local_blocks = (
            _local_worsening_stats(
                orbit_corrected, terrain_corrected, terrain_guard_mask,
                terrain_local_guard_pixels, terrain_local_guard_tolerance,
                terrain_local_guard_min_pixels,
            )
        )
        terrain_accepted = bool(
            np.isfinite(terrain_before_scale) and np.isfinite(terrain_after_scale)
            and terrain_after_scale < terrain_before_scale
            and terrain_local_bad_blocks == 0
        )
        if terrain_accepted:
            terrain_surface = terrain_surface_candidate
            corrected = terrain_corrected
            final_mask = terrain_mask

    surface = orbit_surface + terrain_surface
    corrected[~valid] = np.nan
    surface[~valid] = np.nan
    orbit_surface[~valid] = np.nan
    terrain_surface[~valid] = np.nan
    after_scale = robust_std(corrected[final_mask])
    accepted = bool(orbit_accepted or terrain_accepted)
    correlation_before = correlation_after = float("nan")
    if terrain_enabled:
        correlation_mask = final_mask & np.isfinite(elevation_normalized)
        if np.count_nonzero(correlation_mask) >= 3:
            correlation_before = float(np.corrcoef(
                orbit_corrected[correlation_mask], elevation_normalized[correlation_mask]
            )[0, 1])
            # Correlation is numerically meaningless when the residual is
            # essentially constant at floating-point precision.
            if robust_std(corrected[correlation_mask]) < 1e-10:
                correlation_after = 0.0
            else:
                correlation_after = float(np.corrcoef(
                    corrected[correlation_mask], elevation_normalized[correlation_mask]
                )[0, 1])
    turbulent_surface = np.full_like(image, np.nan, dtype=float)
    if turbulent_qa and int(turbulent_smoothing_pixels) > 0:
        turbulent_surface = lowpass_residual(
            corrected, final_mask, int(turbulent_smoothing_pixels)
        )
    after_dem = corrected.copy()
    after_turbulent_qa = corrected - turbulent_surface
    after_turbulent_qa[~valid] = np.nan
    stats = DetrendStats(
        accepted=accepted, orbit_accepted=orbit_accepted, terrain_accepted=terrain_accepted,
        valid_pixels=int(np.count_nonzero(valid)), initial_fit_pixels=initial_count,
        final_fit_pixels=int(np.count_nonzero(final_mask)),
        orbit_iterations=orbit_iterations, terrain_iterations=terrain_iterations,
        terrain_fit_method=terrain_fit_method,
        terrain_bins_used=terrain_bins_used,
        terrain_local_worst_ratio=terrain_local_worst_ratio,
        terrain_local_bad_blocks=terrain_local_bad_blocks,
        terrain_local_blocks=terrain_local_blocks,
        before_robust_std_rad=before_scale,
        after_robust_std_rad=after_scale,
        correction_robust_range_rad=_robust_range(surface, valid),
        orbit_robust_range_rad=_robust_range(orbit_surface, valid),
        terrain_robust_range_rad=_robust_range(terrain_surface, valid),
        turbulent_qa_robust_range_rad=_robust_range(turbulent_surface, valid),
        terrain_enabled=terrain_enabled, terrain_degree=terrain_degree if terrain_enabled else 0,
        elevation_phase_correlation_before=correlation_before,
        elevation_phase_correlation_after=correlation_after,
    )
    components = {
        "orbit": orbit_surface,
        "terrain": terrain_surface,
        "turbulent_qa": turbulent_surface,
        "after_orbit": after_orbit,
        "after_dem": after_dem,
        "after_turbulent_qa": after_turbulent_qa,
    }
    return corrected, surface, final_mask, stats, components


def build_consensus_stable_mask(stack: np.ndarray, valid_fraction: float = 0.8,
                                dispersion_percentile: float = 70.0) -> np.ndarray:
    """Find broadly stable pixels without external reference observations."""
    stack = np.asarray(stack, dtype=float)
    valid = np.mean(np.isfinite(stack), axis=0) >= valid_fraction
    centered = stack - np.nanmedian(stack, axis=(1, 2), keepdims=True)
    # All-NaN water pixels are intentionally excluded; suppress only NumPy's
    # expected empty-slice warning while preserving their NaN dispersion.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        temporal_center = np.nanmedian(centered, axis=0)
        temporal_mad = np.nanmedian(np.abs(centered - temporal_center), axis=0)
    threshold = np.nanpercentile(temporal_mad[valid], dispersion_percentile) if np.any(valid) else np.nan
    return valid & np.isfinite(temporal_mad) & (temporal_mad <= threshold)
