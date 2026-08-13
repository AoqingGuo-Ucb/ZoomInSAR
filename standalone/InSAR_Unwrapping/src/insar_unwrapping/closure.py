"""Temporal phase-closure diagnostics and integer-cycle repair."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, label, median_filter

TWOPI = 2.0 * np.pi


def build_triangles(pairs: list[tuple[str, str]]) -> list[tuple[int, int, int]]:
    lookup = {pair: index for index, pair in enumerate(pairs)}
    dates = sorted({date for pair in pairs for date in pair})
    triangles: list[tuple[int, int, int]] = []
    for ia, a in enumerate(dates):
        for ib in range(ia + 1, len(dates)):
            b = dates[ib]
            for c in dates[ib + 1:]:
                if (a, b) in lookup and (b, c) in lookup and (a, c) in lookup:
                    triangles.append((lookup[(a, b)], lookup[(b, c)], lookup[(a, c)]))
    return triangles


def closure_cycles(stack: np.ndarray, triangle: tuple[int, int, int]) -> np.ndarray:
    iab, ibc, iac = triangle
    return (stack[iab] + stack[ibc] - stack[iac]) / TWOPI


@dataclass
class ClosureStats:
    triangles: int
    iterations: int
    components_repaired: int
    pixels_repaired: int
    bad_pixels_before: int
    bad_pixels_after: int
    spatially_rejected_components: int = 0
    bad_pixels_history: tuple[int, ...] = ()


def _boundary_median_difference(
    phase: np.ndarray, component: np.ndarray, width: int = 3
) -> float:
    """Return the absolute median phase step across a component boundary.

    Compare multi-pixel inner and outer boundary bands. If a component is too
    narrow for the requested inner width, erosion vanishes naturally and all
    available component pixels are used. NaN/masked pixels (for example water)
    are excluded. A region without a measurable two-sided boundary cannot be
    spatially validated and returns infinity.
    """
    width = max(1, int(width))
    inner = component & ~binary_erosion(component, iterations=width)
    outer = binary_dilation(component, iterations=width) & ~component
    inner_values = phase[inner & np.isfinite(phase)]
    outer_values = phase[outer & np.isfinite(phase)]
    if inner_values.size == 0 or outer_values.size == 0:
        return float("inf")
    return float(abs(np.nanmedian(inner_values) - np.nanmedian(outer_values)))


def _bad_pixel_count(stack: np.ndarray, triangles: list[tuple[int, int, int]], tolerance: float) -> int:
    total = 0
    for triangle in triangles:
        cycles = closure_cycles(stack, triangle)
        total += int(np.count_nonzero(np.isfinite(cycles) & (np.abs(cycles - np.rint(cycles)) < 0.20) & (np.abs(np.rint(cycles)) >= tolerance)))
    return total


def repair_phase_closure(
    stack: np.ndarray,
    pairs: list[tuple[str, str]],
    coherence: np.ndarray | None = None,
    max_iterations: int = 8,
    tolerance_cycles: float = 0.45,
    min_component_pixels: int = 20,
    boundary_width: int = 3,
) -> tuple[np.ndarray, np.ndarray, ClosureStats]:
    """Iteratively repair integer closure errors in interferogram triangles.

    Integer errors are peeled one cycle per iteration using cumulative masks:
    a +2 region participates in the first +1 layer together with surrounding +1
    pixels, then its remaining cycle is reconsidered after closure is recomputed.
    For every layer component, all three triangle edges are tried independently.
    A correction is eligible only if it lowers both the local triangle closure
    magnitude and the median phase step across that component's boundary. The
    eligible edge with the largest spatial-boundary improvement is selected;
    coherence is used only as a secondary tie-breaker.
    """
    corrected = np.asarray(stack, dtype=float).copy()
    triangles = build_triangles(pairs)
    correction_map = np.zeros(stack.shape, dtype=np.int16)
    before = _bad_pixel_count(corrected, triangles, tolerance_cycles)
    bad_pixels_history = [before]
    repaired_components = repaired_pixels = completed_iterations = 0
    spatially_rejected_components = 0
    if coherence is None:
        coherence = np.ones(stack.shape, dtype=float)
    for iteration in range(max_iterations):
        changed = 0
        for iab, ibc, iac in triangles:
            cycles = closure_cycles(corrected, (iab, ibc, iac))
            raw_rounded = np.rint(np.nan_to_num(cycles)).astype(int)
            # Use the median-filtered map only to suppress isolated salt-and-pepper
            # errors, but retain boundary pixels that agree with a supported core.
            core = np.rint(median_filter(np.nan_to_num(cycles), size=3)).astype(int)
            rounded = np.where((raw_rounded == core) | (core == 0), raw_rounded, core)
            finite = np.isfinite(cycles)
            # Peel one signed cycle at a time. Cumulative masks preserve nested
            # regions: +2 pixels belong to the first >=+1 layer, and -2 pixels
            # belong to the first <=-1 layer. Exact-value masks would cut holes
            # around nested errors and force all cycles onto a single edge.
            for layer_sign in (1, -1):
                if layer_sign > 0:
                    layer_mask = finite & (rounded >= max(1, tolerance_cycles))
                else:
                    layer_mask = finite & (rounded <= -max(1, tolerance_cycles))
                labels, count = label(layer_mask)
                for component_id in range(1, count + 1):
                    component = labels == component_id
                    pixels = int(np.count_nonzero(component))
                    if pixels < min_component_pixels:
                        continue
                    edge_indices = (iab, ibc, iac)
                    old_error = float(np.nanmedian(np.abs(cycles[component])))
                    eligible = []
                    any_closure_improved = False
                    # Closure signs are +AB +BC -AC. Try assigning the integer
                    # ambiguity to each edge instead of assuming the lowest-
                    # coherence edge is necessarily wrong.
                    for selected_position, selected_edge in enumerate(edge_indices):
                        sign = (1, 1, -1)[selected_position]
                        applied = -sign * layer_sign
                        old_boundary_error = _boundary_median_difference(
                            corrected[selected_edge], component, boundary_width
                        )
                        trial = corrected[selected_edge].copy()
                        trial[component] += applied * TWOPI
                        new_boundary_error = _boundary_median_difference(
                            trial, component, boundary_width
                        )
                        if selected_edge == iab:
                            trial_closure = trial + corrected[ibc] - corrected[iac]
                        elif selected_edge == ibc:
                            trial_closure = corrected[iab] + trial - corrected[iac]
                        else:
                            trial_closure = corrected[iab] + corrected[ibc] - trial
                        new_error = float(np.nanmedian(
                            np.abs(trial_closure[component] / TWOPI)
                        ))
                        closure_improved = new_error + 0.05 < old_error
                        any_closure_improved |= closure_improved
                        boundary_improvement = old_boundary_error - new_boundary_error
                        if closure_improved and boundary_improvement > 0:
                            quality = float(np.nanmedian(
                                coherence[selected_edge][component]
                            ))
                            if not np.isfinite(quality):
                                quality = 1.0
                            eligible.append((
                                boundary_improvement,
                                old_error - new_error,
                                -quality,
                                selected_edge,
                                applied,
                                trial,
                            ))

                    if eligible:
                        _, _, _, selected_edge, applied, trial = max(
                            eligible, key=lambda item: item[:3]
                        )
                        corrected[selected_edge] = trial
                        correction_map[selected_edge][component] += applied
                        repaired_components += 1; repaired_pixels += pixels; changed += pixels
                    elif any_closure_improved:
                        spatially_rejected_components += 1
        completed_iterations = iteration + 1
        current_bad = _bad_pixel_count(corrected, triangles, tolerance_cycles)
        bad_pixels_history.append(current_bad)
        # Stop at convergence. Accepted local changes can interact through
        # shared triangle edges; continuing is useful only while the global
        # unresolved integer-closure count strictly decreases.
        if not changed or current_bad >= bad_pixels_history[-2]:
            break
    after = _bad_pixel_count(corrected, triangles, tolerance_cycles)
    return corrected, correction_map, ClosureStats(
        triangles=len(triangles), iterations=completed_iterations,
        components_repaired=repaired_components, pixels_repaired=repaired_pixels,
        bad_pixels_before=before, bad_pixels_after=after,
        spatially_rejected_components=spatially_rejected_components,
        bad_pixels_history=tuple(bad_pixels_history),
    )
