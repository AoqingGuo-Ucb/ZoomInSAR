"""Spatial integer-cycle and low-frequency corrections."""

from __future__ import annotations

import heapq
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    binary_closing,
    distance_transform_edt,
    gaussian_filter,
    label,
)

TWOPI = 2.0 * np.pi


def normalized_gaussian(values: np.ndarray, valid: np.ndarray, sigma: float) -> np.ndarray:
    numerator = gaussian_filter(np.where(valid, values, 0.0), sigma=sigma, mode="nearest")
    denominator = gaussian_filter(valid.astype(float), sigma=sigma, mode="nearest")
    return np.divide(numerator, denominator, out=np.full(values.shape, np.nan), where=denominator > 1e-6)


def boundary_score(phase: np.ndarray, component: np.ndarray) -> float:
    ring = binary_dilation(component, iterations=1) & ~component & np.isfinite(phase)
    inside = binary_dilation(ring, iterations=1) & component & np.isfinite(phase)
    if not np.any(ring) or not np.any(inside):
        return np.inf
    return float(abs(np.nanmedian(phase[inside]) - np.nanmedian(phase[ring])))


@dataclass
class SpatialRepairStats:
    components_tested: int = 0
    components_repaired: int = 0
    pixels_repaired: int = 0


@dataclass
class RegionGraphRepairStats:
    regions_detected: int = 0
    regions_tested: int = 0
    regions_repaired: int = 0
    pixels_repaired: int = 0
    iterations: int = 0
    boundary_paths_completed: int = 0
    boundary_pixels_added: int = 0


@dataclass
class EdgeRepairStats:
    components_detected: int = 0
    components_tested: int = 0
    components_repaired: int = 0
    pixels_repaired: int = 0
    rejected_noninteger: int = 0
    rejected_low_support: int = 0
    rejected_low_coherence: int = 0
    rejected_no_improvement: int = 0


def _robust_quadratic_trend(
    phase: np.ndarray, valid: np.ndarray, coherence: np.ndarray | None
) -> np.ndarray:
    """Robust degree-two background used only to propose edge candidates."""
    yy,xx=np.indices(phase.shape,dtype=float)
    x=2.0*xx/max(phase.shape[1]-1,1)-1.0
    y=2.0*yy/max(phase.shape[0]-1,1)-1.0
    design=np.stack((np.ones_like(x),x,y,x*x,x*y,y*y),axis=-1)
    # Edge-connected plateaus must not define their own background. Prefer an
    # interior reference obtained by eroding the valid domain; fall back to all
    # valid pixels only for very narrow scenes.
    interior=binary_erosion(valid,iterations=min(30,max(1,min(phase.shape)//3)))
    mask=interior if np.count_nonzero(interior)>=500 else valid.copy()
    if coherence is not None:
        coh=np.asarray(coherence,float)
        finite_coh=coh[valid & np.isfinite(coh)]
        if finite_coh.size:
            mask &= coh >= max(0.05,float(np.nanpercentile(finite_coh,25)))
    # Evenly subsample large rasters so the fit remains cheap for a full stack.
    indices=np.flatnonzero(mask)
    if indices.size>15000: indices=indices[::int(np.ceil(indices.size/15000))]
    keep=np.ones(indices.size,dtype=bool)
    coefficients=np.zeros(6)
    for _ in range(5):
        use=indices[keep]
        if use.size<30: break
        coefficients=np.linalg.lstsq(design.reshape(-1,6)[use],phase.ravel()[use],rcond=None)[0]
        residual=phase.ravel()[indices]-design.reshape(-1,6)[indices]@coefficients
        center=float(np.median(residual[keep]));mad=float(np.median(np.abs(residual[keep]-center)))
        scale=max(1.4826*mad,0.15*TWOPI)
        new_keep=np.abs(residual-center)<=2.8*scale
        if np.array_equal(new_keep,keep): break
        keep=new_keep
    return np.sum(design*coefficients,axis=-1)


def repair_edge_connected_cycles(
    phase: np.ndarray,
    coherence: np.ndarray | None = None,
    boundary_width: int = 5,
    min_pixels: int = 40,
    max_fraction: float = 0.35,
    integer_tolerance_cycles: float = 0.25,
    minimum_integer_support: float = 0.50,
    minimum_improvement: float = 0.30,
) -> tuple[np.ndarray,np.ndarray,EdgeRepairStats]:
    """Repair integer-cycle plateaus connected to image or mask boundaries."""
    corrected=np.asarray(phase,float).copy();valid=np.isfinite(corrected)
    cycle_map=np.zeros(corrected.shape,dtype=np.int16);stats=EdgeRepairStats()
    valid_count=int(np.count_nonzero(valid))
    if valid_count<max(3*min_pixels,30): return corrected,cycle_map,stats
    trend=_robust_quadratic_trend(corrected,valid,coherence)
    residual_cycles=(corrected-trend)/TWOPI
    nearest=np.zeros(corrected.shape,dtype=np.int16)
    nearest[valid]=np.rint(residual_cycles[valid]).astype(np.int16)
    domain_band=binary_dilation(_domain_boundary(valid),iterations=2)
    for integer in np.unique(nearest[valid]):
        integer=int(integer)
        if integer==0: continue
        close=np.abs(residual_cycles-integer)<=integer_tolerance_cycles
        candidates=binary_closing(valid & close & (nearest==integer),iterations=2)
        labels,count=label(candidates)
        for component_id in range(1,count+1):
            component=labels==component_id
            pixels=int(np.count_nonzero(component))
            if pixels<min_pixels or pixels>max_fraction*valid_count or not np.any(component&domain_band): continue
            stats.components_detected+=1;stats.components_tested+=1
            before,difference,_,support=_paired_boundary_measure(corrected,component,boundary_width)
            observed=difference/TWOPI
            boundary_integer=int(np.rint(observed))
            # The robust trend and the local boundary must agree on direction
            # and integer magnitude; otherwise this may be atmosphere/deformation.
            if boundary_integer!=integer or abs(observed-integer)>integer_tolerance_cycles:
                stats.rejected_noninteger+=1
                continue
            if support<minimum_integer_support:
                stats.rejected_low_support+=1
                continue
            if coherence is not None:
                values=np.asarray(coherence,float)[component]
                if np.any(np.isfinite(values)) and np.nanmedian(values)<0.05:
                    stats.rejected_low_coherence+=1
                    continue
            trial=corrected.copy();trial[component]-=integer*TWOPI
            after,_,_,_=_paired_boundary_measure(trial,component,boundary_width)
            gain=(before-after)/max(before,1e-9)
            if np.isfinite(after) and gain>=minimum_improvement and after<0.55*TWOPI:
                corrected=trial;cycle_map[component]-=integer
                stats.components_repaired+=1;stats.pixels_repaired+=pixels
            else:
                stats.rejected_no_improvement+=1
    return corrected,cycle_map,stats


def _domain_boundary(valid: np.ndarray) -> np.ndarray:
    """Valid pixels adjacent to the image edge or an invalid/water pixel."""
    result = valid & binary_dilation(~valid, iterations=1)
    result[0] |= valid[0]; result[-1] |= valid[-1]
    result[:, 0] |= valid[:, 0]; result[:, -1] |= valid[:, -1]
    return result


def _least_cost_path(
    cost: np.ndarray, start: tuple[int, int], goal: tuple[int, int], margin: int = 4
) -> list[tuple[int, int]]:
    """Small bounded Dijkstra search used only for short boundary gaps."""
    h, w = cost.shape
    y0, x0 = start; y1, x1 = goal
    lo_y=max(0,min(y0,y1)-margin); hi_y=min(h,max(y0,y1)+margin+1)
    lo_x=max(0,min(x0,x1)-margin); hi_x=min(w,max(x0,x1)+margin+1)
    queue=[(0.0,y0,x0)]; distance={(y0,x0):0.0}; previous={}
    while queue:
        value,y,x=heapq.heappop(queue)
        if value != distance.get((y,x)): continue
        if (y,x)==(y1,x1):
            path=[]; node=(y,x)
            while node != (y0,x0): path.append(node); node=previous[node]
            path.append((y0,x0)); return path[::-1]
        for ny,nx in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
            if not (lo_y<=ny<hi_y and lo_x<=nx<hi_x) or not np.isfinite(cost[ny,nx]): continue
            trial=value+float(cost[ny,nx])
            if trial < distance.get((ny,nx),np.inf):
                distance[(ny,nx)]=trial; previous[(ny,nx)]=(y,x)
                heapq.heappush(queue,(trial,ny,nx))
    return []


def _complete_boundary_gaps(
    barrier: np.ndarray,
    valid: np.ndarray,
    phase_gradient: np.ndarray,
    coherence: np.ndarray | None,
    branch_cuts: np.ndarray | None,
    max_gap: int,
    boundary_snap: int,
) -> tuple[np.ndarray, int, int]:
    """Conservatively join short edge gaps and snap them to domain boundaries."""
    completed=barrier.copy(); paths=added=0
    coh=np.nan_to_num(coherence, nan=0.0) if coherence is not None else np.ones(barrier.shape)
    grad=phase_gradient/max(float(np.nanpercentile(phase_gradient[valid],95)),1e-6)
    cuts=np.asarray(branch_cuts,bool) if branch_cuts is not None else np.zeros(barrier.shape,bool)
    # Low coherence, high phase gradient and existing branch cuts make a pixel
    # a cheaper and therefore more plausible continuation of a jump boundary.
    cost=1.0+1.5*np.clip(coh,0,1)+1.5*(1-np.clip(grad,0,1))-1.0*cuts
    cost=np.where(valid,np.maximum(cost,0.15),np.inf)
    domain=_domain_boundary(valid)

    for _ in range(2):
        labels,count=label(completed)
        components=[np.column_stack(np.where(labels==i)) for i in range(1,count+1)]
        components=[p for p in components if len(p)>=4]
        # Snap a nearby edge segment to the valid-domain boundary.
        if np.any(domain):
            dist,idx=distance_transform_edt(~domain,return_indices=True)
            for points in components:
                d=dist[points[:,0],points[:,1]]; j=int(np.argmin(d))
                if 1 < d[j] <= boundary_snap:
                    start=tuple(map(int,points[j])); goal=(int(idx[0][start]),int(idx[1][start]))
                    path=_least_cost_path(cost,start,goal)
                    if path and np.mean([cost[y,x] for y,x in path]) <= 3.2:
                        before=np.count_nonzero(completed); yy,xx=zip(*path);completed[yy,xx]=True
                        delta=np.count_nonzero(completed)-before
                        if delta: paths+=1;added+=delta
        # Join the closest pair of distinct segments when their gap is short.
        labels,count=label(completed)
        comps=[np.column_stack(np.where(labels==i)) for i in range(1,count+1)]
        comps=[p for p in comps if len(p)>=4]
        best=None
        for i in range(len(comps)):
            for j in range(i+1,len(comps)):
                # Small images make this vectorised nearest-pair search cheap;
                # subsampling caps work for dense edge components.
                a=comps[i][::max(1,len(comps[i])//300)]
                b=comps[j][::max(1,len(comps[j])//300)]
                dy=a[:,None,0]-b[None,:,0]; dx=a[:,None,1]-b[None,:,1]
                d2=dy*dy+dx*dx; flat=int(np.argmin(d2)); ia,ib=np.unravel_index(flat,d2.shape)
                distance=float(np.sqrt(d2[ia,ib]))
                if 1 < distance <= max_gap and (best is None or distance<best[0]):
                    best=(distance,tuple(map(int,a[ia])),tuple(map(int,b[ib])))
        if best:
            path=_least_cost_path(cost,best[1],best[2])
            if path and np.mean([cost[y,x] for y,x in path]) <= 3.2:
                before=np.count_nonzero(completed); yy,xx=zip(*path);completed[yy,xx]=True
                delta=np.count_nonzero(completed)-before
                if delta: paths+=1;added+=delta
    return completed,paths,added


def _wide_boundary_score(
    phase: np.ndarray, component: np.ndarray, width: int
) -> tuple[float, float, int]:
    """Return robust inner/outer phase step using multi-pixel boundary bands."""
    finite = np.isfinite(phase)
    inner = component & ~binary_erosion(component, iterations=width) & finite
    outer = binary_dilation(component, iterations=width) & ~component & finite
    if np.count_nonzero(inner) < 8 or np.count_nonzero(outer) < 8:
        # Thin regions may not have ``width`` complete inner rows. In that case
        # use every available region pixel, as requested for small components.
        inner = component & finite
    n = min(np.count_nonzero(inner), np.count_nonzero(outer))
    if n < 8:
        return np.inf, np.nan, int(n)
    difference = float(np.nanmedian(phase[inner]) - np.nanmedian(phase[outer]))
    return abs(difference), difference, int(n)


def _paired_boundary_measure(
    phase: np.ndarray, component: np.ndarray, width: int
) -> tuple[float, float, int, float]:
    """Measure a boundary step from spatially paired inner/outer samples.

    Pairing nearby samples largely cancels an orbital/atmospheric ramp that can
    bias a single median computed over an entire irregular region.  ``support``
    is the fraction of pairs agreeing with the median integer-cycle hypothesis.
    """
    finite = np.isfinite(phase)
    inner = component & ~binary_erosion(component, iterations=width) & finite
    if np.count_nonzero(inner) < 8:
        inner = component & finite
    outside = ~component & finite
    if np.count_nonzero(inner) < 8 or np.count_nonzero(outside) < 8:
        return np.inf, np.nan, 0, 0.0
    # For every inner-band pixel, locate the nearest finite pixel outside the
    # region. Invalid water pixels are never allowed to become references.
    search = component | ~finite
    distances, indices = distance_transform_edt(
        search, return_distances=True, return_indices=True
    )
    ys, xs = np.nonzero(inner & (distances <= max(2, 2 * width + 1)))
    if ys.size < 8:
        return np.inf, np.nan, int(ys.size), 0.0
    oy, ox = indices[0, ys, xs], indices[1, ys, xs]
    differences = phase[ys, xs] - phase[oy, ox]
    differences = differences[np.isfinite(differences)]
    if differences.size < 8:
        return np.inf, np.nan, int(differences.size), 0.0
    difference = float(np.median(differences))
    integer = int(np.rint(difference / TWOPI))
    residual_cycles = np.abs(differences / TWOPI - integer)
    support = float(np.mean(residual_cycles <= 0.35))
    # The robust median absolute pair difference is the continuity objective.
    score = float(np.median(np.abs(differences)))
    return score, difference, int(differences.size), support


def repair_region_graph_cycles(
    phase: np.ndarray,
    coherence: np.ndarray | None = None,
    branch_cuts: np.ndarray | None = None,
    edge_sigma: float = 1.5,
    edge_threshold_cycles: float = 0.25,
    boundary_width: int = 5,
    min_pixels: int = 40,
    max_fraction: float = 0.45,
    max_iterations: int = 3,
    integer_tolerance_cycles: float = 0.25,
    minimum_improvement: float = 0.30,
    minimum_integer_support: float = 0.45,
    open_boundary_max_gap: int = 12,
    open_boundary_snap: int = 10,
) -> tuple[np.ndarray, np.ndarray, RegionGraphRepairStats]:
    """Detect and repair spatial integer-cycle plateaus without closure triangles.

    A lightly smoothed phase image is converted to a region-adjacency graph by
    treating sharp phase gradients as barriers.  Each resulting region is then
    compared with its neighbours using wide inner/outer boundary bands.  A
    region is shifted only when the observed step is close to an integer number
    of 2-pi cycles and the shift materially improves spatial continuity.

    This deliberately analyses the floating-point phase, not rendered PNG
    colours.  It complements (rather than replaces) temporal closure repair.
    """
    corrected = np.asarray(phase, dtype=float).copy()
    valid = np.isfinite(corrected)
    valid_count = int(np.count_nonzero(valid))
    cycle_map = np.zeros(corrected.shape, dtype=np.int16)
    stats = RegionGraphRepairStats()
    if valid_count < max(3 * min_pixels, 24):
        return corrected, cycle_map, stats

    jump = edge_threshold_cycles * TWOPI
    barrier = np.zeros(corrected.shape, dtype=bool)
    # Scale-normalised gradients allow both sharp and blurred jump boundaries
    # to vote. A true step remains strong after multiplication by sigma, while
    # isolated high-frequency noise is less consistent across these scales.
    scales = sorted({max(0.6, edge_sigma * 0.6), edge_sigma, edge_sigma * 2.0})
    votes = np.zeros(corrected.shape, dtype=np.uint8)
    for sigma in scales:
        smooth = normalized_gaussian(corrected, valid, sigma)
        scale_barrier = np.zeros(corrected.shape, dtype=bool)
        horizontal = valid[:, 1:] & valid[:, :-1] & (
            np.abs(np.diff(smooth, axis=1)) * sigma >= jump
        )
        vertical = valid[1:, :] & valid[:-1, :] & (
            np.abs(np.diff(smooth, axis=0)) * sigma >= jump
        )
        scale_barrier[:, 1:] |= horizontal
        scale_barrier[:, :-1] |= horizontal
        scale_barrier[1:, :] |= vertical
        scale_barrier[:-1, :] |= vertical
        votes += scale_barrier
    # One fine-scale or two coarse-scale votes are sufficient. Closing bridges
    # short gaps but does not invent long boundaries across the scene.
    barrier = (votes >= 1) & valid
    raw_boundary_contact = barrier & binary_dilation(_domain_boundary(valid), iterations=2)
    barrier = binary_closing(barrier, iterations=2) & valid
    # scipy's binary closing treats the outside of the array as background and
    # would otherwise retract a legitimate jump edge from the image boundary.
    barrier |= raw_boundary_contact
    smooth_for_cost=normalized_gaussian(corrected,valid,max(0.8,edge_sigma))
    gy,gx=np.gradient(np.nan_to_num(smooth_for_cost,nan=0.0))
    gradient=np.hypot(gx,gy)
    barrier,stats.boundary_paths_completed,stats.boundary_pixels_added = _complete_boundary_gaps(
        barrier,valid,gradient,coherence,branch_cuts,
        open_boundary_max_gap,open_boundary_snap,
    )
    # One-pixel expansion closes small noisy gaps in an otherwise continuous
    # jump boundary. The removed pixels are assigned back to their nearest core.
    barrier = binary_dilation(barrier, iterations=1) & valid
    core_labels, count = label(valid & ~barrier)
    if count <= 1 or not np.any(core_labels):
        return corrected, cycle_map, stats
    nearest = distance_transform_edt(
        core_labels == 0, return_distances=False, return_indices=True
    )
    regions = core_labels[tuple(nearest)]
    regions[~valid] = 0
    sizes = np.bincount(regions.ravel(), minlength=count + 1)
    eligible = [
        region_id for region_id in range(1, count + 1)
        if min_pixels <= sizes[region_id] <= max_fraction * valid_count
    ]
    stats.regions_detected = len(eligible)

    def region_component(region_id: int) -> np.ndarray:
        """Return a region with jump-edge pixels assigned by phase similarity."""
        component = regions == region_id
        nearby = binary_dilation(component, iterations=2) & barrier & valid
        inner_values = corrected[component]
        outer_mask = binary_dilation(component, iterations=boundary_width) & ~component & valid
        if np.any(nearby) and np.any(outer_mask):
            inside_level = float(np.nanmedian(inner_values))
            outside_level = float(np.nanmedian(corrected[outer_mask]))
            choose_inside = (
                np.abs(corrected - inside_level) <= np.abs(corrected - outside_level)
            )
            component = component | (nearby & choose_inside)
        return component

    for iteration in range(max_iterations):
        changed = 0
        # Largest boundary steps are handled first; later regions are evaluated
        # against the already improved phase rather than a stale reference.
        candidates = []
        for region_id in eligible:
            component = region_component(region_id)
            before, difference, samples, support = _paired_boundary_measure(
                corrected, component, boundary_width
            )
            if not np.isfinite(before):
                continue
            candidates.append((before, region_id, difference, samples, support))
        for before, region_id, _, _, _ in sorted(candidates, reverse=True):
            component = region_component(region_id)
            before, difference, _, support = _paired_boundary_measure(
                corrected, component, boundary_width
            )
            if not np.isfinite(before):
                continue
            stats.regions_tested += 1
            observed_cycles = difference / TWOPI
            nearest_cycle = int(np.rint(observed_cycles))
            if nearest_cycle == 0 or abs(observed_cycles - nearest_cycle) > integer_tolerance_cycles:
                continue
            if support < minimum_integer_support:
                continue
            if coherence is not None:
                local_coherence = np.asarray(coherence, dtype=float)[component]
                if np.any(np.isfinite(local_coherence)) and np.nanmedian(local_coherence) < 0.05:
                    continue
            applied = -nearest_cycle
            trial = corrected.copy()
            trial[component] += applied * TWOPI
            after, _, _, _ = _paired_boundary_measure(trial, component, boundary_width)
            relative_gain = (before - after) / max(before, 1e-9)
            if np.isfinite(after) and relative_gain >= minimum_improvement and after < 0.55 * TWOPI:
                corrected = trial
                cycle_map[component] += applied
                stats.regions_repaired += 1
                stats.pixels_repaired += int(np.count_nonzero(component))
                changed += 1
        stats.iterations = iteration + 1
        if changed == 0:
            break
    return corrected, cycle_map, stats


def repair_integer_cycle_regions(
    phase: np.ndarray,
    sigma: float = 12.0,
    threshold_cycles: float = 0.45,
    min_pixels: int = 40,
    max_fraction: float = 0.40,
) -> tuple[np.ndarray, np.ndarray, SpatialRepairStats]:
    """Repair compact integer-cycle plateaus using low-frequency residuals.

    A candidate is accepted only when adding an integer number of cycles lowers
    the phase discontinuity across its boundary by at least 25%.
    """
    corrected = np.asarray(phase, dtype=float).copy()
    valid = np.isfinite(corrected)
    smooth = normalized_gaussian(corrected, valid, sigma)
    residual_cycles = (corrected - smooth) / TWOPI
    nearest = np.rint(residual_cycles)
    candidates = valid & (np.abs(residual_cycles) >= threshold_cycles) & (nearest != 0)
    cycle_map = np.zeros(phase.shape, dtype=np.int16)
    stats = SpatialRepairStats()
    for cycle in np.unique(nearest[candidates]).astype(int):
        labels, count = label(candidates & (nearest == cycle))
        for component_id in range(1, count + 1):
            component = labels == component_id
            pixels = int(np.count_nonzero(component))
            if pixels < min_pixels or pixels > max_fraction * np.count_nonzero(valid):
                continue
            stats.components_tested += 1
            before = boundary_score(corrected, component)
            # Residual +cycle means the component is high, so subtract it.
            applied = -int(cycle)
            trial = corrected.copy(); trial[component] += applied * TWOPI
            after = boundary_score(trial, component)
            if np.isfinite(before) and after <= 0.75 * before:
                corrected = trial
                cycle_map[component] += applied
                stats.components_repaired += 1
                stats.pixels_repaired += pixels
    return corrected, cycle_map, stats
