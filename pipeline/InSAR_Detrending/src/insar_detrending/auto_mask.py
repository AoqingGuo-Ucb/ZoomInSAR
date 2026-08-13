"""Build deformation exclusion masks from preliminary velocity maps."""

from __future__ import annotations

import json
from collections import deque
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


@dataclass
class AutoMaskConfig:
    velocity_percentile: float = 85.0
    min_velocity_m_per_year: float = 0.02
    dilation_pixels: int = 6
    min_component_pixels: int = 80
    max_kml_polygons: int = 20
    min_kml_polygon_area_pixels: float = 100.0


def _load_lon_lat(filtering_dataset: Path, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    geo = filtering_dataset / "GEO"
    lon_files, lat_files = sorted(geo.glob("*.lon")), sorted(geo.glob("*.lat"))
    if len(lon_files) != 1 or len(lat_files) != 1:
        raise ValueError(f"{geo} must contain exactly one *.lon and one *.lat")
    dtype = np.dtype("f4").newbyteorder(">")
    longitude = np.asarray(np.memmap(lon_files[0], dtype=dtype, mode="r", shape=shape), float)
    latitude = np.asarray(np.memmap(lat_files[0], dtype=dtype, mode="r", shape=shape), float)
    return longitude, latitude


def _binary_dilate(mask: np.ndarray, radius: int) -> np.ndarray:
    radius = int(radius)
    if radius <= 0:
        return mask.copy()
    padded = np.pad(mask.astype(bool), radius, mode="constant", constant_values=False)
    out = np.zeros_like(mask, dtype=bool)
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            if dr * dr + dc * dc > radius * radius:
                continue
            out |= padded[
                radius + dr:radius + dr + mask.shape[0],
                radius + dc:radius + dc + mask.shape[1],
            ]
    return out


def _remove_small_components(mask: np.ndarray, min_pixels: int) -> tuple[np.ndarray, list[int]]:
    mask = np.asarray(mask, dtype=bool)
    visited = np.zeros(mask.shape, dtype=bool)
    output = np.zeros(mask.shape, dtype=bool)
    sizes: list[int] = []
    rows, cols = mask.shape
    for start_row, start_col in zip(*np.nonzero(mask & ~visited)):
        if visited[start_row, start_col]:
            continue
        queue: deque[tuple[int, int]] = deque([(int(start_row), int(start_col))])
        visited[start_row, start_col] = True
        pixels: list[tuple[int, int]] = []
        while queue:
            row, col = queue.popleft()
            pixels.append((row, col))
            for rr in range(max(0, row - 1), min(rows, row + 2)):
                for cc in range(max(0, col - 1), min(cols, col + 2)):
                    if visited[rr, cc] or not mask[rr, cc]:
                        continue
                    visited[rr, cc] = True
                    queue.append((rr, cc))
        sizes.append(len(pixels))
        if len(pixels) >= min_pixels:
            rr, cc = zip(*pixels)
            output[np.asarray(rr), np.asarray(cc)] = True
    return output, sizes


def _polygon_area_pixels(vertices: np.ndarray) -> float:
    x = vertices[:, 0]
    y = vertices[:, 1]
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) * 0.5)


def _mask_to_kml(mask: np.ndarray, longitude: np.ndarray, latitude: np.ndarray,
                 path: Path, name: str, config: AutoMaskConfig) -> int:
    fig, axis = plt.subplots(figsize=(4, 4))
    contours = axis.contour(mask.astype(float), levels=[0.5])
    candidates: list[tuple[float, np.ndarray]] = []
    for level_segments in contours.allsegs:
        for vertices in level_segments:
            if vertices.shape[0] < 4:
                continue
            area_pixels = _polygon_area_pixels(vertices)
            if area_pixels < config.min_kml_polygon_area_pixels:
                continue
            candidates.append((area_pixels, vertices))
    candidates.sort(key=lambda item: item[0], reverse=True)
    placemarks = []
    for area_pixels, vertices in candidates[:config.max_kml_polygons]:
        cols = np.clip(np.rint(vertices[:, 0]).astype(int), 0, mask.shape[1] - 1)
        rows = np.clip(np.rint(vertices[:, 1]).astype(int), 0, mask.shape[0] - 1)
        lon = longitude[rows, cols]
        lat = latitude[rows, cols]
        valid = np.isfinite(lon) & np.isfinite(lat)
        if np.count_nonzero(valid) < 4:
            continue
        coordinates = "\n".join(
            f"{lo:.8f},{la:.8f},0" for lo, la in zip(lon[valid], lat[valid])
        )
        first = f"{lon[valid][0]:.8f},{lat[valid][0]:.8f},0"
        if not coordinates.endswith(first):
            coordinates = f"{coordinates}\n{first}"
        placemarks.append(
            "    <Placemark>\n"
            f"      <name>{name} deformation exclude mask; area {area_pixels:.0f} pixels</name>\n"
            "      <Polygon><outerBoundaryIs><LinearRing><coordinates>\n"
            f"{coordinates}\n"
            "      </coordinates></LinearRing></outerBoundaryIs></Polygon>\n"
            "    </Placemark>"
        )
    plt.close(fig)
    kml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
        "  <Document>\n"
        f"    <name>{name} auto deformation exclude mask</name>\n"
        + "\n".join(placemarks)
        + "\n  </Document>\n</kml>\n"
    )
    path.write_text(kml, encoding="utf-8")
    return len(placemarks)


def _plot_mask(path: Path, velocity: np.ndarray, mask: np.ndarray,
               threshold: float, reference_title: str) -> None:
    finite = np.isfinite(velocity)
    vmax = max(float(np.nanpercentile(np.abs(velocity[finite]), 98)), threshold, 1e-6)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    image = axes[0].imshow(velocity, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="upper")
    axes[0].set_title("Preliminary mean LOS velocity")
    axes[0].set_xlabel("Range column")
    axes[0].set_ylabel("Azimuth row")
    fig.colorbar(image, ax=axes[0], shrink=0.82, label="m/year")
    axes[1].imshow(mask, cmap="gray", vmin=0, vmax=1, origin="upper", interpolation="nearest")
    axes[1].set_title(f"Auto deformation exclude mask\n|v| >= {threshold:.3f} m/year")
    axes[1].set_xlabel("Range column")
    axes[1].set_ylabel("Azimuth row")
    fig.suptitle(reference_title)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_mask_from_velocity(timeseries_dataset: str | Path, filtering_dataset: str | Path,
                             output: str | Path, config: AutoMaskConfig | None = None) -> dict:
    timeseries_dataset = Path(timeseries_dataset)
    filtering_dataset = Path(filtering_dataset)
    output = Path(output)
    config = config or AutoMaskConfig()
    velocity_path = timeseries_dataset / "mean_los_velocity_m_per_year.npy"
    if not velocity_path.exists():
        raise FileNotFoundError(f"Missing preliminary velocity map: {velocity_path}")
    velocity = np.load(velocity_path).astype(float)
    finite = np.isfinite(velocity)
    if not np.any(finite):
        raise ValueError(f"No finite velocity pixels in {velocity_path}")
    threshold = max(
        float(np.nanpercentile(np.abs(velocity[finite]), config.velocity_percentile)),
        float(config.min_velocity_m_per_year),
    )
    raw_mask = finite & (np.abs(velocity) >= threshold)
    cleaned, component_sizes = _remove_small_components(raw_mask, config.min_component_pixels)
    mask = _binary_dilate(cleaned, config.dilation_pixels) & finite
    output.mkdir(parents=True, exist_ok=True)
    np.save(output / "auto_deformation_exclude_mask.npy", mask)
    np.save(output / "auto_deformation_exclude_mask_raw.npy", raw_mask)
    longitude, latitude = _load_lon_lat(filtering_dataset, velocity.shape)
    polygon_count = _mask_to_kml(
        mask, longitude, latitude, output / "auto_deformation_exclude_mask.kml",
        timeseries_dataset.name, config,
    )
    _plot_mask(
        output / "auto_deformation_exclude_mask.png",
        velocity, mask, threshold, timeseries_dataset.name,
    )
    summary = {
        "timeseries_dataset": str(timeseries_dataset),
        "filtering_dataset": str(filtering_dataset),
        "output_directory": str(output),
        "config": asdict(config),
        "velocity_threshold_m_per_year": threshold,
        "raw_mask_pixels": int(np.count_nonzero(raw_mask)),
        "final_mask_pixels": int(np.count_nonzero(mask)),
        "finite_velocity_pixels": int(np.count_nonzero(finite)),
        "component_sizes_pixels": component_sizes,
        "kml_polygons": int(polygon_count),
        "mask_path": str(output / "auto_deformation_exclude_mask.npy"),
        "kml_path": str(output / "auto_deformation_exclude_mask.kml"),
    }
    (output / "auto_deformation_exclude_mask_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary
