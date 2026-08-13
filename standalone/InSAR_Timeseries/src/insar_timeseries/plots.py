"""Time-series network and raster figures."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .network import NetworkEdge


def _limits(values: np.ndarray, symmetric: bool = True) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return -1.0, 1.0
    if symmetric:
        maximum = max(float(np.nanpercentile(np.abs(finite), 98)), 1e-9)
        return -maximum, maximum
    return tuple(map(float, np.nanpercentile(finite, [2, 98])))


def plot_network(path: Path, dates: list[datetime], edges: list[NetworkEdge], threshold: int) -> None:
    fig, axis = plt.subplots(figsize=(12, 5), constrained_layout=True)
    index = {date: i for i, date in enumerate(dates)}
    for edge in edges:
        if not edge.selected:
            continue
        x = [edge.first, edge.second]
        y = [index[edge.first], index[edge.second]]
        bridge = "connectivity_bridge" in edge.reason
        axis.plot(x, y, color="crimson" if bridge else "steelblue", linewidth=2.5 if bridge else 1.5,
                  label="Connectivity bridge" if bridge else "Within threshold")
    # Deduplicate legend entries.
    handles, labels = axis.get_legend_handles_labels(); unique = dict(zip(labels, handles))
    if unique:
        axis.legend(unique.values(), unique.keys())
    axis.scatter(dates, range(len(dates)), color="black", zorder=3)
    axis.set_yticks(range(len(dates)), [date.strftime("%Y-%m-%d") for date in dates])
    axis.set_title(f"Selected interferogram network (maximum preferred baseline: {threshold} days)")
    axis.set_xlabel("Acquisition date"); axis.set_ylabel("Acquisition index"); axis.grid(alpha=0.25)
    fig.savefig(path, dpi=180); plt.close(fig)


def _mark_points(axes, points: list[dict]) -> None:
    for axis in np.atleast_1d(axes).flat:
        for point in points:
            axis.plot(point["col"], point["row"], marker="*", markersize=13,
                      markerfacecolor="yellow", markeredgecolor="black", zorder=10)
            axis.annotate(point["id"], (point["col"], point["row"]), xytext=(5, 5),
                          textcoords="offset points", color="black", weight="bold", zorder=11)


def plot_velocity_and_residual(path: Path, velocity: np.ndarray, residual: np.ndarray,
                               points: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    vmin, vmax = _limits(velocity)
    first = axes[0].imshow(velocity * 1000.0, cmap="RdBu_r", vmin=vmin * 1000, vmax=vmax * 1000,
                           origin="upper", interpolation="nearest")
    axes[0].set_title("Mean LOS velocity (mm/year)"); fig.colorbar(first, ax=axes[0], shrink=0.82)
    _, residual_max = _limits(residual, symmetric=False)
    second = axes[1].imshow(residual, cmap="magma", vmin=0, vmax=max(residual_max, 1e-6),
                            origin="upper", interpolation="nearest")
    axes[1].set_title("Phase residual RMS (rad)"); fig.colorbar(second, ax=axes[1], shrink=0.82)
    for axis in axes:
        axis.set_xlabel("Range column"); axis.set_ylabel("Azimuth row")
    _mark_points(axes, points)
    fig.savefig(path, dpi=180); plt.close(fig)


def plot_displacement_maps(path: Path, displacement: np.ndarray, dates: list[datetime],
                           points: list[dict]) -> None:
    columns = 3
    rows = int(np.ceil(len(dates) / columns))
    vmin, vmax = _limits(displacement)
    fig, axes = plt.subplots(rows, columns, figsize=(15, 4.5 * rows), constrained_layout=True, squeeze=False)
    image = None
    for index, axis in enumerate(axes.flat):
        if index >= len(dates):
            axis.set_axis_off(); continue
        image = axis.imshow(displacement[index] * 1000.0, cmap="RdBu_r", vmin=vmin * 1000,
                            vmax=vmax * 1000, origin="upper", interpolation="nearest")
        axis.set_title(dates[index].strftime("%Y-%m-%d")); axis.set_xlabel("Range column"); axis.set_ylabel("Azimuth row")
        _mark_points([axis], points)
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.7, label="Cumulative LOS displacement (mm)")
    fig.savefig(path, dpi=180); plt.close(fig)


def _plot_danger_bridges(axis, dates, series: np.ndarray,
                         danger_bridges: list[tuple[datetime, datetime]]) -> None:
    """Overlay low-quality connectivity-bridge intervals in red."""
    date_index = {date: index for index, date in enumerate(dates)}
    label_used = False
    for first, second in danger_bridges:
        if first not in date_index or second not in date_index:
            continue
        i, j = date_index[first], date_index[second]
        axis.plot(
            [first, second], [series[i], series[j]], color="red", linewidth=3.0,
            marker="o", markersize=5, zorder=8,
            label="Low-quality connectivity bridge (danger)" if not label_used else None,
        )
        axis.axvspan(first, second, color="red", alpha=0.08, zorder=0)
        label_used = True


def plot_representative_timeseries(path: Path, displacement: np.ndarray,
                                   dates: list[datetime], points: list[dict],
                                   danger_bridges: list[tuple[datetime, datetime]] | None = None) -> None:
    if not points:
        return
    fig, axis = plt.subplots(figsize=(9, 4), constrained_layout=True)
    for point in points:
        series = displacement[:, point["row"], point["col"]] * 1000.0
        axis.plot(dates, series, "o-", label=point["id"])
        _plot_danger_bridges(axis, dates, series, danger_bridges or [])
    axis.set_title("Representative time-series points marked on exported maps")
    axis.legend()
    axis.set_ylabel("Cumulative LOS displacement (mm)"); axis.set_xlabel("Date"); axis.grid(alpha=0.3)
    fig.savefig(path, dpi=180); plt.close(fig)


def select_representative_points(velocity: np.ndarray, longitude: np.ndarray,
                                 latitude: np.ndarray) -> list[dict]:
    valid = np.isfinite(velocity) & np.isfinite(longitude) & np.isfinite(latitude)
    if not np.any(valid):
        return []
    masked = np.where(valid, np.abs(velocity), -np.inf)
    row, col = np.unravel_index(int(np.argmax(masked)), velocity.shape)
    return [{
        "id": "P1", "row": int(row), "col": int(col),
        "longitude": float(longitude[row, col]), "latitude": float(latitude[row, col]),
        "velocity_m_per_year": float(velocity[row, col]),
        "velocity_mm_per_year": float(velocity[row, col] * 1000.0),
        "selection": "largest_absolute_velocity",
    }]
