"""Interactive mean-velocity point picker and time-series plotter."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _load_danger_bridges(dataset_output: Path) -> list[tuple[np.datetime64, np.datetime64]]:
    """Load selected low-quality network bridges written by the inversion."""
    network_path = dataset_output / "interferogram_network.csv"
    if not network_path.is_file():
        return []
    bridges = []
    with network_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if (row.get("selected", "").lower() == "true"
                    and row.get("reason") == "low_quality_connectivity_bridge"):
                bridges.append((np.datetime64(row["first"]), np.datetime64(row["second"])))
    return bridges


def _dataset_output(output_root: str | Path, dataset: str) -> Path:
    dataset_path = Path(dataset)
    path = dataset_path if dataset_path.is_dir() else Path(output_root) / dataset
    required = (
        "mean_los_velocity_m_per_year.npy",
        "cumulative_los_displacement_m.npy",
        "dates.json",
    )
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise FileNotFoundError(f"{path} is missing: {', '.join(missing)}")
    return path


def _nearest_valid_pixel(velocity: np.ndarray, x: float, y: float) -> tuple[int, int]:
    row = int(np.clip(round(y), 0, velocity.shape[0] - 1))
    col = int(np.clip(round(x), 0, velocity.shape[1] - 1))
    if np.isfinite(velocity[row, col]):
        return row, col
    rows, cols = np.where(np.isfinite(velocity))
    if rows.size == 0:
        raise ValueError("The mean-velocity raster contains no valid pixels")
    nearest = int(np.argmin((rows - row) ** 2 + (cols - col) ** 2))
    return int(rows[nearest]), int(cols[nearest])


def _load_optional_geolocation(dataset_output: Path, shape: tuple[int, int]):
    summary_path = dataset_output / "timeseries_summary.json"
    if not summary_path.is_file():
        return None, None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    input_dir = Path(summary.get("input_directory", ""))
    try:
        from .geotiff import load_geolocation
        return load_geolocation(input_dir, shape)
    except (FileNotFoundError, ValueError, KeyError, json.JSONDecodeError):
        pass
    return None, None


def interactive_point_timeseries(
    dataset: str,
    output_root: str | Path = "OUTPUT",
    save_name: str = "interactive_point_timeseries.png",
) -> dict:
    """Click a mean-velocity pixel, plot its cumulative LOS time series, and save it.

    The function blocks until the user left-clicks a valid location. Closing the
    window without clicking cancels the selection.
    """
    dataset_output = _dataset_output(output_root, dataset)
    try:
        plt.switch_backend("TkAgg")
    except ImportError as error:
        raise RuntimeError(
            "No interactive Matplotlib GUI backend is available; install or enable Tkinter."
        ) from error
    velocity = np.load(dataset_output / "mean_los_velocity_m_per_year.npy")
    displacement = np.load(dataset_output / "cumulative_los_displacement_m.npy", mmap_mode="r")
    dates = np.asarray(json.loads((dataset_output / "dates.json").read_text(encoding="utf-8")),
                       dtype="datetime64[D]")
    if displacement.shape[0] != dates.size or displacement.shape[1:] != velocity.shape:
        raise ValueError("Velocity, displacement, and date dimensions are inconsistent")
    longitude, latitude = _load_optional_geolocation(dataset_output, velocity.shape)
    danger_bridges = _load_danger_bridges(dataset_output)

    finite = velocity[np.isfinite(velocity)] * 1000.0
    limit = max(float(np.nanpercentile(np.abs(finite), 98)), 1e-6)
    fig, (map_axis, series_axis) = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    image = map_axis.imshow(velocity * 1000.0, origin="upper", interpolation="nearest",
                            cmap="RdBu_r", vmin=-limit, vmax=limit)
    fig.colorbar(image, ax=map_axis, label="Mean LOS velocity (mm/year)")
    map_axis.set_title(f"{dataset_output.name}: click one valid pixel")
    map_axis.set_xlabel("Range column")
    map_axis.set_ylabel("Azimuth row")
    series_axis.set_title("Cumulative LOS displacement")
    series_axis.set_xlabel("Date")
    series_axis.set_ylabel("Displacement (mm)")
    series_axis.grid(alpha=0.3)

    selected: dict = {}
    marker = None

    def on_click(event) -> None:
        nonlocal marker
        if event.inaxes is not map_axis or event.xdata is None or event.ydata is None:
            return
        row, col = _nearest_valid_pixel(velocity, event.xdata, event.ydata)
        if marker is not None:
            marker.remove()
        marker, = map_axis.plot(col, row, marker="*", markersize=16,
                                markerfacecolor="yellow", markeredgecolor="black", zorder=5)
        series_axis.clear()
        values = np.asarray(displacement[:, row, col], float) * 1000.0
        series_axis.plot(dates, values, "o-", color="tab:blue", linewidth=1.6, markersize=4)
        date_index = {date: index for index, date in enumerate(dates)}
        label_used = False
        for first, second in danger_bridges:
            if first not in date_index or second not in date_index:
                continue
            i, j = date_index[first], date_index[second]
            series_axis.plot(
                [first, second], [values[i], values[j]], color="red", linewidth=3.0,
                marker="o", markersize=5, zorder=8,
                label="Low-quality connectivity bridge (danger)" if not label_used else None,
            )
            series_axis.axvspan(first, second, color="red", alpha=0.08, zorder=0)
            label_used = True
        series_axis.set_title(f"row={row}, col={col}; velocity={velocity[row, col] * 1000:.2f} mm/year")
        series_axis.set_xlabel("Date")
        series_axis.set_ylabel("Cumulative LOS displacement (mm)")
        series_axis.grid(alpha=0.3)
        if label_used:
            series_axis.legend()
        selected.clear()
        selected.update({
            "dataset": dataset_output.name,
            "row": row,
            "col": col,
            "longitude": None if longitude is None else float(longitude[row, col]),
            "latitude": None if latitude is None else float(latitude[row, col]),
            "mean_velocity_m_per_year": float(velocity[row, col]),
            "mean_velocity_mm_per_year": float(velocity[row, col] * 1000.0),
        })
        fig.canvas.draw_idle()

    connection = fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show(block=True)
    fig.canvas.mpl_disconnect(connection)
    if not selected:
        plt.close(fig)
        raise RuntimeError("No point was selected")

    figure_path = dataset_output / "figures" / save_name
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    selected["figure"] = str(figure_path)
    (dataset_output / "interactive_point.json").write_text(
        json.dumps(selected, indent=2), encoding="utf-8"
    )
    with (dataset_output / "interactive_point.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=selected.keys())
        writer.writeheader()
        writer.writerow(selected)
    print(json.dumps(selected, indent=2))
    return selected
