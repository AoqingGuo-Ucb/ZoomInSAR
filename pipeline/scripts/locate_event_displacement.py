#!/usr/bin/env python3
"""Locate a short-lived InSAR event from pre/post cumulative displacement.

The script is designed for an existing ZoomInSAR time-series output.  It uses
the median of all acquisitions in a window before an event and the median of
all acquisitions in a window after it. Their difference is a robust event
change map. The largest coherent high-change region is reported as a suggested
ROI; this is a screening product, not a replacement for validation with field
observations or optical flood mapping.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Locate an event from ZoomInSAR pre/post displacement.")
    parser.add_argument("--project-root", type=Path, required=True,
                        help="Project directory containing ZoomInSAR_Results/.")
    parser.add_argument("--dataset", required=True, help="Dataset name, with or without Dataset_.")
    parser.add_argument("--event-date", required=True, help="Event date in YYYY-MM-DD format.")
    parser.add_argument("--before-days", type=int, default=45,
                        help="Days before the event used for the baseline median (default: 45).")
    parser.add_argument("--after-days", type=int, default=60,
                        help="Days after the event used for the post-event median (default: 60).")
    parser.add_argument("--percentile", type=float, default=97.0,
                        help="Absolute-change percentile defining candidate pixels (default: 97).")
    parser.add_argument("--minimum-pixels", type=int, default=25,
                        help="Discard candidate regions smaller than this many pixels (default: 25).")
    parser.add_argument("--roi-padding-pixels", type=int, default=50,
                        help="Padding around the detected region for the suggested ROI (default: 50).")
    parser.add_argument("--output-dir", type=Path,
                        help="Directory for event files; default is ZoomInSAR_Results/event_YYYYMMDD.")
    return parser.parse_args()


def data_paths(root: Path, dataset: str) -> tuple[Path, Path, Path]:
    name = dataset if dataset.startswith("Dataset_") else f"Dataset_{dataset}"
    results = root.expanduser() / "ZoomInSAR_Results"
    timeseries = results / "InSAR_Timeseries" / "OUTPUT" / name
    data = results / "Data" / "ROI" / name
    return results, timeseries, data


def shape_and_geolocation(data_dir: Path) -> tuple[tuple[int, int], np.ndarray, np.ndarray]:
    metadata = data_dir / "crop_metadata.json"
    if not metadata.is_file():
        raise FileNotFoundError(f"Missing crop metadata: {metadata}")
    payload = json.loads(metadata.read_text(encoding="utf-8"))
    shape = payload.get("output_shape_lines_width") or payload.get("output_shape") or payload.get("shape")
    if not shape or len(shape) != 2:
        raise ValueError(f"No output shape in {metadata}")
    shape = (int(shape[0]), int(shape[1]))
    geo = data_dir / "GEO"
    lon_files, lat_files = sorted(geo.glob("*.lon")), sorted(geo.glob("*.lat"))
    if not lon_files or not lat_files:
        raise FileNotFoundError(f"Missing GEO *.lon or *.lat files in {geo}")
    longitude = np.fromfile(lon_files[0], dtype=">f4").reshape(shape)
    latitude = np.fromfile(lat_files[0], dtype=">f4").reshape(shape)
    return shape, longitude, latitude


def select_dates(dates: list[datetime], start: datetime, end: datetime, label: str) -> list[int]:
    selected = [index for index, date in enumerate(dates) if start <= date <= end]
    if not selected:
        raise ValueError(f"No InSAR acquisitions in the {label} window: {start:%F} to {end:%F}.")
    return selected


def main() -> None:
    args = arguments()
    event = datetime.strptime(args.event_date, "%Y-%m-%d")
    if args.before_days <= 0 or args.after_days <= 0:
        raise ValueError("--before-days and --after-days must be positive.")
    results, timeseries, data_dir = data_paths(args.project_root, args.dataset)
    stack_path, dates_path = timeseries / "cumulative_los_displacement_m.npy", timeseries / "dates.json"
    if not stack_path.is_file() or not dates_path.is_file():
        raise FileNotFoundError(f"Missing time-series files in {timeseries}")
    dates = [datetime.fromisoformat(value) for value in json.loads(dates_path.read_text(encoding="utf-8"))]
    stack = np.load(stack_path, mmap_mode="r")
    shape, longitude, latitude = shape_and_geolocation(data_dir)
    if stack.shape[1:] != shape or stack.shape[0] != len(dates):
        raise ValueError("Time-series stack, dates, and GEO shape are inconsistent.")

    before = select_dates(dates, event - timedelta(days=args.before_days), event - timedelta(days=1), "pre-event")
    after = select_dates(dates, event, event + timedelta(days=args.after_days), "post-event")
    pre = np.nanmedian(np.asarray(stack[before], dtype=np.float32), axis=0)
    post = np.nanmedian(np.asarray(stack[after], dtype=np.float32), axis=0)
    change_mm = (post - pre) * 1000.0
    valid = np.isfinite(change_mm) & np.isfinite(longitude) & np.isfinite(latitude) & (longitude != 0) & (latitude != 0)
    if np.count_nonzero(valid) < args.minimum_pixels:
        raise ValueError("Too few valid pixels to locate an event.")
    threshold = float(np.nanpercentile(np.abs(change_mm[valid]), args.percentile))
    candidates = valid & (np.abs(change_mm) >= threshold)
    labels, count = ndimage.label(candidates)
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    label = int(np.argmax(sizes)) if count else 0
    if label == 0 or sizes[label] < args.minimum_pixels:
        raise ValueError(
            f"No coherent event region contains at least {args.minimum_pixels} pixels. "
            "Try lowering --percentile or --minimum-pixels."
        )
    region = labels == label
    rows, cols = np.where(region)
    pad = int(args.roi_padding_pixels)
    r0, r1 = max(0, rows.min() - pad), min(shape[0] - 1, rows.max() + pad)
    c0, c1 = max(0, cols.min() - pad), min(shape[1] - 1, cols.max() + pad)
    weights = np.abs(change_mm[region])
    center_lon = float(np.average(longitude[region], weights=weights))
    center_lat = float(np.average(latitude[region], weights=weights))
    roi_valid = valid[r0:r1 + 1, c0:c1 + 1]
    roi_lon, roi_lat = longitude[r0:r1 + 1, c0:c1 + 1], latitude[r0:r1 + 1, c0:c1 + 1]
    west, east = float(np.nanmin(roi_lon[roi_valid])), float(np.nanmax(roi_lon[roi_valid]))
    south, north = float(np.nanmin(roi_lat[roi_valid])), float(np.nanmax(roi_lat[roi_valid]))

    out = args.output_dir or results / f"event_{event:%Y%m%d}"
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "event_change_mm.npy", change_mm.astype(np.float32))
    report = {
        "event_date": args.event_date,
        "pre_event_dates": [dates[i].strftime("%F") for i in before],
        "post_event_dates": [dates[i].strftime("%F") for i in after],
        "threshold_mm": threshold,
        "event_pixels": int(sizes[label]),
        "center_longitude": center_lon,
        "center_latitude": center_lat,
        "suggested_roi_bbox_west_south_east_north": [west, south, east, north],
        "suggested_roi_pixel_rows_cols": [int(r0), int(r1), int(c0), int(c1)],
    }
    (out / "event_location.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (out / "suggested_roi.wkt").write_text(
        f"POLYGON(({west} {south},{east} {south},{east} {north},{west} {north},{west} {south}))\n",
        encoding="utf-8",
    )
    coordinates = f"{west},{south},0 {east},{south},0 {east},{north},0 {west},{north},0 {west},{south},0"
    (out / "suggested_roi.kml").write_text(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<kml xmlns=\"http://www.opengis.net/kml/2.2\"><Document><Placemark>"
        "<name>Suggested event ROI</name><Polygon><outerBoundaryIs><LinearRing>"
        f"<coordinates>{coordinates}</coordinates>"
        "</LinearRing></outerBoundaryIs></Polygon></Placemark></Document></kml>\n",
        encoding="utf-8",
    )
    limit = max(float(np.nanpercentile(np.abs(change_mm[valid]), 98)), 1.0)
    fig, axis = plt.subplots(figsize=(9, 7))
    image = axis.pcolormesh(longitude, latitude, change_mm, shading="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    axis.contour(longitude, latitude, region.astype(np.uint8), levels=[0.5], colors="black", linewidths=1.0)
    axis.plot(center_lon, center_lat, "y*", markersize=12, markeredgecolor="black", label="Candidate event centre")
    axis.plot([west, east, east, west, west], [south, south, north, north, south], "k--", label="Suggested ROI")
    axis.set(xlabel="Longitude", ylabel="Latitude", title=f"LOS displacement change around {event:%Y-%m-%d}")
    axis.legend(loc="best")
    fig.colorbar(image, ax=axis, label="Post-event minus pre-event LOS displacement (mm)")
    fig.tight_layout()
    fig.savefig(out / "event_change_map.png", dpi=300)
    print(json.dumps(report, indent=2))
    print(f"Wrote {out / 'event_change_map.png'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        raise SystemExit(f"ERROR: {error}") from error
