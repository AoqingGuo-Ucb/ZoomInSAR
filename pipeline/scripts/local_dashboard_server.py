#!/usr/bin/env python3
"""Serve a private local dashboard for many ZoomInSAR results.

Run this on the Linux machine holding the results.  The browser receives only a
downsampled velocity grid and the selected point's time series; large NPY stacks
remain on the Linux server.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np


STATIC_ROOT = Path(__file__).resolve().parents[1] / "local_dashboard"


def valid_geo(longitude: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    return (
        np.isfinite(longitude) & np.isfinite(latitude)
        & (longitude >= -180) & (longitude <= 180)
        & (latitude >= -90) & (latitude <= 90)
        & ~((longitude == 0) & (latitude == 0))
    )


def find_single(directory: Path, suffix: str) -> Path:
    matches = sorted(directory.glob(f"*.{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one *.{suffix} file in {directory}")
    return matches[0]


class Project:
    def __init__(self, result_dir: Path, data_root: Path):
        self.result_dir = result_dir
        self.dataset = result_dir.name
        self.stack_path = result_dir / "cumulative_los_displacement_m.npy"
        self.velocity_path = result_dir / "mean_los_velocity_m_per_year.npy"
        self.dates_path = result_dir / "dates.json"
        self.geo_dir = data_root / "Data" / "ROI" / self.dataset / "GEO"
        self.dates = [datetime.fromisoformat(value) for value in json.loads(self.dates_path.read_text())]
        self.shape = tuple(np.load(self.velocity_path, mmap_mode="r").shape)
        self.lon_path = find_single(self.geo_dir, "lon")
        self.lat_path = find_single(self.geo_dir, "lat")
        stack = np.load(self.stack_path, mmap_mode="r")
        if stack.shape != (len(self.dates), *self.shape):
            raise ValueError("Cumulative displacement stack dimensions do not match dates and velocity shape")

    def longitude_latitude(self) -> tuple[np.ndarray, np.ndarray]:
        dtype = np.dtype(">f4")
        return (
            np.memmap(self.lon_path, dtype=dtype, mode="r", shape=self.shape),
            np.memmap(self.lat_path, dtype=dtype, mode="r", shape=self.shape),
        )

    def summary(self, key: str) -> dict:
        longitude, latitude = self.longitude_latitude()
        mask = valid_geo(longitude, latitude)
        if not np.any(mask):
            raise ValueError("No valid GEO longitude/latitude pixels")
        return {
            "id": key,
            "name": self.dataset.removeprefix("Dataset_"),
            "dataset": self.dataset,
            "shape": list(self.shape),
            "dates": [date.strftime("%Y-%m-%d") for date in self.dates],
            "bounds": [
                float(np.nanmin(longitude[mask])), float(np.nanmin(latitude[mask])),
                float(np.nanmax(longitude[mask])), float(np.nanmax(latitude[mask])),
            ],
        }

    def velocity_for_period(self, start: datetime | None, end: datetime | None) -> np.ndarray:
        if ((start is None or start <= self.dates[0])
                and (end is None or end >= self.dates[-1])):
            return np.asarray(np.load(self.velocity_path, mmap_mode="r"), dtype=np.float32)
        selected = [i for i, date in enumerate(self.dates) if (start is None or date >= start) and (end is None or date <= end)]
        if len(selected) < 2:
            raise ValueError("Choose a period containing at least two acquisitions")
        stack = np.load(self.stack_path, mmap_mode="r")
        years = np.array([(self.dates[i] - self.dates[selected[0]]).days / 365.25 for i in selected], dtype=float)
        centered = years - years.mean()
        denom = float(np.sum(centered * centered))
        velocity = np.full(self.shape, np.nan, dtype=np.float32)
        if denom <= 0:
            return velocity

        # Work in row blocks so a California-scale stack is never loaded into
        # memory all at once. A 128-row block is usually only tens of MB.
        for first_row in range(0, self.shape[0], 128):
            last_row = min(first_row + 128, self.shape[0])
            values = np.asarray(stack[selected, first_row:last_row, :], dtype=np.float32)
            valid = np.all(np.isfinite(values), axis=0)
            if not np.any(valid):
                continue
            series = values[:, valid]
            velocity_block = velocity[first_row:last_row, :]
            velocity_block[valid] = (centered @ (series - series.mean(axis=0))) / denom
        return velocity

    def map_payload(self, start: datetime | None, end: datetime | None, maximum: int = 360) -> dict:
        velocity = self.velocity_for_period(start, end) * 1000.0
        longitude, latitude = self.longitude_latitude()
        valid = valid_geo(longitude, latitude) & np.isfinite(velocity)
        if not np.any(valid):
            raise ValueError("No valid velocity and geographic pixels for this period")
        west, south = float(np.nanmin(longitude[valid])), float(np.nanmin(latitude[valid]))
        east, north = float(np.nanmax(longitude[valid])), float(np.nanmax(latitude[valid]))
        aspect = max((east - west) / max(north - south, 1e-9), 0.2)
        width = maximum
        height = max(100, min(maximum, int(round(maximum / aspect))))
        col = np.clip(((longitude[valid] - west) / max(east - west, 1e-12) * (width - 1)).astype(int), 0, width - 1)
        row = np.clip(((north - latitude[valid]) / max(north - south, 1e-12) * (height - 1)).astype(int), 0, height - 1)
        index = row * width + col
        sums = np.bincount(index, weights=velocity[valid], minlength=width * height)
        counts = np.bincount(index, minlength=width * height)
        grid = np.full(width * height, np.nan, dtype=float)
        grid[counts > 0] = sums[counts > 0] / counts[counts > 0]
        finite = grid[np.isfinite(grid)]
        limit = max(float(np.nanpercentile(np.abs(finite), 98)), 1.0)
        return {
            "bounds": [west, south, east, north], "width": width, "height": height,
            "limit_mm_year": limit,
            "values": [None if not math.isfinite(value) else round(float(value), 4) for value in grid],
            "period": [
                start.strftime("%Y-%m-%d") if start else self.dates[0].strftime("%Y-%m-%d"),
                end.strftime("%Y-%m-%d") if end else self.dates[-1].strftime("%Y-%m-%d"),
            ],
        }

    def point_payload(self, longitude_value: float, latitude_value: float) -> dict:
        longitude, latitude = self.longitude_latitude()
        valid = valid_geo(longitude, latitude)
        if not np.any(valid):
            raise ValueError("No valid geographic pixels")
        distance = np.where(valid, (longitude - longitude_value) ** 2 + (latitude - latitude_value) ** 2, np.inf)
        row, col = np.unravel_index(np.argmin(distance), self.shape)
        stack = np.load(self.stack_path, mmap_mode="r")
        series = np.asarray(stack[:, row, col], dtype=float) * 1000.0
        return {
            "longitude": float(longitude[row, col]), "latitude": float(latitude[row, col]),
            "row": int(row), "col": int(col),
            "dates": [date.strftime("%Y-%m-%d") for date in self.dates],
            "displacement_mm": [None if not math.isfinite(value) else round(float(value), 3) for value in series],
        }


def discover_projects(data_root: Path, recursive: bool = False) -> dict[str, Project]:
    projects: dict[str, Project] = {}
    relative_pattern = "ZoomInSAR_Results/InSAR_Timeseries/OUTPUT/Dataset_*/mean_los_velocity_m_per_year.npy"
    if recursive:
        candidates = data_root.rglob(relative_pattern)
    else:
        # The usual layout is <data-root>/<project>/ZoomInSAR_Results/....
        # Targeted glob patterns avoid walking every interferogram and raster.
        candidates = list(data_root.glob(relative_pattern)) + list(data_root.glob(f"*/{relative_pattern}"))
    for velocity in sorted(candidates):
        result_dir = velocity.parent
        try:
            results_root = result_dir.parents[2]
            project = Project(result_dir, results_root)
            key = str(result_dir.relative_to(data_root))
            projects[key] = project
        except (FileNotFoundError, ValueError, json.JSONDecodeError) as error:
            print(f"Skipping {result_dir}: {error}")
    return projects


class DashboardHandler(SimpleHTTPRequestHandler):
    projects: dict[str, Project] = {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(STATIC_ROOT), **kwargs)

    def json_response(self, payload: dict, status: HTTPStatus = HTTPStatus.OK) -> None:
        encoded = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers(); self.wfile.write(encoded)

    def overview_payload(self) -> dict:
        """Return compact, full-time-span velocity maps for every project."""
        maps = []
        for key, project in self.projects.items():
            try:
                payload = project.map_payload(None, None, maximum=220)
            except ValueError as error:
                print(f"Skipping {project.dataset} in overview: {error}")
                continue
            payload["id"] = key
            payload["name"] = project.dataset.removeprefix("Dataset_")
            maps.append(payload)
        if not maps:
            raise ValueError("No projects have a displayable velocity map")
        values = [abs(value) for payload in maps for value in payload["values"] if value is not None]
        limit = max(float(np.percentile(values, 98)), 1.0) if values else 1.0
        for payload in maps:
            payload["limit_mm_year"] = limit
        return {"projects": maps, "limit_mm_year": limit}

    def do_GET(self):  # noqa: N802
        parsed = urlparse(self.path)
        if not parsed.path.startswith("/api/"):
            return super().do_GET()
        try:
            query = parse_qs(parsed.query)
            if parsed.path == "/api/projects":
                return self.json_response({"projects": [project.summary(key) for key, project in self.projects.items()]})
            if parsed.path == "/api/overview":
                return self.json_response(self.overview_payload())
            project_id = query.get("project", [""])[0]
            project = self.projects.get(project_id)
            if project is None:
                raise ValueError("Unknown project")
            if parsed.path == "/api/velocity":
                start = datetime.fromisoformat(query["start"][0]) if query.get("start", [""])[0] else None
                end = datetime.fromisoformat(query["end"][0]) if query.get("end", [""])[0] else None
                return self.json_response(project.map_payload(start, end))
            if parsed.path == "/api/timeseries":
                return self.json_response(project.point_payload(float(query["lon"][0]), float(query["lat"][0])))
            self.json_response({"error": "Unknown API path"}, HTTPStatus.NOT_FOUND)
        except (KeyError, ValueError) as error:
            self.json_response({"error": str(error)}, HTTPStatus.BAD_REQUEST)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the private local ZoomInSAR dashboard")
    parser.add_argument("--data-root", type=Path, required=True,
                        help="Parent folder containing projects and ZoomInSAR_Results directories")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--recursive", action="store_true",
                        help="Search all nested folders (slower; use only for a non-standard layout)")
    args = parser.parse_args()
    data_root = args.data_root.expanduser().resolve()
    print(f"Scanning {data_root} for ZoomInSAR time-series outputs...", flush=True)
    DashboardHandler.projects = discover_projects(data_root, recursive=args.recursive)
    if not DashboardHandler.projects:
        raise SystemExit("No complete ZoomInSAR time-series outputs were found below --data-root")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), DashboardHandler)
    print(f"Found {len(DashboardHandler.projects)} project(s). Open http://127.0.0.1:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nDashboard stopped.")


if __name__ == "__main__":
    main()
