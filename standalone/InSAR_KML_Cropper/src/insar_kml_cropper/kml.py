"""Small, dependency-free KML coordinate parser."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path


def kml_coordinate_rings(path: Path) -> list[list[tuple[float, float]]]:
    """Return all coordinate sequences as (longitude, latitude) rings/paths."""
    root = ET.parse(path).getroot()
    rings: list[list[tuple[float, float]]] = []
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] != "coordinates" or not element.text:
            continue
        points: list[tuple[float, float]] = []
        for token in element.text.split():
            fields = token.split(",")
            if len(fields) >= 2:
                points.append((float(fields[0]), float(fields[1])))
        if points:
            rings.append(points)
    return rings


def kml_bounds(path: Path) -> tuple[float, float, float, float]:
    """Return (west, east, south, north) across all KML coordinate elements."""
    rings = kml_coordinate_rings(path)
    points = [point for ring in rings for point in ring]
    if not points:
        raise ValueError(f"No coordinates found in {path}")
    longitude, latitude = zip(*points)
    return min(longitude), max(longitude), min(latitude), max(latitude)


def expand_bounds(
    bounds: tuple[float, float, float, float], fraction: float
) -> tuple[float, float, float, float]:
    """Expand every side by fraction of the original width/height."""
    west, east, south, north = bounds
    dx, dy = (east - west) * fraction, (north - south) * fraction
    return west - dx, east + dx, south - dy, north + dy
