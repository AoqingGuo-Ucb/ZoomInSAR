#!/usr/bin/env python3
"""Create reviewable landfill profiles from the public CalRecycle SWIS service.

This script matches each dashboard project to the nearest SWIS point.  It is a
coordinate-based preliminary match, not a regulatory determination.  Entries
farther than the configured review distance are explicitly flagged for review.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


SWIS_QUERY = "https://gis.cnra.ca.gov/arcgis/rest/services/Ocean/CSMW_Solid_Waste_Facilities/MapServer/0/query"
FIELDS = "SiteName,SwisNo,County,Operator,Owner,EnforAgent,Category,Activity,RegStatus,OpStatus,Latitude,Longitude"


def fetch_json(url: str) -> dict:
    with urlopen(url, timeout=90) as response:  # nosec B310: fixed public endpoint
        return json.load(response)


def fetch_swis() -> list[dict]:
    features: list[dict] = []
    offset = 0
    while True:
        query = urlencode({
            "where": "1=1", "outFields": FIELDS, "returnGeometry": "false",
            "resultOffset": offset, "resultRecordCount": 2000, "f": "json",
        })
        page = fetch_json(f"{SWIS_QUERY}?{query}")
        if "error" in page:
            raise RuntimeError(page["error"])
        batch = [item["attributes"] for item in page.get("features", [])]
        features.extend(batch)
        if not page.get("exceededTransferLimit"):
            return features
        offset += len(batch)


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    radians = math.pi / 180
    a = (math.sin((lat2 - lat1) * radians / 2) ** 2
         + math.cos(lat1 * radians) * math.cos(lat2 * radians)
         * math.sin((lon2 - lon1) * radians / 2) ** 2)
    return 2 * 6371.0088 * math.asin(math.sqrt(a))


def confidence(distance_km: float, review_distance_km: float, in_california: bool) -> str:
    if not in_california:
        return "outside_california"
    if distance_km <= 1:
        return "high"
    if distance_km <= review_distance_km:
        return "medium"
    return "review_required"


def project_records(payload: dict) -> list[dict]:
    records = []
    for project in payload["projects"]:
        west, south, east, north = project["bounds"]
        records.append({
            "project_id": project["id"], "project_name": project["name"],
            "center_longitude": (west + east) / 2, "center_latitude": (south + north) / 2,
            "bounds": project["bounds"],
        })
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description="Build reviewable CalRecycle SWIS dashboard profiles")
    parser.add_argument("--projects-url", default="http://127.0.0.1:8765/api/projects")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory, normally <data-root>/landfill_profiles")
    parser.add_argument("--review-distance-km", type=float, default=2.0)
    args = parser.parse_args()

    projects = project_records(fetch_json(args.projects_url))
    facilities = [item for item in fetch_swis() if item.get("Latitude") is not None and item.get("Longitude") is not None]
    profiles = []
    for project in projects:
        in_california = -125 < project["center_longitude"] < -114 and 32 < project["center_latitude"] < 43
        nearest = min(facilities, key=lambda site: haversine_km(
            project["center_longitude"], project["center_latitude"], site["Longitude"], site["Latitude"]
        ))
        distance = haversine_km(project["center_longitude"], project["center_latitude"], nearest["Longitude"], nearest["Latitude"])
        profiles.append({
            **project,
            "match": {"method": "nearest_swis_coordinate", "distance_km": round(distance, 3),
                      "confidence": confidence(distance, args.review_distance_km, in_california)},
            "facility": {
                "name": nearest.get("SiteName"), "swis_number": nearest.get("SwisNo"),
                "county": nearest.get("County"), "operator": nearest.get("Operator"),
                "owner": nearest.get("Owner"), "enforcement_agency": nearest.get("EnforAgent"),
                "category": nearest.get("Category"), "activity": nearest.get("Activity"),
                "regulatory_status": nearest.get("RegStatus"), "operational_status": nearest.get("OpStatus"),
            },
            "deformation_notes": [],
            "sources": [{"title": "CalRecycle Solid Waste Information System (SWIS)",
                         "url": "https://www2.calrecycle.ca.gov/solidwaste/site/search",
                         "accessed": str(date.today())}],
        })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"generated": str(date.today()), "source": "CalRecycle SWIS public feature service",
                "review_distance_km": args.review_distance_km, "profiles": profiles}
    (args.output_dir / "index.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for profile in profiles:
        safe_name = profile["project_name"].replace("/", "_").replace(" ", "_")
        (args.output_dir / f"{safe_name}.json").write_text(json.dumps(profile, indent=2) + "\n")
    counts = {level: sum(item["match"]["confidence"] == level for item in profiles)
              for level in ("high", "medium", "review_required", "outside_california")}
    print(f"Wrote {len(profiles)} profiles to {args.output_dir}")
    print("Match confidence:", counts)


if __name__ == "__main__":
    main()
