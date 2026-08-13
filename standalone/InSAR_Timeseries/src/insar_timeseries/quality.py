"""Propagate unwrapping and detrending quality into time-series inversion."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class EdgeQuality:
    tier: int = 0
    weight: float = 1.0
    qc_status: str = "quality_metadata_unavailable"
    detrending_accepted: bool | None = None


def _csv_by_pair(path: Path, key: str) -> dict[str, dict]:
    if not path.exists():
        return {}
    result = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            pair = row.get(key, "").replace("_", "-")
            if pair:
                result[pair] = row
    return result


def load_edge_quality(dataset: Path) -> dict[str, EdgeQuality]:
    """Read QC metadata for either Detrending or Unwrapping Dataset_* input."""
    summary_path = dataset / "run_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    source = summary.get("source_unwrapping_summary")
    unwrapping_dataset = Path(source["output_directory"]) if source else dataset
    qc = _csv_by_pair(unwrapping_dataset / "interferogram_exclusion_report.csv", "pair")
    detrending = _csv_by_pair(dataset / "detrending_quality.csv", "interferogram")
    output: dict[str, EdgeQuality] = {}
    for pair in set(qc) | set(detrending):
        qrow = qc.get(pair, {})
        status = qrow.get("status", "quality_metadata_unavailable")
        protected = status == "network_protected"
        tier = 2 if protected else (3 if status.startswith("excluded") else 0)
        fraction = float(qrow.get("spatial_residual_fraction") or 0.0)
        phase_range = float(qrow.get("robust_phase_range_cycles") or 0.0)
        weight = math.exp(-8.0 * fraction - 0.4 * max(phase_range - 1.0, 0.0))
        accepted: bool | None = None
        if pair in detrending:
            accepted = detrending[pair].get("accepted", "").lower() == "true"
            if not accepted:
                tier = max(tier, 3)
                weight *= 0.20
        if protected:
            weight *= 0.25
        output[pair] = EdgeQuality(tier, max(0.03, min(1.0, weight)), status, accepted)
    return output
