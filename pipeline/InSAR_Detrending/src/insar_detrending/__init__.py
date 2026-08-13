"""Robust spatial trend removal for unwrapped interferograms."""

from .detrend import robust_polynomial_detrend
from .pipeline import run_all_datasets, run_dataset

__all__ = ["robust_polynomial_detrend", "run_all_datasets", "run_dataset"]
__version__ = "1.0.0"
