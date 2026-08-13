import numpy as np

from insar_detrending.detrend import robust_polynomial_detrend


def test_removes_quadratic_trend_while_retaining_local_signal():
    y, x = np.mgrid[-1:1:80j, -1:1:90j]
    trend = 1.2 + 0.8*x - 0.5*y + 0.35*x*x + 0.2*x*y
    blob = 3.0 * np.exp(-((x-0.25)**2 + (y+0.15)**2) / 0.015)
    image = trend + blob
    corrected, surface, mask, _ = robust_polynomial_detrend(image, degree=2, mad_scale=2.5)
    background = blob < 0.05
    assert np.sqrt(np.mean((surface[background] - trend[background])**2)) < 0.03
    assert np.nanmax(corrected) > 2.8
    assert np.count_nonzero(mask) < image.size


def test_preserves_nan_mask():
    image = np.arange(400, dtype=float).reshape(20, 20)
    image[:3] = np.nan
    corrected, surface, _, _ = robust_polynomial_detrend(image, degree=1)
    assert np.all(np.isnan(corrected[:3])) and np.all(np.isnan(surface[:3]))


def test_jointly_removes_orbit_and_elevation_correlated_terms():
    y, x = np.mgrid[-1:1:80j, -1:1:90j]
    elevation = 500 + 180 * np.sin(2.5 * x) + 120 * np.cos(2.0 * y)
    orbit = 0.7 + 0.5 * x - 0.3 * y + 0.2 * x * y
    image = orbit + 0.006 * elevation
    corrected, surface, _, stats = robust_polynomial_detrend(
        image, degree=2, elevation=elevation, terrain_degree=1,
    )
    assert np.sqrt(np.mean(corrected ** 2)) < 1e-8
    assert np.sqrt(np.mean((surface - image) ** 2)) < 1e-8
    assert stats.terrain_enabled
    assert abs(stats.elevation_phase_correlation_after) < 1e-6
