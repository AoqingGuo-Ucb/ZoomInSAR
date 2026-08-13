import numpy as np

from insar_filtering.filters import nonlocal_phase_filter


def circular_rmse(estimate, truth):
    return np.sqrt(np.mean(np.angle(estimate * np.conj(truth)) ** 2))


def test_nonlocal_filter_reduces_phase_noise_and_preserves_shape():
    rng = np.random.default_rng(7)
    y, x = np.mgrid[:48, :56]
    truth = np.exp(1j * (0.08 * x + 0.05 * y))
    noisy = truth * np.exp(1j * rng.normal(0, 0.65, truth.shape))
    filtered, concentration, support = nonlocal_phase_filter(
        noisy, np.full(truth.shape, 0.45), search_radius=3, patch_radius=1
    )
    assert filtered.shape == concentration.shape == support.shape == truth.shape
    assert circular_rmse(filtered, truth) < circular_rmse(noisy, truth)
    assert np.all((concentration >= 0) & (concentration <= 1))
