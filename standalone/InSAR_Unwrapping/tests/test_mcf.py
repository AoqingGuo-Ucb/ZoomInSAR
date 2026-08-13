import numpy as np

from insar_unwrapping.mcf import unwrap_mcf
from insar_unwrapping.residues import wrap


def test_unwraps_smooth_ramp_modulo_constant():
    y, x = np.mgrid[:40, :50]
    truth = 0.18 * x + 0.12 * y
    result, offsets, cuts, stats = unwrap_mcf(wrap(truth), np.ones(truth.shape), 0.0)
    difference = result - truth
    difference -= np.nanmedian(difference)
    assert np.nanmax(np.abs(difference)) < 1e-10
    assert stats.unwrapped_pixels == truth.size
    assert offsets.shape == cuts.shape == truth.shape
