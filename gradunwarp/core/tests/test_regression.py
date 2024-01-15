import os
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

import numpy as np
from numpy.testing import assert_equal, assert_array_equal, \
    assert_array_almost_equal, assert_almost_equal

from gradunwarp.core import coeffs, utils
from gradunwarp.core.unwarp_resample import siemens_B, cart2sph

def test_siemens_B():
    gradfile = os.path.join(DATA_DIR, 'gradunwarp_coeffs.grad')

    siemens_coeffs = coeffs.get_coefficients('siemens', gradfile)
    R0 = siemens_coeffs.R0_m  * 1000

    vec = np.linspace(-300, 300, 60, dtype=np.float32)
    x, y ,z = np.meshgrid(vec, vec, vec)
    r, cos_theta, theta, phi = cart2sph(x, y, z)

    ref_b = np.load(os.path.join(DATA_DIR, 'siemens_B_output.npz'))

    for d in 'xyz':
        alpha_d = getattr(siemens_coeffs, f"alpha_{d}")
        beta_d = getattr(siemens_coeffs, f"beta_{d}")
        bd = siemens_B(alpha_d, beta_d, r, cos_theta, theta, phi, R0)

        # changes in legendre function is causing differences at 6th decimal
        assert_array_almost_equal(ref_b[f"b{d}"], bd, decimal=4)
