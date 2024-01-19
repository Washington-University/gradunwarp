import os

import numpy as np
from numpy.testing import assert_array_almost_equal

from gradunwarp.core.unwarp_resample import cart2sph, siemens_B
from gradunwarp.core import coeffs

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def test_siemens_B():
    gradfile = os.path.join(DATA_DIR, 'gradunwarp_coeffs.grad')

    siemens_coeffs = coeffs.get_coefficients('siemens', gradfile)
    R0 = siemens_coeffs.R0_m  * 1000

    vec = np.linspace(-300, 300, 60, dtype=np.float32)
    x, y ,z = np.meshgrid(vec, vec, vec)
    r, cos_theta, theta, phi = cart2sph(x, y, z)

    ref_b = np.load(os.path.join(DATA_DIR, 'siemens_B_output.npz'))

    for d in 'xyz':
        alpha_d = getattr(siemens_coeffs, "alpha_%s"%d)
        beta_d = getattr(siemens_coeffs, "beta_%s"%d)
        bd = siemens_B(alpha_d, beta_d, r, cos_theta, theta, phi, R0)

        # changes in legendre function is causing differences at 5th decimal
        assert_array_almost_equal(ref_b["b%s"%d], bd, decimal=5)
