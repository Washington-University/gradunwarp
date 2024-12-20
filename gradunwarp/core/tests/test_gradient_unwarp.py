import os

import numpy as np
import nibabel as nb
from nibabel.tmpdirs import InTemporaryDirectory

from ..gradient_unwarp import GradientUnwarpRunner


class Arguments:
    """Just something to dump args into"""


def test_trivial_unwarp():
    with InTemporaryDirectory():
        # Siemens Allegra coefficient arrays are 15x15, keeping things small and fast
        coef_file = "allegra.coef"
        open(coef_file, 'wb').close()  # touch
        assert os.path.exists(coef_file)

        orig_arr = np.arange(24).reshape(2, 3, 4)

        # Use centered LAS affine for simplicity. Easiest way to get it is
        # creating the image and asking nibabel to make it for us.
        img = nb.Nifti1Image(orig_arr.astype("float32"), None)
        img.set_sform(img.header.get_base_affine(), 1)
        img.set_qform(img.header.get_base_affine(), 1)
        img.to_filename("test.nii")

        args = Arguments()
        args.infile = "test.nii"
        args.outfile = "out.nii"
        args.vendor = "siemens"
        args.coeffile = coef_file

        unwarper = GradientUnwarpRunner(args)
        unwarper.run()
        unwarper.write()

        # No change
        unwarped_img = nb.load("out.nii")
        assert np.allclose(unwarped_img.affine, img.affine)
        assert np.allclose(unwarped_img.get_fdata(), orig_arr)

        # Rerun with right-handed image
        ras_img = nb.as_closest_canonical(img)
        ras_img.to_filename("test_ras.nii")
        # confirm we do have affine/data flip
        assert np.allclose(ras_img.affine, np.diag([-1, 1, 1, 1]).dot(img.affine))
        assert np.allclose(ras_img.get_fdata(), np.flip(orig_arr, axis=0))

        args.infile = "test_ras.nii"
        args.outfile = "out_ras.nii"

        unwarper = GradientUnwarpRunner(args)
        unwarper.run()
        unwarper.write()

        # Output matches input
        unwarped_img = nb.load("out_ras.nii")
        assert np.allclose(unwarped_img.affine, ras_img.affine)
        assert np.allclose(unwarped_img.get_fdata(), np.flip(orig_arr, axis=0))

        # Make sure we don't have open filehandles that would prevent Windows
        # from cleaning up tmpdir
        del unwarper
        del img, unwarped_img, ras_img
