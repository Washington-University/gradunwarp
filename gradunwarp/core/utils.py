### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the gradunwarp package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
from collections import namedtuple

import numpy as np
import nibabel as nib
from nibabel.affines import apply_affine

# This is a container class that has 3 np.arrays which contain
# the x, y and z coordinates respectively. For example, the output
# of a meshgrid belongs to this
# x, y, z = meshgrid(np.arange(5), np.arange(6), np.arange(7))
# cv = CoordsVector(x=x, y=y, z=z)
CoordsVector = namedtuple('CoordsVector', 'x, y, z')

def transform_coordinates(A, M):
    transformed = apply_affine(M, np.stack(A).T).T
    return CoordsVector(*transformed)

def get_vol_affine(infile):
    nibimage = nib.load(infile)
    return np.asanyarray(nibimage.dataobj), nibimage.affine
