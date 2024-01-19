### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the gradunwarp package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
from __future__ import print_function
import numpy as np
from collections import namedtuple
import math


# This is a container class that has 3 np.arrays which contain
# the x, y and z coordinates respectively. For example, the output
# of a meshgrid belongs to this
# x, y, z = meshgrid(np.arange(5), np.arange(6), np.arange(7))
# cv = CoordsVector(x=x, y=y, z=z)
CoordsVector = namedtuple('CoordsVector', 'x, y, z')

from nibabel.affines import apply_affine

def transform_coordinates(A, M):
    transformed = apply_affine(M, np.stack(A).T).T
    return CoordsVector(*transformed)

def get_vol_affine(infile):
    try:
        import nibabel as nib
    except ImportError:
        raise ImportError('gradunwarp needs nibabel for I/O of mgz/nifti files.'
                          ' Please install')
    nibimage = nib.load(infile)
    return np.asanyarray(nibimage.dataobj), nibimage.affine

factorial = math.factorial
