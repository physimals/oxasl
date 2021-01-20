"""This module provides wrapper functions for the FSL `FNIRTFILEUTILS
<https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/>`_ tool

.. autosummary::
   :nosignatures:

   fnirtfileutils
"""

import fsl.utils.assertions as asrt
from fsl.wrappers import wrapperutils  as wutils

@wutils.fileOrImage('src', 'ref', 'out', 'jac')
@wutils.fslwrapper
def fnirtfileutils(src, **kwargs):
    """Wrapper for the ``fnirt`` command.
    
    Compulsory arguments (You MUST set one or more of):
	    -i,--in		filename of input coefficient volume (to be converted)

    Optional arguments (You may optionally specify one or more of):
        -r,--ref	filename for reference volume
        -o,--out	filename for output (field/coef) volume - uses relative warp convention
        -f,--outformat	Output format [field spline], default=field
        -w,--warpres	Warp resolution (mm), only relevant when --outformat=spline
        -k,--knotspace	Knot-spacing (voxels), only relevant when --outformat=spline
        -j,--jac	filename for output (jacobian map) volume
        -a,--withaff	If set, the affine transform is included in the field/jacobian
        -v,--verbose	switch on diagnostic messages
        -h,--help	display this message
    """

    asrt.assertIsNifti(src)

    cmd = ['fnirtfileutils', '--in={}'.format(src)]
    cmd += wutils.applyArgStyle('--=', **kwargs)

    return cmd

