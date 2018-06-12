import os
import numpy as np

import oxasl.fslwrap as fsl

def test_numpy_image():
    d = np.random.rand(5, 5, 5, 6)
    img = fsl.Image("asldata", data=d)
    cwd = os.getcwd()
    assert img.iname == "asldata"
    assert img.ipath == os.path.join(cwd, "asldata")
    assert img.fname == "asldata.nii.gz"
    assert img.fpath == os.path.join(cwd, "asldata.nii.gz")
