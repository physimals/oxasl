
"""
FSL-compatible wrappers for Fabber tools
"""
from fsl.wrappers import wrapperutils  as wutils
import fsl.utils.assertions as asrt

def fast(img, out_base="fast", n_classes=3, **kwargs):
    """
    Wrapper for FAST.

    :param img: Input image
    :param out_base: Base name for output data
    :param n_classes: Number of tissue classes

    Keyword arguments as per command line
    """
    fast_wrapper = _get_fast_wrapper(out_base, n_classes)
    return fast_wrapper(img, **kwargs)

def _get_fast_wrapper(out_base, n_classes):
    """
    Construct a FAST wrapper for a given base output name and number of classes

    This is required because the output images are dependent on these parameters
    """
    file_or_image = ['img']
    file_or_image += ["%s_%s" % (out_base, imtype) for imtype in ("seg", "pveseg", "mixeltype", "restore", "bias")]
    for idx in range(n_classes):
        for imtype in ("pve", "seg", "pveseg"):
            file_or_image.append('%s_%s_%i' % (out_base, imtype, idx))

    @wutils.fileOrImage(*file_or_image)
    @wutils.fslwrapper
    def _wrapper(img, **kwargs):
        asrt.assertIsNifti(img)

        # Output images are not command line parameters themselves, so 
        # remove them from kwargs
        for output_img in file_or_image[1:]:
            kwargs.pop(output_img, None)

        cmd = ['fast', '-v', '--out=%s' % out_base, '--class=%i' % n_classes] + wutils.applyArgStyle('--=', **kwargs)
        cmd.append(img)

        return cmd
    return _wrapper
