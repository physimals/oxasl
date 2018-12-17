"""
Wrapper for avscale command
"""

import fsl.utils.assertions as asrt
from fsl.wrappers import wrapperutils  as wutils

def extract_avscale_output(fn):
    def _wrapper(*args, **kwargs):
        capture = six.StringIO()
        log = kwargs.pop("log", {})
        old_stdout = log.get("stdout", None)
        log["stdout"] = capture
        result = fn(*args, log=log, **kwargs)
        print(result)
        if old_stdout is not None:
            old_stdout.write(capture.getvalue())
        for line in capture.getvalue().splitlines():
            print(line)
        return capture.getvalue()
        
    return _wrapper

@extract_avscale_output
@wutils.fileOrArray('matrixfile')
@wutils.fslwrapper
def avscale(matrixfile):
    """Wrapper for the ``avscale`` command.
    
    Required options:
    
    :arg matrixfile:    FLIRT transformation matrix 
    """
    cmd = ['avscale', '--allparams', matrixfile]
    return cmd
