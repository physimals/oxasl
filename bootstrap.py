import os
import subprocess
import re

def git_version():
    # Full version includes the Git commit hash
    full_version = subprocess.check_output('git describe --dirty', shell=True).strip(" \n")

    # Standardized version in form major.minor.patch-build
    p = re.compile("v?(\d+\.\d+\.\d+(-\d+)?).*")
    m = p.match(full_version)
    if m is not None:
        std_version = m.group(1)
    else:
        raise RuntimeError("Failed to parse version string %s" % full_version)

    return full_version, std_version

def set_python_version(rootdir, version):
    vfile = open(os.path.join(rootdir, "oxasl", "__init__.py"), "w")
    vfile.write("__version__='%s'" % version)
    vfile.close()

def set_conda_version(rootdir, version):
    cfile = os.path.join(rootdir, "meta.yaml.in")
    with open(cfile, "r") as f:
        content = f.read()
    content = content.replace("@GIT_VERSION@", version)
    ofile = os.path.join(rootdir, "meta.yaml") 
    with open(ofile, "w") as f:
        f.write(content)

rootdir = os.path.join(os.path.abspath(os.path.dirname(__file__)))
fullv, stdv = git_version()
set_python_version(rootdir, stdv)
set_conda_version(rootdir, stdv)
