#!/usr/bin/env python
import os
import subprocess
import re

from setuptools import setup
from setuptools import find_packages

def git_version():
    # Full version includes the Git commit hash
    full_version = subprocess.check_output('git describe --dirty', shell=True).decode("utf-8").strip(" \n")

    # Python standardized version in form major.minor.patch.dev<build>
    p = re.compile(r"v?(\d+\.\d+\.\d+(-\d+)?).*")
    m = p.match(full_version)
    if m is not None:
        std_version = m.group(1).replace("-", ".dev")
    else:
        raise RuntimeError("Failed to parse version string %s" % full_version)

    return full_version, std_version

def set_python_version(rootdir, version):
    vfile = open(os.path.join(rootdir, "oxasl", "_version.py"), "w")
    vfile.write("__version__ = '%s'" % version)
    vfile.close()

# Read in requirements from the requirements.txt file.
with open('requirements.txt', 'rt') as f:
    requirements = [l.strip() for l in f.readlines()]

# Generate a list of all of the packages that are in your project.
packages = find_packages()

rootdir = os.path.join(os.path.abspath(os.path.dirname(__file__)))
fullv, stdv = git_version()
set_python_version(rootdir, stdv)

setup(
    name='oxasl',
    description='Python library for manipulating and modelling ASL data',
    url='',
    author='Martin Craig',
    author_email='martin.craig@eng.ox.ac.uk',
    license='',
    packages=packages,
    version=stdv,
    install_requires=requirements,
    entry_points = {
        'console_scripts' : [
            "asl_preproc=oxasl.preproc:main", 
            "asl_basil=oxasl.basil:main",
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules'],
)
