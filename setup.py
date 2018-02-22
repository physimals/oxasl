#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages

from oxasl import __version__

# Read in requirements from the requirements.txt file.
with open('requirements.txt', 'rt') as f:
    requirements = [l.strip() for l in f.readlines()]

# Generate a list of all of the packages that are in your project.
packages = find_packages()

setup(
    name='oxasl',
    description='Python library for manipulating and modelling ASL data',
    url='',
    author='Martin Craig',
    author_email='martin.craig@eng.ox.ac.uk',
    license='',
    packages=packages,
    version=__version__,
    install_requires=requirements,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules'],
)
