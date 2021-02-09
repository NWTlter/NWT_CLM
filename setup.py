#!/usr/bin/env python

"""The setup script."""

import sys
from setuptools import setup, find_packages
# import versioneer


# with open("requirements.txt") as f:
#     INSTALL_REQUIRES = f.read().strip().split("\n")

with open("README.md") as f:
    LONG_DESCRIPTION = f.read()

PYTHON_REQUIRES = '>=3.6'

description = ("Utilities for working with ctsm data")
setup(
    name="NWT_CLM",
    description=description,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    maintainer="Will Wieder",
    maintainer_email="wwieder@ucar.edu",
    url="https://github.com/NWTlter/NWT_CLM",
    py_modules=['NWT_CLM'],
    packages=find_packages(),
    python_requires=PYTHON_REQUIRES,
    license="Apache",
    keywords="NWT_CLM",
    version='0.0.1',
)
