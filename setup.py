#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Standalone command line tool to visualize coverage from a BAM file
"""

from setuptools import setup, find_packages
import codecs

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

name = "tinycov"

MAJOR = 0
MINOR = 3
MAINTENANCE = 1
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MAINTENANCE)

LICENSE = "MIT"
URL = "https://github.com/cmdoret/tinycov"

DESCRIPTION = __doc__.strip("\n")

with codecs.open("README.md", encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

with open("requirements.txt", "r") as f:
    REQUIREMENTS = f.read().splitlines()

with open("tinycov/version.py", "w") as f:
    f.write("__version__ = '{}'\n".format(VERSION))


setup(
    name=name,
    author="cmdoret",
    author_email="cyril.matthey-doret@pasteur.fr",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    version=VERSION,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    url=URL,
    packages=find_packages(),
    python_requires=">=3.6",
    include_package_data=True,
    install_requires=REQUIREMENTS,
    long_description_content_type="text/markdown",
    entry_points={"console_scripts": ["tinycov=tinycov.__main__:cli"]},
)
