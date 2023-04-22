#!/usr/bin/env python
from glob import glob
import os
from setuptools import setup
import sys


def get_alntools_version():
    sys.path.insert(0, "alntools")
    import version
    return version.__version__


with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = []
test_requirements = []


on_rtd = os.environ.get("READTHEDOCS", None)

if not on_rtd:
    with open("requirements.txt") as requirements_file:
        requirements_lines = requirements_file.readlines()
        for line in requirements_lines:
            requirements.append(line)

setup(
    name="alntools",
    version=get_alntools_version(),
    description="Alignment tools",
    long_description=readme + "\n\n" + history,
    author="Matthew Vincent",
    author_email="matt.vincent@jax.org",
    url="https://github.com/churchill-lab/alntools",
    packages=[
        "alntools",
    ],
    package_dir={"alntools": "alntools"},
    include_package_data=True,
    scripts=glob("bin/*"),
    install_requires=requirements,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords="alntools",
    classifiers=[
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.10",
    ],
    test_suite="tests",
    tests_require=test_requirements,
)
