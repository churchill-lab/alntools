#!/usr/bin/env python
# -*- coding: utf-8 -*-
from glob import glob
import os
from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = []
test_requirements = []


on_rtd = os.environ.get('READTHEDOCS', None)

if not on_rtd:
    requirements.append('Click>=6.0')
    requirements.append('numexpr>=2.3.1')
    requirements.append('numpy==1.8.2')
    requirements.append('scipy==0.13.3')
    requirements.append('pysam>=0.6')
    requirements.append('cython>=0.13')
    requirements.append('tables==3.1.0')
    requirements.append('biopython>=1.63')
    requirements.append('future>=0.16')
    requirements.append('flask>=1.0.2')

setup(
    name='alntools',
    version='0.1.0',
    description="Alignment tools",
    long_description=readme + '\n\n' + history,
    author="Matthew Vincent",
    author_email='matt.vincent@jax.org',
    url='https://github.com/churchill-lab/alntools',
    packages=[
        'alntools',
    ],
    package_dir={'alntools':
                 'alntools'},
    include_package_data=True,
    scripts=glob("bin/*"),
    #setup_requires=requirements,
    install_requires=requirements,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='alntools',
    classifiers=[
        'Development cous :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
