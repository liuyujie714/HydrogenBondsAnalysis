#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

# Read the contents of README file if it exists
try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = ''

setup(
    name='mdaHbonds',
    version='1.0.0',
    author='Yujie Liu',
    author_email='',
    description='Hydrogen bond analysis using MDAnalysis with GROMACS compatibility',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/liuyujie714/HydrogenBondsAnalysis',
    py_modules=['mdaHbonds'],
    classifiers=[
        'Development Status :: 1 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.8',
    install_requires=[
        'MDAnalysis>=2.0.0',
        'numpy',
        'tqdm',
    ],
    keywords='hydrogen-bonds, molecular-dynamics, GROMACS, MDAnalysis, trajectory-analysis',
    project_urls={
        'Bug Reports': 'https://github.com/liuyujie714/HydrogenBondsAnalysis/issues',
        'Source': 'https://github.com/liuyujie714/HydrogenBondsAnalysis',
        'Documentation': 'https://github.com/liuyujie714/HydrogenBondsAnalysis',
    }
)
