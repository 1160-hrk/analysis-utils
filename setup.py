#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Setup script for analysis-utils package
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements from requirements.txt
def read_requirements():
    with open('requirements.txt', 'r', encoding='utf-8') as f:
        requirements = []
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Remove version constraints for setup.py
                package = line.split('>=')[0].split('==')[0].split('<')[0]
                requirements.append(package)
        return requirements

setup(
    name="analysis-utils",
    version="0.1.0",
    author="hirokitsusaka",
    author_email="",
    description="Python-based analysis toolkit for scientific data processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourname/analysis-utils",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.12",
    install_requires=[
        "numpy>=1.24.0",
        "scipy>=1.10.0",
        "matplotlib>=3.6.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "pylint>=2.17.0",
            "mypy>=1.0.0",
            "isort>=5.12.0",
        ],
        "jupyter": [
            "jupyter>=1.0.0",
            "jupyterlab>=4.0.0",
            "notebook>=6.5.0",
            "ipywidgets>=8.0.0",
        ],
        "docs": [
            "sphinx>=6.0.0",
            "sphinx-rtd-theme>=1.2.0",
        ],
        "all": [
            "pandas>=2.0.0",
            "seaborn>=0.12.0",
            "tqdm>=4.65.0",
        ]
    },
    include_package_data=True,
    zip_safe=False,
) 