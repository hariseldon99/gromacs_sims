#!/usr/bin/env python3
"""
Setup script for gmx_protlig.
"""

from setuptools import setup, find_packages
import os

here = os.path.abspath(os.path.dirname(__file__))

long_description = ""
readme_path = os.path.join(here, "README.md")
if os.path.exists(readme_path):
    with open(readme_path, "r", encoding="utf-8") as fh:
        long_description = fh.read()

setup(
    name="gmx_protlig",
    version="1.0.0",
    description="Modular Python package & CLI suite for GROMACS Protein-Ligand MD simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(where=".", include=["gmx_protlig*", "protlig_api*"]),
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "gmx-protlig = gmx_protlig.cli:main",
            "gmx-protlig-step1 = gmx_protlig.steps:step1_cli",
            "gmx-protlig-step2 = gmx_protlig.steps:step2_cli",
            "gmx-protlig-step3 = gmx_protlig.steps:step3_cli",
            "gmx-protlig-step4 = gmx_protlig.steps:step4_cli",
            "gmx-protlig-step5 = gmx_protlig.steps:step5_cli",
            "gmx-protlig-step6 = gmx_protlig.steps:step6_cli",
            "gmx-protlig-step7 = gmx_protlig.steps:step7_cli",
            "gmx-protlig-step8 = gmx_protlig.steps:step8_cli",
            "gmx-protlig-batch = gmx_protlig.batch:batch_cli",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
