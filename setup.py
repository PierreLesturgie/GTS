#!/usr/bin/env python


from setuptools import setup, find_packages

about = {}
with open("gts/_version.py") as f:
    exec(f.read(), about)

setup(
    name="gts",
    version=about["__version__"],
    description="Gene Tree Statistics Tools",
    author="Pierre Lesturgie",
    packages=find_packages(),
    install_requires=[
        # Add required packages here
        "numpy",
        "pandas",
        "msprime",
        "tqdm",
        "numpy",
        "pyslim",
        "tskit",
        "argparse"
    ],
    entry_points={
        "console_scripts": [
            "gts = gts.GTS:main",  # if GTS.py has a main()
        ]
    },
    python_requires=">=3.12",
)
