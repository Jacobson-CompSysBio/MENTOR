#!/usr/bin/env python3

'''
Based on template from:
https://packaging.python.org/tutorials/packaging-projects/
'''

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="functional_partitioning",
    version="0.1.1",
    author='J. I. Miller',
    author_email='millerji@ornl.gov',
    description='Functional partitioning.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/izaakm/jail-functional-partitioning',
    project_urls={
        "Bug Tracker": "https://github.com/izaakm/jail-functional-partitioning/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.9",
)

# END.
