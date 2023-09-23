#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name='Mana',
    version='1.0',
    author='Ted Wong',
    author_email='t.wong@garvan.org.au',
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/scchess/Mana',
    license='LICENSE.txt',
    description='Command-line tool that analyses ONT alignments (.BAM) to report quality control statistics.',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    zip_safe=True
)
