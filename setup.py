#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''
This program is a basic python conversion of Mick Watson's Ideel.
It reads one or more input FASTA files and for each file it will use
prodigal for rapid annotation, then run diamond blast, then compare the
query length to hit length.

It was built with the help of 'Bionitio'
'''


setup(
    name='pydeel',
    version='0.1.0.0',
    author='Alistair Legione',
    author_email='legionea@unimelb.edu.au',
    packages=['pydeel'],
    package_dir={'pydeel': 'pydeel'},
    entry_points={
        'console_scripts': ['pydeel = pydeel.pydeel:main']
    },
    url='https://github.com/alegione/pydeel',
    license='LICENSE',
    description=('Assembly completion by annotation assessment'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["argparse", "pandas", "altair", "seaborn", "selenium", "datetime"],
)
