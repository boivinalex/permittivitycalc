# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 19:10:34 2018

@author: alex

A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='permittivitycalc',  

    version='0.6.0',

    description='Scripts to calculate and plot the complex permittivity from S-parameter data',  

    long_description=long_description,  

    url='https://github.com/boivinalex/permittivitycalc',
    
    author='Alexandre Boivin',  

    author_email='alex.boivin@utoronto.ca',
    
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='python dielectric permittivity s-parameter data analysis',  

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    
    install_requires=['numpy',
                      'uncertainties',
                      'matplotlib',
                      'seaborn',
                      'cycler',
                      'scipy',
                      'lmfit==0.9.14',
                      'emcee==3.0.2',
                      'corner',
                      'tqdm',
                      ],
                      
    python_requires='>=3.6',

    package_data={  # Optional
        'permittivitycalc': ['data/rexolite_PAL.txt',
                 'data/serpentine_dry.txt',
                 ],
    },
)
