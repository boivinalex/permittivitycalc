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
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='permittivity-calc',  

    version='0.0.4dev0',

    description='Scripts to calculate and plot the complex permittivity from S-parameter data',  

    long_description=long_description,  

    url='https://github.com/boivinalex/permittivity-calc',
    
    author='Alexandre Boivin',  

    author_email='alex.boivin@mail.utoronto.ca',
    
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='python dielectric permittivity s-parameter data analysis',  

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    
    install_requires=['tkinter',
                      'numpy',
                      'uncertainties',
                      'matplotlib',
                      'seaborn',
                      'cycler',
                      'scipy',
                      ],
                      
    python_requires='>=3',

#    extras_require={
#        'dev': ['check-manifest'],
#        'test': ['coverage'],
#    },

    package_data={  # Optional
        'data': ['2.5hrs.txt',
                 'atm.txt',
                 ],
    },

#    data_files=[('my_data', ['data/data_file'])],

#    entry_points={
#        'console_scripts': [
#            'sample=sample:main',
#        ],
#    },
)