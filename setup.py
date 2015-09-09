#! /usr/bin/python
"""Setuptools-based setup script for MDSynthesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='mdsynthesis',
      version='0.6.0-dev',
      maintainer='David Dotson', 
      maintainer_email='dotsdl@gmail.com',
      packages=['mdsynthesis', 'mdsynthesis.tests'],
      license='GPL 2',
      long_description=open('README.rst').read(),
      dependency_links=[
      'http://github.com/dotsdl/datreant/tarball/develop#egg=datreant-0.6.0-dev',
      ],
      install_requires=[
                'numpy',
                'numexpr',
                'Cython',
                'datreant>=0.6.0-dev',
                'pandas>=0.16.1',
                'tables>=3.2.0',
                'h5py>=2.5.0',
                'MDAnalysis>=0.11.0',
                'scandir>=1.0',
                'PyYAML>=3.11'
                ],
     )
