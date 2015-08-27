#! /usr/bin/python
"""Setuptools-based setup script for MDSynthesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='mdsynthesis',
      version='0.5.0',
      maintainer='David Dotson', 
      maintainer_email='dotsdl@gmail.com',
      packages=['mdsynthesis', 'mdsynthesis.tests'],
      license='GPL 2',
      long_description=open('README.rst').read(),
      requires=['datreant',
                'pandas',
                'tables',
                'h5py',
                'MDAnalysis',
                'scandir']
     )
