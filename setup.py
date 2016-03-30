#! /usr/bin/python
"""Setuptools-based setup script for MDSynthesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages

setup(name='mdsynthesis',
      version='0.6.0-dev',
      maintainer='David Dotson', 
      maintainer_email='dotsdl@gmail.com',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      license='GPL 2',
      long_description=open('README.rst').read(),
      install_requires=[
                'datreant.core>=0.6.0',
                'datreant.data>=0.6.0',
                'MDAnalysis>=0.11.0',
                ],
     )
