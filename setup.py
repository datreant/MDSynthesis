#! /usr/bin/python
"""Setuptools-based setup script for MDSynthesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages

setup(name='mdsynthesis',
      version='0.6.2-dev',
      description='a persistence engine for molecular dynamics data',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      url='http://mdsynthesis.readthedocs.org/',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
      packages=find_packages('src'),
      package_dir={'': 'src'},
      license='GPL 2',
      long_description=open('README.rst').read(),
      entry_points={'console_scripts':
                  ['mds_06to1=mdsynthesis.scripts.mds_06to1:main']},
      install_requires=[
                'datreant',
                'MDAnalysis>=0.16.0',
                'tables', 'h5py', 'numpy', 'pandas'
                ],
      )
