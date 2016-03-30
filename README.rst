=============================================================
MDSynthesis: a persistence engine for molecular dynamics data
=============================================================

|zen| |docs| |build| |cov|

Although the raw data for any study involving molecular dynamics simulations are
the full trajectories themselves, often we are most interested in
lower-dimensional measures of what is happening. These measures may be as simple
as the distance between two specific atoms, or as complex as the percentage of
contacts relative to some native structure. In any case, it may be time-consuming
to obtain these lower-dimensional intermediate data, and so it is useful to store
them.

Stay organized
==============
MDSynthesis is designed to perform the logistics of medium-to-large-scale
analysis of many trajectories, individually or as entire groups. It should
allow the scientist to operate at a high level when working with the data,
while MDSynthesis handles the details of storing and recalling this data.

In other words, MDSynthesis lets the computer do the boring work of keeping
track of where things are and how they are stored.

Efficiently store intermediate data from individual simulations for easy recall
-------------------------------------------------------------------------------
The MDSynthesis **Sim** object gives an interface to raw simulation data
through `MDAnalysis`_. Data structures generated from raw trajectories (pandas
objects, numpy arrays, or any pure python structure) can then be stored and
easily recalled later. Under the hood, datasets are stored in the efficient
HDF5 format when possible.

.. _MDAnalysis: http://www.mdanalysis.org

Powered by ``datreant`` under the hood
--------------------------------------
MDSynthesis is built on top of the general-purpose `datreant`_ library.  The
Sim is a :ref:`~datreant.core.Treant` with special features for working with
molecular dynamics data, but every feature of datreant applies to MDSynthesis.

.. _`datreant`: http://datreant.org/

Documentation
=============
A brief user guide is available on `Read the Docs
<http://mdsynthesis.readthedocs.org/>`__.

Contributing
============
This project is still under heavy development, and there are certainly rough
edges and bugs. Issues and pull requests welcome!

MDSynthesis follows the development model of `datreant`_; see the
`contributor's guide`_ to learn how to get started with contributing back.

.. _`datreant`: http://datreant.readthedocs.org/
.. _`contributor's guide`: http://datreant.readthedocs.org/en/latest/contributing.html

.. |docs| image:: https://readthedocs.org/projects/mdsynthesis/badge/?version=develop
    :alt: Documentation Status
    :scale: 100%
    :target: https://readthedocs.org/projects/mdsynthesis

.. |build| image:: https://travis-ci.org/datreant/MDSynthesis.svg?branch=develop
    :alt: Build Status
    :target: https://travis-ci.org/datreant/MDSynthesis

.. |cov| image:: http://codecov.io/github/datreant/MDSynthesis/coverage.svg?branch=develop
    :alt: Code Coverage
    :scale: 100%
    :target: http://codecov.io/github/datreant/MDSynthesis?branch=develop

.. |zen| image:: https://zenodo.org/badge/doi/10.5281/zenodo.18851.svg   
    :alt: Citation
    :target: http://dx.doi.org/10.5281/zenodo.18851
