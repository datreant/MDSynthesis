=============================================================
MDSynthesis: a persistence engine for molecular dynamics data
=============================================================

|zen| |docs| |build| |cov|

As computing power increases, it is now possible to produce hundreds of
molecular dynamics simulation trajectories that vary widely in length,
system size, composition, starting conditions, and other parameters. Managing
this complexity in ways that allow use of the data to answer scientific
questions has itself become a bottleneck. MDSynthesis is an answer to this
problem.

Built on top of `datreant`_, MDSynthesis gives a Pythonic interface to
molecular dynamics trajectories using `MDAnalysis`_, giving the ability to work
with the data from many simulations scattered throughout the filesystem
with ease. It makes it possible to write analysis code that can work across
many varieties of simulation, but even more importantly, MDSynthesis allows
interactive work with the results from hundreds of simulations at once without
much effort. 

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
Sim is a `Treant`_ with special features for working with molecular dynamics
data, but every feature of datreant applies to MDSynthesis.

.. _Treant: http://datreant.readthedocs.org/en/latest/treants.html

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

.. _datreant: http://datreant.readthedocs.org/
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

.. |zen| image:: https://zenodo.org/badge/13742/datreant/MDSynthesis.svg
    :alt: Citation
    :target: https://zenodo.org/badge/latestdoi/13742/datreant/MDSynthesis
