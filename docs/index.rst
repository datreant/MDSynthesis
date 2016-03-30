=============================================================
MDSynthesis: a persistence engine for molecular dynamics data
=============================================================
Although the raw data for any study involving molecular dynamics simulations are
the full trajectories themselves, often we are most interested in
lower-dimensional measures of what is happening. These measures may be as simple
as the distance between two specific atoms, or as complex as the percentage of
contacts relative to some native structure. Some measures may even be
comparisons of two or more trajectories against each other. In any case, it may
be time-consuming to obtain these lower-dimensional intermediate data, and so
it is useful to store them.

.. warning:: This package is **experimental**. It is not API stable, and has
             many rough edges and limitations. It is, however, usable.

Stay organized
==============
MDSynthesis is designed to perform the logistics of medium-to-large-scale
analysis of many trajectories, individually or as entire groups. It is intended
to allow the scientist to operate at a high level when working with the data,
while letting MDSynthesis handle the details of storing and recalling this
data. 

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
Sim is a :class:`~datreant.core.Treant` with special features for working with
molecular dynamics data, but every feature of datreant applies to MDSynthesis.

.. _`datreant`: http://datreant.org/

Getting MDSynthesis
===================
See the `installation instructions`_ for installation details. The package
itself is pure Python, but many of its dependencies are not.

If you want to work on the code, either for yourself or to contribute back to
the project, clone the repository to your local machine with::

    git clone https://github.com/datreant/MDSynthesis.git

.. _`installation instructions`: http://mdsynthesis.readthedocs.org/en/develop/install.html

Contributing
============
This project is still under heavy development, and there are certainly rough
edges and bugs. Issues and pull requests welcome!

MDSynthesis follows the development model of `datreant`_; see the
`contributor's guide`_ to learn how to get started with contributing back.

.. _`contributor's guide`: http://datreant.readthedocs.org/en/latest/contributing.html

--------------------------------------------------------------------------------

.. toctree::
    :maxdepth: 1
    :caption: User Documentation

    install
    sim
    api 
