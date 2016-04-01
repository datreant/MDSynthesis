=============================================================
MDSynthesis: a persistence engine for molecular dynamics data
=============================================================
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

.. warning:: This package is **experimental**. It is not API stable, and has
             many rough edges and limitations. It is, however, usable.

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
    sims
    datreant
    api 
