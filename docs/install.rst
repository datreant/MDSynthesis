======================
Installing MDSynthesis
======================
Since MDSynthesis uses HDF5 as the file format of choice for persistence of
intermediate datasets, you will first need to install the HDF5 libraries either
using your package manager or manually.

On **Ubuntu 14.04** this will be ::

    apt-get install libhdf5-serial-1.8.4 libhdf5-serial-dev

and on **Arch Linux** ::

    pacman -S hdf5

You can then install MDSynthesis from `PyPI <https://pypi.python.org/>`_
using pip::

    pip install mdsynthesis

It is also possible to use ``--user`` to install into your user's site-packages
directory::

    pip install --user mdsynthesis

Be aware that some dependencies may require ``numpy`` and/or Cython to
be installed beforehand!

Installing using Conda
======================

We also provide conda packages. Using the conda packages has the advantage that
they contain the hdf5 library, so you can install MDSynthesis easily on all unix
systems including macOS. To install MDSynthesis with conda type::

    conda install -c datreant mdsynthesis

Dependencies
============
The dependencies of MDSynthesis are:

- `MDAnalysis`_: 0.18 or higher
- `datreant`_: 1.0.0 or higher
- `pandas`_: 0.16.1 or higher
- `PyTables`_: 3.2.0 or higher
- `h5py`_: 2.5.0 or higher

.. _`MDAnalysis`: http://www.mdanalysis.org
.. _`datreant`: http://datreant.readthedocs.org/
.. _`pandas`: http://pandas.pydata.org/
.. _`PyTables`: http://www.pytables.org/
.. _`h5py`: http://www.h5py.org/


Installing from source
======================
To install from source, clone the repository and switch to the master branch ::

    git clone git@github.com:datreant/MDSynthesis.git
    cd MDSynthesis
    git checkout master

Installation of the packages is as simple as ::

    pip install .

which installs ``mdsynthesis`` in the system wide python directory (this may
require administrative privileges).

It is also possible to use ``--user`` to install into your user's site-packages
directory::

    pip install --user .
