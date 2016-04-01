======================
Installing MDSynthesis
======================
Since MDSynthesis requires ``datreant.data``, which uses HDF5 as the file
format of choice for persistence, you will first need to install the HDF5
libraries either using your package manager or manually. 

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

Dependencies
============
The dependencies of MDSynthesis are:

- `MDAnalysis`_: 0.15 or higher
- `datreant.core`_: 0.6.0 or higher
- `datreant.data`_: 0.6.0 or higher

.. _`MDAnalysis`: http://www.mdanalysis.org
.. _`datreant.core`: http://datreant.readthedocs.org/
.. _`datreant.data`: http://datreantdata.readthedocs.org/

Installing from source
======================
To install from source, clone the repository and switch to the master branch ::

    git clone git@github.com:datreant/MDSynthesis.git
    cd MDSynthesis
    git checkout master

Installation of the packages is as simple as

    pip install .

which installs ``mdsynthesis`` in the system wide python directory (this may
require administrative privileges).

It is also possible to use ``--user`` to install into your user's site-packages
directory::

    pip install --user .
