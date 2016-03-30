Individual simulations
======================
The Sim is the core unit of functionality of ``mdsynthesis``. They function as
markers for simulation data, giving an interface to raw topologies and
trajectories by way of MDAnalysis.

The components documented here are those included within ``mdsynthesis``.
However, the API elements of `:mod:datreant.core` and `:mod:datreant.data`
are also available for use with Sims.

.. _Sim_api:

Sim
---
The class :class:`mdsynthesis.Sim` is the central object of ``mdsynthesis``. 

.. autoclass:: mdsynthesis.Sim
    :members:
    :inherited-members:

.. _Selections_api:

Udef
````
The class :class:`mdsynthesis.limbs.Udef` is the interface used by a Sim to
define its :class:`MDAnalysis.Universe`.

.. autoclass:: mdsynthesis.limbs.Udef
    :members:
    :inherited-members:

AtomSelections
``````````````
The class :class:`mdsynthesis.limbs.AtomSelections` is the interface used by Sims to
get :class:`MDAnalysis.AtomGroup` objects from stored selection definitions. 

.. autoclass:: mdsynthesis.limbs.AtomSelections
    :members:
    :inherited-members:
