============================================
Leveraging molecular dynamics data with Sims
============================================
A :class:`~mdsynthesis.Sim` is a :class:`~datreant.core.Treant` with
specialized components for working with molecular dynamics data. In particular,
it can store a definition for an MDAnalysis
:class:`~MDAnalysis.core.AtomGroup.Universe` for painless recall, as well as
custom atom selections.

.. note:: Since Sims are Treants, everything that applies to Treants applies
          to Sims as well. See the `datreant documentation
          <http://datreant.readthedocs.org/>`_ for how Treants work and
          effectively used.

As with a normal Treant, to generate a Sim from scratch, we need only give it a
name ::

    >>> from mdsynthesis import Sim
    >>> s = Sim('adk')
    >>> s
    <Sim: 'adk'>

And we can immediately give the Sim characteristics like :ref:`tags
<datreantcore:Tags_guide>`::

    >>> s.tags.add('biased', 'closed-to-open', 'mep')
    >>> s.tags
    <Tags(['biased', 'closed-to-open', 'MEP'])>

and :ref:`categories <datreantcore:Categories_guide>`::

    >>> s.categories['sampling method'] = 'DIMS'
    >>> s.categories['sampling '] = 'heay atom'
    <Categories({'sampling ': 'heay atom', 'sampling method': 'DIMS'})>

These can be used later to filter and group Sims when we have many to work
with. They can also be used as switches for analysis code, since we may need to
do different things depending on e.g. the type of sampling method used to
produce the trajectory.

Defining the Universe
=====================
What makes a Sim different from a basic Treant is that it can store a
:class:`MDAnalysis.Universe` definition. We can access the Sim's Universe
directly with::

    >>> s.universe

but doing so now gives back ``None``. However, if we define a topology for the
Universe with::

    >>> s.udef.topology = 'path/to/topology.psf'

then we get back a Universe built with this topology instead::

    >>> s.universe
    <Universe with 3341 atoms and 3365 bonds>

we can also add a trajectory::

    >>> s.udef.trajectory = 'path/to/trajectory.dcd'

and our Universe is re-initialized with both the defined topology and trajectory::

    >>> s.universe.trajectory
    <DCDReader /home/bob/research/path/to/trajectory.dcd with 98 frames of 3341 atoms>

We can define our Universe as having multiple trajectories by giving a list of
paths instead, and this will work as well. Internally, the Universe generated
will use the :class:`MDAnalysis.coordinates.base.ChainReader` for treating the
trajectories as a contiguous whole.

The Universe definition is persistent, so we can get back an identical Universe
later from another Python session with our Sim::

    >>> import mdsynthesis as mds
    >>> s = mds.Sim('adk')
    >>> s.universe
    <Universe with 3341 atoms and 3365 bonds>

.. note:: Changing the topology or trajectory definition will reload the
          Universe automatically. This means that any AtomGroups you are
          working with will not point to the new Universe, but perhaps the old
          one, so it's generally best to regenerate them manually.

Storing keyword arguments
-------------------------
If the Universe needed requires keyword arguments on initialization, these can
be stored as well. For example, if our topology was a PDB file and we wanted
bonds to be guessed upfront, we could make this happen every time::

    >>> s.udef.kwargs = {'guess_bonds': True}

Reinitializing the Universe
---------------------------
If you make modifications to the Universe but you want to restore the original
from its definition, you can force it to reload with::

    >>> s.udef.reload()

API Reference: UniverseDefinition
---------------------------------
See the :ref:`UniverseDefinition_api` API reference for more details.


Storing custom atom selections
==============================
MDAnalysis includes its own selection language for extracting
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` objects, which function as an
ordered list of atoms from the system. The selection strings needed to specify
these can be long and complex, and sometimes multiple selection strings are
required in a particular order to extract a given AtomGroup from all the atoms
in the Universe. What's more, for different simulation systems the same
selection of atoms, e.g. the "solvent", might require a different set of
strings.

To make this more manageable, we can store custom atom selections within our
Sim. Say we want to select the lid, core, and NMP domains of adenylate
kinase, the protein we simulated. We can store these immediately::

    >>> s.atomselections['lid'] = 'resid 122:159'
    >>> s.atomselections['nmp'] = 'resid 30:59'
    >>> s.atomselections['core'] = ['resid 1:29', 'resid 60:121', 'resid 160:214']

We can now get new AtomGroups back for each selection at any time ::

    >>> s.atomselections['lid']
    <AtomGroup with 598 atoms>
    
    >>> s.atomselections['core']
    <AtomGroup with 2306 atoms>

and we don't have to remember or know how 'lid' or 'core' are defined for this
particular system. If we have other simulations of adenylate kinase performed
with other molecular dynamics engines or with different forcefields, we can
store the atom selection strings required for those systems in the same way,
perhaps using the same names 'lid', 'core', etc. This abstraction allows us to
work with many variants of a simulation system without getting stuck on the
details every time.

.. note:: Storing a list of strings as a selection will apply them in order,
          producing an AtomGroup concatenated from each one in the same way
          as providing multiple strings to
          :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` does. This
          is especially useful when storing selections used for structural
          alignments.

Want just the selection strings back? We can use
:meth:`~mdsynthesis.limbs.AtomSelections.define`::

    >>> s.atomselections.define('lid')
    ['resid 122:159']

Note that selections are always stored as lists, even if only a single
selection string was given.

Atom selections from atom indices 
---------------------------------
Already have an AtomGroup, and just want to define it according to its atom
indices instead of as a selection string? We can do that, too::

    >>> lid = s.universe.select_atoms('resid 122:159')
    >>> s.atomselections['lid'] = lid
    >>> s.atomselections['lid']
    <AtomGroup with 598 atoms>

API Reference: AtomSelections
-----------------------------
See the :ref:`AtomSelections_api` API reference for more details.

API Reference: Sim
==================
See the :ref:`Sim_api` API reference for more details.
