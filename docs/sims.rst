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
    >>> s.categories['sampling '] = 'heavy atom'
    <Categories({'sampling ': 'heavy atom', 'sampling method': 'DIMS'})>

These can be used later to filter and aggregate Sims when we have many to work
with. They can also be used as switches for analysis code, since we may need to
do different things depending on, for example, the type of sampling method used to
produce the trajectory.

Defining the Universe
=====================
What makes a Sim different from a basic Treant is that it can store an
MDAnalysis :class:`~MDAnalysis.core.AtomGroup.Universe` definition. We can
access the Sim's Universe directly with::

    >>> s.universe

At this point, we get back ``None``. However, we can define the Sim's Universe
by giving it one direclty::

    >>> import MDAnalysis as mda
    >>> s.universe = mda.Universe('path/to/topology.psf',
                                  'path/to/trajectory.dcd')
    >>> s.universe
    <Universe with 3341 atoms and 3365 bonds>

The Universe definition is persistent, so we can get back an identical Universe
later from another Python session with our Sim::

    >>> import mdsynthesis as mds
    >>> s = mds.Sim('adk')
    >>> s.universe
    <Universe with 3341 atoms and 3365 bonds>


Changing the Universe definition
--------------------------------
.. warning:: This interface may be removed in a future release, but remains for
             now due to limitations in MDAnalysis. It is encouraged to set the
             Universe definition directly as shown above.

We can directly change the topology used for the Sim's Universe with :: 

    >>> s.universedef.topology = 'path/to/another/topology.psf'

then we get back a Universe built with this topology instead::

    >>> s.universe.filename
    '/home/bob/research/path/to/another/topology.psf'

We can also change the trajectory::

    >>> s.universedef.trajectory = 'path/to/another/trajectory.dcd'

which re-initializes the Universe with both the defined topology and trajectory::

    >>> s.universe.trajectory
    <DCDReader /home/bob/research/path/to/another/trajectory.dcd with 98 frames of 3341 atoms>

We can also define our Universe as having multiple trajectories by giving a
list of filepaths instead. Internally, the Universe generated will use the
MDAnalysis :class:`~MDAnalysis.coordinates.base.ChainReader` for treating the
trajectories as a contiguous whole.

.. note:: Changing the topology or trajectory definition will reload the
          Universe automatically. This means that any AtomGroups you are
          working with will not point to the new Universe, but perhaps the old
          one, so it's generally best to regenerate them manually.

Storing keyword arguments
-------------------------
If the Universe needed requires keyword arguments on initialization, these can
be stored as well. For example, if our topology was a PDB file and we wanted
bonds to be guessed upfront, we could make this happen every time::

    >>> s.universedef.kwargs = {'guess_bonds': True}

Reinitializing the Universe
---------------------------
If you make modifications to the Universe but you want to restore the original
from its definition, you can force it to reload with::

    >>> s.universedef.reload()

API Reference: UniverseDefinition
---------------------------------
See the :ref:`UniverseDefinition_api` API reference for more details.


Storing custom atom selections
==============================
MDAnalysis includes its own selection language for extracting
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` objects, which function as
ordered lists of (selected) atoms from the system. The selection strings needed to specify
these can be long and complex, and sometimes multiple selection strings are
required in a particular order to extract a given AtomGroup from all the atoms
in the Universe. Moreover, given different simulation systems, the same
selection of atoms (e.g. the "solvent") might require a different set of
selection strings.

Fortunately, Sims provide a mechanism for storing (many) atom selections.
Say we want to select the LID, CORE, and NMP domains of adenylate
kinase, the protein we simulated. We can store these immediately::

    >>> s.atomselections['lid'] = 'resid 122:159'
    >>> s.atomselections['nmp'] = 'resid 30:59'
    >>> s.atomselections['core'] = ('resid 1:29', 'resid 60:121', 'resid 160:214')

We can now get new AtomGroups back for each selection at any time with the 
:meth:`~mdsynthesis.limbs.AtomSelections.create` method::

    >>> s.atomselections.create('lid')
    <AtomGroup with 598 atoms>
    
    >>> s.atomselections.create('core')
    <AtomGroup with 2306 atoms>

and we don't have to remember or know how 'lid' or 'core' are defined for this
particular system. If we have other simulations of adenylate kinase performed
with other molecular dynamics engines or with different forcefields, we can
store the atom selection strings required for those systems in the same way,
perhaps using the same names 'lid', 'core', etc. This abstraction allows us to
work with many variants of a simulation system without having to micromanage.

.. note:: Storing a list of strings as a selection will apply them in order,
          producing an AtomGroup concatenated from each one in the same way
          as providing multiple strings to
          :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` does. This
          is especially useful when storing selections used for structural
          alignments.

Want just the selection strings back? We can use
:meth:`~mdsynthesis.limbs.AtomSelections.get`::

    >>> s.atomselections.get('lid')
    'resid 122:159'

    # or using getitem syntax
    >>> s.atomselections['lid']
    'resid 122:159'

Atom selections from atom indices 
---------------------------------
Do you already have an AtomGroup and prefer to define it according to its atom
indices instead of as a selection string? That can be done, too::

    >>> lid = s.universe.select_atoms('resid 122:159')
    >>> s.atomselections['lid'] = lid.indices
    >>> s.atomselections.create('lid')
    <AtomGroup with 598 atoms>

Lists/tuples of selection strings or atom indices can be stored in any combination
as a selection. These are applied in order to yield the AtomGroup when calling the 
:meth:`~mdsynthesis.limbs.AtomSelections.create` method.

API Reference: AtomSelections
-----------------------------
See the :ref:`AtomSelections_api` API reference for more details.

API Reference: Sim
==================
See the :ref:`Sim_api` API reference for more details.
