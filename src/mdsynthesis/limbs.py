"""
Limbs are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Sim) in
ways that are user-friendly.

In short, an Limb is designed to be user friendly on its own, but are
often used as components of a Treant.

"""
import os
from six import string_types
import numpy as np
from numpy.lib.utils import deprecate
import warnings

from datreant.core import Leaf
from datreant.core.metadata import Metadata
import MDAnalysis as mda


class UniverseDefinition(Metadata):
    """The defined universe of the Sim.

    Universes are interfaces to raw simulation data, with stored selections for
    this universe directly available via ``Sim.atomselections``.

    """
    # _name = 'universedef'
    # _filepaths = ['abs', 'rel']
    _statefilename = 'universedef.json'

    @staticmethod
    def _init_state(jsonfile):
        """Used solely for initializing JSONFile state for storing tag
        information.

        """
        jsonfile._state = {
            'topology': dict(),
            'trajectory': list(),
            'kwargs': dict()
        }

    @property
    def topology(self):
        """The topology file for this Sim's universe.

        To change the topology file for this Sim's universe, set the path to
        the file. Setting to ``None`` will disable the Sim's universe,
        but the trajectory paths will be retained.

        """
        with self._read:
            mdsdict = self._statefile._state
            topstate = mdsdict['topology']
            if not topstate:
                return None
            else:
                return topstate['abspath']

    @topology.setter
    def topology(self, path):
        if isinstance(path, string_types):
            path = path
        elif isinstance(path, Leaf):
            path = path.abspath
        elif path is None:
            pass
        else:
            raise TypeError("Path to topology must be a string or Leaf")

        self._set_topology(path)

        # Move into Sim
        # reset universe, if present
        # if self._universe:
        #     self._activate()

    def _set_topology(self, path):
        with self._write:
            mdsdict = self._statefile._state
            topstate = mdsdict['topology']

            if path is None:
                mdsdict['topology'] = dict()
            else:
                topstate['abspath'] = os.path.abspath(path)

    @property
    def trajectory(self):
        """The trajectory file for this Sim's universe.

        To change the trajectory file(s) for this Sim's Universe, set the path
        to the file. A single path will use a single trajectory file, while a
        list or tuple of paths will use each trajectory file in order for
        building the universe. ``None`` indicates that the Universe will not
        load a trajectory file at all (but the topology may have coordinates).

        """
        with self._read:
            mdsdict = self._statefile._state
            traj = mdsdict['trajectory']
            if not traj:
                return None
            elif len(traj) == 1:
                return traj[0][0]
            else:
                return tuple([t[0] for t in traj])

    @trajectory.setter
    def trajectory(self, path):
        if isinstance(path, string_types):
            trajs = [path]
        elif isinstance(path, Leaf):
            trajs = [path.abspath]
        elif isinstance(path, (list, tuple)):
            trajs = list(path)
        elif path is None:
            trajs = []
        else:
            raise TypeError("Path to topology must be a string, Leaf, or"
                            " list/tuple of paths.")

        self._set_trajectory(trajs)

    def _set_trajectory(self, trajs):
        with self._write:
            mdsdict = self._statefile._state
            mdsdict['trajectory'] = []
            trajstate = mdsdict['trajectory']

            for traj in trajs:
                trajstate.append([
                    os.path.abspath(traj),
                ])

    @property
    def kwargs(self):
        """The keyword arguments applied to the Sim's universe when building
        it.

        Set these with a dictionary of keywords and values to change them.
        Keywords must be strings and values must be strings, ints, floats,
        bools, or ``None``.

        """
        with self._read:
            return self._statefile._state['kwargs']

    @kwargs.setter
    def kwargs(self, kwargs):
        if kwargs is None:
            pass
        elif isinstance(kwargs, dict):
            # check that values are serializable
            for key, value in kwargs.items():
                if not (isinstance(value, (string_types, bool, int, float))
                        or value is None):
                    raise ValueError("Cannot store keyword '{}' for Universe; "
                                     "value must be a string, bool, int, "
                                     "float, or ``None``, "
                                     "not '{}'".format(key, type(value)))
        else:
            raise TypeError("Must be a dictionary or ``None``")

        with self._write:
            self._statefile._state['kwargs'] = kwargs

    @property
    def _args(self):
        """dict to generate a universe"""
        if self.topology is None:
            return None
        args = [self.topology, ]
        if self.trajectory is not None:
            args.append(self.trajectory)
        return args

    def _clear(self):
        self._set_topology(None)
        self._set_trajectory([])
        self.kwargs = None

    def update(self, universe):
        if universe is None:
            self._clear()
        elif not isinstance(universe, mda.Universe):
            raise TypeError("Cannot set to {}; must be Universe".format(
                type(universe)))
        else:
            self.topology = universe.filename
            try:  # ChainReader?
                traj = universe.trajectory.filenames
            except AttributeError:
                try:  # Reader?
                    traj = [universe.trajectory.filename]
                except AttributeError:  # Only topology
                    traj = []
            self.trajectory = traj
            self.kwargs = universe.kwargs


class AtomSelections(Metadata):
    """Stored atom selections for the universe.

    Useful atom selections can be stored for the universe and recalled later.

    """
    _statefilename = 'atomselections.json'

    @staticmethod
    def _init_state(jsonfile):
        """Used solely for initializing JSONFile state for storing tag
        information.

        """
        jsonfile._state = {}

    def __repr__(self):
        return "<AtomSelections({})>".format(
            {x: self.get(x)
             for x in self.keys()})

    def __getitem__(self, handle):
        """Get selection for given handle.

        Parameters
        ----------
        handle : str
            Name of selection to return.

        Returns
        -------
        selection : str, list, array_like
            The named selection definition.

        """
        return self.get(handle)

    def __setitem__(self, handle, selection):
        """Selection for the given handle.

        """
        if isinstance(selection, (string_types, np.ndarray)):
            selection = [selection]
        self.add(handle, *selection)

    def __iter__(self):
        return self.keys().__iter__()

    def __delitem__(self, handle):
        """Remove stored selection for given handle and the active universe.

        """
        self.remove(handle)

    def add(self, handle, *selection):
        """Add an atom selection for the attached universe.

        AtomGroups are needed to obtain useful information from raw coordinate
        data. It is useful to store AtomGroup selections for later use, since
        they can be complex and atom order may matter.

        If a selection with the given `handle` already exists, it is replaced.

        Parameters
        ----------
        handle : str
            Name to use for the selection.
        selection : str, AtomGroup
            Selection string or AtomGroup indices; multiple selections may be
            given and their order will be preserved, which is useful for e.g.
            structural alignments.

        """
        if len(selection) == 1:
            sel = selection[0]
            if isinstance(sel, np.ndarray):
                outsel = sel.tolist()
            elif isinstance(sel, string_types):
                outsel = sel
        else:
            outsel = list()
            for sel in selection:
                if isinstance(sel, np.ndarray):
                    outsel.append(sel.tolist())
                elif isinstance(sel, string_types):
                    outsel.append(sel)
                else:
                    raise ValueError("Selections must be strings, arrays of "
                                     "atom indices, or tuples/lists of these.")

        with self._write:
            seldict = self._statefile._state
            seldict[handle] = outsel

    def remove(self, *handle):
        """Remove an atom selection for the universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection(s) to remove.

        """
        with self._write:
            seldict = self._statefile._state
            for item in handle:
                try:
                    del seldict[item]
                except KeyError:
                    raise KeyError("No such selection '{}'".format(item))

    def keys(self):
        """Return a list of all selection handles.

        """
        with self._read:
            seldict = self._statefile._state
            return seldict.keys()

    def create(self, handle):
        """Generate AtomGroup from universe from the given named selection.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection to return as an AtomGroup.

        Returns
        -------
        AtomGroup
            The named selection as an AtomGroup of the universe.

        """
        pass

    #     sel = self.get(handle)

    #     # Selections might be either
    #     # - a single string
    #     # - a single list of indices
    #     # - a list of strings
    #     # - a list of indices

    #     if isinstance(sel, string_types):
    #         # if we have a single string
    #         ag = self._treant.universe.select_atoms(sel)
    #     elif all([isinstance(i, int) for i in sel]):
    #         # if we have a single array_like of indices
    #         ag = self._treant.universe.atoms[sel]
    #     else:
    #         ag = None
    #         for item in sel:
    #             if isinstance(item, string_types):
    #                 if ag:
    #                     ag += self._treant.universe.select_atoms(item)
    #                 else:
    #                     ag = self._treant.universe.select_atoms(item)
    #             else:
    #                 if ag:
    #                     ag += self._treant.universe.atoms[item]
    #                 else:
    #                     ag = self._treant.universe.atoms[item]

    #     return ag

    def get(self, handle):
        """Get selection definition for given handle.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection to get definition of.

        Returns
        -------
        definition : list
            list of strings defining the atom selection

        """
        pass
        # with self._treant._read:
        #     seldict = self._treant._state['mdsynthesis']['atomselections']

        # try:
        #     seldef = seldict[handle]
        # except KeyError:
        #     raise KeyError("No such selection '{}'".format(handle))

        # if isinstance(seldef, string_types):
        #     # if we have a single string
        #     out = seldef
        # elif all([isinstance(i, int) for i in seldef]):
        #     # if we have a single list of indices
        #     out = np.array(seldef)
        # else:
        #     out = []
        #     for item in seldef:
        #         if isinstance(item, string_types):
        #             out.append(item)
        #         else:
        #             out.append(np.array(item))

        #     out = tuple(out)

        # return out
