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
from datreant.core.limbs import Limb
from MDAnalysis import Universe
from MDAnalysis.core.AtomGroup import AtomGroup

from .filesystem import Universehound


class UniverseDefinition(Limb):
    """The defined universe of the Sim.

    Universes are interfaces to raw simulation data, with stored selections for
    this universe directly available via ``Sim.atomselections``.

    """
    _name = 'udef'
    _filepaths = ['abs', 'rel']

    def __init__(self, treant):
        super(UniverseDefinition, self).__init__(treant)

        subitems = {'topology': dict,
                    'trajectory': list,
                    'universe_kwargs': dict}

        # init state if for udef not already there;
        # if read-only, check that it is there,
        # and raise exception if it is not
        try:
            with self._treant._write:
                try:
                    self._treant._state['mdsynthesis'][self._name]
                except KeyError:
                    self._treant._state['mdsynthesis'][self._name] = dict()

                for key, value in subitems.items():
                    try:
                        self._treant._state['mdsynthesis'][self._name][key]
                    except KeyError:
                        self._treant._state['mdsynthesis'][self._name][key] = value()

        except (IOError, OSError):
            with self._treant._read:
                try:
                    self._treant._state['mdsynthesis'][self._name]
                except KeyError:
                    raise KeyError(
                            ("Missing '{}' data, and cannot write to "
                             "Treant '{}'".format(self._name,
                                                  self._treant.filepath)))

    def _activate(self):
        """Make the universe and attach it.

        """
        if not self.topology:
            self._treant._universe = None
            return

        uh = Universehound(self)
        paths = uh.fetch()
        topology = paths['top']
        trajectory = paths['traj']

        if not trajectory:
            self._treant._universe = Universe(topology, **self.universe_kwargs)
        else:
            self._treant._universe = Universe(topology, *trajectory,
                                              **self.universe_kwargs)

        self._apply_resnums()

        # update the universe definition; will automatically build current
        # path variants for each file
        # if read-only, move on
        try:
            self._set_topology(topology)
            self._set_trajectory(trajectory)
        except OSError:
            warnings.warn(
                "Cannot update paths for universe; "
                " state file is read-only.")

    def reload_universe(self):
        """Re-load the universe from its stored definition.

        """
        self._activate()

    @property
    def topology(self):
        """The topology file for this Sim's universe.

        To change the topology file for this Sim's universe, set the path to
        the file. Setting to ``None`` will disable the Sim's universe,
        but the trajectory paths will be retained.

        """
        with self._treant._read:
            topstate = self._treant._state['mdsynthesis']['udef']['topology']

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

        # reset universe, if present
        if self._treant._universe:
            self._activate()

    def _set_topology(self, path):
        with self._treant._write:
            topstate = self._treant._state['mdsynthesis']['udef']['topology']

            if path is None:
                self._treant._state['mdsynthesis']['udef']['topology'] = dict()
            else:
                topstate['abspath'] = os.path.abspath(path)
                topstate['relpath'] = os.path.relpath(path,
                                                      self._treant.abspath)

    @property
    def trajectory(self):
        """The trajectory file for this Sim's universe.

        To change the trajectory file(s) for this Sim's Universe, set the path
        to the file. A single path will use a single trajectory file, while a
        list or tuple of paths will use each trajectory file in order for
        building the universe. ``None`` indicates that the Universe will not
        load a trajectory file at all (but the topology may have coordinates).

        """
        with self._treant._read:
            traj = self._treant._state['mdsynthesis']['udef']['trajectory']
            if not traj:
                return None
            elif len(traj) == 1:
                return traj[0][0]
            else:
                return [t[0] for t in traj]

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

        # reset universe, if present
        if self._treant._universe:
            self._activate()

    def _set_trajectory(self, trajs):
        with self._treant._write:
            self._treant._state['mdsynthesis']['udef']['trajectory'] = []
            trajstate = self._treant._state['mdsynthesis']['udef']['trajectory']

            for traj in trajs:
                trajstate.append(
                        [os.path.abspath(traj),
                         os.path.relpath(traj, self._treant.abspath)])

    def _apply_resnums(self):
        """Apply resnum definition to universe.

        """
        with self._treant._read:
            simdict = self._treant._state['mdsynthesis']['udef']
            try:
                resnums = simdict['resnums']
            except KeyError:
                resnums = None

        if resnums:
            self._treant._universe.residues.set_resnums(np.array(resnums))

    @deprecate(message="resnum storage is deprecated")
    def _set_resnums(self, resnums):
        """Define resnums for the universe.

        Resnums are useful for referring to residues by their canonical resid,
        for instance that stored in the PDB. By giving a resnum definition
        for the universe, this definition will be applied to the universe.

        Will overwrite existing resnum definition if it exists.

        Parameters
        ----------
        resnums : array_like
                array or list giving the resnum for each residue in the
                topology, in atom index order; giving ``None`` will delete
                resnum definition
        """
        with self._treant._write:
            simdict = self._treant._state['mdsynthesis']['udef']
            if resnums is None:
                simdict['resnums'] = None
            else:
                simdict['resnums'] = list(resnums)

            if self._treant._universe:
                self._apply_resnums()

    def _define(self, pathtype='abs'):
        """Get the stored path to the topology and trajectory used for the
        universe.

        .. note:: Does no checking as to whether these paths are valid. To
                  check this, try using the universe.

        Parameters
        ----------
        pathtype : {'abs', 'rel'}
            type of path to return; 'abs' gives an absolute path, 'rel' gives a
            path relative to the Sim's state file

        Returns
        -------
        topology : str
            path to the topology file
        trajectory : list
            list of paths to trajectory files

        """
        with self._treant._read:
            top = self._treant._state['mdsynthesis']['udef']['topology']
            outtop = {'abs': top['abspath'],
                      'rel': top['relpath']}

            trajs = self._treant._state['mdsynthesis']['udef']['trajectory']
            outtraj = {key: [] for key in self._filepaths}

        for traj in trajs:
            for i, key in enumerate(self._filepaths):
                outtraj[key].append(traj[i])

        return outtop[pathtype], outtraj[pathtype]

    @property
    def universe_kwargs(self):
        """The keyword arguments applied to the Sim's universe when building
        it.

        Set these with a dictionary of keywords and values to change them.
        Keywords must be strings and values must be strings, ints, floats,
        bools, or ``None``.

        """
        with self._treant._read:
            return self._treant._state['mdsynthesis']['udef']['universe_kwargs']

    @universe_kwargs.setter
    def universe_kwargs(self, kwargs):
        if not isinstance(kwargs, dict):
            raise TypeError("Must be a dictionary")

        with self._treant._write:
            simdict = self._treant._state['mdsynthesis']['udef']
            simdict['universe_kwargs'] = kwargs


class AtomSelections(Limb):
    """Stored atom selections for the universe.

    Useful atom selections can be stored for the universe and recalled later.

    """
    _name = 'atomselections'

    def __init__(self, treant):
        super(AtomSelections, self).__init__(treant)

        # init state if for selections not already there;
        # if read-only, check that it is there,
        # and raise exception if it is not
        try:
            with self._treant._write:
                try:
                    self._treant._state['mdsynthesis'][self._name]
                except KeyError:
                    self._treant._state['mdsynthesis'][self._name] = dict()
        except (IOError, OSError):
            with self._treant._read:
                try:
                    self._treant._state['mdsynthesis'][self._name]
                except KeyError:
                    raise KeyError(
                            ("Missing '{}' data, and cannot write to "
                             "Treant '{}'".format(self._name,
                                                  self._treant.filepath)))

    def __repr__(self):
        return "<Selections({})>".format(
                {x: self.define(x) for x in self.keys()})

    def __str__(self):
        selections = self.keys()
        agg = "Selections"
        majsep = "="
        minsep = "-"
        subsep = "| "
        seplength = len(agg)

        if not self._treant._uname:
            out = "No universe attached; no Selections to show"
        elif not selections:
            out = "No selections for universe '{}'".format(
                self._treant._uname)
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for selection in selections:
                out = out + "'{}'\n".format(selection)
                for item in self.define(selection):
                    out = out + subsep + "'{}'\n".format(item)
                out = out + minsep * seplength + '\n'

        return out

    def __getitem__(self, handle):
        """Get selection as an AtomGroup for given handle and the universe.

        Parameters
        ----------
        handle : str
            Name of selection to return as an AtomGroup.

        Returns
        -------
        AtomGroup
            The named selection as an AtomGroup of the universe.

        """
        return self._asAtomGroup(handle)

    def __setitem__(self, handle, selection):
        """Selection for the given handle and the active universe.

        """
        if isinstance(selection, (string_types, AtomGroup)):
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

        If a selection with the given *handle* already exists, it is replaced.

        Parameters
        ----------
        handle : str
            Name to use for the selection.
        selection : str, AtomGroup
            Selection string or AtomGroup; multiple selections may be given and
            their order will be preserved, which is useful for e.g. structural
            alignments.
        """
        # Conversion function, leave strings alone,
        # turn AtomGroups into their indices
        def conv(x):
            return x if isinstance(x, string_types) else x.indices

        selection = map(conv, selection)

        with self._treant._write:
            outsel = list()
            for sel in selection:
                if isinstance(sel, np.ndarray):
                    outsel.append(sel.tolist())
                elif isinstance(sel, string_types):
                    outsel.append(sel)

            seldict = self._treant._state['mdsynthesis']['atomselections']
            seldict[handle] = outsel

    def remove(self, *handle):
        """Remove an atom selection for the universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection(s) to remove.

        """
        with self._treant._write:
            seldict = self._treant._state['mdsynthesis']['atomselections']
            for item in handle:
                try:
                    del seldict[item]
                except KeyError:
                    raise KeyError("No such selection '{}'".format(item))

    def keys(self):
        """Return a list of all selection handles.

        """
        with self._treant._read:
            seldict = self._treant._state['mdsynthesis']['atomselections']
            return seldict.keys()

    def _asAtomGroup(self, handle):
        """Get AtomGroup from universe from the given named selection.

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
        sel = self.define(handle)

        # Selections might be either
        # - a list of strings
        # - a list of indices

        ag = None
        for item in sel:
            if isinstance(item, string_types):
                if ag:
                    ag += self._treant.universe.select_atoms(item)
                else:
                    ag = self._treant.universe.select_atoms(item)
            else:
                if ag:
                    ag += self._treant.universe.atoms[item]
                else:
                    ag = self._treant.universe.atoms[item]

        return ag

    def define(self, handle):
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
        with self._treant._read:
            seldict = self._treant._state['mdsynthesis']['atomselections']

        try:
            selstring = seldict[handle]
        except KeyError:
            raise KeyError("No such selection '{}'".format(handle))

        return selstring
