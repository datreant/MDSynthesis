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
import warnings

from datreant.core.limbs import Limb
from MDAnalysis import Universe
from MDAnalysis.core.AtomGroup import AtomGroup

from . import filesystem


class Universes(Limb):
    """Manage the defined universes of the Sim.

    Universes are interfaces to raw simulation data. The Sim can store
    multiple universe definitions corresponding to different versions
    of the same simulation output (e.g. post-processed trajectories derived
    from the same raw trajectory). The Sim has at most one universe
    definition that is "active" at one time, with stored selections for
    this universe directly available via ``Sim.selections``.

    The Sim can also store a preference for a "default" universe, which is
    activated on a call to ``Sim.universe`` when no other universe is
    active.

    """
    _name = 'universes'
    _filepaths = ['abs', 'rel']

    def __init__(self, sim):
        super(Universes, self).__init__(sim)

        # init state if tags not already there;
        # if read-only, check that they are there,
        # and raise exception if they are not
        try:
            with self._treant._write:
                try:
                    self._treant._state['mds']['universes']
                except KeyError:
                    self._treant._state['mds']['universes'] = dict()
        except (IOError, OSError):
            with self._treant._read:
                try:
                    self._treant._state['mds']['universes']
                except KeyError:
                    raise KeyError(
                            ("Missing 'universes' data, and cannot write to "
                             "Sim '{}'".format(self._treant.filepath)))

    def __repr__(self):
        return "<Universes({})>".format(self.keys())

    def __str__(self):
        universes = self.keys()
        agg = "Universes"
        majsep = "="
        seplength = len(agg)

        if not universes:
            out = "No universes"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for universe in universes:
                out = out + "'{}'".format(universe)
                if self.default == universe:
                    out = out + ' (default)'
                if self._treant._uname == universe:
                    out = out + ' (active)'
                out = out + '\n'
        return out

    def __contains__(self, item):
        return item in self.keys()

    def __iter__(self):
        return self.keys().__iter__()

    def __getitem__(self, handle):
        """Attach universe and return a reference to it.

        :Arguments:
            *handle*
                given name for selecting the universe

        :Returns:
            *universe*
                a reference to the newly attached universe
        """
        self.activate(handle)

        return self._treant.universe

    def __setitem__(self, handle, universe):
        """Add a Universe object directly as a Universe definition.

        Does not at the moment store any special init keywords used to create
        the Universe originally.

        """
        top = universe.filename
        try:
            traj = universe.trajectory.filename
        except AttributeError:
            try:
                traj = universe.trajectory.filenames
            except AttributeError:
                traj = []

        self.add(handle, top, traj)

    def add(self, handle, topology, *trajectory):
        """Add a universe definition to the Sim object.

        A universe is an MDAnalysis object that gives access to the details
        of a simulation trajectory. A Sim object can contain multiple universe
        definitions (topology and trajectory pairs), since it is often
        convenient to have different post-processed versions of the same
        raw trajectory.

        Using an existing universe handle will replace the topology and
        trajectory for that definition; selections for that universe will be
        retained.

        If there is no current default universe, then the added universe will
        become the default.

        :Arguments:
            *handle*
                given name for selecting the universe
            *topology*
                path to the topology file
            *trajectory*
                path to the trajectory file; multiple files may be given
                and these will be used in order as frames for the trajectory

        """
        outtraj = []
        for traj in trajectory:
            if isinstance(traj, list):
                for t in traj:
                    outtraj.append(t)
            else:
                outtraj.append(traj)

        with self._treant._write:
            # if universe schema already exists, don't overwrite it
            if handle not in self._treant._state['mds']['universes']:
                self._treant._state['mds']['universes'][handle] = dict()

            udict = self._treant._state['mds']['universes'][handle]

            # add topology paths
            udict['top'] = [os.path.abspath(topology),
                            os.path.relpath(topology, self._treant.abspath)]

            # add trajectory paths
            udict['traj'] = list()
            for segment in trajectory:
                udict['traj'].append(
                        [os.path.abspath(segment),
                         os.path.relpath(segment, self._treant.abspath)])

            # add selections schema
            if 'sels' not in udict:
                udict['sels'] = dict()

            # add resnums schema
            if 'resnums' not in udict:
                udict['resnums'] = None

            # if no default universe, make this default
            if not self._treant._state['mds']['default']:
                self._treant._state['mds']['default'] = handle

            if not self.default:
                self.default = handle

    def remove(self, *handle):
        """Remove a universe definition.

        Also removes any selections associated with the universe.

        :Arguments:
            *handle*
                name of universe(s) to delete
        """
        with self._treant._write:
            for item in handle:
                try:
                    del self._treant._state['mds']['universes'][item]
                except KeyError:
                    raise KeyError(
                            "No such universe '{}';".format(handle) +
                            " nothing to remove.")

                if self._treant._uname == item:
                    self._treant._universe = None
                    self._treant._uname = None

                if self.default == item:
                    self.default = None

    def rename(self, handle, newname):
        """Rename a universe definition.

        :Arguments:
            *handle*
                name of universe to rename
            *newname*
                new name of universe
        """
        with self._treant._write:
            if newname in self._treant._state['mds']['universes']:
                raise ValueError("Universe '{}' already exixts; "
                                 "remove it first".format(newname))

            udicts = self._treant._state['mds']['universes']
            try:
                udicts[newname] = udicts.pop(handle)
            except KeyError:
                raise KeyError("No such universe '{}'; "
                               "nothing to rename.".format(handle))

        if self._treant._uname == handle:
            self._treant._uname = newname

            if self.default == handle:
                self.default = newname

    def keys(self):
        """Get handles for all universe definitions as a list.

        :Returns:
            *handles*
                list of all universe handles
        """
        with self._treant._read:
            return self._treant._state['mds']['universes'].keys()

    def activate(self, handle=None):
        """Make the selected universe active.

        Only one universe definition can be active in a Sim at one time. The
        active universe can be accessed from ``Sim.universe``. Stored
        selections for the active universe can be accessed as items in
        ``Sim.selections``.

        If no handle given, the default universe is loaded.

        If a resnum definition exists for the universe, it is applied.

        Raises :exc:`ValueError` if no handle given and no universe is set
        as the default.

        :Arguments:
            *handle*
                given name for selecting the universe; if ``None``, default
                universe selected
        """
        if not handle:
            handle = self.default

        if handle:
            uh = filesystem.Universehound(self, handle)
            paths = uh.fetch()
            topology = paths['top']
            trajectory = paths['traj']

            self._treant._universe = Universe(topology, *trajectory)
            self._treant._uname = handle
            self._apply_resnums()

            # update the universe definition; will automatically build current
            # path variants for each file
            # if read-only, move on
            try:
                self.add(handle, topology, *trajectory)
            except OSError:
                warnings.warn(
                    "Cannot update paths for universe '{}';".format(handle) +
                    " state file is read-only.")
        else:
            raise ValueError("No handle given, and no default universe set;" +
                             " no universe activated.")

    @property
    def current(self):
        """Name of the currently active universe.

        ``None`` indicates no universe is active.

        """
        return self._treant._uname

    def deactivate(self):
        """Deactivate the current universe.

        Deactivating the current universe may be necessary to conserve
        memory, since the universe can then be garbage collected.

        """
        self._treant._universe = None
        self._treant._uname = None

    def _apply_resnums(self):
        """Apply resnum definition to active universe.

        """
        with self._treant._read:
            udict = self._treant._state['mds']['universes'][self.current]
            resnums = udict['resnums']

        if resnums:
            self._treant._universe.residues.set_resnums(np.array(resnums))

    def resnums(self, handle, resnums):
        """Define resnums for the given universe.

        Resnums are useful for referring to residues by their canonical resid,
        for instance that stored in the PDB. By giving a resnum definition
        for the universe, this definition will be applied to the universe
        on activation.

        Will overwrite existing resnum definition if it exists.

        :Arguments:
            *handle*
                name of universe to apply resnums to
            *resnums*
                array or list giving the resnum for each residue in the
                topology, in atom index order; giving ``None`` will delete
                resnum definition
        """
        with self._treant._write:
            udict = self._treant._state['mds']['universes'][handle]
            if resnums is None:
                udict['resnums'] = None
            else:
                udict['resnums'] = list(resnums)

            if handle == self._treant._uname:
                self._apply_resnums()

    @property
    def default(self):
        """The default universe.

        Set to an existing universe handle to change.

        The default universe is loaded on calls to ``Sim.universe`` or
        ``Sim.selections`` when no other universe is attached.

        """
        with self._treant._read:
            return self._treant._state['mds']['default']

    @default.setter
    def default(self, value):
        with self._treant._write:
            if value is None or value in self:
                self._treant._state['mds']['default'] = value
            else:
                ValueError("'{}' not an existing "
                           "universe handle".format(value))

    def define(self, handle, pathtype='abs'):
        """Get the stored path to the topology and trajectory used for the
        specified universe.

        *Note*: Does no checking as to whether these paths are valid. To
                check this, try activating the universe.

        :Arguments:
            *handle*
                name of universe to get definition for

        :Keywords:
            *pathtype*
                type of path to return; 'abs' gives an absolute path,
                'rel' gives a path relative to the Sim's state file

        :Returns:
            *topology*
                path to the topology file
            *trajectory*
                list of paths to trajectory files
        """
        with self._treant._read:
            top = self._treant._state['mds']['universes'][handle]['top']
            outtop = {key: value for key, value in zip(self._filepaths, top)}

            trajs = self._treant._state['mds']['universes'][handle]['traj']
            outtraj = {key: [] for key in self._filepaths}

        for traj in trajs:
            for i, key in enumerate(self._filepaths):
                outtraj[key].append(traj[i])

        return outtop[pathtype], outtraj[pathtype]


class Selections(Limb):
    """Stored atom selections for the active universe.

    Useful atom selections can be stored for the active universe and
    recalled later. Selections are stored separately for each defined
    universe, since the same selection may require a different selection
    string for different universes.

    """
    _name = 'selections'

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
        """Get selection as an AtomGroup for given handle and the active universe.

        :Arguments:
            *handle*
                name of selection to return as an AtomGroup

        :Returns:
            *AtomGroup*
                the named selection as an AtomGroup of the active universe

        """
        return self.asAtomGroup(handle)

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

        :Arguments:
            *handle*
                name to use for the selection
            *selection*
                selection string or AtomGroup; multiple selections may be given
                and their order will be preserved, which is useful for e.g.
                structural alignments
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

            udicts = self._treant._state['mds']['universes']
            udicts[self._treant._uname]['sels'][handle] = outsel

    def remove(self, *handle):
        """Remove an atom selection for the attached universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        :Arguments:
            *handle*
                name of selection(s) to remove
        """
        with self._treant._write:
            udicts = self._treant._state['mds']['universes']
            for item in handle:
                try:
                    del udicts[self._treant._uname]['sels'][item]
                except KeyError:
                    raise KeyError("No such selection '{}'".format(item))

    def keys(self, universe=None):
        """Return a list of all selection handles.

        :Keywords:
            *universe*
                universe definition to list selection keys for; ``None``
                will give keys for currently active universe

        """
        with self._treant._read:
            udicts = self._treant._state['mds']['universes']

            if universe is None and self._treant._uname:
                return udicts[self._treant._uname]['sels'].keys()
            else:
                try:
                    return udict[universe]['sels'].keys()
                except KeyError:
                    raise KeyError("No such universe '{}'".format(universe))

    def asAtomGroup(self, handle):
        """Get AtomGroup from active universe from the given named selection.

        If named selection doesn't exist, :exc:`KeyError` raised.

        :Arguments:
            *handle*
                name of selection to return as an AtomGroup

        :Returns:
            *AtomGroup*
                the named selection as an AtomGroup of the active universe
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

    def define(self, handle, universe=None):
        """Get selection definition for given handle and the active universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        :Arguments:
            *handle*
                name of selection to get definition of

        :Keywords:
            *universe*
                universe definition to list selection definitions for; ``None``
                will give definitions for currently active universe

        :Returns:
            *definition*
                list of strings defining the atom selection
        """
        with self._treant._read:
            udicts = self._treant._state['mds']['universes']

            if universe is None and self._treant._uname:
                sels = udicts[self._treant._uname]['sels']
            else:
                try:
                    sels = udict[universe]['sels']
                except KeyError:
                    raise KeyError("No such universe '{}'".format(universe))

        try:
            selstring = sels[handle]
        except KeyError:
            raise KeyError("No such selection '{}'".format(handle))

        return selstring

    def copy(self, universe):
        """Copy defined selections of another universe to the active universe.

        :Arguments:
            *universe*
                name of universe definition to copy selections from
        """
        with self._treant._write:
            if self._treant._uname:
                try:
                    selections = self.keys(universe=universe)
                except KeyError:
                    raise KeyError("No such universe '{}';".format(universe) +
                                   " cannot copy selections.")

                for sel in selections:
                    seldef = udicts[universe]['sels'][sel]
                    self.add(sel, *seldef)
