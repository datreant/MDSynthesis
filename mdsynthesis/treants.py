"""
Basic Treant objects: the organizational units for :mod:`mdsynthesis`.

"""
import warnings
import os
from six import string_types

import numpy as np
from MDAnalysis import Universe
from numpy.lib.utils import deprecate

from datreant.core import Treant, Leaf
from . import limbs
from .backends import statefiles
from .filesystem import Universehound


class Sim(Treant):
    """The Sim object is an interface to data for a single simulation.

    `sim` should be a base directory of a new or existing Sim. An existing
    Sim will be regenerated if a state file is found.  If no state file is
    found, a new Sim will be created.

    A Tree object may also be used in the same way as a directory string.

    If multiple Treant state files are in the given directory, a
    :exc:`MultipleTreantsError` will be raised; specify the full path to the
    desired state file to regenerate the desired Treant in this case. It is
    generally better to avoid having multiple state files in the same
    directory.

    Use the `new` keyword to force generation of a new Sim at the given path.

    Parameters
    ----------
    sim : str or Tree
        Base directory of a new or existing Sim; will regenerate a Sim if a
        state file is found, but will genereate a new one otherwise; may also
        be a Tree object
    new : bool
        Generate a new Sim even if one already exists at the given location
    categories : dict
        dictionary with user-defined keys and values; used to give Sims
        distinguishing characteristics
    tags : list
        list with user-defined values; like categories, but useful for adding
        many distinguishing descriptors
    """
    _treanttype = 'Sim'
    _backendclass = statefiles.SimFile

    def __init__(self, sim, new=False, categories=None, tags=None):
        super(Sim, self).__init__(sim,
                                  new=new,
                                  categories=categories,
                                  tags=tags)

        self._selections = None
        self._universe = None     # universe 'dock'

    def __repr__(self):
        return "<{}: '{}'>".format(self._treanttype, self.name)

    @property
    def universe(self):
        """The universe of the Sim.

        Universes are interfaces to raw simulation data, with stored selections
        for this universe directly available via ``Sim.selections``.

        Setting this to a :class:`MDAnalysis.Universe` will set that as the
        universe definition for this Sim. It will not preserve any keyword
        arguments used to initialize it, however; you will need to add
        these to :attr:`universe_kwargs` as a dictionary if you want these
        to apply next time the universe is loaded.

        """
        # TODO: include check for changes to universe definition, not just
        # definition absence
        if self._universe:
            return self._universe
        else:
            self._activate()
            return self._universe

    @universe.setter
    def universe(self, universe):
        if not isinstance(universe, Universe):
            raise TypeError("Cannot set to {}; must be Universe".format(
                                type(universe)))

        self.topology = universe.filename
        try:
            traj = universe.trajectory.filename
        except AttributeError:
            try:
                traj = universe.trajectory.filenames
            except AttributeError:
                traj = None

        self.trajectory = traj

        # finally, just use this instance
        self._universe = universe

    def _activate(self):
        """Make the universe and attach it.

        """
        if not self.topology:
            self._universe = None
            return

        uh = Universehound(self)
        paths = uh.fetch()
        topology = paths['top']
        trajectory = paths['traj']

        if not trajectory:
            self._universe = Universe(topology, **self.universe_kwargs)
        else:
            self._universe = Universe(topology, *trajectory,
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
        with self._read:
            topstate = self._state['mdsynthesis']['sim']['topology']

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
        if self._universe:
            self._activate()

    def _set_topology(self, path):
        with self._write:
            topstate = self._state['mdsynthesis']['sim']['topology']

            if path is None:
                topstate = dict()
            else:
                topstate['abspath'] = os.path.abspath(path)
                topstate['relpath'] = os.path.relpath(path,
                                                      self.abspath)

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
            traj = self._state['mdsynthesis']['sim']['trajectory']
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
        if self._universe:
            self._activate()

    def _set_trajectory(self, trajs):
        with self._write:
            self._state['mdsynthesis']['sim']['trajectory'] = []
            trajstate = self._state['mdsynthesis']['sim']['trajectory']

            for traj in trajs:
                trajstate.append(
                        [os.path.abspath(traj),
                         os.path.relpath(traj, self.abspath)])

    def _apply_resnums(self):
        """Apply resnum definition to universe.

        """
        with self._read:
            simdict = self._state['mdsynthesis']['sim']
            try:
                resnums = simdict['resnums']
            except KeyError:
                resnums = None

        if resnums:
            self._universe.residues.set_resnums(np.array(resnums))

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
        with self._write:
            simdict = self._state['mdsynthesis']['sim']
            if resnums is None:
                simdict['resnums'] = None
            else:
                simdict['resnums'] = list(resnums)

            if self._universe:
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
        filepaths = ('abs', 'rel')
        with self._read:
            top = self._state['mdsynthesis']['sim']['topology']
            outtop = {'abs': top['abspath'],
                      'rel': top['relpath']}

            trajs = list(self._state['mdsynthesis']['sim']['trajectory'])
            outtraj = {key: [] for key in filepaths}

        for traj in trajs:
            for i, key in enumerate(filepaths):
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
        with self._read:
            return self._state['mdsynthesis']['sim']['universe_kwargs']

    @universe_kwargs.setter
    def universe_kwargs(self, kwargs):
        if not isinstance(kwargs, dict):
            raise TypeError("Must be a dictionary")

        with self._write:
            self._state['mdsynthesis']['sim']['universe_kwargs'] = kwargs

    @property
    def selections(self):
        """Stored atom selections for the universe.

        Useful atom selections can be stored for the universe and
        recalled later. Selections are stored separately for each defined
        universe, since the same selection may require a different selection
        string for different universes.

        """
        # attach universe if not attached, and only give results if a
        # universe is present thereafter
        if self.universe:
            if not self._selections:
                self._selections = limbs.Selections(self)

            return self._selections
