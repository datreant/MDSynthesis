"""
Basic Treant objects: the organizational units for :mod:`mdsynthesis`.

"""
import warnings
import os
from six import string_types

from MDAnalysis import Universe

from datreant.core import Treant, Leaf
from . import limbs
from .backends import statefiles


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

        self._universedef = None
        self._atomselections = None
        self._universe = None     # universe 'dock'

    def __repr__(self):
        return "<{}: '{}'>".format(self._treanttype, self.name)

    @property
    def universe(self):
        """The universe of the Sim.

        Universes are interfaces to raw simulation data, with stored selections
        for this universe directly available via ``Sim.selections``.

        Setting this to a :class:`MDAnalysis.Universe` will set that as the
        universe definition for this Sim. Setting to ``None`` will remove
        the universe definition entirely.

        """
        # TODO: include check for changes to universe definition, not just
        # definition absence
        if self._universe:
            return self._universe
        else:
            self.universedef._activate()
            return self._universe

    @universe.setter
    def universe(self, universe):
        if universe is None:
            self.universedef._set_topology(None)
            self.universedef._set_trajectory([])
            self.universedef.kwargs = None
            self.universedef.activate()

        elif not isinstance(universe, Universe):
            raise TypeError("Cannot set to {}; must be Universe".format(
                                type(universe)))
        else:
            self.universedef._set_topology(universe.filename)
            try:  # ChainReader?
                traj = universe.trajectory.filenames
            except AttributeError:
                try:  # Reader?
                    traj = [universe.trajectory.filename]
                except AttributeError:  # Only topology
                    traj = []

            self.universedef._set_trajectory(traj)

            # try and store keyword arguments
            try:
                self.universedef.kwargs = universe.kwargs
            except AttributeError:
                warnings.warn("Universe did not keep keyword arguments; "
                              "cannot store keyword arguments for Universe.")

            # finally, just use this instance
            self._universe = universe

    @universe.deleter
    def universe(self):
        self.universedef._set_topology(None)
        self.universedef._set_trajectory([])
        self.universedef.kwargs = None
        self.universedef.activate()

    @property
    def universedef(self):
        """The universe definition for this Sim.

        """
        # attach universe if not attached, and only give results if a
        # universe is present thereafter
        if not self._universedef:
            self._universedef = limbs.UniverseDefinition(self)

        return self._universedef

    @property
    def atomselections(self):
        """Stored atom selections for the universe.

        Useful atom selections can be stored for the universe and
        recalled later.
        """
        if not self._atomselections:
            self._atomselections = limbs.AtomSelections(self)

        return self._atomselections
