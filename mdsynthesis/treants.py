"""
Basic Treant objects: the organizational units for :mod:`mdsynthesis`.

"""
import os

from datreant.core.treants import Treant
from . import limbs
from .backends import statefiles
from MDAnalysis import Universe


class Sim(Treant):
    """The Sim object is an interface to data for single simulations.

    """
    _treanttype = 'Sim'
    _backendclass = statefiles.SimFile

    def __init__(self, treant, new=False, categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) Treant object.

        :Required arguments:
            *treant*
                base directory of a new or existing Treant; will regenerate
                a Treant if a state file is found, but will genereate a new
                one otherwise

                if multiple Treant state files are in the given directory,
                will raise :exception:`MultipleTreantsError`; specify
                the full path to the desired state file to regenerate the
                desired Treant in this case

                use the *new* keyword to force generation of a new Treant
                at the given path

        :Optional arguments:
            *new*
                generate a new Treant even if one already exists at the given
                location *treant*
            *categories*
                dictionary with user-defined keys and values; used to give
                Treants distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors
        """
        super(Sim, self).__init__(treant,
                                  new=new,
                                  categories=categories,
                                  tags=tags)

        self._universes = None
        self._selections = None
        self._universe = None     # universe 'dock'
        self._uname = None        # attached universe name

    def __repr__(self):
        if not self._uname:
            out = "<Sim: '{}'>".format(self.name)
        else:
            out = "<Sim: '{}' | active universe: '{}'>".format(self.name,
                                                               self._uname)

        return out

    @property
    def universe(self):
        """The active universe of the Sim.

        Universes are interfaces to raw simulation data. The Sim can store
        multiple universe definitions corresponding to different versions
        of the same simulation output (e.g. post-processed trajectories derived
        from the same raw trajectory). The Sim has at most one universe
        definition that is "active" at one time, with stored selections for
        this universe directly available via ``Sim.selections``.

        To have more than one universe available as "active" at the same time,
        generate as many instances of the Sim object from the same statefile on
        disk as needed, and make a universe active for each one.

        """
        # TODO: include check for changes to universe definition, not just
        # definition absence
        if self._uname in self._backend.list_universes():
            return self._universe
        elif not self._universe:
            self.universes.activate()
            return self._universe
        else:
            self._universe = None
            self._logger.info('This universe is no longer defined. '
                              'It has been detached.')

    @property
    def universes(self):
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
        if not self._universes:
            self._universes = limbs.Universes(self)

        return self._universes

    @property
    def selections(self):
        """Stored atom selections for the active universe.

        Useful atom selections can be stored for the active universe and
        recalled later. Selections are stored separately for each defined
        universe, since the same selection may require a different selection
        string for different universes.

        """
        # attach default universe if not attached, and only give results if a
        # universe is present thereafter
        if self.universe:
            if not self._selections:
                self._selections = limbs.Selections(self)

            return self._selections
