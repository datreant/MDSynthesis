"""
Basic Treant objects: the organizational units for :mod:`mdsynthesis`.

"""
import os

from datreant.treants import Treant
from mdsynthesis import aggregators
from MDAnalysis import Universe


class Sim(Treant):
    """The Sim object is an interface to data for single simulations.

    """
    _treanttype = 'Sim'

    def __init__(self, sim, universe=None, uname='main', location='.',
                 coordinator=None, categories=None, tags=None):
        """Generate a new or regenerate an existing (on disk) Sim object.

        :Required arguments:
            *sim*
                if generating a new Sim, the desired name to give it;
                if regenerating an existing Sim, string giving the path
                to the directory containing the Sim object's state file

        :Optional arguments when generating a new Sim:
            *uname*
                desired name to associate with universe; this universe
                will be made the default (can always be changed later)
            *universe*
                arguments usually given to an MDAnalysis Universe
                that defines the topology and trajectory of the atoms
            *location*
                directory to place Sim object; default is the current directory
            *coordinator*
                directory of the Coordinator to associate with the Sim; if the
                Coordinator does not exist, it is created; if ``None``, the Sim
                will not associate with any Coordinator
            *categories*
                dictionary with user-defined keys and values; used to give Sims
                distinguishing characteristics
            *tags*
                list with user-defined values; like categories, but useful for
                adding many distinguishing descriptors

        *Note*: optional arguments are ignored when regenerating an existing
                Sim

        """
        if os.path.exists(sim):
            self._regenerate('Sim', sim)
        else:
            self._generate('Sim', sim, universe=universe, uname=uname,
                           location=location, coordinator=coordinator,
                           categories=categories, tags=tags)

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
            self._universes = aggregators.Universes(
                    self, self._backend, self._logger)

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
                self._selections = aggregators.Selections(
                        self, self._backend, self._logger)

            return self._selections

    def _generate(self, treanttype, treant, universe=None, uname='main',
                  location='.', coordinator=None, categories=None, tags=None):
        """Generate new Sim object.

        """
        super(Sim, self)._generate(treanttype, treant, location=location,
                                   coordinator=coordinator,
                                   categories=categories, tags=tags)

        # add universe
        if (uname and universe):
            self.universes.add(uname, *universe)

    def _placeholders(self):
        """Necessary placeholders for aggregator instances.

        """
        super(Sim, self)._placeholders()

        self._universes = None
        self._selections = None
        self._universe = None     # universe 'dock'
        self._uname = None        # attached universe name
