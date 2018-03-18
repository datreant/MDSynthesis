"""
Basic Treant objects: the organizational units for :mod:`mdsynthesis`.

"""
from __future__ import absolute_import

import warnings
import os
from functools import wraps

import MDAnalysis as mda

from datreant import Treant
from datreant.names import TREANTDIR_NAME
from datreant.util import makedirs
from .names import SIMDIR_NAME
from . import metadata
from .data import Data


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
    categories : dict
        dictionary with user-defined keys and values; used to give Sims
        distinguishing characteristics
    tags : list
        list with user-defined values; like categories, but useful for adding
        many distinguishing descriptors
    """
    _treanttype = 'Sim'

    def __init__(self, sim, categories=None, tags=None):
        super(Sim, self).__init__(sim,
                                  categories=categories,
                                  tags=tags)

        self._universedef = metadata.UniverseDefinition(self)
        self._universe = None
        self._args = None
        self._atomselections = metadata.AtomSelections(self, parent=self)

        # make simdir
        self._make_simdir()

        # make simdir
        self._make_simdir()
        self._data = Data(self)

    def __repr__(self):
        return "<{}: '{}'>".format(self._treanttype, self.name)

    def _make_simdir(self):
        abspath = self._path.absolute()
        simdir = abspath / os.path.join(TREANTDIR_NAME, SIMDIR_NAME)

        if not simdir.exists():
            # build mdsynthesis dir; stop if we hit a permissions error
            try:
                makedirs(simdir, exist_ok=True)
            except OSError as e:
                if e.errno == 13:
                    raise OSError(13, "Permission denied; " +
                                  "cannot create '{}'".format(simdir))
                else:
                    raise

    @property
    def _simdir(self):
        return os.path.join(self.abspath, TREANTDIR_NAME, SIMDIR_NAME)

    @property
    def universe(self):
        """The universe of the Sim.

        Universes are interfaces to raw simulation data, with stored selections
        for this universe directly available via ``Sim.selections``.

        Setting this to a :class:`MDAnalysis.Universe` will set that as the
        universe definition for this Sim. Setting to ``None`` will remove
        the universe definition entirely.

        """
        _args = self.universedef._args
        if _args != self._args:
            self._args = _args
            kwargs = self.universedef.kwargs
            if _args is None:
                self._universe = None
            else:
                self._universe = mda.Universe(*_args, **kwargs)
        return self._universe

    @universe.setter
    def universe(self, universe):
        self.universedef.update(universe)
        self._universe = universe

    @property
    def universedef(self):
        """The universe definition for this Sim.

        """
        return self._universedef

    @property
    def atomselections(self):
        """Stored atom selections for the universe.

        Useful atom selections can be stored for the universe and
        recalled later.
        """
        return self._atomselections

    @property
    def data(self):
        return self._data
