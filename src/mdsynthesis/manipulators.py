"""User-level functions for manipulating Sims.

"""
import os

from datreant import discover as _discover
from datreant import Bundle
from datreant.names import TREANTDIR_NAME

from .treants import Sim
from .names import SIMDIR_NAME


def _is_sim(treant):
    return os.path.exists(os.path.join(treant.abspath,
                                       TREANTDIR_NAME,
                                       SIMDIR_NAME))


def discover(dirpath='.', depth=None, treantdepth=None):
    treants = _discover(dirpath=dirpath,
                        depth=depth,
                        treantdepth=treantdepth)

    return Bundle([Sim(treant) for treant in treants if _is_sim(treant)])
