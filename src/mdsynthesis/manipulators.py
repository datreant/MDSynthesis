"""User-level functions for manipulating Sims.

"""
from datreant.core import discover as _discover

from .collections import Bundle

def discover(dirpath='.', depth=None, treantdepth=None):
    return Bundle(*_discover(dirpath=dirpath,
                            depth=depth,
                            treantdepth=treantdepth))
