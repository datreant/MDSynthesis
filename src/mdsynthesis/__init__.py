# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# mdsynthesis

"""
MDSynthesis --- a persistence engine for molecular dynamics data
================================================================
"""
# Bring some often used objects into the current namespace
from datreant import Tree, Leaf, View, Bundle
from datreant import discover

from .treants import Sim
from .manipulators import discover

__all__ = ['Sim', 'Bundle', 'discover', 'Tree', 'Leaf', 'View']
__version__ = "0.6.2-dev"  # NOTE: keep in sync with RELEASE in setup.py
