# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# mdsynthesis

"""
MDSynthesis --- a persistence engine for molecular dynamics data
================================================================
"""
# Bring some often used objects into the current namespace
from datreant.core import Tree, Leaf, View
from datreant.core import discover

from .treants import Sim
from .collections import Bundle
from .manipulators import discover

__all__ = ['Sim', 'Bundle', 'discover', 'Tree', 'Leaf', 'View']
__version__ = "0.6.2-dev"  # NOTE: keep in sync with RELEASE in setup.py
