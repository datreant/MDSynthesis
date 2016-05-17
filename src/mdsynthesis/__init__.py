# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# mdsynthesis

"""
MDSynthesis --- a persistence engine for molecular dynamics data
================================================================
"""
# Bring some often used objects into the current namespace
from datreant.core import Treant, Group, Bundle, Tree, Leaf, View
from datreant.core import discover

from .treants import Sim
from . import attach

__all__ = ['Sim', 'Group', 'Bundle']
__version__ = "0.6.1"  # NOTE: keep in sync with RELEASE in setup.py
