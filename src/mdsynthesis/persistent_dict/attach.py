"""
Modifications made to :mod:`datreant` classes on import of module.

"""

from datreant.core import Treant, Tree, Bundle, View, Group
from . import limbs
from . import agglimbs

Tree._attach_limb_class(limbs.Data)
Bundle._attach_agglimb_class(agglimbs.AggData)
View._attach_aggtreelimb_class(agglimbs.AggData)
