"""
Limb attachments for :mod:`mdsynthesis` classes.

"""
from datreant.core.treants import Group
from datreant.core.collections import CollectionBase
from datreant.data.limbs import Data
from datreant.data.agglimbs import MemberData

from .treants import Sim

Sim._attach_limb_class(Data)
Group._attach_limb_class(Data)
CollectionBase._attach_agglimb_class(MemberData)
