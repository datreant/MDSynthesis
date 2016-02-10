"""
Limb attachments for :mod:`mdsynthesis` classes.

"""

from .treants import Sim
from datreant.core.collections import CollectionBase
from datreant.data.limbs import Data
from datreant.data.agglimbs import MemberData

Sim._attach_limb_class(Data)
CollectionBase._attach_agglimb_class(MemberData)
