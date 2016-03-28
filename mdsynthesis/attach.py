"""
Limb attachments for :mod:`mdsynthesis` classes.

"""
from datreant.core.treants import Group
from datreant.core.collections import Bundle, View
from datreant.data.limbs import Data
from datreant.data.agglimbs import AggData

from .treants import Sim

Sim._attach_limb_class(Data)
Group._attach_limb_class(Data)
Bundle._attach_agglimb_class(AggData)
View._attach_aggtreelimb_class(AggData)
