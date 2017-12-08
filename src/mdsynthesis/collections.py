import os

from datreant.core import Bundle as _Bundle
from datreant.core.names import TREANTDIR_NAME

from .names import SIMDIR_NAME
from .treants import Sim

class Bundle(_Bundle):

    @staticmethod
    def _is_sim(treant):
        return os.path.exists(os.path.join(treant, TREANTDIR_NAME, SIMDIR_NAME))

    def _add(self, *sims):

        # we only add Treants that are actually Sims
        filtered = [Sim(sim) for sim in sims if self._is_sim(sim)]
        super(Bundle, self)._add(*filtered)
