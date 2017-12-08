import os

from datreant.core import Bundle as _Bundle
from datreant.core import View
from datreant.core.names import TREANTDIR_NAME

from .names import SIMDIR_NAME
from .treants import Sim

class Bundle(_Bundle):

    @staticmethod
    def _is_sim(treant):
        return os.path.exists(os.path.join(treant, TREANTDIR_NAME, SIMDIR_NAME))

    def _add(self, *sims):

        filtered = list()
        for sim in sims:
            if sim is None:
                pass
            elif isinstance(sim, (list, tuple, View, _Bundle)):
                self._add(*sim)
            elif self._is_sim(sim):
                filtered.append(Sim(os.path.abspath(sim)))

        super(Bundle, self)._add(*filtered)
