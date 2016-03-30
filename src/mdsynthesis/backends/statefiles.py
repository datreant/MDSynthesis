"""
Interface classes for state files.

"""
from datreant.core.backends.statefiles import TreantFile


class SimFile(TreantFile):

    def _init_state(self):
        super(SimFile, self)._init_state()
        self._state['mdsynthesis'] = dict()
