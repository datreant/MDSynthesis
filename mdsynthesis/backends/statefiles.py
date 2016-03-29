"""
Interface classes for state files.

"""
from datreant.core.backends.statefiles import TreantFile


class SimFile(TreantFile):

    def _init_state(self):
        super(SimFile, self)._init_state()
        self._state['mdsynthesis'] = dict()
        self._state['mdsynthesis']['sim'] = dict()
        self._state['mdsynthesis']['sim']['topology'] = dict()
        self._state['mdsynthesis']['sim']['trajectory'] = list()
        self._state['mdsynthesis']['sim']['selections'] = dict()
        self._state['mdsynthesis']['sim']['universe_kwargs'] = dict()
