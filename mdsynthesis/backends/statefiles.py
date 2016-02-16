"""
Interface classes for state files.

"""

import os
import sys
import fcntl
import logging
import warnings
import json
from functools import wraps

from six import string_types
import numpy as np

from datreant.core.backends.core import FileSerial
from datreant.core.backends.statefiles import TreantFile


class SimFile(TreantFile):
    filepaths = ['abs', 'rel']

    def _init_state(self):
        super(SimFile, self)._init_state()
        self._state['mds'] = dict()
        self._state['mds']['universes'] = dict()
        self._state['mds']['default'] = None

    @FileSerial._read
    def get_mds_version(self):
        """Get Sim mdsynthesis version.

        :Returns:
            *version*
                mdsynthesis version of Sim

        """
        return self._state['mds']['version']

    # TODO: need a proper schema update mechanism
    @FileSerial._write
    def update_mds_schema(self):
        """Update mdsynthesis-specific schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            version = self._state['version']
        except KeyError:
            version = mdsynthesis.__version__

        return version

    @FileSerial._write
    def update_mds_version(self, version):
        """Update mdsynthesis version of Sim.

        :Arugments:
            *version*
                new mdsynthesis version of Treant
        """
        self._state['mds']['version'] = version
