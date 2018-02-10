"""
Abstract interface components for reading and writing datasets.

"""

import os

import numpy as np
import pandas as pd

from . import pydata
from . import npdata
from . import pddata


class DataFile(object):
    """Interface to data files.

    This is an abstraction layer to the pddata.pdDataFile, npdata.npDataFile,
    and pydata.pyDataFile objects. This can be used by higher level objects
    without worrying about whether to use pandas storers, numpy storers, or
    pickle.

    """

    def __init__(self, datadir, datafiletype=None, **kwargs):
        """Initialize data interface.

        :Arguments:
           *datadir*
              path to data directory
           *datafiletype*
              If known, either pddata.pddatafile or npdata.npdatafile

        """
        self.datadir = datadir
        self.datafile = None

        # if given, can get data
        self.datafiletype = datafiletype

    def add_data(self, key, data):
        """Add a pandas data object (Series, DataFrame, Panel), numpy array,
        or pickleable python object to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                name given to the data; used as the index for retrieving
                the data later
            *data*
                the data object to store; should be either a pandas Series,
                DataFrame, Panel, or a numpy array
        """
        if isinstance(data, np.ndarray):
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile))
        elif isinstance(data, (pd.Series, pd.DataFrame, pd.Panel, pd.Panel4D)):
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile))
        else:
            self.datafile = pydata.pyDataFile(
                os.path.join(self.datadir, pydata.pydatafile))

        self.datafile.add_data(key, data)

        # dereference
        self.datafile = None

    def append_data(self, key, data):
        """Append rows to an existing pandas data object stored in the data file.

        Note that column names of new data must match those of the existing
        data. Columns cannot be appended due to the technical details of the
        HDF5 standard. To add new columns, store as a new dataset.

        :Arguments:
            *key*
                name of existing data object to append to
            *data*
                the data object whose rows are to be appended to the existing
                stored data; must have same columns (with names) as existing
                data

        """
        # TODO: add exceptions where appending isn't possible
        if isinstance(data, np.ndarray):
            raise TypeError('Cannot append numpy arrays.')
        elif isinstance(data, (pd.Series, pd.DataFrame, pd.Panel, pd.Panel4D)):
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile))
            self.datafile.append_data(key, data)
        else:
            raise TypeError('Cannot append python object.')

            # dereference
            self.datafile = None

    def get_data(self, key, **kwargs):
        """Retrieve data object stored in file.

        :Arguments:
            *key*
                name of data to retrieve

        :Keywords:
            *where*
                for pandas objects, conditions for what rows/columns to return
            *start*
                for pandas objects, row number to start selection
            *stop*
                for pandas objects, row number to stop selection
            *columns*
                for pandas objects, list of columns to return; all columns
                returned by default
            *iterator*
                for pandas objects, if True, return an iterator [``False``]
            *chunksize*
                for pandas objects, number of rows to include in iteration;
                implies ``iterator=True``

        :Returns:
            *data*
                the selected data
        """
        if self.datafiletype == npdata.npdatafile:
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile))
            out = self.datafile.get_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddata.pddatafile:
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile))
            out = self.datafile.get_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydata.pydatafile:
            self.datafile = pydata.pyDataFile(
                os.path.join(self.datadir, pydata.pydatafile))
            out = self.datafile.get_data(key)
            self.datafile = None
        else:
            raise TypeError('Cannot return data without knowing datatype.')
            out = None

        return out

    def del_data(self, key, **kwargs):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete

        :Keywords:
            *where*
                for pandas objects, conditions for what rows/columns to remove
            *start*
                for pandas objecs, row number to start selection
            *stop*
                for pandas objects, row number to stop selection

        """
        if self.datafiletype == npdata.npdatafile:
            self.datafile = npdata.npDataFile(
                os.path.join(self.datadir, npdata.npdatafile))
            out = self.datafile.del_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pddata.pddatafile:
            self.datafile = pddata.pdDataFile(
                os.path.join(self.datadir, pddata.pddatafile))
            out = self.datafile.del_data(key, **kwargs)
            self.datafile = None
        elif self.datafiletype == pydata.pydatafile:
            pass
        else:
            raise TypeError('Cannot return data without knowing datatype.')
            out = None

        return out
