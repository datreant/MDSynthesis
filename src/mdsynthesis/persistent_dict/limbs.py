"""
Limbs for convenient Treant data storage and retrieval.

"""
import os
import six
from functools import wraps

from datreant.core.limbs import TreeLimb
from . import npdata, pddata, pydata
from .core import DataFile


class Data(TreeLimb):
    """Interface to stored data.

    """
    _name = 'data'

    def __repr__(self):
        return "<Data({})>".format(self.keys())

    def __str__(self):
        data = self.keys()
        agg = "Data"
        majsep = "="
        seplength = len(agg)

        if not data:
            out = "No Data"
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for datum in data:
                out = out + "'{}'\n".format(datum)
        return out

    def __iter__(self):
        return self.keys().__iter__()

    def _makedirs(self, p):
        """Make directories and all parents necessary.

        :Arguments:
            *p*
                directory path to make
        """
        try:
            os.makedirs(p)
        except OSError:
            pass

    def _get_datafile(self, handle):
        """Return path to datafile corresponding to given handle.

        :Arguments:
            *handle*
                name of dataset whose datafile path to return

        :Returns:
            *datafile*
                datafile path; None if does not exist
            *proxyfile*
                proxyfile path; None if does not exist
            *datafiletype*
                datafile type; either ``persistence.pddatafile`` or
                ``persistence.npdatafile``

        """
        datafile = None
        datafiletype = None
        for dfiletype in (pddata.pddatafile, npdata.npdatafile,
                          pydata.pydatafile):
            dfile = os.path.join(self._tree.abspath,
                                 handle, dfiletype)
            if os.path.exists(dfile):
                datafile = dfile
                proxyfile = os.path.join(self._tree.abspath,
                                         handle, ".{}.proxy".format(dfiletype))
                datafiletype = dfiletype

        if datafile is None and datafiletype is None:
            raise KeyError("No data for '{}'".format(handle))

        return (datafile, proxyfile, datafiletype)

    def _read_datafile(func):
        """Decorator for generating DataFile instance for reading data.

        DataFile instance is generated and mounted at self._datafile. It
        is then dereferenced after the method call. Since data files can be
        deleted in the filesystem, this should handle cleanly the scenarios
        in which data appears, goes missing, etc. while a Treant is loaded.

        .. note:: Methods wrapped with this decorator need to have *handle* as
                  the first argument.

        """
        @wraps(func)
        def inner(self, handle, *args, **kwargs):
            filename, proxy, filetype = self._get_datafile(handle)

            if filename:
                self._datafile = DataFile(
                        os.path.join(self._tree.abspath,
                                     handle),
                        datafiletype=filetype)
                try:
                    out = func(self, handle, *args, **kwargs)
                finally:
                    del self._datafile
            else:
                out = None

            return out

        return inner

    def _write_datafile(func):
        """Decorator for generating DataFile instance for writing data.

        DataFile instance is generated and mounted at self._datafile. It
        is then dereferenced after the method call. Since data files can be
        deleted in the filesystem, this should handle cleanly the scenarios
        in which data appears, goes missing, etc. while a Treant is loaded.

        .. note:: Methods wrapped with this decorator need to have *handle* as
                  the first argument.

        """
        @wraps(func)
        def inner(self, handle, *args, **kwargs):
            dirname = os.path.join(self._tree.abspath, handle)

            self._makedirs(dirname)
            self._datafile = DataFile(dirname)

            try:
                out = func(self, handle, *args, **kwargs)
            finally:
                del self._datafile

            return out

        return inner

    def __getitem__(self, handle):
        """Get dataset corresponding to given handle(s).

        If dataset doesn't exist, ``None`` is returned.

        :Arguments:
            *handle*
                name of data to retrieve; may also be a list of names

        :Returns:
            *data*
                stored data; if *handle* was a list, will be a list
                of equal length with the stored data as members; will yield
                ``None`` if requested data is nonexistent

        """
        if isinstance(handle, list):
            out = list()
            for item in handle:
                out.append(self.retrieve(item))
        elif isinstance(handle, six.string_types):
            out = self.retrieve(handle)

        return out

    def __setitem__(self, handle, data):
        """Set dataset corresponding to given handle.

        A data instance must be either a pandas Series, DataFrame, or Panel
        object. If dataset doesn't exist, it is added. If a dataset already
        exists for the given handle, it is replaced.

        :Arguments:
            *handle*
                name given to data; needed for retrieval
            *data*
                data to store; must be a pandas Series, DataFrame, or Panel

        """
        self.add(handle, data)

    def __delitem__(self, handle):
        """Remove a dataset.

        Note: the directory containing the dataset file (``Data.h5``) will NOT
        be removed if it still contains file after the removal of the dataset
        file.

        :Arguments:
            *handle*
                name of dataset to delete

        """
        self.remove(handle)

    @_write_datafile
    def add(self, handle, data):
        """Store data in Treant.

        A data instance can be a pandas object (Series, DataFrame, Panel),
        a numpy array, or a pickleable python object. If the dataset doesn't
        exist, it is added. If a dataset already exists for the given handle,
        it is replaced.

        :Arguments:
            *handle*
                name given to data; needed for retrieval
            *data*
                data structure to store

        """
        self._datafile.add_data('main', data)

    def remove(self, handle, **kwargs):
        """Remove a dataset, or some subset of a dataset.

        Note: in the case the whole dataset is removed, the directory
        containing the dataset file (``Data.h5``) will NOT be removed if it
        still contains file(s) after the removal of the dataset file.

        For pandas objects (Series, DataFrame, or Panel) subsets of the whole
        dataset can be removed using keywords such as *start* and *stop* for
        ranges of rows, and *columns* for selected columns.

        :Arguments:
            *handle*
                name of dataset to delete

        :Keywords:
            *where*
                conditions for what rows/columns to remove
            *start*
                row number to start selection
            *stop*
                row number to stop selection
            *columns*
                columns to remove

        """
        datafile, proxy, datafiletype = self._get_datafile(handle)

        if kwargs and datafiletype == pddata.pddatafile:
            self._delete_data(handle, **kwargs)
        elif datafile:
            os.remove(datafile)
            os.remove(proxy)
            top = self._tree.abspath
            directory = os.path.dirname(datafile)
            while directory != top:
                try:
                    os.rmdir(directory)
                    directory = os.path.dirname(directory)
                except OSError:
                    break

    @_write_datafile
    def _delete_data(self, handle, **kwargs):
        """Remove a dataset, or some subset of a dataset.

        This method loads the given data instance before doing anything,
        and should generally not be used to remove data since it will create
        a datafile object if one is not already present, which could have
        side-effects for other instances of this Treant.

        Note: in the case the whole dataset is removed, the directory
        containing the dataset file (``Data.h5``) will NOT be removed if it
        still contains file(s) after the removal of the dataset file.

        :Arguments:
            *handle*
                name of dataset to delete

        :Keywords:
            *where*
                conditions for what rows/columns to remove
            *start*
                row number to start selection
            *stop*
                row number to stop selection
            *columns*
                columns to remove

        """
        # only called for pandas objects at the moment
        filename, proxy, filetype = self._get_datafile(handle)
        self._datafile.datafiletype = filetype
        try:
            self._datafile.del_data('main', **kwargs)
        except NotImplementedError:
            datafile = self._get_datafile(handle)[0]

            os.remove(datafile)
            os.remove(proxy)
            top = self._tree.abspath
            directory = os.path.dirname(datafile)
            while directory != top:
                try:
                    os.rmdir(directory)
                    directory = os.path.dirname(directory)
                except OSError:
                    break

    @_read_datafile
    def retrieve(self, handle, **kwargs):
        """Retrieve stored data.

        The stored data structure is read from disk and returned.

        If dataset doesn't exist, ``None`` is returned.

        For pandas objects (Series, DataFrame, or Panel) subsets of the whole
        dataset can be returned using keywords such as *start* and *stop* for
        ranges of rows, and *columns* for selected columns.

        Also for pandas objects, the *where* keyword takes a string as input
        and can be used to filter out rows and columns without loading the full
        object into memory. For example, given a DataFrame with handle 'mydata'
        with columns (A, B, C, D), one could return all rows for columns A and
        C for which column D is greater than .3 with::

            retrieve('mydata', where='columns=[A,C] & D > .3')

        Or, if we wanted all rows with index = 3 (there could be more than
        one)::

            retrieve('mydata', where='index = 3')

        See :meth:`pandas.HDFStore.select` for more information.

        :Arguments:
            *handle*
                name of data to retrieve

        :Keywords:
            *where*
                conditions for what rows/columns to return
            *start*
                row number to start selection
            *stop*
                row number to stop selection
            *columns*
                list of columns to return; all columns returned by default
            *iterator*
                if True, return an iterator [``False``]
            *chunksize*
                number of rows to include in iteration; implies
                ``iterator=True``

        :Returns:
            *data*
                stored data; ``None`` if nonexistent

        """
        return self._datafile.get_data('main', **kwargs)

    @_write_datafile
    def append(self, handle, data):
        """Append rows to an existing dataset.

        The object must be of the same pandas class (Series, DataFrame, Panel)
        as the existing dataset, and it must have exactly the same columns
        (names included).

        :Arguments:
            *handle*
                name of data to append to
            *data*
                data to append

        """
        self._datafile.append_data('main', data)

    def keys(self):
        """List available datasets.

        :Returns:
            *handles*
                list of handles to available datasets

        """
        datasets = list()
        top = self._tree.abspath
        for root, dirs, files in os.walk(top):
            if ((pddata.pddatafile in files) or
                    (npdata.npdatafile in files) or
                    (pydata.pydatafile in files)):
                datasets.append(os.path.relpath(root, start=top))

        datasets.sort()

        return datasets
