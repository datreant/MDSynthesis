"""
File backends for storing numpy arrays.

"""

import h5py

from datreant.core.backends.core import File


npdatafile = 'npData.h5'


class npDataFile(File):
    """Interface to numpy object data files.

    Data is stored as numpy arrays in the HDF5 format. This class gives the
    needed components for storing and retrieving stored data. It uses h5py as
    its backend.

    """
    def _open_file_r(self):
        return h5py.File(self.filename, 'r')

    def _open_file_w(self):
        return h5py.File(self.filename, 'w')

    def add_data(self, key, data):
        """Add a numpy array to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                name given to the data; used as the index for retrieving
                the data later
            *data*
                the numpy array to store
        """
        with self.write():
            try:
                self.handle.create_dataset(key, data=data)
            except RuntimeError:
                del self.handle[key]
                self.handle.create_dataset(key, data=data)

    def get_data(self, key, **kwargs):
        """Retrieve numpy array stored in file.

        :Arguments:
            *key*
                name of data to retrieve

        :Returns:
            *data*
                the selected data
        """
        with self.read():
            return self.handle[key].value

    def del_data(self, key, **kwargs):
        """Delete a stored data object.

        :Arguments:
            *key*
                name of data to delete

        """
        with self.write():
            del self.handle[key]

    def list_data(self):
        """List names of all stored datasets.

        Although the true names start with '\' indicating the root of the
        HDF5 data tree, we wish to abstract this away. We remove the leading
        '\' from the output. This shouldn't cause any problems since the
        leading '\' can be neglected when referring to stored objects by name
        using all of h5py.File's methods anyway.

        """
        with self.read():
            keys = self.handle.keys()
            return keys
