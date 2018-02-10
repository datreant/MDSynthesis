"""
File backends for storing general python objects.

"""

from six.moves import cPickle as pickle

from datreant.core.backends.core import File


pydatafile = 'pyData.pkl'


class pyDataFile(File):
    """Interface to python object data files.

    Arbitrary python objects are stored as pickled objects on disk. This class
    gives the needed components for storing and retrieving stored data in the
    same basic way as for pandas and numpy objects. It uses pickle files for
    serialization.

    """
    def _open_file_r(self):
        return open(self.filename, 'rb')

    def _open_file_w(self):
        return open(self.filename, 'wb+')

    def add_data(self, key, data):
        """Add a numpy array to the data file.

        If data already exists for the given key, then it is overwritten.

        :Arguments:
            *key*
                not used, but needed to give consistent interface
            *data*
                the numpy array to store
        """
        # use highest python 2 pickle protocol. This allows efficient storage
        # across python versions
        with self.write():
            pickle.dump(data, self.handle, 2)

    def get_data(self, key, **kwargs):
        """Retrieve numpy array stored in file.

        :Arguments:
            *key*
                not used, but needed to give consistent interface

        :Returns:
            *data*
                the selected data
        """
        with self.read():
            # load in bytes with python3 to ensure successful read EVERYTIME.
            # Using the ascii encoding it can happen that python3 can read a
            # pickle written with python 2.
            try:
                return pickle.load(self.handle, encoding='bytes')
            except TypeError:
                return pickle.load(self.handle)
