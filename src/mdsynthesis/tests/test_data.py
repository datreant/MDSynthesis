import pandas as pd
import numpy as np
import pytest
import os
import py

import mdsynthesis as mds
from mdsynthesis.tests import data


class TestTreant:
    treantname = 'testtreant'

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            c = mds.Sim(TestTreant.treantname)
        return c

    class TestData:
        """Test data storage and retrieval"""

        class DataMixin:
            """Mixin class for data storage tests.

            Contains general tests to be used for all storable data formats.

            """
            handle = 'testdata'

            def test_add_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                assert os.path.exists(os.path.join(treant.abspath,
                                                   self.handle,
                                                   self.datafile))

            def test_remove_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                assert os.path.exists(os.path.join(treant.abspath,
                                                   self.handle,
                                                   self.datafile))

                treant.data.remove('testdata')
                assert not os.path.exists(os.path.join(treant.abspath,
                                                       self.handle,
                                                       self.datafile))

                # check that directory got deleted, too
                assert not os.path.exists(os.path.join(treant.abspath,
                                                       self.handle))

            def test_retrieve_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                np.testing.assert_equal(treant.data.retrieve(self.handle),
                                        datastruct)
                np.testing.assert_equal(treant.data[self.handle],
                                        datastruct)

        class PandasMixin(DataMixin):
            """Mixin class for pandas tests"""
            datafile = mds.persistent_dict.pddata.pddatafile

            def test_retrieve_data(self, treant, datastruct):
                treant.data.add(self.handle, datastruct)
                np.testing.assert_equal(
                        treant.data.retrieve(self.handle).values,
                        datastruct.values)
                np.testing.assert_equal(
                        treant.data[self.handle].values,
                        datastruct.values)

        # TODO
        class AppendablesMixin:
            """Mixin class for pandas objects that we expect should append"""

            def test_append_data(self, treant, datastruct):
                for i in range(5):
                    treant.data.append(self.handle, datastruct)

                stored = treant.data.retrieve(self.handle)
                equiv = pd.concat([datastruct]*5)

                np.testing.assert_equal(stored, equiv)

        class Test_Series(data.Series, PandasMixin):
            pass

        class Test_DataFrame(data.DataFrame, PandasMixin):
            pass

        class Test_Blank_DataFrame(data.Blank_DataFrame, PandasMixin):
            pass

        class Test_Wide_Blank_DataFrame(data.Wide_Blank_DataFrame,
                                        PandasMixin):
            pass

        class Test_Thin_Blank_DataFrame(data.Thin_Blank_DataFrame,
                                        PandasMixin):
            pass

        class NumpyMixin(DataMixin):
            """Test numpy datastructure storage and retrieval"""
            datafile = mds.persistent_dict.npdata.npdatafile

        class Test_NumpyScalar(data.NumpyScalar, NumpyMixin):
            pass

        class Test_Numpy1D(data.Numpy1D, NumpyMixin):
            pass

        class Test_Numpy2D(data.Numpy2D, NumpyMixin):
            pass

        class Test_Wide_Numpy2D(data.Wide_Numpy2D, NumpyMixin):
            pass

        class Test_Thin_Numpy2D(data.Thin_Numpy2D, NumpyMixin):
            pass

        class Test_Numpy3D(data.Numpy3D, NumpyMixin):
            pass

        class Test_Numpy4D(data.Numpy4D, NumpyMixin):
            pass

        class PythonMixin(DataMixin):
            """Test pandas datastructure storage and retrieval"""
            datafile = mds.persistent_dict.pydata.pydatafile

            def test_overwrite_data(self, treant, datastruct):
                treant.data[self.handle] = datastruct

                # overwrite the data with a scalar
                treant.data[self.handle] = 23
                assert treant.data[self.handle] == 23

        class Test_List(data.List, PythonMixin):
            pass

        class Test_Dict(data.Dict, PythonMixin):
            pass

        class Test_Tuple(data.Tuple, PythonMixin):
            pass

        class Test_Set(data.Set, PythonMixin):
            pass

        class Test_Dict_Mix(data.Dict_Mix, PythonMixin):
            pass
