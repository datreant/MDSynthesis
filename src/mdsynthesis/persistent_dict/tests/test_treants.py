
import datreant.core as dtr
import datreant.data.attach
import pandas as pd
import numpy as np
import pytest
import os
import py

from datreant.data.tests import test_data


class TestTreant:
    treantname = 'testtreant'
    treanttype = 'Treant'
    treantclass = dtr.treants.Treant

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            c = dtr.treants.Treant(TestTreant.treantname)
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
            datafile = datreant.data.pddata.pddatafile

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
                index = datastruct.index
                for i in range(5):
                    treant.data.append(self.handle, datastruct)

                stored = treant.data.retrieve(self.handle)
                equiv = pd.concat([datastruct]*5)

                np.testing.assert_equal(stored, equiv)

        class Test_Series(test_data.Series, PandasMixin):
            pass

        class Test_DataFrame(test_data.DataFrame, PandasMixin):
            pass

        class Test_Blank_DataFrame(test_data.Blank_DataFrame, PandasMixin):
            pass

        class Test_Wide_Blank_DataFrame(test_data.Wide_Blank_DataFrame,
                                        PandasMixin):
            pass

        class Test_Thin_Blank_DataFrame(test_data.Thin_Blank_DataFrame,
                                        PandasMixin):
            pass

        class Test_Panel(test_data.Panel, PandasMixin):
            pass

        class Test_Panel4D(test_data.Panel4D, PandasMixin):
            pass

        class NumpyMixin(DataMixin):
            """Test numpy datastructure storage and retrieval"""
            datafile = datreant.data.npdata.npdatafile

        class Test_NumpyScalar(test_data.NumpyScalar, NumpyMixin):
            pass

        class Test_Numpy1D(test_data.Numpy1D, NumpyMixin):
            pass

        class Test_Numpy2D(test_data.Numpy2D, NumpyMixin):
            pass

        class Test_Wide_Numpy2D(test_data.Wide_Numpy2D, NumpyMixin):
            pass

        class Test_Thin_Numpy2D(test_data.Thin_Numpy2D, NumpyMixin):
            pass

        class Test_Numpy3D(test_data.Numpy3D, NumpyMixin):
            pass

        class Test_Numpy4D(test_data.Numpy4D, NumpyMixin):
            pass

        class PythonMixin(DataMixin):
            """Test pandas datastructure storage and retrieval"""
            datafile = datreant.data.pydata.pydatafile

            def test_overwrite_data(self, treant, datastruct):
                treant.data[self.handle] = datastruct

                # overwrite the data with a scalar
                treant.data[self.handle] = 23
                assert treant.data[self.handle] == 23

        class Test_List(test_data.List, PythonMixin):
            pass

        class Test_Dict(test_data.Dict, PythonMixin):
            pass

        class Test_Tuple(test_data.Tuple, PythonMixin):
            pass

        class Test_Set(test_data.Set, PythonMixin):
            pass

        class Test_Dict_Mix(test_data.Dict_Mix, PythonMixin):
            pass
