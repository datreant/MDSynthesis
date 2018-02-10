import pandas as pd
import numpy as np
import pytest

import datreant.core as dtr
import datreant.data.attach
from . import test_data


class CollectionTests():
    """Tests for common elements of Group.members and Bundle"""

    @pytest.fixture
    def testtreant(self, tmpdir, request):
        with tmpdir.as_cwd():
            t = dtr.Treant('dummytreant')
        return t

    @pytest.fixture
    def testgroup(self, tmpdir, request):
        with tmpdir.as_cwd():
            g = dtr.Group('dummygroup')
            g.members.add(dtr.Treant('bark'), dtr.Treant('leaf'))
        return g

    class TestAggData:
        """Test member data functionality"""

        @pytest.fixture
        def collection(self, collection, tmpdir):
            with tmpdir.as_cwd():
                s1 = dtr.Treant('lark')
                s2 = dtr.Treant('hark')
                g3 = dtr.Group('linus')

            collection.add(s1, [g3, s2])
            return collection

        class DataMixin:
            """Mixin class for data storage tests.

            Contains general tests to be used for all storable data formats.

            """
            handle = 'testdata'

            def test_retrieve_data(self, collection, datastruct):
                agg = dict()
                for member in collection:
                    member.data.add(self.handle, datastruct)
                    agg[member.abspath] = datastruct

                np.testing.assert_equal(collection.data.retrieve(self.handle),
                                        agg)
                np.testing.assert_equal(collection.data[self.handle],
                                        agg)

        class PanelMixin:
            """Mixin class for pandas structures that don't support
            MultiIndexing.

            """
            handle = 'testdata'

            def test_retrieve_data(self, collection, datastruct):
                agg = dict()
                for member in collection:
                    member.data.add(self.handle, datastruct)
                    agg[member.abspath] = datastruct

                stored = collection.data.retrieve(self.handle)
                for item in agg:
                    np.testing.assert_equal(stored[item].values,
                                            agg[item].values)

                stored = collection.data[self.handle]
                for item in agg:
                    np.testing.assert_equal(stored[item].values,
                                            agg[item].values)

        class MultiIndexMixin(DataMixin):
            """Mixin class for pandas structures that support MultiIndexes"""

            def test_retrieve_data(self, collection, datastruct):
                agg = dict()
                for member in collection:
                    member.data.add(self.handle, datastruct)
                    agg[member.abspath] = datastruct

                def dict2multiindex(agg):
                    agg_mi = None
                    for member in agg:
                        d = agg[member].copy(deep=True)
                        label = len(d.index)*[member]
                        index = pd.MultiIndex.from_arrays([label, d.index])
                        d.index = index

                        if agg_mi is not None:
                            agg_mi = agg_mi.append(d)
                        else:
                            agg_mi = d

                    return agg_mi

                np.testing.assert_equal(
                        collection.data.retrieve(self.handle).values,
                        dict2multiindex(agg).values)
                np.testing.assert_equal(
                        collection.data[self.handle].values,
                        dict2multiindex(agg).values)

        class Test_Series(test_data.Series, MultiIndexMixin):
            pass

        class Test_DataFrame(test_data.DataFrame, MultiIndexMixin):
            pass

        class Test_Blank_DataFrame(test_data.Blank_DataFrame, MultiIndexMixin):
            pass

        class Test_Wide_Blank_DataFrame(test_data.Wide_Blank_DataFrame,
                                        MultiIndexMixin):
            pass

        class Test_Thin_Blank_DataFrame(test_data.Thin_Blank_DataFrame,
                                        MultiIndexMixin):
            pass

        class Test_Panel(test_data.Panel, PanelMixin):
            pass

        class Test_Panel4D(test_data.Panel4D, PanelMixin):
            pass

        class Test_NumpyScalar(test_data.NumpyScalar, DataMixin):
            pass

        class Test_Numpy1D(test_data.Numpy1D, DataMixin):
            pass

        class Test_Numpy2D(test_data.Numpy2D, DataMixin):
            pass

        class Test_Wide_Numpy2D(test_data.Wide_Numpy2D, DataMixin):
            pass

        class Test_Thin_Numpy2D(test_data.Thin_Numpy2D, DataMixin):
            pass

        class Test_Numpy3D(test_data.Numpy3D, DataMixin):
            pass

        class Test_Numpy4D(test_data.Numpy4D, DataMixin):
            pass

        class Test_List(test_data.List, DataMixin):
            pass

        class Test_Dict(test_data.Dict, DataMixin):
            pass

        class Test_Tuple(test_data.Tuple, DataMixin):
            pass

        class Test_Set(test_data.Set, DataMixin):
            pass

        class Test_Dict_Mix(test_data.Dict_Mix, DataMixin):
            pass


class TestBundle(CollectionTests):
    """Test Bundle features"""

    @pytest.fixture
    def collection(self):
        return dtr.Bundle()
