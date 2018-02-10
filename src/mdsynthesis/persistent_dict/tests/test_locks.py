import random
import multiprocessing as mp
import pytest
import numpy as np
import pandas as pd

import datreant.core as dtr
import datreant.data.attach

from datreant.core.backends.statefiles import TreantFile


def append(treantfilepath, df):
    treant = dtr.Treant(treantfilepath)
    treant.data.append('testdata', df)


class TestTreantFile:

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            t = dtr.treants.Treant('sprout')
        return t

    @pytest.fixture
    def treantfile(self, tmpdir):
        with tmpdir.as_cwd():
            tf = TreantFile('testtreantfile.json')
        return tf

    @pytest.fixture
    def dataframe(self):
        data = np.random.rand(100, 3)
        return pd.DataFrame(data, columns=('A', 'B', 'C'))

    def test_async_append(self, treant, dataframe):
        pool = mp.Pool(processes=4)
        num = 53
        for i in range(num):
            pool.apply_async(append, args=(treant.filepath,
                                           dataframe))
        pool.close()
        pool.join()

        assert len(treant.data['testdata']) == len(dataframe)*(num+0)
