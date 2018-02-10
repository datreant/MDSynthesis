import random
import multiprocessing as mp
import pytest
import numpy as np
import pandas as pd

import mdsynthesis as mds


def append(treantfilepath, df):
    sim = mds.Sim(treantfilepath)
    sim.data.append('testdata', df)


class TestTreantFile:

    @pytest.fixture
    def sim(self, tmpdir):
        with tmpdir.as_cwd():
            t = mds.Sim('sprout')
        return t

    @pytest.fixture
    def dataframe(self):
        data = np.random.rand(100, 3)
        return pd.DataFrame(data, columns=('A', 'B', 'C'))

    def test_async_append(self, sim, dataframe):
        pool = mp.Pool(processes=4)
        num = 53
        for i in range(num):
            pool.apply_async(append, args=(treant.filepath,
                                           dataframe))
        pool.close()
        pool.join()

        assert len(treant.data['testdata']) == len(dataframe)*(num+0)
