"""Interface tests for Treants.

"""

import mdsynthesis as mds
import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py
from datreant.core.tests.test_treants import TestTreant

import MDAnalysis
from MDAnalysisTests.datafiles import GRO, XTC


class TestSim(TestTreant):
    """Test Sim-specific features"""
    treantname = 'testsim'
    treanttype = 'Sim'

    @pytest.fixture
    def treantclass(self):
        return mds.Sim

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            s = mds.Sim(TestSim.treantname)
        return s

    class TestUniverse:
        """Test universe functionality"""

        def test_add_universe(self, treant):
            """Test adding a new unverse definition"""
            treant.topology = GRO
            treant.trajectory = XTC

            assert isinstance(treant.universe, MDAnalysis.Universe)

            assert treant.universe.filename == GRO
            assert treant.universe.trajectory.filename == XTC

        def test_remove_universe(self, treant):
            """Test universe removal"""
            treant.topology = None

            assert treant.topology is None

    class TestSelections:
        """Test stored selections functionality"""
        @pytest.fixture
        def treant(self, tmpdir):
            with tmpdir.as_cwd():
                s = mds.Sim(TestSim.treantname)

                s.topology = GRO
                s.trajectory = XTC
            return s

        def test_add_selection(self, treant):
            """Test adding new selection definitions"""

            treant.selections.add('CA', 'protein and name CA')
            treant.selections.add('someres', 'resid 12')

            CA = treant.universe.select_atoms('protein and name CA')
            someres = treant.universe.select_atoms('resid 12')

            assert (CA.indices == treant.selections['CA'].indices).all()
            assert (someres.indices ==
                    treant.selections['someres'].indices).all()

        def test_remove_selection(self, treant):
            """Test universe removal"""

            treant.selections.add('CA', 'protein and name CA')
            treant.selections.add('someres', 'resid 12')

            assert 'CA' in treant.selections
            assert 'someres' in treant.selections

            treant.selections.remove('CA')

            assert 'CA' not in treant.selections
            assert 'someres' in treant.selections

            del treant.selections['someres']
            treant.selections.add('moreres', 'resid 12:20')

            assert 'CA' not in treant.selections
            assert 'someres' not in treant.selections
            assert 'moreres' in treant.selections

            with pytest.raises(KeyError):
                treant.selections.remove('someres')

            with pytest.raises(KeyError):
                del treant.selections['CA']

        def test_selection_keys(self, treant):
            treant.selections.add('CA', 'protein and name CA')
            treant.selections.add('someres', 'resid 12')

            assert set(('CA', 'someres')) == set(treant.selections.keys())

        def test_selection_define(self, treant):
            CA = 'protein and name CA'
            treant.selections.add('CA', CA)

            assert treant.selections.define('CA')[0] == CA

        def test_selection_get(self, treant):
            with pytest.raises(KeyError):
                treant.selections['CA']

        def test_add_selections_multiple_strings_via_add(self, treant):
            """Add a selection that has multiple selection strings"""
            treant.selections.add('funky town', 'name N', 'name CA')
            assert 'funky town' in treant.selections

            ref = treant.universe.select_atoms('name N', 'name CA')
            sel = treant.selections['funky town']
            assert (ref.indices == sel.indices).all()

        def test_add_selections_multiple_strings_via_setitem(self, treant):
            """Add a selection that has multiple selection strings"""
            treant.selections['funky town 2'] = 'name N', 'name CA'
            assert 'funky town 2' in treant.selections

            ref = treant.universe.select_atoms('name N', 'name CA')
            sel = treant.selections['funky town 2']
            assert (ref.indices == sel.indices).all()

        def test_add_selection_as_atomgroup_via_add(self, treant):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = treant.universe.atoms[:10:2]

            treant.selections.add('ag sel', ag)
            assert 'ag sel' in treant.selections

            ag2 = treant.selections['ag sel']
            assert (ag.indices == ag2.indices).all()

        def test_add_selection_as_atomgroup_via_setitem(self, treant):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = treant.universe.atoms[25:50:3]

            treant.selections['ag sel 2'] = ag
            assert 'ag sel 2' in treant.selections

            ag2 = treant.selections['ag sel 2']
            assert (ag.indices == ag2.indices).all()


class TestReadOnly:
    """Test Sim functionality when read-only"""

    GRO = 'md.gro'
    XTC = 'md.xtc'

    @pytest.fixture
    def sim(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = mds.Sim('testsim')

            # copy universe files to within the Sim's tree
            sub = py.path.local(c.abspath).mkdir('sub')
            GRO_t = sub.join(self.GRO)
            XTC_t = sub.join(self.XTC)
            py.path.local(GRO).copy(GRO_t)
            py.path.local(XTC).copy(XTC_t)

            c.topology = GRO_t.strpath
            c.trajectory = XTC_t.strpath

            py.path.local(c.abspath).chmod(0550, rec=True)

        def fin():
            py.path.local(c.abspath).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    def test_sim_universe_access(self, sim):
        """Test that Sim can access Universe when read-only.
        """
        assert isinstance(sim.universe, MDAnalysis.Universe)

    def test_sim_moved_universe_access(self, sim, tmpdir):
        """Test that Sim can access Universe when read-only, especially
        when universe files have moved with it (stale paths).
        """
        py.path.local(sim.abspath).chmod(0770, rec=True)
        sim.location = tmpdir.mkdir('test').strpath
        py.path.local(sim.abspath).chmod(0550, rec=True)

        assert isinstance(sim.universe, MDAnalysis.Universe)

        py.path.local(sim.abspath).chmod(0770, rec=True)
