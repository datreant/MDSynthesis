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

import MDAnalysis as mda
from MDAnalysisTests.datafiles import PDB, GRO, XTC


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
            treant.udef.topology = GRO
            treant.udef.trajectory = XTC

            assert isinstance(treant.universe, mda.Universe)

            assert treant.universe.filename == GRO
            assert treant.universe.trajectory.filename == XTC

        def test_set_universe(self, treant):
            """Test setting the universe with a Universe object"""

            u = mda.Universe(GRO, XTC)

            # set the universe directly
            treant.universe = u

            assert isinstance(treant.universe, mda.Universe)
            assert treant.udef.topology == GRO
            assert treant.udef.trajectory == XTC

        def test_add_univese_typeerror(self, treant):
            """Test checking of what is passed to setter"""
            with pytest.raises(TypeError):
                treant.universe = 72

        def test_remove_universe(self, treant):
            """Test universe removal"""
            treant.udef.topology = GRO
            treant.udef.trajectory = XTC

            assert isinstance(treant.universe, mda.Universe)

            treant.udef.topology = None
            assert treant.udef.topology is None
            assert treant.universe is None

            # trajectory definition should still be there
            assert treant.udef.trajectory

            # this should give us a universe again
            treant.udef.topology = PDB

            assert isinstance(treant.universe, mda.Universe)

            # removing trajectories should keep universe, but with PDB as
            # coordinates
            treant.udef.trajectory = None

            assert isinstance(treant.universe, mda.Universe)
            assert treant.universe.trajectory.n_frames == 1

        def test_set_resnums(self, treant):
            """Test that we can add resnums to a universe."""
            treant.udef.topology = GRO
            treant.udef.trajectory = XTC

            protein = treant.universe.select_atoms('protein')
            resids = protein.residues.resids
            protein.residues.set_resnum(resids + 3)

            treant.udef._set_resnums(treant.universe.residues.resnums)

            treant.udef.reload()

            protein = treant.universe.select_atoms('protein')
            assert (resids + 3 == protein.residues.resnums).all()

            # test resetting of resnums
            protein.residues.set_resnum(resids + 6)

            assert (protein.residues.resnums == resids + 6).all()
            treant.udef._set_resnums(treant.universe.residues.resnums)

            treant.udef.reload()

            protein = treant.universe.select_atoms('protein')
            assert (resids + 6 == protein.residues.resnums).all()

    class TestSelections:
        """Test stored atomselections functionality"""
        @pytest.fixture
        def treant(self, tmpdir):
            with tmpdir.as_cwd():
                s = mds.Sim(TestSim.treantname)

                s.udef.topology = GRO
                s.udef.trajectory = XTC
            return s

        def test_add_selection(self, treant):
            """Test adding new selection definitions"""

            treant.atomselections.add('CA', 'protein and name CA')
            treant.atomselections.add('someres', 'resid 12')

            CA = treant.universe.select_atoms('protein and name CA')
            someres = treant.universe.select_atoms('resid 12')

            assert (CA.indices == treant.atomselections['CA'].indices).all()
            assert (someres.indices ==
                    treant.atomselections['someres'].indices).all()

        def test_remove_selection(self, treant):
            """Test universe removal"""

            treant.atomselections.add('CA', 'protein and name CA')
            treant.atomselections.add('someres', 'resid 12')

            assert 'CA' in treant.atomselections
            assert 'someres' in treant.atomselections

            treant.atomselections.remove('CA')

            assert 'CA' not in treant.atomselections
            assert 'someres' in treant.atomselections

            del treant.atomselections['someres']
            treant.atomselections.add('moreres', 'resid 12:20')

            assert 'CA' not in treant.atomselections
            assert 'someres' not in treant.atomselections
            assert 'moreres' in treant.atomselections

            with pytest.raises(KeyError):
                treant.atomselections.remove('someres')

            with pytest.raises(KeyError):
                del treant.atomselections['CA']

        def test_selection_keys(self, treant):
            treant.atomselections.add('CA', 'protein and name CA')
            treant.atomselections.add('someres', 'resid 12')

            assert set(('CA', 'someres')) == set(treant.atomselections.keys())

        def test_selection_define(self, treant):
            CA = 'protein and name CA'
            treant.atomselections.add('CA', CA)

            assert treant.atomselections.define('CA')[0] == CA

        def test_selection_get(self, treant):
            with pytest.raises(KeyError):
                treant.atomselections['CA']

        def test_add_atomselections_multiple_strings_via_add(self, treant):
            """Add a selection that has multiple selection strings"""
            treant.atomselections.add('funky town', 'name N', 'name CA')
            assert 'funky town' in treant.atomselections

            ref = treant.universe.select_atoms('name N', 'name CA')
            sel = treant.atomselections['funky town']
            assert (ref.indices == sel.indices).all()

        def test_add_atomselections_multiple_strings_via_setitem(self, treant):
            """Add a selection that has multiple selection strings"""
            treant.atomselections['funky town 2'] = 'name N', 'name CA'
            assert 'funky town 2' in treant.atomselections

            ref = treant.universe.select_atoms('name N', 'name CA')
            sel = treant.atomselections['funky town 2']
            assert (ref.indices == sel.indices).all()

        def test_add_selection_as_atomgroup_via_add(self, treant):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = treant.universe.atoms[:10:2]

            treant.atomselections.add('ag sel', ag)
            assert 'ag sel' in treant.atomselections

            ag2 = treant.atomselections['ag sel']
            assert (ag.indices == ag2.indices).all()

        def test_add_selection_as_atomgroup_via_setitem(self, treant):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = treant.universe.atoms[25:50:3]

            treant.atomselections['ag sel 2'] = ag
            assert 'ag sel 2' in treant.atomselections

            ag2 = treant.atomselections['ag sel 2']
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

            c.udef.topology = GRO_t.strpath
            c.udef.trajectory = XTC_t.strpath

            py.path.local(c.abspath).chmod(0550, rec=True)

        def fin():
            py.path.local(c.abspath).chmod(0770, rec=True)

        request.addfinalizer(fin)

        return c

    def test_sim_universe_access(self, sim):
        """Test that Sim can access Universe when read-only.
        """
        assert isinstance(sim.universe, mda.Universe)

    def test_sim_moved_universe_access(self, sim, tmpdir):
        """Test that Sim can access Universe when read-only, especially
        when universe files have moved with it (stale paths).
        """
        py.path.local(sim.abspath).chmod(0770, rec=True)
        sim.location = tmpdir.mkdir('test').strpath
        py.path.local(sim.abspath).chmod(0550, rec=True)

        assert isinstance(sim.universe, mda.Universe)

        py.path.local(sim.abspath).chmod(0770, rec=True)
