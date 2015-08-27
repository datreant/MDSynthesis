"""Interface tests for Treants.

"""

import mdsynthesis as mds
import pandas as pd
import numpy as np
import pytest
import os
import shutil
import py
from datreant.tests.test_treants import TestTreant

import MDAnalysis
from MDAnalysisTests.datafiles import GRO, XTC


class TestSim(TestTreant):
    """Test Sim-specific features"""
    treantname = 'testsim'
    treanttype = 'Sim'
    treantclass = mds.Sim

    @pytest.fixture
    def treant(self, tmpdir):
        with tmpdir.as_cwd():
            s = mds.Sim(TestSim.treantname)
        return s

    class TestUniverses:
        """Test universe functionality"""

        def test_add_universe(self, treant):
            """Test adding new unverse definitions"""
            treant.universes.add('spam', GRO, XTC)

            assert 'spam' in treant.universes
            assert isinstance(treant.universes['spam'], MDAnalysis.Universe)

            assert treant.universe.filename == GRO
            assert treant.universe.trajectory.filename == XTC

        def test_remove_universe(self, treant):
            """Test universe removal"""
            treant.universes.add('spam', GRO, XTC)
            treant.universes.add('eggs', GRO, XTC)

            assert 'spam' in treant.universes

            treant.universes.remove('spam')

            assert 'spam' not in treant.universes
            assert 'eggs' in treant.universes

            with pytest.raises(KeyError):
                treant.universes.remove('ham')

        def test_rename_universe(self, treant):
            """Test universe renaming"""
            treant.universes.add('spam', GRO, XTC)
            treant.universes.add('eggs', GRO, XTC)

            assert 'spam' in treant.universes

            treant.universes.rename('spam', 'boots')

            assert 'spam' not in treant.universes
            assert 'boots' in treant.universes
            assert 'eggs' in treant.universes

            with pytest.raises(KeyError):
                treant.universes.rename('ham', 'lark')

            with pytest.raises(ValueError):
                treant.universes.rename('boots', 'eggs')

        def test_set_default_universe(self, treant):
            """Test that a default universe exists, and that it's settable"""
            treant.universes.add('lolcats', GRO, XTC)

            assert treant.universe.filename == GRO
            assert treant.universe.trajectory.filename == XTC
            assert treant._uname == 'lolcats'
            treant.universes.deactivate()

            treant.universes.add('megaman', GRO)
            treant.universes.default('megaman')

            assert treant.universes.default() == 'megaman'

            assert treant.universe.filename == GRO
            assert treant.universe.trajectory.filename == GRO
            assert treant._uname == 'megaman'

            treant.universes.remove('megaman')

            assert treant.universes.default() == None

        def test_set_resnums(self, treant):
            """Test that we can add resnums to a universe."""
            treant.universes.add('lolcats', GRO, XTC)

            protein = treant.universe.selectAtoms('protein')
            resids = protein.residues.resids()
            protein.residues.set_resnum(resids + 3)

            treant.universes.resnums('lolcats',
                                     treant.universe.atoms.resnums())

            treant.universes['lolcats']

            protein = treant.universe.selectAtoms('protein')
            assert (resids + 3 == protein.residues.resnums()).all()

            # BUG IN MDANALYSIS PREVENTS RESETTING OF RESNUMS
            # protein.residues.set_resnum(resids + 6)

            # assert (protein.residues.resnums == resids + 6).all()
            # treant.universes.resnums('lolcats',
            #                            treant.universe.atoms.resnums())

            # treant.universes['lolcats']

            # protein = treant.universe.selectAtoms('protein')
            # assert (resids + 6 == protein.residues.resnums()).all()

        def test_KeyError(self, treant):
            """Test that a KeyError raised when trying to activate a Universe
            that doesn't exist.
            """
            with pytest.raises(KeyError):
                treant.universes.activate('ham')

            treant.universes.add('lolcats', GRO, XTC)

            with pytest.raises(KeyError):
                treant.universes.activate('eggs')

    class TestSelections:
        """Test stored selections functionality"""
        @pytest.fixture
        def treant(self, tmpdir):
            with tmpdir.as_cwd():
                s = mds.Sim(TestSim.treantname)
                s.universes.add('spam', GRO, XTC)
            return s

        def test_add_selection(self, treant):
            """Test adding new selection definitions"""

            treant.selections.add('CA', 'protein and name CA')
            treant.selections.add('someres', 'resid 12')

            CA = treant.universe.selectAtoms('protein and name CA')
            someres = treant.universe.selectAtoms('resid 12')

            assert (CA.indices() == treant.selections['CA'].indices()).all()
            assert (someres.indices() ==
                    treant.selections['someres'].indices()).all()

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

            ref = treant.universe.selectAtoms('name N', 'name CA')
            sel = treant.selections['funky town']
            assert (ref.indices() == sel.indices()).all()

        def test_add_selections_multiple_strings_via_setitem(self, treant):
            """Add a selection that has multiple selection strings"""
            treant.selections['funky town 2'] = 'name N', 'name CA'
            assert 'funky town 2' in treant.selections

            ref = treant.universe.selectAtoms('name N', 'name CA')
            sel = treant.selections['funky town 2']
            assert (ref.indices() == sel.indices()).all()

        def test_add_selection_as_atomgroup_via_add(self, treant):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = treant.universe.atoms[:10:2]

            treant.selections.add('ag sel', ag)
            assert 'ag sel' in treant.selections

            ag2 = treant.selections['ag sel']
            assert (ag.indices() == ag2.indices()).all()

        def test_add_selection_as_atomgroup_via_setitem(self, treant):
            """Make an arbitrary AtomGroup then save selection as AG"""
            ag = treant.universe.atoms[25:50:3]

            treant.selections['ag sel 2'] = ag
            assert 'ag sel 2' in treant.selections

            ag2 = treant.selections['ag sel 2']
            assert (ag.indices() == ag2.indices()).all()


class TestReadOnly:
    """Test Sim functionality when read-only"""

    GRO = 'md.gro'
    XTC = 'md.xtc'

    @pytest.fixture
    def sim(self, tmpdir, request):
        with tmpdir.as_cwd():
            c = mds.Sim('testsim')

            # copy universe files to within the Sim's tree
            sub = py.path.local(c.basedir).mkdir('sub')
            GRO_t = sub.join(self.GRO)
            XTC_t = sub.join(self.XTC)
            py.path.local(GRO).copy(GRO_t)
            py.path.local(XTC).copy(XTC_t)

            c.universes.add('main', GRO_t.strpath, XTC_t.strpath)

            py.path.local(c.basedir).chmod(0550, rec=True)

        def fin():
            py.path.local(c.basedir).chmod(0770, rec=True)

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
        py.path.local(sim.basedir).chmod(0770, rec=True)
        sim.location = tmpdir.mkdir('test').strpath
        py.path.local(sim.basedir).chmod(0550, rec=True)

        assert isinstance(sim.universe, MDAnalysis.Universe)

        py.path.local(sim.basedir).chmod(0770, rec=True)
