import pytest
from glob import glob
from os import path
import json

from MDAnalysisTests.datafiles import PDB, XTC

import mdsynthesis as mds
from mdsynthesis.scripts.mds_06to1 import convert


@pytest.fixture
def old_sim():
    return {
        'mdsynthesis': {
            'atomselections': {
                'aa': 'protein',
                'bb': 'not protein'
            },
            'universedef': {
                'trajectory': [[
                    XTC,
                    path.basename(XTC)
                ]],
                'kwargs': {
                    'a': 1
                },
                'topology': {
                    'abspath': PDB,
                    'relpath': path.basename(PDB)
                }
            }
        },
        'categories': {
            'timestep': '1 ps',
            'protein': '1duf',
            'md-engine': 'gromacs 4.x.',
        },
        'tags': ['1duf', 'dna']
    }


def test_convert(old_sim, tmpdir):
    sim_folder = str(tmpdir)
    with open(path.join(sim_folder, 'Sim-uuid.json'), 'w') as fh:
        json.dump(old_sim, fh)

    convert(sim_folder)
    sim = mds.Sim(sim_folder)

    assert sim.tags == old_sim['tags']
    assert sim.categories == old_sim['categories']
    assert sim.atomselections['aa'] == old_sim['mdsynthesis'][
        'atomselections']['aa']
    assert sim.atomselections['bb'] == old_sim['mdsynthesis'][
        'atomselections']['bb']
    assert sim.universedef.trajectory == old_sim['mdsynthesis']['universedef'][
        'trajectory'][0][0]
    assert sim.universedef.topology == old_sim['mdsynthesis']['universedef'][
        'topology']['abspath']
    assert sim.universedef.kwargs == old_sim['mdsynthesis']['universedef'][
        'kwargs']


def test_convert_no_sim(tmpdir):
    sim_folder = str(tmpdir)
    with pytest.warns(UserWarning, match="No sim found.*"):
        convert(sim_folder)


def test_convert_to_many_sims(tmpdir, old_sim):
    sim_folder = str(tmpdir)
    with open(path.join(sim_folder, 'Sim-uuid.json'), 'w') as fh:
        json.dump(old_sim, fh)

    with open(path.join(sim_folder, 'Sim-uuid2.json'), 'w') as fh:
        json.dump(old_sim, fh)

    with pytest.warns(UserWarning, match="Multiple sims found.*"):
        convert(sim_folder)


def test_convert_no_atomselection(old_sim, tmpdir):
    del old_sim['mdsynthesis']['atomselections']

    sim_folder = str(tmpdir)
    with open(path.join(sim_folder, 'Sim-uuid.json'), 'w') as fh:
        json.dump(old_sim, fh)

    convert(sim_folder)
    sim = mds.Sim(sim_folder)

    assert sim.tags == old_sim['tags']
    assert sim.categories == old_sim['categories']
    assert sim.universedef.trajectory == old_sim['mdsynthesis']['universedef'][
        'trajectory'][0][0]
    assert sim.universedef.topology == old_sim['mdsynthesis']['universedef'][
        'topology']['abspath']
    assert sim.universedef.kwargs == old_sim['mdsynthesis']['universedef'][
        'kwargs']
