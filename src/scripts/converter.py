#! /usr/bin/env python
import mdsynthesis as mds
import MDAnalysis as mda
import json
from os import path
from glob import glob
import six


def convert(folder):
    json_file = glob(path.join(folder, 'Sim*'))[0]
    with open(json_file) as fh:
        old = json.load(fh)
    sim = mds.Sim(folder, categories=old['categories'], tags=old['tags'])

    old_sim = old['mdsynthesis']
    # update universe definition
    udef = old_sim['universedef']
    if udef['topology']:
        args = [udef['topology']['abspath'],]
        if udef['trajectory']:
            args.append([abspath for abspath, relpath in udef['trajectory']])
        u = mda.Universe(*args, **udef['kwargs'])
        sim.universe = u

    # update atom selection
    try:
        atom_sel = old_sim['atomselections']
        for key, val in six.iteritems(atom_sel):
            sim.atomselections[key] = val
    except KeyError:
        pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="Convert old sims to new sims.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "directories",
        metavar="DIRECTORY",
        nargs="+",
        help="one or more directories that are being processed")
    args = parser.parse_args()

    for dir in args.directories:
        convert(dir)
