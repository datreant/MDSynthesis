#! /usr/bin/env python
import mdsynthesis as mds
import json
from os import path
from glob import glob

def convert(folder):
    json_file = glob(path.join(folder, 'Sim*'))[0]
    with open(json_file) as fh:
        old = json.load(fh)
    mds.Sim(folder, categories=old['categories'], tags=old['tags'])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Convert old sims to new sims.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("directories", metavar="DIRECTORY", nargs="+", help="one or more directories that are being processed")
    args = parser.parse_args()

for dir in args.directories:
    convert(dir)
