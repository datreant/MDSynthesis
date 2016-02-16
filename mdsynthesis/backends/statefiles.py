"""
Interface classes for state files.

"""

import os
import sys
import fcntl
import logging
import warnings
import json
from functools import wraps

from six import string_types
import numpy as np

from datreant.core.backends.core import FileSerial
from datreant.core.backends.statefiles import TreantFile


class SimFile(TreantFile):
    filepaths = ['abs', 'rel']

    def _init_record(self):
        super(SimFile, self)._init_record()
        self._record['mds'] = dict()
        self._record['mds']['universes'] = dict()
        self._record['mds']['default'] = None

    @FileSerial._read
    def get_mds_version(self):
        """Get Sim mdsynthesis version.

        :Returns:
            *version*
                mdsynthesis version of Sim

        """
        return self._record['mds']['version']

    # TODO: need a proper schema update mechanism
    @FileSerial._write
    def update_mds_schema(self):
        """Update mdsynthesis-specific schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            version = self._record['version']
        except KeyError:
            version = mdsynthesis.__version__

        return version

    @FileSerial._write
    def update_mds_version(self, version):
        """Update mdsynthesis version of Sim.

        :Arugments:
            *version*
                new mdsynthesis version of Treant
        """
        self._record['mds']['version'] = version

    @FileSerial._write
    def update_default(self, universe=None):
        """Mark the given universe as the default.

        :Arguments:
            *universe*
                name of universe to mark as default; if ``None``, remove
                default preference
        """
        self._record['mds']['default'] = universe

    @FileSerial._read
    def get_default(self):
        """Get default universe.

        :Returns:
            *default*
                name of default universe; if no default universe, returns
                ``None``

        """
        return self._record['mds']['default']

    @FileSerial._read
    def list_universes(self):
        """List universe names.

        :Returns:
            *universes*
                list giving names of all defined universes

        """
        return self._record['mds']['universes'].keys()

    @FileSerial._read
    def get_universe(self, universe):
        """Get topology and trajectory paths for the desired universe.

        Returns multiple path types, including absolute paths (abspath)
        and paths relative to the Sim object (relCont).

        :Arguments:
            *universe*
                given name for selecting the universe

        :Returns:
            *topology*
                dictionary containing all paths to topology
            *trajectory*
                dictionary containing all paths to trajectories

        """
        top = self._record['mds']['universes'][universe]['top']
        outtop = {key: value for key, value in zip(self.filepaths, top)}

        trajs = self._record['mds']['universes'][universe]['traj']
        outtraj = {key: [] for key in self.filepaths}

        for traj in trajs:
            for i, key in enumerate(self.filepaths):
                outtraj[key].append(traj[i])

        return outtop, outtraj

    @FileSerial._write
    def add_universe(self, universe, topology, *trajectory):
        """Add a universe definition to the Sim object.

        A Universe is an MDAnalysis object that gives access to the details
        of a simulation trajectory. A Sim object can contain multiple universe
        definitions (topology and trajectory pairs), since it is often
        convenient to have different post-processed versions of the same
        raw trajectory.

        :Arguments:
            *universe*
                given name for selecting the universe
            *topology*
                path to the topology file
            *trajectory*
                path to the trajectory file; multiple files may be given
                and these will be used in order as frames for the trajectory

        """
        # if universe schema already exists, don't overwrite it
        if universe not in self._record['mds']['universes']:
            self._record['mds']['universes'][universe] = dict()

        udict = self._record['mds']['universes'][universe]

        # add topology paths
        udict['top'] = [os.path.abspath(topology),
                        os.path.relpath(topology, self.get_location())]

        # add trajectory paths
        udict['traj'] = list()
        for segment in trajectory:
            udict['traj'].append(
                    [os.path.abspath(segment),
                     os.path.relpath(segment, self.get_location())])

        # add selections schema
        if 'sels' not in udict:
            udict['sels'] = dict()

        # add resnums schema
        if 'resnums' not in udict:
            udict['resnums'] = None

        # if no default universe, make this default
        if not self._record['mds']['default']:
            self._record['mds']['default'] = universe

    @FileSerial._write
    def del_universe(self, universe):
        """Delete a universe definition.

        Deletes any selections associated with the universe.

        :Arguments:
            *universe*
                name of universe to delete
        """
        del self._record['mds']['universes'][universe]

    @FileSerial._write
    def rename_universe(self, universe, newname):
        """Rename a universe definition.

        :Arguments:
            *universe*
                name of universe to rename
            *newname*
                new name of universe
        """
        if newname in self._record['mds']['universes']:
            raise ValueError(
                    "Universe '{}' already exixts;".format(newname) +
                    " remove it first")
        udicts = self._record['mds']['universes']
        udicts[newname] = udicts.pop(universe)

    @FileSerial._read
    def list_selections(self, universe):
        """List selection names.

        :Arguments:
            *universe*
                name of universe the selections apply to

        :Returns:
            *selections*
                list giving names of all defined selections for the given
                universe

        """
        return self._record['mds']['universes'][universe]['sels'].keys()

    @FileSerial._read
    def get_selection(self, universe, handle):
        """Get a stored atom selection for the given universe.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name to use for the selection

        :Returns:
            *selection*
                list of the selection strings making up the atom selection; may
                also be a list of atom indices

        """
        selections = self._record['mds']['universes'][universe]['sels']

        return selections[handle]

    @FileSerial._write
    def add_selection(self, universe, handle, *selection):
        """Add an atom selection definition for the named Universe definition.

        AtomGroups are needed to obtain useful information from raw coordinate
        data. It is useful to store AtomGroup selections for later use, since
        they can be complex and atom order may matter.

        Will overwrite existing definition if it exists.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name to use for the selection
            *selection*
                selection string or numpy array of indices; multiple selections
                may be given and their order will be preserved, which is
                useful for e.g. structural alignments

        """
        outsel = list()
        for sel in selection:
            if isinstance(sel, np.ndarray):
                outsel.append(sel.tolist())
            elif isinstance(sel, string_types):
                outsel.append(sel)

        self._record['mds']['universes'][universe]['sels'][handle] = outsel

    @FileSerial._write
    def del_selection(self, universe, handle):
        """Delete an atom selection from the specified universe.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name of the selection

        """
        del self._record['mds']['universes'][universe]['sels'][handle]

    @FileSerial._write
    def update_resnums(self, universe, resnums):
        """Update resnum definition for the given universe.

        Resnums are useful for referring to residues by their canonical resid,
        for instance that stored in the PDB. By giving a resnum definition
        for the universe, this definition can be applied to the universe
        on activation.

        Will overwrite existing definition if it exists.

        :Arguments:
            *universe*
                name of universe to associate resnums with
            *resnums*
                list giving the resnum for each residue in the topology, in
                index order
        """
        try:
            udict = self._record['mds']['universes'][universe]
        except KeyError:
            self.logger.info(
                "Universe definition '{}'".format(universe) +
                " does not exist. Add it first.")
            return

        # add resnums
        udict['resnums'] = resnums

    @FileSerial._read
    def get_resnums(self, universe):
        """Get the resnum definition for the given universe.

        :Arguments:
            *universe*
                name of universe the resnum definition applies to

        :Returns:
            *resnums*
                list of the resnums for each atom in topology; None if
                no resnums defined
        """
        return self._record['mds']['universes'][universe]['resnums']

    @FileSerial._write
    def del_resnums(self, universe):
        """Delete resnum definition from specified universe.

        :Arguments:
            *universe*
                name of universe to remove resnum definition from
        """
        self._record['mds']['universes'][universe]['resnums'] = None
