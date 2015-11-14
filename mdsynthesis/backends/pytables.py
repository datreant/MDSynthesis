"""
Interface classes for state files and data files.

"""

import tables
import numpy as np
import os
from datreant.backends.core import File
from datreant.backends.pytables import TreantFile
import mdsynthesis

# max length in characters for all paths
pathlength = 511

# max character length of strings used for handles, tags, categories
namelength = 55


class SimFile(TreantFile):
    """Main Sim state file.

    This file contains all the information needed to store the state of a
    Sim object. It includes accessors, setters, and modifiers for all
    elements of the data structure, as well as the data structure definition.

    """
    class _MDSversion(tables.IsDescription):
        """Table definition for storing version number of file schema.

        All strings limited to hardcoded size for now.

        """
        # version of MDS file schema corresponds to allows future-proofing
        # of old objects so that formats of new releases can be automatically
        # built from old ones
        version = tables.StringCol(15)

    class _Default(tables.IsDescription):
        """Table definition for storing default universe preference.

        Stores which universe is marked as default.

        """
        default = tables.StringCol(namelength)

    class _Topology(tables.IsDescription):
        """Table definition for storing universe topology paths.

        Two versions of the path to a topology are stored: the absolute path
        (abspath) and the relative path from the Sim object's directory
        (relCont). This allows the Sim object to use some heuristically good
        starting points when trying to find missing files using Finder.

        """
        abspath = tables.StringCol(pathlength)
        relCont = tables.StringCol(pathlength)

    class _Trajectory(tables.IsDescription):
        """Table definition for storing universe trajectory paths.

        The paths to trajectories used for generating the Universe
        are stored in this table.

        See UniverseTopology for path storage descriptions.

        """
        abspath = tables.StringCol(255)
        relCont = tables.StringCol(255)

    class _Resnums(tables.IsDescription):
        """Table definition for storing resnums.

        """
        resnum = tables.UInt32Col()

    def __init__(self, filename, logger=None, **kwargs):
        """Initialize Sim state file.

        :Arguments:
            *filename*
                path to file
            *logger*
                logger to send warnings and errors to

        :Keywords:
            *name*
                user-given name of Treant object
            *coordinator*
                directory in which coordinator state file can be found [None]
            *categories*
                user-given dictionary with custom keys and values; used to
                give distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search

        """
        super(SimFile, self).__init__(filename, logger=logger, **kwargs)

    def create(self, **kwargs):
        """Build Sim data structure.

        :Keywords:
            *name*
                user-given name of Sim object
            *coordinator*
                directory in which Coordinator state file can be found
                [``None``]
            *categories*
                user-given dictionary with custom keys and values; used to give
                distinguishing characteristics to object for search
            *tags*
                user-given list with custom elements; used to give
                distinguishing characteristics to object for search

        .. Note:: kwargs passed to :meth:`create`

        """
        super(SimFile, self).create(treanttype='Sim', **kwargs)

        self._make_universegroup()
        try:
            self.get_default()
        except tables.NoSuchNodeError:
            self.update_default()

    @File._write
    def _make_universegroup(self):
        """Make universes and universe groups.

        Intended for file initialization.

        """
        try:
            group = self.handle.get_node('/', 'universes')
        except tables.NoSuchNodeError:
            group = self.handle.create_group('/', 'universes', 'universes')

    @File._write
    def _make_default_table(self):
        """Make table for storing default universe.

        Used only on file creation.

        """
        try:
            table = self.handle.get_node('/', 'default')
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'default', self._Default, 'default')

    @File._read
    def get_MDS_version(self):
        """Get Sim MDS version.

        :Returns:
            *version*
                MDS version of Treant

        """
        table = self.handle.get_node('/', 'mds_version')
        return table.cols.version[0]

    # TODO: need a proper schema update mechanism
    @File._write
    def update_MDS_schema(self):
        """Update MDS schema of file.

        :Returns:
            *version*
                version number of file's new schema
        """
        try:
            table = self.handle.get_node('/', 'mds_version')
            version = table.cols.version[0]
        except tables.NoSuchNodeError:
            version = mdsynthesis.__version__

        return version

    @File._write
    def update_MDS_version(self, version):
        """Update MDS version of Sim.

        :Arugments:
            *version*
                new MDS version of Treant
        """
        try:
            table = self.handle.get_node('/', 'mds_version')
            table.cols.version[0] = version
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'mds_version', self._MDSversion, 'mds_version')
            table.row['version'] = version
            table.row.append()

    @File._write
    def update_default(self, universe=None):
        """Mark the given universe as the default.

        :Arguments:
            *universe*
                name of universe to mark as default; if ``None``,
                remove default preference
        """
        try:
            table = self.handle.get_node('/', 'default')
            table.cols.default[0] = universe
        except tables.NoSuchNodeError:
            table = self.handle.create_table(
                '/', 'default', self._Default, 'default')
            table.row['default'] = universe
            table.row.append()

    @File._read
    def get_default(self):
        """Get default universe.

        :Returns:
            *default*
                name of default universe; if no default
                universe, returns ``None``

        """
        table = self.handle.get_node('/', 'default')
        default = table.cols.default[0]

        if default == 'None':
            default = None

        return default

    @File._read
    def list_universes(self):
        """List universe names.

        :Returns:
            *universes*
                list giving names of all defined universes

        """
        group = self.handle.get_node('/', 'universes')

        return group.__members__

    @File._read
    def get_universe(self, universe):
        """Get topology and trajectory paths for the desired universe.

        Returns multiple path types, including absolute paths (abspath)
        and paths relative to the Sim object (relCont).

        :Arguments:
            *universe*
                given name for selecting the universe

        :Returns:
            *topology*
                structured array containing all paths to topology
            *trajectory*
                structured array containing all paths to trajectory(s)

        """
        try:
            # get topology file
            table = self.handle.get_node('/universes/{}'.format(universe),
                                         'topology')
            topology = table.read()

            # get trajectory files
            table = self.handle.get_node('/universes/{}'.format(universe),
                                         'trajectory')
            trajectory = table.read()

        except tables.NoSuchNodeError:
            raise KeyError(
                    "No such universe '{}'; add it first.".format(universe))

        return (topology, trajectory)

    @File._write
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

        # build this universe's group; if it exists, do nothing
        try:
            group = self.handle.create_group(
                '/universes', universe, universe, createparents=True)
        except tables.NodeError:
            self.handle.remove_node(
                '/universes/{}'.format(universe), 'topology')
            self.handle.remove_node(
                '/universes/{}'.format(universe), 'trajectory')

        # construct topology table
        table = self.handle.create_table(
            '/universes/{}'.format(universe), 'topology', self._Topology,
            'topology')

        # add topology paths to table
        table.row['abspath'] = os.path.abspath(topology)
        table.row['relCont'] = os.path.relpath(topology, self.get_location())
        table.row.append()

        # construct trajectory table
        table = self.handle.create_table(
            '/universes/{}'.format(universe), 'trajectory', self._Trajectory,
            'trajectory')

        # add trajectory paths to table
        for segment in trajectory:
            table.row['abspath'] = os.path.abspath(segment)
            table.row['relCont'] = os.path.relpath(segment,
                                                   self.get_location())
            table.row.append()

        # construct selection group; necessary to catch NodError
        # exception when a Universe is re-added because selections are
        # maintained
        try:
            group = self.handle.create_group(
                '/universes/{}'.format(universe), 'selections', 'selections')
        except tables.NodeError:
            pass

    @File._write
    def del_universe(self, universe):
        """Delete a universe definition.

        Deletes any selections associated with the universe.

        :Arguments:
            *universe*
                name of universe to delete
        """
        try:
            self.handle.remove_node('/universes', universe, recursive=True)
        except tables.NoSuchNodeError:
            raise KeyError(
                    "No such universe '{}';".format(universe) +
                    " nothing to remove.")

    @File._write
    def rename_universe(self, universe, newname):
        """Rename a universe definition.

        :Arguments:
            *universe*
                name of universe to rename
            *newname*
                new name of universe
        """
        try:
            self.handle.rename_node('/universes', newname, name=universe)
        except tables.NoSuchNodeError:
            raise KeyError(
                    "No such universe '{}';".format(universe) +
                    " nothing to rename.")
        except tables.NodeError:
            raise ValueError(
                    "A universe '{}' already exists;".format(universe) +
                    " remove or rename it first.")

    @File._write
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
                list giving the resnum for each atom in the topology, in index
                order
        """
        try:
            table = self.handle.create_table(
                '/universes/{}'.format(universe), 'resnums', self._Resnums,
                'resnums')
        except tables.NoSuchNodeError:
            self.logger.info(
                "Universe definition '{}'".format(universe) +
                " does not exist. Add it first.")
            return
        except tables.NodeError:
            self.logger.info(
                "Replacing existing resnums for '{}'.".format(universe))
            self.handle.remove_node(
                '/universes/{}'.format(universe), 'resnums')
            table = self.handle.create_table(
                '/universes/{}'.format(universe), 'resnums', self._Resnums,
                'resnums')

        # add resnums to table
        for item in resnums:
            table.row['resnum'] = item
            table.row.append()

    @File._read
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
        try:
            table = self.handle.get_node(
                '/universes/{}'.format(universe), 'resnums')
            resnums = [x['resnum'] for x in table.iterrows()]
        except tables.NoSuchNodeError:
            resnums = None

        return resnums

    @File._write
    def del_resnums(self, universe):
        """Delete resnum definition from specified universe.

        :Arguments:
            *universe*
                name of universe to remove resnum definition from
        """
        self.handle.remove_node('/universes/{}'.format(universe), 'resnums')

    @File._read
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
        try:
            group = self.handle.get_node(
                '/universes/{}'.format(universe), 'selections')
        except tables.NoSuchNodeError:
            raise KeyError("No such universe '{}';".format(universe) +
                           " cannot copy selections.")

        return group.__members__

    @File._read
    def get_selection(self, universe, handle):
        """Get a stored atom selection for the given universe.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name to use for the selection

        :Returns:
            *selection*
                list of the selection strings making up the atom selection
        """
        try:
            table = self.handle.get_node(
                '/universes/{}/selections'.format(universe), handle)
            selection = [x for x in table.read()]
        except tables.NoSuchNodeError:
            raise KeyError(
                    "No such selection '{}'; add it first.".format(handle))

        return selection

    @File._write
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
        # construct selection table
        if isinstance(selection[0], np.ndarray):
            selection = selection[0]

        try:
            array = self.handle.create_array(
                '/universes/{}/selections'.format(universe), handle, selection,
                handle)
        except tables.NodeError:
            self.logger.info(
                "Replacing existing selection '{}'.".format(handle))
            self.handle.remove_node(
                '/universes/{}/selections'.format(universe), handle)
            table = self.handle.create_array(
                '/universes/{}/selections'.format(universe), handle, selection,
                handle)

    @File._write
    def del_selection(self, universe, handle):
        """Delete an atom selection from the specified universe.

        :Arguments:
            *universe*
                name of universe the selection applies to
            *handle*
                name of the selection

        """
        try:
            self.handle.remove_node(
                '/universes/{}/selections'.format(universe), handle)
        except tables.NoSuchNodeError:
            raise KeyError(
                    "No such selection '{}';".format(handle) +
                    " nothing to remove.")
