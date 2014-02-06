"""
Base classes for :mod:`MDSynthesis` objects.

"""
import os, sys
import time
import yaml
import cPickle
import logging
import glob
from uuid import uuid4
from multiprocessing import Process

metafile = 'metadata.yaml'
logfile = 'logfile.log'
datafile = 'data.pkl'
dbfile = 'MDSdatabase.yaml'

class ObjectCore(object):
    """Lowest-level mixin; functionality common to all MDSynthesis objects.
    
    """
    def __init__(self):
        """Low-level attribute initialization.

        """
        self.util = Utilities()

class ContainerCore(ObjectCore):
    """Mixin class for all Containers.

    The ContainerCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Container objects.

    """
    _metafile = metafile
    _logfile = logfile
    _datafile = datafile
    _dbfile = dbfile

    def __init__(self):
        """
        
        """
        super(ContainerCore, self).__init__()
        self.metadata = dict()              # information about object; defines base object
        self.data = dict()              # analysis data 'modular dock'
    
    def save(self, *args):
        """Save Container metadata and, if desired, analysis data instances.

        By providing names of loaded datasets as arguments, you can save the
        loaded versions of the data to their source files. This is useful if
        you have need to manually edit the data in a source file, as you can
        load it, make changes, then write it out.

        If no arguments are given, then no datasets are saved to their source
        files. Only metadata is saved.

        Metadata is also pushed to Database.

        If 'all' is in argument list, every dataset that is loaded is written
        to its file.

        :Arguments:
            *args*
                datasets to save

        """
        self.util.makedirs(self.metadata['basedir'])

        with self.util.open(os.path.join(self.metadata['basedir'], self._metafile), 'w') as f:
            yaml.dump(self.metadata, f)

        # update database
        self.push()

        if args:
            if 'all' in args:
                self._logger.info("Saving all loaded data into source files for '{}'...".format(self.metadata['name']))
                savelist = self.data
            else:
                self._logger.info("Saving selected data into source files for '{}'...".format(self.metadata['name']))
                savelist = args

            for i in savelist:
                self._logger.info("Saving {}...".format(i))
                with self.util.open(os.path.join(self.metadata['basedir'], '{}/{}'.format(i, self._datafile)), 'wb') as f:
                    cPickle.dump(self.data[i], f)
            self._logger.info("All selected data saved.")

    def refresh(self):
        """Reloads metadata from file.

        """
        metafile = os.path.join(self.metadata['basedir'], self._metafile)
        with self.util.open(metafile, 'r') as f:
            self.metadata = yaml.load(f)

    def load(self, *args, **kwargs):
        """Load data instances into object.

        If 'all' is in argument list, every available dataset is loaded.

        :Arguments:
            *args*
                datasets to load
            
        :Keywords:
            *force*
                if True, reload data even if already loaded; default False
        """

        force = kwargs.pop('force', False)

        if 'all' in args:
            self._logger.info("Loading all known data into object '{}'...".format(self.metadata['name']))
            loadlist = self.metadata['data']
        else:
            self._logger.info("Loading selected data into object '{}'...".format(self.metadata['name']))
            loadlist = args

        for i in loadlist:
            if (i not in self.data) or (force == True):
                self._logger.info("Loading {}...".format(i))
                with self.util.open(os.path.join(self.metadata['basedir'], '{}/{}'.format(i, self._datafile)), 'rb') as f:
                    self.data[i] = cPickle.load(f)
            else:
                self._logger.info("Skipping reload of {}...".format(i))
        self._logger.info("Object '{}' loaded with selected data.".format(self.metadata['name']))

    def unload(self, *args):
        """Unload data instances from object.

        If 'all' is in argument list, every loaded dataset is unloaded.

        :Arguments:
            *args*
                datasets to unload
        """
        if 'all' in args:
            self.data.clear()
            self._logger.info("Object '{}' unloaded of all data.".format(self.metadata['name']))
        else:
            self._logger.info("Unloading selected data from object {}...".format(self.metadata['name']))
            for i in args:
                self._logger.info("Unloading {}...".format(i))
                self.data.pop(i, None)
            self._logger.info("Object '{}' unloaded of selected data.".format(self.metadata['name']))

    def _makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)
    
    def _build_metadata(self, **kwargs):
        """Build metadata. Runs each time object is generated.
        
        Only adds keys; never modifies existing ones.

        :Keywords:
            *name*
                desired name of object, used for logging and referring to
                object in some analyses; default None
        """
        # building core items
        uuid = self._generate_uuid()
        attributes = {'uuid': uuid,
                      'name': kwargs.pop('name', None),
                      'data': list(),
                      'class': self.__class__.__name__,
                      'category': kwargs.pop('category', dict()),
                      'tag': kwargs.pop('tag', list()),
                      }

        for key in attributes:
            if not key in self.metadata:
                self.metadata[key] = attributes[key]
    
    def _start_logger(self):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.metadata['name']))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)
    
            # file handler
            if 'basedir' in self.metadata:
                logfile = os.path.join(self.metadata['basedir'], self._logfile)
                fh = logging.FileHandler(logfile)
                ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
                fh.setFormatter(ff)
                self._logger.addHandler(fh)
    
            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)
    
    def _generate_uuid(self):
        """Generate a 'unique' identifier.

        """
        return str(uuid4())

    def push(self):
        """Update metadata stored in Database for this Container.
    
        """
        # use restore_database if database file does not exist

    def pull(self):
        """Update metadata from information stored in Database for this Container.

        Note: This will overwrite metadata file with Database version!

        """

    def _init_database(self, database, **kwargs):
        """On generation of Container, perform standard Database interactions.

        If the Database specified already exists, this Container will interface
        with it. If it does not exist, a new Database will be created where the
        specified one should have been.

        If the *locate* keyword is set to True, then in the event it cannot find
        the specified Database, the Container will search directories upward
        from its location for a Database. If it finds one, it will interface
        with that Database. If it does not find one, then a new Database will
        be created where the specified one should have been.

        :Arguments:
            *database*
                directory of database

        :Keywords:
            *locate*
                if True, automatically try to find a database if not found ``[False]``
        
        """
        locate = kwargs.pop('locate', False)

        if database:
            database = os.path.abspath(database)
            self._logger.info("Attempting to connect to database in {}".format(database))
            dbname = self._connect_database(database)
        else:
            self._logger.info("No database specified.".format(database))
            dbname = None

        if not dbname:
            if locate:
                new_database = self._locate_database(startdir='.')
                if new_database:
                    database = new_database
            dbname = self._connect_database(database, new=True)

        if dbname:
            self._logger.info("Successfully connected to database {}.".format(dbname))
        else:
            self._logger.info("Could not connect to a database.".format(dbname))
            
    def _restore_database(self, database, **kwargs):
        """When Database can't be reached, find or make a new one.

        If the *locate* keyword is set to True, then in the event it cannot find
        the specified Database, the Container will search directories upward
        from its location for a Database. If it finds one, it will interface
        with that Database. If it does not find one, then a new Database will
        be created where the specified one should have been.

        :Arguments:
            *database*
                directory of database
        
        :Keywords:
            *locate*
                if True, automatically try to find a database if not found ``[False]``
        """
        new_database = self._locate_database(startdir='.')
        if new_database:
            database = new_database
        self._connect_database(database, new=True)
    
    def _connect_database(self, database, **kwargs):
        """Connect Container to a Database.

        If the Database doesn't exist, it can be created with the *new* keyword.

        :Arguments:
            *database*
                directory containing a Database 

        :Keywords:
            *new*
                if True, will make a new Database if an existing one is not
                found  ``[False]``

        :Returns:
            *dbname*
                name of database if connection success; otherwise None
        """
        new = kwargs.pop('new', False)

        if not database:
            return False

        # attempt to open existing database
        database = os.path.abspath(database)
        if os.path.exists(os.path.join(database, self._dbfile)):
            db = Database(database)
        
            if db._handshake():
                self._logger.info("Handshake success; database in {}".format(db.database['basedir']))
                self.metadata['database'] = db.database['basedir']
                db.add(self)
                dbname = db.database['name']
            else:
                self._logger.warning("Specified database failed handshake; not a real database?")
                dbname = None
        # make a new database in location
        elif new:
            db = Database(database)
            self.metadata['database'] = db.database['basedir']
            db.add(self)
            self._logger.info("Created new database in {}".format(db.database['basedir']))
            dbname = db.database['name']
    
        return dbname

    def _locate_database(self, **kwargs):
        """Find database; to be used if it can't be found.

        The Container looks upward from its location on the filesystem through
        the file heirarchy, looking for a Database file. The directory containing
        the first such file found will be returned. None is returned if no such
        files found.

        :Keywords:
            *startdir*
                directory from which to begin upward search; default is
                Container basedir

        :Returns:
            *database*
                directory of located Database; if no Database found, is None
        
        """
        startdir = kwargs.pop('startdir', None)
        
        if not startdir:
            startdir = self.metadata['basedir']

        # search upward for a database
        startdir = os.path.abspath(startdir)
        directory = startdir
        found = False
        
        self._logger.info("Beginning search for database from {}".format(directory))

        while (directory != '/') and (not found):
            directory, tail = os.path.split(directory)
            candidates = glob.glob(os.path.join(directory, self._dbfile))
            
            if candidates:
                self._logger.info("Database candidate located: {}".format(candidates[0]))
                basedir = os.path.dirname(candidates[0])
                db = Database(basedir)
                found = db._handshake()
        
        if not found:
            self._logger.warning("No database found!")
            basedir = None

        return basedir

class OperatorCore(ObjectCore):
    """Mixin class for all Operators.

    The OperatorCore object is not intended to be useful on its own, but
    instead contains methods and attributes common to all Operator objects.

    """
    _datafile = datafile

    def __init__(self, *args, **kwargs):
        """
        
        """
        super(OperatorCore, self).__init__()
        self.containers = list(args)

    def run(self, **kwargs):
        """Obtain compute-intensive data, usually timeseries.

        :Keywords:
            *force*
                If True, force recollection of data; default False

            **kwargs passed to `:meth:self._run_container()`
        """
        joblist = []
        force = kwargs.pop('force', False)

        # run analysis on each container as a separate process
        for container in self.containers:
            if (not self._datacheck(system)) or force:
                p = (Process(target=self._run_container, args=(container,), kwargs=kwargs))
                p.start()
                joblist.append(p)
            else:
                system._logger.info('{} data already present; skipping data collection.'.format(self.__class__.__name__))

        for p in joblist:
            p.join()
    
        # finish up
        for container in self.containers:
            # update analysis list in each object
            if not self.__class__.__name__ in container.metadata['data']:
                container.metadata['data'].append(self.__class__.__name__)
                container.save()

    def analyze(self, **kwargs):
        """Perform analysis of compute-intensive data.

        Does not require stepping through any trajectories.

        """
        # make sure data loaded into each container; should use try/catch here
        self._load()

    def _class(self):
        """Return class name.

        """
        return self.__class__.__name__

    def _save(self, container, cont_results):
        """Save results to main data file.

        :Arguments:
            *container*
                container to save data for
            *cont_results*
                results for container
        """
        outputdir = self._make_savedir(container)
        datafile = self._datafile(container)

        with self.util.open(datafile, 'wb') as f:
            cPickle.dump(cont_results, f)

    def _make_savedir(self, container):
        """Make directory where all output files are placed.

        :Arguments:
            *container*
                Container for which to save data 

        :Returns:
            *outputdir*
                full path to output directory

        """
        outputdir = self._outputdir(container)
        self.util.makedirs(outputdir)
        return outputdir
    
    def _load(self, **kwargs):
        """Load data for each container if not already loaded.

        :Keywords:
            *force*
                If True, force reload of data; default False
        """
        force = kwargs.pop('force', False)

        # make sure data loaded into each system; should use try/catch here
        for container in self.containers:
            if (not self.__class__.__name__ in container.data.keys()) or force:
                container.load(self.__class__.__name__)
    
    def _datacheck(self, container):
        """Check if data file already present.

        :Arguments:
            *container*
                Container object to check
                
        :Returns:
            *present*
                True if data is already present; False otherwise
        """
        return os.path.isfile(self._datafile(container))

    def _outputdir(self, container):
        """Return path to output directory for a particular Container.

        :Arguments:
            *container*
                Container object

        :Returns:
            *analysis_path*
                path to output directory
        """
        return os.path.join(container.metadata['basedir'], self.__class__.__name__)
    
    def _datafile(self, container):
        """Return path to main output datafile for a particular Container.

        :Arguments:
            *container*
                Container object

        :Returns:
            *datafile_path*
                path to datafile
        """
        return os.path.join(self._outputdir(container), self._datafile)

class Database(ObjectCore):
    """Database object for tracking and coordinating Containers.
    
    """
    _metafile = metafile
    _dbfile = dbfile
    _logfile = logfile

    def __init__(self, database, **kwargs):
        """Generate Database object for the first time, or interface with an existing one.

        :Arguments:
            *database*
                directory containing a Database file; if no Database file is
                found, a new one is generated

        """
        super(Database, self).__init__()
        self.database = dict()              # the database data itself

        database = os.path.abspath(database)
        dbfile = os.path.join(database, self._dbfile)
        if os.path.exists(dbfile):
            self._regenerate(database, **kwargs)
        else:
            self._generate(database, **kwargs)
    
    def _generate(self, database):
        """Generate a new database.
        
        """
        self.database['basedir'] = database
        self._build_metadata()
        self._build_attributes()

        # write to database file
        self.commit()
        self._start_logger()

    def _regenerate(self, database):
        """Re-generate existing database.
    
        """
        self.database['basedir'] = database
        self.refresh()
        self._start_logger(database)
        
        self._check_location(database)

        # rebuild missing parts
        self._build_metadata()
        self._build_attributes()

    def _handshake(self):
        """Run check to ensure that database is fine.

        """
        #TODO: add various checks to ensure things are in working order
        return ('basedir' in self.database)

    def search(self, searchstring):
        """Search the Database for Containers that match certain criteria.

        Results are printed in digest form to the screen. To print full
        metadata for all matching containers, use print='full'

        :Arguments:
            *searchstring*
                string giving the search query

        :Keywords:
            *print*
                format of results printed to ouptut

        :Returns:
            *locations*
                filesystem paths to Containers that match criteria

        """
        #TODO: Add in selection system similar to that implemented in
        # MDAnalysis for atom selections. This one, however, will parse
        # metadata elements, and shouldn't be quite so complex

        return

    def add(self, *containers, **kwargs):
        """Add Container to Database.

        :Arguments:
            *containers*
                Containers to add, each given as a path to a Container directory
                or as a generated Container object
            
        """
        for container in containers:
            if isinstance(container, basestring) and os.path.isdir(container):
                with self.util.open(os.path.join(container, self._metafile), 'r') as f:
                    meta = yaml.load(f)
                uuid = meta['uuid']
                meta['basedir'] = os.path.abspath(container)
                self.database['container'][uuid] = meta
                with self.util.open(os.path.join(container, self._metafile), 'w') as f:
                    yaml.dump(meta, f)
            else:
                uuid = container.metadata['uuid']
                self.database['container'][uuid] = container.metadata
    
            self.database['container'][uuid]['database'] = self.database['basedir']
            self.push(uuid)
            self._logger.info("Added {} container '{}' to database.".format(self.database['container'][uuid]['class'], self.database['container'][uuid]['name']))

    def remove(self, *containers, **kwargs):
        """Remove Container from Database.

        Note: if Container name is used to specify removal and more than one
        Container has that name, then both will be removed.

        :Arguments:
            *containers*
                Containers to remove, each given as a path to a Container directory,
                a Container UUID, or a Container's given name

        :Keywords:
            *hard*
                if True, delete Container object from filesystem too ``[False]``
            *all*
                if True, will remove all entries ``[False]``
        """
        all_conts = kwargs.pop('all', False)

        if all_conts:
            containers = [ x for x in self.database['container'] ]

        for container in containers:
            if os.path.isdir(container):
                basedir = os.path.abspath(container)
                contype = ['basedir']
            else:
                contype = ['uuid', 'name']

            matches = []
            for entry in self.database['container'].values():
                for criteria in contype:
                    if entry[criteria] == container:
                        matches.append(entry['uuid'])
    
            for match in matches:
                self.database['container'].pop(match, None)

    def clean(self):
        """Clear entries from Database corresponding to Containers that can't be found.

        """

    def commit(self):
        """Save the current state of the database to its file.
        
        """
        self.util.makedirs(self.database['basedir'])
        with self.util.open(os.path.join(self.database['basedir'], self._dbfile), 'w') as f:
            yaml.dump(self.database, f)

    def refresh(self):
        """Reload contents of database file.

        """
        dbfile = os.path.join(self.database['basedir'], self._dbfile)
        with self.util.open(dbfile, 'r') as f:
            self.database = yaml.load(f)

    def pull(self, *containers, **kwargs):
        """Update information stored in Database from Container metadata.

        Note: if Container name is used to specify the update, all Containers
        with that name will be updated in the Database.

        :Arguments:
            *args*
                Containers to update, each given as a path to a Container directory,
                a Container UUID, or a Container's given name

        :Keywords:
            *all*
                if True, will update entries for all known Containers from metadata files
        """
        all_conts = kwargs.pop('all', False)

        if all_conts:
            containers = [ x for x in self.database['container'] ]
    
        for container in containers:
            if os.path.isdir(container):
                basedir = os.path.abspath(container)
                contype = ['basedir']
            else:
                contype = ['uuid', 'name']

            matches = []
            for entry in self.database['container'].values():
                for criteria in contype:
                    if entry[criteria] == container:
                        matches.append(entry['uuid'])

            for match in matches:
                # ensure we are finding the right Container
                basedir = self._get_container(match)
                if not basedir:
                    self._logger.warning("Not found: {}".format(match))
                else:
                    with self.util.open(os.path.join(basedir, self._metafile), 'r') as f:
                        meta = yaml.load(f)

    def push(self, *containers, **kwargs):
        """Update Container metadata with information stored in Database.

        This is the opposite of `:meth:self.pull()`

        Note: if Container name is used to specify the update, all Containers
        with that name will have metadata updated.

        :Arguments:
            *containers*
                Containers to update; either a path to a Container directory,
                Container UUID, or a Container's given name
        :Keywords:
            *all*
                if True, will update all known Container metadata files from entries
        """
        all_conts = kwargs.pop('all', False)

        if all_conts:
            containers = [ x for x in self.database['container'] ]

        for container in containers:
            if os.path.isdir(container):
                basedir = os.path.abspath(container)
                contype = ['basedir']
            else:
                contype = ['uuid', 'name']

            matches = []
            for entry in self.database['container'].values():
                for criteria in contype:
                    if entry[criteria] == container:
                        matches.append(entry['uuid'])
    
            for match in matches:
                # since this method is used for Container init, basedir may not
                # be defined in metadata yet
                try:
                    # ensure we are finding the right Container
                    basedir = self._get_container(match)
                    if not basedir:
                        self._logger.warning("Not found: {}".format(match))
                    else:
                        with self.util.open(os.path.join(basedir, self._metafile), 'w') as f:
                            yaml.dump(self.database['container'][match], f)
                except KeyError:
                    self.database['container'][match]['basedir'] = self._build_basedir(match)
                    with self.util.open(os.path.join(self.database['container'][match]['basedir'], self._metafile), 'w') as f:
                        yaml.dump(self.database['container'][match], f)
            self.commit()

    def _get_container(self, uuid):
        """Get path to a Container.

        Will perform checks to ensure the Container returned matches the uuid given.
        It will go looking for the Container if not found at last known location.

        :Arguments:
            *uuid*
                unique id for Container to return

        :Returns:
            *container*
                path to Container
        """
        if os.path.exists(os.path.join(self.database['container'][uuid]['basedir'], self._metafile)):
            with self.util.open(os.path.join(self.database['container'][uuid]['basedir'], self._metafile), 'r') as f:
                meta = yaml.load(f)
            if meta['uuid'] == uuid:
                container = self.database['container'][uuid]['basedir']
            else:
                container = self._locate_container(uuid)
        else:
            container = self._locate_container(uuid)

        return container

    def discover(self):
        """Traverse filesystem downward from Database directory and add all new Containers found.
        
        """
        for root, dirs, files in os.walk(self.database['basedir']):
            if self._metafile in files:
                dirs = []
                self.add(root)
    
    def merge(self, database):
        """Merge another database's contents into this one.
        
        :Arguments:
            *database*
                path to database or Database object

        """

    def split(self, database):
        """Split selected Containers off of database into another.
        
        :Arguments:
            *database*
                path to destination database or Database object
        """

    def _check_location(self, database, **kwargs):
        """Check Database location; if changed, send new location to all Containers.

        :Keywords:
            *force*
                if True, new location sent to all Containers even if unchanged;
                default False
        """
        force = kwargs.pop('force', False)

        if (database != self.database['basedir']) or force:
            self.database['basedir'] = database

            # update entries first
            self.pull(all=True)

            for entry in self.database['container'].values():
                entry['database'] = self.database['basedir']
                            
            self.commit()
            self.push(all=True)

    def _build_metadata(self, **kwargs):
        """Build metadata. Runs each time object is generated.
        
        Only adds keys; never modifies existing ones.

        """
        attributes = {'class': self.__class__.__name__,
                      'name': kwargs.pop('name', os.path.basename(self.database['basedir'])),
                      'container': dict(),
                      }
    
        for key in attributes:
            if not key in self.database:
                self.database[key] = attributes[key]

    def _build_attributes(self):
        """Build attributes of Database. Runs each time object is generated.

        """

    def _locate_container(self, uuid):
        """Find a Container by traversing downward through the filesystem. 

        Looks in each directory below the Database. If found, the basedir for the
        Container is updated in both metadata and the Database.

        :Arguments:
            *uuid*
                unique id for Container to return
        """
        self._logger.info("Searching for Container: {} ({})".format(self.database['container'][uuid]['name'], uuid))
        container = None
        for root, dirs, files in os.walk(self.database['basedir']):
            if self._metafile in files:
                dirs = []
                with self.util.open(os.path.join(root, self._metafile), 'r') as f:
                    meta = yaml.load(f)
                if meta['uuid'] == uuid:
                    container = os.path.abspath(root)
                    meta['basedir'] = container
                    self._logger.info("Found: {}".format(container))

                    # update basedir in Container metadata and in Database
                    with self.util.open(os.path.join(root, self._metafile), 'w') as f:
                        yaml.dump(meta, f)
                    self.database['container'][uuid]['basedir'] = container

        if not container:
            self._logger.warning("Could not find Container.")
                
        return container

    def _build_basedir(self, uuid):
        """Build basedir location based on database location, Container class, and Container name.

        :Arguments:
            *database*
                directory where database resides
            *name*
        """
        database = self.database['basedir']
        container = self.database['container'][uuid]

        # if name given and directory with name doesn't already exist, make named basedir
        if container['name'] and not os.path.exists(os.path.join(database, container['class'], container['name'])):
            dest = container['name']
        # if basedir already exists, use UUID instead
        else:
            dest = container['uuid']

        dest = os.path.join(database, container['class'], dest)
        self.util.makedirs(dest)

        return dest

    def _start_logger(self, basedir):
        """Start up the logger.

        """
        # set up logging
        self._logger = logging.getLogger('{}.{}'.format(self.__class__.__name__, self.database['name']))

        if not self._logger.handlers:
            self._logger.setLevel(logging.INFO)

            # file handler
            logfile = os.path.join(basedir, self._logfile)
            fh = logging.FileHandler(logfile)
            ff = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
            fh.setFormatter(ff)
            self._logger.addHandler(fh)

            # output handler
            ch = logging.StreamHandler(sys.stdout)
            cf = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            ch.setFormatter(cf)
            self._logger.addHandler(ch)

class Utilities(object):
    """Lowest level utilities; contains all methods that are common to every 
       MDSynthesis object.

    """
    def open(self, *args, **kwargs):
        """Open a file for i/o and apply an exclusive lock.
    
        Arguments and keywords are passed directly to the open() builtin.
    
        """
        F = File(*args, **kwargs)
        return F
    
    def makedirs(self, p):
        if not os.path.exists(p):
            os.makedirs(p)
    
class File(object):
    """File object class. Implements needed file locking.

    """
    def __init__(self, *args, **kwargs):
        """Create File instance for loading files and consistently handling locking.

        """
        self.file = open(*args, **kwargs)
        self.lockname = "{}.lock".format(self.file.name) 
        self.lock()

    def __enter__(self):
        return self.file

    def __exit__(self, *args):
        self.close()

    def lock(self):
        """Get lock.
    
        """
        # if lockfile already present, wait until it disappears
        while os.path.exists(self.lockname):
            time.sleep(5)

        # open lockfile (will appear in filesystem)
        self.lockfile = open(self.lockname, 'w')

    def close(self):
        """Close file stream.
        
        """
        self.file.close()
        self.lockfile.close()
        os.remove(self.lockname)

class Attributes(object):
    """Class for user-defined attributes.

    """
    def __init__(self):
        pass