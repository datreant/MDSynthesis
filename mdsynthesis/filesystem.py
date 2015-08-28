"""
Functions and classes for finding Treants in the filesystem.

"""
import os
import glob

from datreant.treants import Treant, Group


class Universehound(object):
    """Locator for Universe files.

    This object is used by Sims to find their Universe files, even when
    they go missing.

    """
    def __init__(self, caller, uname):
        """Generate a Universehound to track down Universe files.

        :Arguments:
            *caller*
                object that summoned the Univershound
            *uname*
                universe handle to find files for

        """
        self.caller = caller
        self.uname = uname

        # once found: uuids as keys, absolute paths as values
        self.ufiles = list()

    def fetch(self):
        """Find the Universe files.

        For finding Universe files, the Foxhound begins by looking for
        the files with the same basename among the paths it was given.
        If the files can't be found in those places, it stops looking.

        Raises :exc:`IOError` if a file could not be found.

        :Returns:
            *results*
                dictionary giving 'top' and 'traj' as keys, with
                lists of topology and trajectory file(s) as
                values, respectively; ``None`` as a member of
                a list indicates that a particular file could not
                be found

        """
        # search last-known locations
        results = self._check_basedirs()

        if (None in results['top']) or (None in results['traj']):
            raise IOError("At least one file for" +
                          " universe '{}' could not".format(self.uname) +
                          " be found from stored absolute and relative" +
                          " paths. Re-add the universe with correct paths.")

        # TODO: include hash check; will require stored hashes for each
        # universe file

        return results

    def _check_basedirs(self):
        """Check last-known locations for Universe files.

        :Returns:
            *results*
                dictionary giving 'top' and 'traj' as keys, with
                lists of topology and trajectory file(s) as
                values, respectively; ``None`` as a member of
                a list indicates that a particular file could not
                be found
        """
        paths = self.caller._backend.get_universe(self.uname)
        paths = {'top': paths[0], 'traj': paths[1]}

        # initialize output list with None
        outpaths = dict()

        # iterate through topology and trajectory files
        for filetype in paths:
            outpaths[filetype] = [None]*len(paths[filetype])

            for i, entry in enumerate(paths[filetype]):

                # check absolute path first
                if 'abspath' in entry.dtype.names:
                    if os.path.exists(entry['abspath']):
                        outpaths[filetype][i] = entry['abspath']
                if 'relCont' in entry.dtype.names:
                    candidate = os.path.join(
                        self.caller._backend.get_location(), entry['relCont'])

                    # if both abspath and relCont exist, check that they point
                    # to the same file. if so, accept; if not, choose abspath
                    # and log a warning
                    if os.path.exists(candidate):
                        if outpaths[filetype][i]:
                            if not os.path.samefile(candidate,
                                                    outpaths[filetype][i]):
                                outpaths[filetype][i] = entry['abspath']
                                raise IOError(
                                    "Absolute and relative paths for a file" +
                                    " in universe '{}'".format(self.uname) +
                                    " point to different files; update paths" +
                                    " by re-adding this universe")
                        # otherwise, accept relCont
                        else:
                            outpaths[filetype][i] = candidate

        return outpaths
