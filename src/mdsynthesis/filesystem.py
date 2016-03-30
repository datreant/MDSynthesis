"""
Functions and classes for finding Treants in the filesystem.

"""
import os
import glob

from datreant.core import Treant, Group


class Universehound(object):
    """Locator for Universe files.

    This object is used by Sims to find their Universe files, even when
    they go missing.

    """
    def __init__(self, caller):
        """Generate a Universehound to track down Universe files.

        :Arguments:
            *caller*
                object that summoned the Univershound

        """
        self.caller = caller
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

        if (results['top'] is None) or (None in results['traj']):
            raise IOError("At least one file for"
                          " universe could not"
                          " be found from stored absolute and relative"
                          " paths. Set correct paths manually.")

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
        top = dict()
        traj = dict()
        for pathtype in ('abs', 'rel'):
            top[pathtype], traj[pathtype] = self.caller._define(
                    pathtype=pathtype)
        outpaths = dict()

        def check(paths):

            out = {'abs': None,
                   'rel': None}

            chosen = None

            # check absolute path first
            if os.path.exists(paths['abs']):
                out['abs'] = paths['abs']

            candidate = os.path.join(self.caller._treant.abspath, paths['rel'])
            if os.path.exists(candidate):
                out['rel'] = candidate

            # if both abs and rel exist, check that they point to the same
            # file. if so, accept; if not, throw exception
            if ((out['abs'] and out['rel']) and not
                    os.path.samefile(out['abs'], out['rel'])):
                raise IOError("Absolute and relative paths for file '{}'"
                              " point to different files.".format(out['abs']))
            elif out['rel']:
                # otherwise, accept rel
                chosen = os.path.abspath(out['rel'])
            elif out['abs']:
                # if rel path gives no file, accept abs
                chosen = out['abs']
            else:
                # if none of the paths resolve, raise exception
                raise IOError(
                        "File '{}' not found for universe".format(out['abs']))

            return chosen

        outpaths['top'] = check(top)

        outpaths['traj'] = [None]*len(traj['abs'])
        for i, entry in enumerate(traj['abs']):
            outpaths['traj'][i] = check({'abs': traj['abs'][i],
                                         'rel': traj['rel'][i]})

        return outpaths
