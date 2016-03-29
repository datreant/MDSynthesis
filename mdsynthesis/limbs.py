"""
Limbs are user interfaces for accessing stored data, as well as querying
the state of an object (data loaded, universe attached, etc.). They are also
used to aggregate the functionality of higher level objects (such as Sim) in
ways that are user-friendly.

In short, an Limb is designed to be user friendly on its own, but are
often used as components of a Treant.

"""
import os
from six import string_types
import numpy as np
import warnings

from datreant.core.limbs import Limb
from MDAnalysis import Universe
from MDAnalysis.core.AtomGroup import AtomGroup

from . import filesystem


class Selections(Limb):
    """Stored atom selections for the universe.

    Useful atom selections can be stored for the universe and recalled later.

    """
    _name = 'selections'

    def __repr__(self):
        return "<Selections({})>".format(
                {x: self.define(x) for x in self.keys()})

    def __str__(self):
        selections = self.keys()
        agg = "Selections"
        majsep = "="
        minsep = "-"
        subsep = "| "
        seplength = len(agg)

        if not self._treant._uname:
            out = "No universe attached; no Selections to show"
        elif not selections:
            out = "No selections for universe '{}'".format(
                self._treant._uname)
        else:
            out = agg + '\n'
            out = out + majsep * seplength + '\n'
            for selection in selections:
                out = out + "'{}'\n".format(selection)
                for item in self.define(selection):
                    out = out + subsep + "'{}'\n".format(item)
                out = out + minsep * seplength + '\n'

        return out

    def __getitem__(self, handle):
        """Get selection as an AtomGroup for given handle and the universe.

        Parameters
        ----------
        handle : str
            Name of selection to return as an AtomGroup.

        Returns
        -------
        AtomGroup
            The named selection as an AtomGroup of the universe.

        """
        return self.asAtomGroup(handle)

    def __setitem__(self, handle, selection):
        """Selection for the given handle and the active universe.

        """
        if isinstance(selection, (string_types, AtomGroup)):
            selection = [selection]
        self.add(handle, *selection)

    def __iter__(self):
        return self.keys().__iter__()

    def __delitem__(self, handle):
        """Remove stored selection for given handle and the active universe.

        """
        self.remove(handle)

    def add(self, handle, *selection):
        """Add an atom selection for the attached universe.

        AtomGroups are needed to obtain useful information from raw coordinate
        data. It is useful to store AtomGroup selections for later use, since
        they can be complex and atom order may matter.

        If a selection with the given *handle* already exists, it is replaced.

        Parameters
        ----------
        handle : str
            Name to use for the selection.
        selection : str, AtomGroup
            Selection string or AtomGroup; multiple selections may be given and
            their order will be preserved, which is useful for e.g. structural
            alignments.
        """
        # Conversion function, leave strings alone,
        # turn AtomGroups into their indices
        def conv(x):
            return x if isinstance(x, string_types) else x.indices

        selection = map(conv, selection)

        with self._treant._write:
            outsel = list()
            for sel in selection:
                if isinstance(sel, np.ndarray):
                    outsel.append(sel.tolist())
                elif isinstance(sel, string_types):
                    outsel.append(sel)

            seldict = self._treant._state['mdsynthesis']['sim']['selections']
            seldict[handle] = outsel

    def remove(self, *handle):
        """Remove an atom selection for the universe.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection(s) to remove.

        """
        with self._treant._write:
            seldict = self._treant._state['mdsynthesis']['sim']['selections']
            for item in handle:
                try:
                    del seldict[item]
                except KeyError:
                    raise KeyError("No such selection '{}'".format(item))

    def keys(self):
        """Return a list of all selection handles.

        """
        with self._treant._read:
            seldict = self._treant._state['mdsynthesis']['sim']['selections']
            return seldict.keys()

    def asAtomGroup(self, handle):
        """Get AtomGroup from universe from the given named selection.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection to return as an AtomGroup.

        Returns
        -------
        AtomGroup
            The named selection as an AtomGroup of the universe.
        """
        sel = self.define(handle)

        # Selections might be either
        # - a list of strings
        # - a list of indices

        ag = None
        for item in sel:
            if isinstance(item, string_types):
                if ag:
                    ag += self._treant.universe.select_atoms(item)
                else:
                    ag = self._treant.universe.select_atoms(item)
            else:
                if ag:
                    ag += self._treant.universe.atoms[item]
                else:
                    ag = self._treant.universe.atoms[item]

        return ag

    def define(self, handle):
        """Get selection definition for given handle.

        If named selection doesn't exist, :exc:`KeyError` raised.

        Parameters
        ----------
        handle : str
            Name of selection to get definition of.

        Returns
        -------
        definition : list
            list of strings defining the atom selection
        """
        with self._treant._read:
            seldict = self._treant._state['mdsynthesis']['sim']['selections']

        try:
            selstring = seldict[handle]
        except KeyError:
            raise KeyError("No such selection '{}'".format(handle))

        return selstring
