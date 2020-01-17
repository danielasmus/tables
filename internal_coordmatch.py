#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.0"

"""
HISTORY:
    - 2020-01-17: created by Daniel Asmus


NOTES:
    -

TO-DO:
    -
"""


import numpy as np
import astropy.table as T
from astropy import units as u
from astropy import coordinates

from .identify_groups import identify_groups as _identify_groups


def internal_coordmatch(d, racol, deccol, conrad, verbose=False):
    """
    Perform an internal crossmatch by sky coordinates on an astropy table d,
    whereas the column names for the coordinates have to be provided as well
    as the cone radius (conrad) in arcsec. Match rows are marked by a
    "groupID" column  of groupsize equal to the possible matches (not involved
    members).
    """

    ntot = len(d)

    coords = coordinates.SkyCoord(ra=d[racol]*u.degree, dec=d[deccol]*u.degree)

    idxr, idxl, sep, _ = coords.search_around_sky(coords, conrad*u.arcsec)

    sep = sep.value * 3600


    groupids, _, _, _ = _identify_groups(idxl, idxr, ntot, ntot, verbose=verbose)

    if "groupID" not in d.colnames:
        d.add_column(T.MaskedColumn(np.ma.zeros(ntot, dtype=int), name="groupID"))
    else:
        print("INTERNAL_COORDMATCH: WARNING: column 'groupID' already present. Will overwrite...")

    if "groupsize" not in d.colnames:
        d.add_column(T.MaskedColumn(np.ma.zeros(ntot, dtype=int), name="groupsize"))
    else:
        print("INTERNAL_COORDMATCH: WARNING: column 'groupsize' already present. Will overwrite...")


    d["groupID"][:] = np.ma.masked
    d["groupsize"][:] = np.ma.masked


    d["groupID"][idxl] = groupids[idxl]


    # --- the group size needs to be calculated a new because the one returned from
    #     the identify_groups routine is not applicable for an internal match
    ngroups = np.ma.max(d["groupID"])

    for i in range(ngroups):
        ids = np.where(d["groupID"] == i+1)[0]

        d["groupsize"][ids] = len(ids)

    if verbose:
        print("INTERNAL_COORDMATCH: Real largest groupsize: ", np.ma.max(d["groupsize"]))


