#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.0"

"""
HISTORY:
    - 2020-01-15: created by Daniel Asmus


NOTES:
    -

TO-DO:
    -
"""


import numpy as np

from tqdm import tqdm


def select_closest(dsel, sepcolnames, namecheckcols=None):

    """
    go through the groups (marked after a coordinate crossmatch)
    in an astropy table, dsel, and mark those objects that are closest in terms of
    coordinates to the counterparts in the provided priority order.
    sepcolnames gives the column names for the counterpart separations on sky
    namecheckcols gives the column names for the object names of the
    counterparts to avoid the same object to be selected in case the crossmatch
    was an internal one.
    """

    ncp = len(sepcolnames)

    ngroups = np.max(dsel["groupID"])

    #dsel["review"][:] = True


    if not isinstance(sepcolnames, list):
        sepcolnames = [sepcolnames]


    # --- loop over all the groups
    for i in tqdm(range(ngroups)):

        idg = np.where(dsel["groupID"] == i+1)[0]

    #    print(i+1, len(idg), NEDnames_dupl[i])


        # --- in case there are more than one column given for the separation,
        #     loop over all of them
        for j in range(ncp):

            colname = sepcolnames[j]

            ids = np.where(np.invert(dsel[colname].mask[idg]))[0]

            # --- check whether a counterpart is available at all
            if len(ids) == 0:
                continue

            # --- if yes, then preselect the object that is closest to that
            #     counterpart
            idm = np.argmin(dsel[colname][idg])

            dsel["select"][idg[idm]] = True

            # --- check if it is the same counterpart though
            if namecheckcols is not None:

                colname = namecheckcols[j]

                for x in ids:
                    if dsel[colname][idg[x]] != dsel[colname][idg[idm]]:
                        break

            # --- if they are all the same, the situation seems clear
            if "review" in dsel.colnames:
                dsel["review"][idg] = False

            break


