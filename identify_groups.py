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
import datetime
import time


#%%
# --- Helper routine
def timestamp():
    """
    return a handy string with the current time
    """

    return(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))


#%%
# --- Helper routine
def identify_groups(idxl, idxr, nd1, nd2, match_sel="all", d2IDs=None,
                    verbose=False):
    """
    Internal helper routine to identify groups after an astropy coordinate
    crossmatch providing the indices idxl and idxr to the left and right tables
    of the match. This routine is used by the crossmatch_tables and
    internal_coordmatch routines.
    """

    n = len(idxl)

    lgroupids = np.ma.zeros(nd1, dtype=int)
    rgroupids = np.ma.zeros(nd2, dtype=int)

    rgroupids[:] = np.ma.masked
    lgroupids[:] = np.ma.masked


    lgroupsizes = np.ma.zeros(nd1, dtype=int)
    lgroupsizes[:] = np.ma.masked

    rgroupsizes = np.ma.zeros(nd2, dtype=int)
    rgroupsizes[:] = np.ma.masked


    # --- first look for all groups from the left side
    if match_sel != "best for 1":

        # --- get again all left side duplicates
        id_uni_l, lcounts = np.unique(idxl, return_counts=True)
        id_dupl_l = np.array(id_uni_l[lcounts > 1])
        n_dupl_l = len(id_dupl_l)

        print(timestamp() + "IDENTIFY_GROUPS: Left side:")

        for i in tqdm(range(n_dupl_l)):

            ids = np.where(idxl == id_dupl_l[i])[0]

        #    groupsize = len(ids)
            groupid = i+1

            lgroupids[idxl[ids]] = groupid
        #    dN["groupsize"][idxl[ids]] = groupsize

            rgroupids[idxr[ids]] = groupid
        #    dC["groupsize"][idxr[ids]] = groupsize

        #    print(i, dN["NEDname"][idxl[ids[0]]], groupid, " --> ", dC["CDSname"][idxr[ids]].data)   # , groupsize)


    # --- now go through all groups from the right side
        #print("\n\n ---- Right side: ")
    if match_sel != "best for 2":

        # --- get all right side duplicates
        if d2IDs is not None:
            id_uni_r, rcounts = np.unique(d2IDs, return_counts=True)
            id_dupl_r = np.array(id_uni_r[rcounts > 1])
            n_dupl_r = len(id_dupl_r)

        else:
            id_uni_r, rcounts = np.unique(idxr, return_counts=True)
            id_dupl_r = id_uni_r[rcounts > 1]
            n_dupl_r = len(id_dupl_r)

#            return(dC, idxr, select)
        print(timestamp() + "IDENTIFY_GROUPS: Right side:")

        for i in tqdm(range(n_dupl_r)):

            if d2IDs is not None:
                ids = np.where(d2IDs == id_dupl_r[i])[0]

            else:
                ids = np.where(idxr == id_dupl_r[i])[0]

            nr = len(ids)

            groupid = rgroupids[idxr[ids[0]]]

            # --- check whether right-side object is in multiple groups
            for j in range(nr):
                if rgroupids[idxr[ids[j]]] != groupid:

                    # --- previous element(s) not in a group so far
                    if groupid == 0:
                        groupid = rgroupids[idxr[ids[j]]]

                    # --- different group! --> merge
                    else:
                       if ((rgroupids[idxr[ids[j]]] < groupid)
                           & (rgroupids[idxr[ids[j]]] > 0)):
                           groupid = rgroupids[idxr[ids[j]]]

                       # --- we also need to set all other groups members to the new id
                       idso = np.where(rgroupids[idxr] == rgroupids[idxr[ids[j]]])[0]
                       lgroupids[idxl[idso]] = groupid
                       rgroupids[idxr[idso]] = groupid


            # --- if no existing group, new number
            if groupid ==0:
                groupid = n_dupl_l + i

#            print(i, len(lgroupids), len(idxl), ids)
            lgroupids[idxl[ids]] = groupid
        #    dN["groupsize"][idxl[ids]] = groupsize

            rgroupids[idxr[ids]] = groupid
        #    dC["groupsize"][idxr[ids]] = groupsize

        #    print(i, dC["CDSname"][idxr[ids[0]]], groupid)  # , groupsize)

        #    print(i, dC["CDSname"][idxr[ids[0]]], groupid, " --> ", dN["NEDname"][idxl[ids]].data)   # , groupsize)


    # --- resolve remaining group conflicts:
    if match_sel != "best for 2" and  match_sel != "best for 1":
        if verbose:
            print(timestamp() + ": IDENTIFY_GROUPS: Resolving group conflicts...")

        idb = np.where(lgroupids[idxl]!= rgroupids[idxr])[0]

        nconfl = len(idb)

        #for i in range(nmatch):
        #    if dN["groupID"][idxl[i]] != dC["groupID"][idxr[i]]:
        #        print(i, idxl[i], dN["groupID"][idxl[i]], dN["NEDname"][idxl[i]],
        #              "<-->", idxr[i], dC["groupID"][idxr[i]], dC["CDSname"][idxr[i]])

        maxiter = 10
        n_iter = 0

        while ((nconfl > 0) & (n_iter < maxiter)):

            for i in range(nconfl):
                # --- select the larger ID of the two conflicting
                group1 = lgroupids[idxl[idb[i]]]
                group2 = rgroupids[idxr[idb[i]]]

                selgroup = np.max([group1, group2])

                # --- change all groupIDs on the left side to the selected
                idsl = np.where((lgroupids[idxl] == group1)
                                | (lgroupids[idxl] == group2))[0]

                lgroupids[idxl[idsl]] = selgroup
                rgroupids[idxr[idsl]] = selgroup

                # --- change all groupIDs on the right side to the selected
                idsr = np.where((rgroupids[idxr] == group1)
                                | (rgroupids[idxr] == group2))[0]

                rgroupids[idxr[idsr]] = selgroup
                lgroupids[idxl[idsr]] = selgroup

            # --- check for any remaining conflicts
            idb = np.where(lgroupids[idxl] != rgroupids[idxr])[0]

            nconfl = len(idb)
            n_iter += 1


        if n_iter >= maxiter & nconfl > 0:

            print(timestamp() + ": IDENTIFY_GROUPS: ERROR: not all group conflicts could be resolved: ")

            for i in range(n):
                if lgroupids[idxl[i]] != rgroupids[idxr[i]]:
                    print(i, idxl[i], lgroupids[idxl[i]],
                          "<-->", idxr[i], rgroupids[idxr[i]])


    # --- clean up the groupIDs to run from 1 to max and determine their sizes
    if verbose:
        print(timestamp() + ": IDENTIFY_GROUPS: Cleaning groupIDs and determine groupsizes...")

    uniquegroups = np.unique(lgroupids[lgroupids > 0])
    ngroups = len(uniquegroups)


    for i in range(ngroups):
        ids = np.where(lgroupids[idxl] == uniquegroups[i])[0]

        rgroupids[idxr[ids]] = i+1
        lgroupids[idxl[ids]] = i+1

        # --- determine the group sizes
    #    groupsize = len(np.unique(idxl[ids])) + len(np.unique(idxr[ids]))
        lgroupsizes[idxl[ids]] = len(ids)  # --- we want the number of possible matches (not the number of members)
        rgroupsizes[idxr[ids]] = len(ids)

    if verbose:
        print(timestamp() + ": IDENTIFY_GROUPS: Groups...")
        print("                      - Number of groups: ", ngroups)
        print("                      - Number of objects involved: ", np.sum(np.invert(lgroupids.mask)))
        print("                      - Largest group size: ", np.ma.max(lgroupsizes))


    return(lgroupids, lgroupsizes, rgroupids, rgroupsizes)

