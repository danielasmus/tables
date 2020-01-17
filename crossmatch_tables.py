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
import sys
from tqdm import tqdm
from astropy.io import ascii
import astropy.table as T
from astropy import units as u
from astropy import coordinates
from astroquery.xmatch import XMatch
import datetime
import time

from .identify_groups import identify_groups as _identify_groups

#%%
# --- Helper routine
def timestamp():
    """
    return a handy string with the current time
    """

    return(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))



#%%

def crossmatch_tables(d1, d2, fout=None, join_type="1 and 2", ra1='RA_deg',
                      dec1="DEC_deg", ra2='RA_deg', vizier=False,
                      dec2="DEC_deg", conrad=3.0, verbose=False,
                      d2_unicol=None, match_sel="best for both",
                      mark_groups=True):

    """
    Crossmatch two astropy tables by coordinates. At the moment all matches
    are considered not only the closest.
    If the vizier flag is set, then d1 has to be the local
    table, while d2 has to be the vizier ID of the online table (in this case,
    obviously, the options "1 or 2" and "2 not 1" do not work).
    """


    # --- make sure the tables support masking
    dN = T.Table(d1, masked=True)
    nN = len(dN)

    ncolsN = len(dN.columns)

#    print(dN.columns)

    join_type = join_type.lower()
    match_sel = match_sel.lower()

    joins = ["1 and 2", "1 not 2", "2 not 1", "1 or 2", "all from 1", "all from 2"]
    matches = ["best for both", "best for 1", "best for 2", "all"]

    if join_type not in joins:
        print("CROSSMATCH_TABLES: ERROR: selected join type not valid: ", join_type)
        print("  - available types are: ", joins)
        return(-1)

    if match_sel not in matches:
        print("CROSSMATCH_TABLES: ERROR: selected match selection not valid: ", match_sel)
        print("  - available selections are: ", matches)
        return(-1)

    if not vizier:

        dC = T.Table(d2, masked=True)

        nC = len(dC)
        ncolsC = len(dC.columns)

        if np.sum(dC[ra2].mask) + np.sum(dC[dec2].mask) >0:
            print("CROSSMATCH_TABLES: ERROR: Table 2 contains empty coordinates! Aborting...")
            return(-1)

        coordsN = coordinates.SkyCoord(ra=dN[ra1]*u.degree, dec=dN[dec1]*u.degree)
        coordsC = coordinates.SkyCoord(ra=dC[ra2]*u.degree, dec=dC[dec2]*u.degree)


    if verbose:
         print(timestamp() + ": CROSSMATCH_TABLES: Input parameters: ")
         print("                       - input rows (left): ", nN)
         if not vizier:
             print("                       - input rows (right): ", nC)
         print("                       - join type: ", join_type)
         print("                       - match selection: ", match_sel)
         print("                       - cone radius: ", conrad)
         print("                       - vizier XMATCH: ", vizier)
         print("                       - RA 1: ", ra1)
         print("                       - DEC 1: ", dec1)
         print("                       - RA 2: ", ra2)
         print("                       - DEC 2: ", dec2)
         print("                       - output file: ", fout)


    # --- the crossmatch can not handle masked coordinates!
    if np.sum(dN[ra1].mask) + np.sum(dN[dec1].mask) >0:
        print("CROSSMATCH_TABLES: ERROR: Table 1 contains empty coordinates! Aborting...")
        return(-1)


    if verbose:
        print(timestamp() + ": CROSSMATCH_TABLES: Doing the crossmatch...")
        sys.stdout.flush()


    # --- in case of the online XMATCH with vizier table:
    if vizier:

        if join_type == "1 or 2" or join_type == "2 not 1":
            print("CROSSMATCH_TABLES: ERROR: For XMATCH with VISIR '1 or 2' and '2 not 1' are not available! Aborting...")
            return(-1)

        # --- crerate small table for upload (the columns can not be of masked
        #     type either)
        dtemp = T.Table({"UniID": range(nN), ra1: np.array(dN[ra1]),
                      dec1: np.array(dN[dec1])})

        dC = XMatch.query(cat1=dtemp, cat2=d2, max_distance=conrad * u.arcsec,
                           colRA1=ra1, colDec1=dec1)

        idxl = np.array(dC["UniID"])
        idxr = np.array(range(len(dC)))
        sep = dC["angDist"]

        nC = len(dC)

        del dC[ra1]
        del dC[dec1]
        del dC["UniID"]
        del dC["angDist"]

        ncolsC = len(dC.columns)

        # --- now we need to reconstruct the right side table
        if d2_unicol is not None:

            d2IDs = np.array(dC[d2_unicol])
            d2unique = np.unique(d2IDs)
            nrunique = len(d2unique)

        else:
            d2IDs = None


    else:
        idxr, idxl, sep, _ = coordsN.search_around_sky(coordsC, conrad*u.arcsec)

        nrunique = len(np.unique(idxr))
        sep = sep.value * 3600
        d2IDs = None


    nmatch = len(idxr)

    #nmatch_vol = np.sum(dN['NEDredshift'][idxr] < zlim)

    lunique = np.unique(idxl)


    nlunique = len(lunique)


    if verbose:
        print(timestamp() + ": CROSSMATCH_TABLES: Total number of found matches: ", nmatch)
        print("                     CROSSMATCH_TABLES: Number of unique left sources match: ", nlunique)
        print("                     CROSSMATCH_TABLES: Number of unique right sources match: ", nrunique)

    # --- if we are only interested in the non-matches we are done here:
    if join_type == "1 not 2":

        # --- ids of those not in the match
        idNnomatch = [x for x in range(nN) if x not in idxl]
        dout = dN[idNnomatch]

        if fout is not None:
            dout.write(fout, delimiter=',', format='ascii',
                       fill_values=[(ascii.masked, '')], overwrite=True)

        if verbose:
            print("CROSSMATCH_TABLES: Number of non-matches from 1: ", len(dout))

        return(dout)

    elif join_type == "2 not 1":

         # --- ids of those not in the match
        idCnomatch = [x for x in range(nC) if x not in idxr]
        dout = dC[idCnomatch]

        if fout is not None:
            dout.write(fout, delimiter=',', format='ascii',
                       fill_values=[(ascii.masked, '')], overwrite=True)

        if verbose:
            print("CROSSMATCH_TABLES: Number of non-matches from 2: ", len(dout))

        return(dout)


    # --- Groups and Duplicates

    # --- get allleft side duplicates
    id_uni_l, lcounts = np.unique(idxl, return_counts=True)
    id_dupl_l = np.array(id_uni_l[lcounts > 1])
    n_dupl_l = len(id_dupl_l)

    # --- get all right side duplicates
    if vizier and d2_unicol is not None:

        id_uni_r, rcounts = np.unique(d2IDs, return_counts=True)
        id_dupl_r = np.array(id_uni_r[rcounts > 1])
        n_dupl_r = len(id_dupl_r)

    else:
        id_uni_r, rcounts = np.unique(idxr, return_counts=True)
        id_dupl_r = id_uni_r[rcounts > 1]
        n_dupl_r = len(id_dupl_r)



    # --- In case of best matches only select correspondingly

    select = np.ones(nmatch, dtype=bool)

    if match_sel == "best for both" or match_sel == "best for 1":

        if verbose:
            print(timestamp() + "CROSSMATCH_TABLES: Finding best match on the left side...")

        for i in tqdm(range(n_dupl_l)):
            ids = np.where(idxl == id_dupl_l[i])[0]

            exclude = np.where(sep[ids] > np.nanmin(sep[ids]))[0]

            keep = np.where(sep[ids] == np.nanmin(sep[ids]))[0]

            if len(keep) > 1:
                if verbose:
                    print(" - WARNING: more than 1 object at minimum separation: ", len(keep))

            select[ids[exclude]] = False

#            print(i, idxl[ids], np.nanmin(sep[idxl[ids]]), exclude, idxl[ids[exclude]])


    if match_sel == "best for both" or match_sel == "best for 2":

        if verbose:
            print(timestamp() + "CROSSMATCH_TABLES: Finding best match on the right side...")

        for i in tqdm(range(n_dupl_r)):
            ids = np.where(idxr == id_dupl_r[i])[0]

            exclude = np.where(sep[ids] > np.nanmin(sep[ids]))[0]

            select[ids[exclude]] = False


    if verbose:
        print("CROSSMATCH_TABLES: Number of excluded matches: ", np.sum(np.invert(select)))

    idxl = idxl[select]
    idxr = idxr[select]
    sep = sep[select]

    if vizier:
#        idxr = np.array(range(len(idxr)))
        d2IDs = d2IDs[select]

    # --- Group business is only requierd if multiple matches are allowed
    if match_sel != "best for both" and mark_groups:

        # --- Add the columns of the tables to each other
        if verbose:
            print(timestamp() + ": CROSSMATCH_TABLES: Adding columns...")

        for i in range(ncolsN):
            if dN.columns[i].name not in dC.columns:
                dC.add_column(T.MaskedColumn(np.ma.zeros(nC, dtype=dN.columns[i].dtype),
                                             name=dN.columns[i].name))

                dC[dN.columns[i].name] = np.ma.masked

    #            print(" Added column: ", dN.columns[i].name)

        dC.add_column(T.MaskedColumn(np.ma.zeros(nC, dtype=float), name="separation_as"))
        dC["separation_as"] = np.ma.masked

        dC.add_column(T.MaskedColumn(np.ma.zeros(nC, dtype=int), name="groupID"))
        dC.add_column(T.MaskedColumn(np.ma.zeros(nC, dtype=int), name="groupsize"))

        for i in range(ncolsC):
            if dC.columns[i].name not in dN.columns:
                dN.add_column(T.MaskedColumn(np.ma.zeros(nN, dtype=dC.columns[i].dtype),
                                             name=dC.columns[i].name))

                dN[dC.columns[i].name] = np.ma.masked

#        print(dN.columns)
        dN.add_column(T.MaskedColumn(np.ma.zeros(nN, dtype=float), name="separation_as"))
        dN["separation_as"] = np.ma.masked

        dN.add_column(T.MaskedColumn(np.ma.zeros(nN, dtype=int), name="groupID"))
        dN.add_column(T.MaskedColumn(np.ma.zeros(nN, dtype=int), name="groupsize"))


        # --- now identify the groups
        if verbose:
            print(timestamp() + ": CROSSMATCH_TABLES: Identifying groups...")


        lgroupids, lgroupsizes, rgroupids, rgroupsizes = _identify_groups(idxl,
                                               idxr, nN, nC, match_sel=match_sel,
                                               d2IDs=d2IDs, verbose=verbose)


        dN['groupID'][idxl] = lgroupids[idxl]
        dC['groupID'][idxr] = rgroupids[idxr]

        dN['groupsize'][idxl] = lgroupsizes[idxl]
        dC['groupsize'][idxr] = rgroupsizes[idxr]


    # --- then we build a table with the matches
    dmatch = dN[idxl]

    for i in range(ncolsC):
        dmatch[dC.columns[i].name] = dC[dC.columns[i].name][idxr]

    # --- fill in the separations
    dmatch['separation_as'] = sep


    #for i in range(len(dmatch)): print(dmatch['NEDname'][i], " <--> ", dmatch['CDSname'][i])


    if verbose:
        print(timestamp() + ": CROSSMATCH_TABLES: Preparing output...")

    # --- now the output options:
    if join_type == "1 or 2":

        # --- ids of those not in the match
        idNnomatch = [x for x in range(nN) if x not in idxl]
        idCnomatch = [x for x in range(nC) if x not in idxr]
        dout = T.vstack([dN[idNnomatch], dmatch, dC[idCnomatch]])

    elif join_type == "all from 1":

         # --- ids of those not in the match
        idNnomatch = [x for x in range(nN) if x not in idxl]
        dout = T.vstack([dN[idNnomatch], dmatch])

    elif join_type == "all from 2":

        # --- ids of those not in the match
        idCnomatch = [x for x in range(nC) if x not in idxr]
        dout = T.vstack([dC[idCnomatch], dmatch])

    else:  #  out == "1 and 2"

        dout = dmatch


    # --- mask the empy values
    if 'groupID' in dout.colnames:
        idnull = dout['groupID'] == 0
        dout['groupID'][idnull] = np.ma.masked
        dout['groupsize'][idnull] = np.ma.masked


    if fout is not None:
        dout.write(fout, delimiter=',', format='ascii',
                   fill_values=[(ascii.masked, '')], overwrite=True)

    if verbose:

        print(timestamp() + ": CROSSMATCH_TABLES: Number of output lines: ", len(dout))

    return(dout, idxl, idxr)
