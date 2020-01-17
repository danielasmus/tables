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
import astropy.table as T


def groups_by_value(d, colname):

    """
    find all the groups in a given astropy table by having the same value in a
    given colume. Group ID and size are added as new columns to the subset
    table with all the rows that beleong to groups that is returned.
    """
    uni, ncounts = np.unique(d[colname], return_counts=True)
    dupl = np.array(uni[ncounts > 1])

    ntot = len(d)

    # --- for unknown reasons, a "0" gets inserted into the name dupl list in
    #     the above command
    dupl = dupl[dupl != "0"]

    ngroups = len(dupl)

    select = np.zeros(ntot, dtype=bool)
    groupid = np.zeros(ntot, dtype=int)
    groupsize = np.zeros(ntot, dtype=int)


    for i in tqdm(range(ngroups)):

        ids = np.where(d[colname] == dupl[i])[0]

        select[ids] = True
        groupid[ids] = i+1
        groupsize[ids] = len(ids)


    dsel = d[select]

    dsel.add_column(T.Column(groupid[select], name="groupID"))
    dsel.add_column(T.Column(groupsize[select], name="groupsize"))

    return(dsel, ngroups)


