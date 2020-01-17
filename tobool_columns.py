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


def tobool_column(table, colname, trueval='true'):

    """
    Little helper function to convert columns (back) to bool type, e.g., when
    written by Topcat (which uses 'true' instead of 'True')
    """

    # --- test if the column is already in bool:
    if table[colname].dtype != bool:
        n = len(table)
        newcol = T.MaskedColumn(np.zeros(n, dtype=bool), name=colname)
        idt = table[colname] == trueval
        newcol[idt] = True
        table.replace_column(colname, newcol)


def tobool_columns(table, trueval='true', falseval='false'):

    """
    Little helper function to convert all columns (back) to bool type, e.g., when
    written by Topcat (which uses 'true' instead of 'True')
    """

    n = len(table)

    # --- go over all the columns and check which ones are to be converted
    for c in table.colnames:
        if (table[c].dtype == '<U5') | (table[c].dtype == '<U4'):
            if trueval in table[c] or falseval in table[c]:

                # --- replace the column with a bool column without masked values
                newcol = T.Column(np.zeros(n, dtype=bool), name=c)

                # --- old column masked or not?
                if np.ma.is_masked(table[c]):
                    idt = (np.invert(table[c].mask)) & (table[c] == trueval)
                else:
                    idt = (table[c] == trueval)

                # --- set trie values to true
                newcol[idt] = True

                # --- replace the old column with the new
                table.replace_column(c, newcol)
