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


import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'


from .crossmatch_tables import crossmatch_tables
from .duplicates_in_column import duplicates_in_column
from .groups_by_value import groups_by_value
from .identify_groups import identify_groups
from .internal_coordmatch import internal_coordmatch
from .select_closest import select_closest
from .tobool_columns import tobool_column, tobool_columns
