#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exports results in Tecplot (TM) ASCII format.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def particles_tecplot(file_name, pt, n):
# -----------------------------------------------------------------------------
    file_id = open(file_name + ".plt", "w")
    col = 6
    # VisIt can't read Tecplot (TM) files which contain comments
    verbatim = False

    # --------------------------
    # Write the file header out
    # --------------------------
    if verbatim:
        file_id.write("# File header \n")
    file_id.write("title=\"PyNS Output\"\n")
    file_id.write("variables=\"xp\" \"yp\" \"zp\" \"up\" \"vp\" \"wp\" ")
    file_id.write("\n")
    file_id.write("zone i=%d"   % (n) +   \
                      " j=%d"   % (1) +   \
                      " k=%d\n" % (1))
    file_id.write("datapacking = block\n")
    file_id.write("varlocation=([1-6]=nodal)\n")

    # -------------------------------------------------------------------
    # Write the coordinates out (remember - those are nodal coordinates)
    # -------------------------------------------------------------------
    
    if verbatim:
        file_id.write("\n# X coordinates\n")

    c = 0 
    for i in range(0,n):
        file_id.write("%12.5e " % pt[i].x)
        c = c + 1
        if c % col == 0:
           file_id.write("\n")
           

    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    if verbatim:
        file_id.write("# Y coordinates\n")
    c = 0                                # column counter
           
    for i in range(0,n):
        file_id.write("%12.5e " % pt[i].y)
        c = c + 1
        if c % col == 0:
           file_id.write("\n")
    
    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    if verbatim:
        file_id.write("# Z coordinates\n")
    c = 0                                # column counter
           
    for i in range(0,n):
        file_id.write("%12.5e " % pt[i].z)
        c = c + 1
        if c % col == 0:
           file_id.write("\n")

    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    # ------------------------
    #
    # Write the variables out
    #
    # ------------------------
           
    c = 0 
    for i in range(0,n):
        file_id.write("%12.5e " % pt[i].u)
        c = c + 1
        if c % col == 0:
           file_id.write("\n")
           
    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    c = 0                                # column counter
           
    for i in range(0,n):
        file_id.write("%12.5e " % pt[i].v)
        c = c + 1
        if c % col == 0:
           file_id.write("\n")
    
    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")
    c = 0                                # column counter
           
    for i in range(0,n):
        file_id.write("%12.5e " % pt[i].w)
        c = c + 1
        if c % col == 0:
           file_id.write("\n")

    file_id.close()

    return  # end of function

