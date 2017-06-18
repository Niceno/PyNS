"""
Exports results in Tecplot (TM) ASCII format.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def tecplot_nodal_uvw(file_name, xyzn, uvwn):
# -----------------------------------------------------------------------------
    """
    Args:
      file_name: String containing name of the file to be created.
      xyzn:      Tuple containing one-dimensional arrays with "x", "y"
                 and "z" coordinates.
      vars:      Tuple containing "Unknowns" to be exported to Tecplot (TM).
                 Individual unknowns can be either collocated or staggered.

    Returns:
      none!
    """

    # Unpack tuples
    xn, yn, zn = xyzn
    un, vn, wn = uvwn

    # Compute cell resolutions (remember: cell resolutions)
    nx = len(xn)-1
    ny = len(yn)-1
    nz = len(zn)-1

    file_id = open(file_name + ".plt", "w")

    # VisIt can't read Tecplot (TM) files which contain comments
    verbatim = False

    # --------------------------
    # Write the file header out
    # --------------------------
    if verbatim:
        file_id.write("# File header \n")
    file_id.write("title=\"PyNS Output\"\n")
    file_id.write("variables=\"x\" \"y\" \"z\" \"u\" \"v\" \"w\" ")
    file_id.write("\n")
    file_id.write("zone i=%d"   % (nx+1) +   \
                      " j=%d"   % (ny+1) +   \
                      " k=%d\n" % (nz+1))
    file_id.write("datapacking = block\n")
    file_id.write("varlocation=([1-6]=nodal)\n")

    # -------------------------------------------------------------------
    # Write the coordinates out (remember - those are nodal coordinates)
    # -------------------------------------------------------------------
    if verbatim:
        file_id.write("\n# X coordinates\n")
    c = 0                                  # column counter
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e " % xn[i])
                c = c + 1
                if c % 4 == 0:             # go to new line after 4th column
                    file_id.write("\n")

    if verbatim:
        file_id.write("\n# Y coordinates\n")
    c = 0                                  # column counter
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e " % yn[j])
                c = c + 1
                if c % 4 == 0:             # go to new line after 4th column
                    file_id.write("\n")

    if verbatim:
        file_id.write("\n# Z coordinates\n")
    c = 0                                  # column counter
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e " % zn[k])
                c = c + 1
                if c % 4 == 0:             # go to new line after 4th column
                    file_id.write("\n")

    # ------------------------
    # Write the variables out
    # ------------------------

    # Average values to be written for staggered variables
    for vel in (un, vn, wn):
        c = 0
        for k in range(0, nz+1):
            for j in range(0, ny+1):
                for i in range(0, nx+1):
                    file_id.write("%12.5e " % vel[i,j,k])
                    c = c + 1
                    if c % 4 == 0:
                        file_id.write("\n")

    file_id.close()

    return  # end of function
