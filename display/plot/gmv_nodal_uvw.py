"""
Exports results in GMV (TM) ASCII format.

Source:
  http://www.generalmeshviewer.com/doc.color.pdf
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def gmv_nodal_uvw(file_name, xyzn, uvwn):
# -----------------------------------------------------------------------------
    """
    Args:
      file_name: String containing name of the file to be created.
      xyzn:      Tuple containing one-dimensional arrays with "x", "y"
                 and "z" coordinates.
      vars:      Tuple containing "Unknowns" to be exported to GMV (TM).
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

    file_id = open(file_name + ".gmv", "w")

    # --------------------------
    # Write the file header out
    # --------------------------
    file_id.write("gmvinput ascii\n")

    # -------------------------------------------------------------------
    # Write the coordinates out (remember - those are nodal coordinates)
    # -------------------------------------------------------------------
    file_id.write("nodev %d\n" % ((nx+1)*(ny+1)*(nz+1)))
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e %12.5e %12.5e\n" % (xn[i], yn[j], zn[k]))

    # --------------------
    # Write the cells out
    # --------------------
    file_id.write("cells %d\n" % (nx*ny*nz))

    #   cells' local numbering
    #
    #     ^ z
    #     |
    #     | n011-------n111
    #     | /|         /|
    #     |/ |        / |
    #    n001-------n101|
    #     |  |       |  |
    #     |  | / y   |  |
    #     |  |/      |  |
    #     | n010-----|-n110
    #     | /        | /
    #     |/         |/
    #    n000-------n100 ---->
    #                          x
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
               n000 = 1 + i + j * (nx+1) + k * ((nx+1)*(ny+1))
               n100 = n000 + 1
               n010 = n000 + (nx+1)
               n110 = n000 + (nx+1) + 1
               n001 = n000 + (nx+1)*(ny+1)
               n101 = n001 + 1
               n011 = n001 + (nx+1)
               n111 = n001 + (nx+1) + 1

               file_id.write("  hex 8\n")
               file_id.write("  %d %d %d %d %d %d %d %d\n"  \
                             % (n000,n100,n110,n010, n001,n101,n111,n011))

    # ------------------------
    # Write the variables out
    # ------------------------
    file_id.write("variables\n")

    # Average values to be written for staggered variables
    c = 0
    for vel in (un, vn, wn):
        if c == 0:
            file_id.write("u 1\n")
        if c == 1:
            file_id.write("v 1\n")
        if c == 2:
            file_id.write("w 1\n")
        for k in range(0, nz+1):
            for j in range(0, ny+1):
                for i in range(0, nx+1):
                    file_id.write("%12.5e\n" % vel[i,j,k])
        c += 1            

    file_id.write("endvars\n")

    # --------------------------
    # Write the file footer out
    # --------------------------
    file_id.write("endgmv\n")

    file_id.close()

    return  # end of function
