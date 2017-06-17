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
def gmv(file_name, xyzn, vars):
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
    for v in vars:
        if v.pos == C:
            val = v.val
        elif v.pos == X:
            val = avg_x(cat_x((v.bnd[W].val[:1,:,:],   \
                               v.val,                  \
                               v.bnd[E].val[:1,:,:])))
        elif v.pos == Y:
            val = avg_y(cat_y((v.bnd[S].val[:,:1,:],   \
                               v.val,                  \
                               v.bnd[N].val[:,:1,:])))
        elif v.pos == Z:
            val = avg_z(cat_z((v.bnd[B].val[:,:,:1],   \
                               v.val,                  \
                               v.bnd[T].val[:,:,:1])))

        file_id.write("%s 0\n" % v.name)
        c = 0
        for k in range(0, nz):
            for j in range(0, ny):
                for i in range(0, nx):
                    file_id.write("%12.5e\n" % val[i,j,k])

    file_id.write("endvars\n")

    # --------------------------
    # Write the file footer out
    # --------------------------
    file_id.write("endgmv\n")

    file_id.close()

    return  # end of function

# =============================================================================
def gmv_one(file_name, xyzn, var, pos):
# -----------------------------------------------------------------------------
    """
    Args:
      file_name: String containing name of the file to be created.
      xyzn:      Tuple containing one-dimensional arrays with "x", "y"
                 and "z" coordinates.
      var:       Variable to be exported to GMV (TM); collocated or staggered.
      pos:       Position of the variable.
    Returns:
      none!
    """

    # Unpack tuples
    xn, yn, zn = xyzn

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
    # file_id.write("cells %d\n" % (nx * ny * nz))
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

    if pos == C:
        val = var
    elif v.pos == X:
        val = avg_x(cat_x((var[ :1,:,:], var, var[-1:,:,:])))
    elif v.pos == Y:
        val = avg_y(cat_y((var[:, :1,:], var, var[:,-1:,:])))
    elif v.pos == Z:
        val = avg_z(cat_z((var[:,:, :1], var, var[:,:,-1:])))

    file_id.write("unnamed 0\n")
    c = 0
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                file_id.write("%12.5e\n" % val[i,j,k])

    file_id.write("endvars\n")

    # --------------------------
    # Write the file footer out
    # --------------------------
    file_id.write("endgmv\n")

    file_id.close()

    return  # end of function

    return
