"""
Exports results in Tecplot (TM) ASCII format.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def tecplot(file_name, xyzn, vars):
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

    # Compute cell resolutions (remember: cell resolutions)
    nx = len(xn)-1
    ny = len(yn)-1
    nz = len(zn)-1

    file_id = open(file_name, 'w')

    # VisIt can't read Tecplot (TM) files which contain comments 
    verbatim = False

    # --------------------------
    # Write the file header out
    # --------------------------
    if verbatim:
        file_id.write("# File header \n")
    file_id.write("title=\"PyNS Output\"\n")
    file_id.write("variables=\"x\" \"y\" \"z\" ")
    for v in variables:
        file_id.write("\"%s\" " % v.name)
    file_id.write("\n")
    file_id.write("zone i=%d"   % (nx+1) +   \
                      " j=%d"   % (ny+1) +   \
                      " k=%d\n" % (nz+1))
    file_id.write("datapacking = block\n")
    file_id.write("varlocation=([1-3]=nodal ")
    file_id.write("[4-%d]=cellcentered)" % (3+len(variables)) )

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
    for v in variables:
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

        if verbatim:
            file_id.write("\n# %s \n" % v.name)
        c = 0
        for k in range(0, nz):
            for j in range(0, ny):
                for i in range(0, nx):
                    file_id.write("%12.5e " % val[i,j,k])
                    c = c + 1
                    if c % 4 == 0:
                        file_id.write("\n")

    file_id.close()

    return  # end of function

# =============================================================================
def tecplot(file_name, xyzn, var, pos):
# -----------------------------------------------------------------------------
    """
    Args:
      file_name: String containing name of the file to be created.
      xyzn:      Tuple containing one-dimensional arrays with "x", "y"
                 and "z" coordinates.
      var:       Variable to be exported to Tecplot (TM)
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

    file_id = open(file_name, 'w')

    # VisIt can't read Tecplot (TM) files which contain comments 
    verbatim = False

    # --------------------------
    # Write the file header out
    # --------------------------
    if verbatim:
        file_id.write("# File header \n")
    file_id.write("title=\"PyNS Output\"\n")
    file_id.write("variables=\"x\" \"y\" \"z\" \"unnamed\"")
    file_id.write("\n")
    file_id.write("zone i=%d"   % (nx+1) +   \
                      " j=%d"   % (ny+1) +   \
                      " k=%d\n" % (nz+1))
    file_id.write("datapacking = block\n")
    file_id.write("varlocation=([1-3]=nodal [4]=cellcentered)")

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

    # -----------------------
    # Write the variable out
    # -----------------------

    # Average values to be written for staggered variables
    if pos == C:
        val = var
    elif v.pos == X:
        val = avg_x(cat_x((var[ :1,:,:], var, var[-1:,:,:])))
    elif v.pos == Y:
        val = avg_y(cat_y((var[:, :1,:], var, var[:,-1:,:])))
    elif v.pos == Z:
        val = avg_z(cat_z((var[:,:, :1], var, var[:,:,-1:])))

    if verbatim:
        file_id.write("\n# %s \n" % v.name)
    c = 0
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                file_id.write("%12.5e " % var[i,j,k])
                c = c + 1
                if c % 4 == 0:
                    file_id.write("\n")

    file_id.close()

    return  # end of function
