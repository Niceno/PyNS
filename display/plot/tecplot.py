"""
Exports results in Tecplot (TM) ASCII format.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def tecplot(file_name, xyzn, unknowns = (), arrays = ()):
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


    # ------------------------------------
    # Unpack tuple with nodal coordinates
    # ------------------------------------
    xn, yn, zn = xyzn

    # Compute cell resolutions (remember: cell resolutions)
    nx = len(xn)-1
    ny = len(yn)-1
    nz = len(zn)-1

    # VisIt can't read Tecplot (TM) files with
    # comments, so keep verbatim "False"
    verbatim = True
    
    # Number of columns in the output file
    col = 10

    # ------------------
    #
    # Analyze the input
    #
    # ------------------
    
    # Browse through input arguments to find their names, positions and values
    key = namedtuple('key', 'name pos val bnd')
    vars = ()

    # First through unknowns (they have their names and positions)
    for u in unknowns:
        unk = key(u.name, u.pos, u.val, u.bnd)
        vars = vars + (unk,)

    # Then through arrays (have no name, and position is implocitly defined)
    c = 0
    for a in arrays:
      
        # Find the position of the array
        pos = C               # first assume central, then prove otherwise ;-)
        ax, ay, az = a.shape
        if ax==nx-1 and ay==ny and az==nz:
            pos = X
        if ax==nx and ay==ny-1 and az==nz:
            pos = X
        if ax==nx and ay==ny and az==nz-1:
            pos = Z
        if ax==nx+1 and ay==ny+1 and az==nz+1:
            pos = N           # here I use N, usually "north" for nodal  
        
        arr = key("pyns-var-%2.2d" % c, pos, a, None)        
        vars = vars + (arr,)
        c += 1

    # Count nodal and cell-centered variables
    c_nod = 0
    c_cel = 0
    for v in range(0, len(vars)):
        if vars[v].pos == N:
            c_nod += 1
        else:
            c_cel += 1
            
    # --------------
    #
    # Open the file
    #
    # --------------
    file_id = open(file_name + ".plt", "w")
            
    # --------------------------
    #
    # Write the file header out
    #
    # --------------------------
    if verbatim:
        file_id.write("# File header \n")
    file_id.write("title=\"PyNS Output\"\n")
    file_id.write("variables=\"x\" \"y\" \"z\" ")
    # First names of the nodal variables
    for v in range(0, len(vars)):
        if vars[v].pos == N:
            file_id.write("\"%s\" " % vars[v].name)
    # Then names of the nodal variables
    for v in range(0, len(vars)):
        if vars[v].pos != N:
            file_id.write("\"%s\" " % vars[v].name)
    file_id.write("\n")
    file_id.write("zone i=%d"   % (nx+1) +   \
                      " j=%d"   % (ny+1) +   \
                      " k=%d\n" % (nz+1))
    file_id.write("datapacking = block\n")
    # In the lines which follow, ghost number "3" is for coordinates
    file_id.write("varlocation=([1-%d]=nodal " % (3 + c_nod))
    file_id.write("[%d-%d]=cellcentered)\n" % (4 + c_nod, 3 + c_nod+c_cel) )

    # -------------------------------------------------------------------
    #
    # Write the coordinates out (remember - those are nodal coordinates)
    #
    # -------------------------------------------------------------------
    if verbatim:
        file_id.write("# X coordinates\n")
    c = 0                                # column counter
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e " % xn[i])
                c = c + 1
                if c % col == 0:         # go to new line after col columns
                    file_id.write("\n")
    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    if verbatim:
        file_id.write("# Y coordinates\n")
    c = 0                                # column counter
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e " % yn[j])
                c = c + 1
                if c % col == 0:         # go to new line after col columns
                    file_id.write("\n")
    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    if verbatim:
        file_id.write("# Z coordinates\n")
    c = 0                                # column counter
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e " % zn[k])
                c = c + 1
                if c % col == 0:         # go to new line after col columns
                    file_id.write("\n")
    if c % col != 0:                     # finish the line if necessary
        file_id.write("\n")

    # ------------------------
    #
    # Write the variables out
    #
    # ------------------------

    # ----------------------
    # First nodal variables
    # ----------------------
    for v in range(0, len(vars)):
        if vars[v].pos == N:
            c = 0
            for k in range(0, nz+1):
                for j in range(0, ny+1):
                    for i in range(0, nx+1):
                        file_id.write("%12.5e " % vars[v].val[i,j,k])
                        c = c + 1
                        if c % col == 0:
                            file_id.write("\n")
            if c % col != 0:             # finish the line if necessary
                file_id.write("\n")

    # ---------------------------------
    # Then the cell-centered variables
    # ---------------------------------
    
    # Average values to be written for staggered variables
    for v in range(0, len(vars)):
        if vars[v].pos != N:
            if vars[v].pos == C:
                val = vars[v].val
            elif vars[v].pos == X:
                val = avg_x(cat_x((vars[v].bnd[W].val[:1,:,:],
                                   vars[v].val,
                                   vars[v].bnd[E].val[:1,:,:])))
            elif vars[v].pos == Y:
                val = avg_y(cat_y((vars[v].bnd[S].val[:,:1,:],
                                   vars[v].val,
                                   vars[v].bnd[N].val[:,:1,:])))
            elif vars[v].pos == Z:
                val = avg_z(cat_z((vars[v].bnd[B].val[:,:,:1],
                                   vars[v].val,
                                   vars[v].bnd[T].val[:,:,:1])))

            if verbatim:
                file_id.write("# %s \n" % vars[v].name)
            c = 0
            for k in range(0, nz):
                for j in range(0, ny):
                    for i in range(0, nx):
                        file_id.write("%12.5e " % val[i,j,k])
                        c = c + 1
                        if c % col == 0:
                            file_id.write("\n")
            if c % col != 0:             # finish the line if necessary
                file_id.write("\n")

    file_id.close()

    return  # end of function
