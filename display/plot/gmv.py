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
def gmv(file_name, xyzn, unknowns = (), arrays = ()):
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

    # ------------------------------------
    # Unpack tuple with nodal coordinates
    # ------------------------------------
    xn, yn, zn = xyzn

    # Compute cell resolutions (remember: cell resolutions)
    nx = len(xn)-1
    ny = len(yn)-1
    nz = len(zn)-1

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
        pos = C               # first assume central, then proove otherwise ;-)
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

    # --------------
    #
    # Open the file
    #
    # --------------
    file_id = open(file_name + ".gmv", "w")

    # --------------------------
    #
    # Write the file header out
    #
    # --------------------------
    file_id.write("gmvinput ascii\n")

    # -------------------------------------------------------------------
    #
    # Write the coordinates out (remember - those are nodal coordinates)
    #
    # -------------------------------------------------------------------
    file_id.write("nodev %d\n" % ((nx+1)*(ny+1)*(nz+1)))
    for k in range(0, nz+1):
        for j in range(0, ny+1):
            for i in range(0, nx+1):
                file_id.write("%12.5e %12.5e %12.5e\n" % (xn[i], yn[j], zn[k]))

    # --------------------
    #
    # Write the cells out
    #
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
    #
    # Write the variables out
    #
    # ------------------------
    file_id.write("variables\n")

    # Average values to be written for staggered variables
    for v in range(0, len(vars)):
        if vars[v].pos != N:
            if vars[v].pos == C:
                val = vars[v].val
            elif vars[v].pos == X:
                val = avg_x(cat_x((vars[v].bnd[W].val[:1,:,:],   \
                                   vars[v].val,                  \
                                   vars[v].bnd[E].val[:1,:,:])))
            elif vars[v].pos == Y:
                val = avg_y(cat_y((vars[v].bnd[S].val[:,:1,:],   \
                                   vars[v].val,                  \
                                   vars[v].bnd[N].val[:,:1,:])))
            elif vars[v].pos == Z:
                val = avg_z(cat_z((vars[v].bnd[B].val[:,:,:1],   \
                                   vars[v].val,                  \
                                   vars[v].bnd[T].val[:,:,:1])))
        else:        
            val = vars[v].val
                      
        # Write out nodal variables
        if vars[v].pos == N:
            file_id.write("%s 1\n" % vars[v].name)
            c = 0
            for k in range(0, nz+1):
                for j in range(0, ny+1):
                    for i in range(0, nx+1):
                        file_id.write("%12.5e\n" % val[i,j,k])
                        
        # Write out cell variables                
        else:
            file_id.write("%s 0\n" % vars[v].name)
            c = 0
            for k in range(0, nz):
                for j in range(0, ny):
                    for i in range(0, nx):
                        file_id.write("%12.5e\n" % val[i,j,k])

    file_id.write("endvars\n")

    # --------------------------
    #
    # Write the file footer out
    #
    # --------------------------
    file_id.write("endgmv\n")

    file_id.close()

    return  # end of function
