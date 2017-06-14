"""
Defines class of the type "Unknown".  "Unknown" is the object which holds the
values inside the computational domain and on the boundaries.  Value inside
the domain is stored in a full three-dimensional matrix.  Boundary values are 
also formally stored in three-dimensional matrices, but depending on their 
position (W, E, S, N, B or T), one dimension is set to one. For boundaries, 
also the type of boundary conditions are specified.

To access an unknown, say its name is "phi" later in the program, the
following syntax is used:

Value inside the domain, at the cell with coordinates (12, 34, 56):

  phi.val[12, 34, 56] = 0.0

Value on the east boundary, at coordinates (34, 56).

  phi.bnd[E].val[1, 34, 56] = 0.0

Type of boundary condition at the same boundary cell from above:

  phi.bnd[E].typ[1, 34, 56] == NEUMANN
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

class Unknown:

    # =========================================================================
    def __init__(self, name, pos, res, def_bc, per=(False, False, False)):
    # -------------------------------------------------------------------------  
        """
        Args:
          name:   String holding the name of the variable.  It is intended to
                  be used for post-processing.
          pos:    Integer specifying if the variable is cell centered 
                  (value C), staggered in "x" (value X), "y" (value Y) 
                  or in "z" direction (value Z).
          res:    Vector specifying resolutions in "x", "y" and "z" directions.
          def_bc: Integer specifying if the default boundary condition is 
                  of Dirichlet, Neumann or Outlet type.
          per:    Tuple holding three Boolean type variables which specify 
                  if the unknown is periodic in "x", "y" or "z" direction.
                  
        Returns:
          Oh well, its own self, isn't it?
        """

        # Store name, position and resolution
        self.name = name
        self.pos  = pos
        self.nx, self.ny, self.nz = res

        # Allocate memory space for new and old values
        self.val = zeros((self.nx, self.ny, self.nz))
        self.old = zeros((self.nx, self.ny, self.nz))

        # Create boundary tuple
        nx, ny, nz = res
        key = namedtuple('key', 'typ val')
        self.bnd = (key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),
                    key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),
                    key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),
                    key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),
                    key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))),
                    key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))))
        
        # Prescribe default boundary conditions
        self.bnd[W].typ[0,:,:] = def_bc
        self.bnd[E].typ[0,:,:] = def_bc
        self.bnd[S].typ[:,0,:] = def_bc
        self.bnd[N].typ[:,0,:] = def_bc
        self.bnd[B].typ[:,:,0] = def_bc
        self.bnd[T].typ[:,:,0] = def_bc

        self.bnd[W].val[0,:,:] = 0
        self.bnd[E].val[0,:,:] = 0
        self.bnd[S].val[:,0,:] = 0
        self.bnd[N].val[:,0,:] = 0
        self.bnd[B].val[:,:,0] = 0
        self.bnd[T].val[:,:,0] = 0

        print("Created unknown ", self.name)
