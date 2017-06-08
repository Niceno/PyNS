"""
Creates a new unkown.  An unknown has a form of a structure, and holds the 
values inside the computational domain and on the boundaries.  Value inside
the domain is stored in a full three-dimensional matrix.  Boundary values are
also formally stored in three-dimensional matrices, but depending on their 
position (W, E, S, N, B or T), one dimension is set to one. For boundaries,
also the type of boundary conditions are specified. 

To access an unknown, say its name is "phi" later in the program, the following
syntax is used:
    
Value inside the domain, at the cell with coordinates (12, 34, 56):

  phi.val[12, 34, 56] = 0.0

Value on the east boundary, at coordinates (34, 56).    

  phi.bnd[E].val[1, 34, 56] = 0.0 

Type of boundary condition at the same boundary cell from above:
                    
  phi.bnd[E].typ[1, 34, 56] == NEUMANN    
"""

# Standard Python modules
from standard import *

# PyNS modules
from constants.all      import *
from operators.all      import *

# =============================================================================
def create_unknown(name, pos, res, def_bc):
# -----------------------------------------------------------------------------
  """
  Args: 
    name:   String holding the name of the variable.  It is intended to be 
            used for post-processing.
    pos:    Integer specifying if the variable is cell centered (value C), 
            staggered in "x" (value X), "y" (value Y) or in "z" direction (Z).
    res:    Vector specifying resolutions in "x", "y" and "z" directions.
    def_bc: Integer specifying if the default boundary condition is of 
            Dirichlet, Neumann or Outlet type.

   Returns:
     phi: Formed unknown.  
   """ 
# -----------------------------------------------------------------------------

  # Fetch resolutions 
  nx, ny, nz = res
  
  # Create boundary tuple
  key = namedtuple('key', 'typ val')  
  bnd = (key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),  \
         key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),  \
         key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),  \
         key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),  \
         key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))),  \
         key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))))
  
  # Create the unknown tuple
  unk = namedtuple('unk', 'name pos val old bnd')
  phi = unk(name, pos, zeros(res), zeros(res), bnd)
  
  # -----------------
  # Set the position
  # -----------------
  if pos != C and pos != X and pos != Y and pos != Z: 
    print('Variable must be defined at positions C, X, Y or Z')

  # --------------------------------------------------
  # Set default boundary conditions, types and values
  # --------------------------------------------------
  if def_bc != DIRICHLET and def_bc != NEUMANN and def_bc != OUTLET: 
    print('Variable must be defined with boundary conditions ',  \
          'DIRICHLET, NEUMANN or OUTLET!')

  phi.bnd[W].typ[0,:,:] = def_bc
  phi.bnd[E].typ[0,:,:] = def_bc
  phi.bnd[S].typ[:,0,:] = def_bc
  phi.bnd[N].typ[:,0,:] = def_bc
  phi.bnd[B].typ[:,:,0] = def_bc
  phi.bnd[T].typ[:,:,0] = def_bc

  phi.bnd[W].val[0,:,:] = 0
  phi.bnd[E].val[0,:,:] = 0 
  phi.bnd[S].val[:,0,:] = 0
  phi.bnd[N].val[:,0,:] = 0
  phi.bnd[B].val[:,:,0] = 0
  phi.bnd[T].val[:,:,0] = 0

  print("Created variable ", name)

  return phi  # end of function
