"""
Defines class of the type "Unknown".  "Unknown" is the object which holds the
values inside the computational domain and on the boundaries.  Value inside
the domain is stored in a full three-dimensional array.  Boundary values are
also formally stored in three-dimensional arrays, but depending on their
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

Unknowns in PyNS are shown below for a two-dimensional grid (for simplicity):

  Here, collocated resolution is:                  Legend:

    nx = 6                                            o   ... scalars
    ny = 4                                           ---  ... u-velocities
                             [N]                      |   ... v-velocities

      +-------+-------+-------+-------+-------+-------+
      |       |       |       |       |       |       |
      |   o  ---  o  ---  o  ---  o  ---  o  ---  o   | j=ny-1
      |       |       |       |       |       |       |
      +---|---+---|---+---|---+---|---+---|---+---|---+     j=ny-2
      |       |       |       |       |       |       |
      |   o  ---  o  ---  o  ---  o  ---  o  ---  o   | ...
      |       |       |       |       |       |       |
 [W]  +---|---+---|---+---|---+---|---+---|---+---|---+     j=1     [E]
      |       |       |       |       |       |       |
      |   o  ---  o  ---  o  ---  o  ---  o  ---  o   | j=1
      |       |       |       |       |       |       |
      +---|---+---|---+---|---+---|---+---|---+---|---+     j=0  (v-velocity)
      |       |       |       |       |       |       |
      |   o  ---  o  ---  o  ---  o  ---  o  ---  o   | j=0      (scalar)
      |       |       |       |       |       |       |
      +-------+-------+-------+-------+-------+-------+
      .  i=0     i=1     ...     ...   i=nx-2  i=nx-1 .    (scalar)
      :      i=0      i=1    ...    i=nx-3  i=nx-2    :    (u-velocity)
      |                      [S]                      |
  y   |                                               |
 ^    |<--------------------------------------------->|
 |                      domain lenght
 +---> x

"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *
from pyns.display   import write

class Unknown:

    # =========================================================================
    def __init__(self, name, pos, res, default_bc, 
                 per = (False, False, False), 
                 verbose = False):
    # -------------------------------------------------------------------------
        """
        Args:
          name: ..... String holding the name of the variable.  
                      It is intended to be used for post-processing.
          pos: ...... Integer specifying if the variable is cell centered
                      (value C), staggered in "x" (value X), "y" (value Y)
                      or in "z" direction (value Z).
          res: ...... Vector specifying resolutions in "x", "y" and "z" 
                      directions.
          default_bc: Integer specifying if the default boundary condition is
                      of Dirichlet, Neumann or Outlet type.
          per: ...... Tuple holding three Boolean type variables which specify
                      if the unknown is periodic in "x", "y" or "z" direction.

        Returns:
          Oh well, its own self, isn't it?
        """

        if verbose is True:
            write.at(__name__)

        # Store name, position and resolution
        self.name = name
        self.pos  = pos
        self.per  = per
        self.nx, self.ny, self.nz = res

        # Allocate memory space for new and old values
        self.val = zeros((self.nx, self.ny, self.nz))
        self.old = zeros((self.nx, self.ny, self.nz))

        # Create boundary tuple
        nx, ny, nz = res
        key = namedtuple("key", "typ val")
        self.bnd = (key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),
                    key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),
                    key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),
                    key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),
                    key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))),
                    key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))))

        # Prescribe default boundary conditions
        if self.per[X] == False:
            self.bnd[W].typ[0,:,:] = default_bc
            self.bnd[E].typ[0,:,:] = default_bc
        else:
            if verbose is True:
                print("  ", self.name, "is periodic in x direction;", end="")
                print("  skipping default boundary condition for it!")
        if self.per[Y] == False:
            self.bnd[S].typ[:,0,:] = default_bc
            self.bnd[N].typ[:,0,:] = default_bc
        else:
            if verbose is True:
                print("  ", self.name, "is periodic in y direction;", end="")
                print("  skipping default boundary condition for it!")
        if self.per[Z] == False:
            self.bnd[B].typ[:,:,0] = default_bc
            self.bnd[T].typ[:,:,0] = default_bc
        else:
            if verbose is True:
                print("  ", self.name, "is periodic in z direction;", end="")
                print("  skipping default boundary condition for it!")

        self.bnd[W].val[0,:,:] = 0
        self.bnd[E].val[0,:,:] = 0
        self.bnd[S].val[:,0,:] = 0
        self.bnd[N].val[:,0,:] = 0
        self.bnd[B].val[:,:,0] = 0
        self.bnd[T].val[:,:,0] = 0

        if verbose is True:
            print("  Created unknown:", self.name)

        return  # end of function

    # =========================================================================
    def exchange(self):
    # -------------------------------------------------------------------------
        """
        Function to refresh buffers.  For the time being it only takes
        care of periodic boundary conditions, but in the future it may
        also refresh buffers used in parallel execution.

        Periodicity for scalar cells:

          - value in cell "0" is identical to value in "nx-1"

          .---<-------<-------<-------<-------<---.
          |                    send to buffers    |
          |               +--->------->------->---)--->------->---.
          |               |                       |               |
          v --+-------+---|---+-------+-------+---|---+-------+-- v
              |       |   ^   |       |       |   ^   |       |
          o   |   o   |   o   |   o   |   o   |   o   |   o   |   o
              |       |       |       |       |       |       |
          - --+-------+-------+-------+-------+-------+-------+-- -
         [W]     i=0     i=1     ...     ...   i=nx-2  i=nx-1    [E]
          =       =                                       =       =
         nx-2    nx-1                                     0       1
                  |        effective domain lenght        |
                  |<------------------------------------->|

        Periodicity for vector cells:

          - value in vector "0" is east from value in "nx-2"
          - value in vector "

                      .--->------->------->------->------->---.
                      |        send to buffers                |
              .---<---)---<-------<-------<-------<---.       |
              |       |                               |       |
              |-------|-------+-------+-------+-------|-------|
              v       ^       |       |       |       ^       v
             ---     ---     ---     ---     ---     ---     ---
              |       |       |       |       |       |       |
              +-------+-------+-------+-------+-------+-------+
             [W]     i=0     i=1     ...     ...   i=nx-2    [E]
              =                                               =
            nx-2                                              0
                  |        effective domain lenght        |
                  |<------------------------------------->|

        """

        if self.per[X] == True:
            if self.pos == X:
                self.bnd[W].val[:] = self.val[-1:,:,:]
                self.bnd[E].val[:] = self.val[ :1,:,:]
            else:
                self.bnd[W].val[:] = self.val[-2:-1,:,:]
                self.bnd[E].val[:] = self.val[ 1: 2,:,:]
        if self.per[Y] == True:
            if self.pos == Y:
                self.bnd[S].val[:] = self.val[:,-1:,:]
                self.bnd[N].val[:] = self.val[:, :1,:]
            else:
                self.bnd[S].val[:] = self.val[:,-2:-1,:]
                self.bnd[N].val[:] = self.val[:, 1: 2,:]
        if self.per[Z] == True:
            if self.pos == Z:
                self.bnd[B].val[:] = self.val[:,:,-1:]
                self.bnd[T].val[:] = self.val[:,:, :1]
            else:
                self.bnd[B].val[:] = self.val[:,:,-2:-1]
                self.bnd[T].val[:] = self.val[:,:, 1: 2]

        return  # end of function
