"""
Program to test implementation of advection schemes in the code, using
one-dimensoinal transport of a step function in X, Y or Z
direction, either in positive or negative sense.

The coordinate direction is specified with the local variable "TEST",
which can be either X, Y or Z.

Sense is specified with the variable "FLOW", which can assume values 'p'
for positive, and 'n' for negative sense.
"""

#!/usr/bin/python

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.discretization import *
from pyns.display        import plot, write
from pyns.physical       import properties

from matplotlib  import pyplot as plt

def main(show_plot=True, time_steps=180, plot_freq=180):

# =============================================================================
#
# Define problem
#
# =============================================================================

    sens = -1, +1
    test = X, Y, Z
    SENS = sens[ floor( random() * len(sens) ) ]
    TEST = test[ floor( random() * len(test) ) ]
    MIN = 10
    MAX = 20
    
    if TEST == X:
        uvw_bulk = SENS, 0.0, 0.0
        xn  = nodes(0, 1,   200)
        yn  = nodes(0, 0.1,   4)
        zn  = nodes(0, 0.1,   4)
        per = (True, False, False)
    if TEST == Y:
        uvw_bulk = 0.0, SENS, 0.0
        xn  = nodes(0, 0.1,   4)
        yn  = nodes(0, 1,   200)
        zn  = nodes(0, 0.1,   4)
        per = (False, True, False)
    elif TEST == Z:
        uvw_bulk = 0.0, 0.0, SENS
        xn  = nodes(0, 0.1,   4)
        yn  = nodes(0, 0.1,   4)
        zn  = nodes(0, 1,   200)
        per = (False, False, True)
    
    # Set domain dimensions and grid resolution
    
    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)
    
    # Set physical properties
    rho   = zeros(rc)
    mu    = zeros(rc)
    kappa = zeros(rc)
    cap   = zeros(rc)
    rho  [:,:,:] = 1.0
    mu   [:,:,:] = 0.0
    kappa[:,:,:] = 0.0
    cap  [:,:,:] = 1.0
    
    # Time-stepping parameters
    dt  = 0.0025      # time step
    ndt = time_steps  # number of time steps
    
    # Create unknowns; names, positions and sizes
    uf = Unknown('face-u-vel',  X, ru, DIRICHLET)
    vf = Unknown('face-v-vel',  Y, rv, DIRICHLET)
    wf = Unknown('face-w-vel',  Z, rw, DIRICHLET)
    t  = Unknown('temperature', Y, rv, NEUMANN, per)
    
    # Make a tuple with all velocity components
    uvwf = uf, vf, wf
    
    # Specify inital and boundary conditions
    t.val[:] = MIN
    if TEST == X:
        t.val[nx//2-25:nx//2+25,:,:] = MAX
    elif TEST == Y:
        t.val[:,ny//2-25:ny//2+25,:] = MAX
    elif TEST == Z:
        t.val[:,:,nz//2-25:nz//2+25] = MAX
    
    adj_n_bnds(t)
    
    # This is from the advection-1D case
    for i in (X,Y,Z):
        uvwf[i].val[:,:,:] = uvw_bulk[i]
        for j in (W,E,S,N,B,T):
            uvwf[i].bnd[j].val[:,:,:] = uvw_bulk[i]
    
    obst = zeros(rc)

# =============================================================================
#
# Solution algorithm
#
# =============================================================================

    # ----------
    #
    # Time loop
    #
    # ----------
    for ts in range(1,ndt+1):
    
        write.time_step(ts)
    
        # -----------------
        # Store old values
        # -----------------
        t.old[:] = t.val[:]
    
        # -----------------------
        # Temperature (enthalpy)
        # -----------------------
        calc_t(t, uvwf, (rho*cap), kappa, dt, (dx,dy,dz), obst)

# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:            
            if ts % plot_freq == 0:
                st = 'ro'
                xc = avg(xn)
                vc = t.val[:,ny//2,nz//2]
                if TEST == Y:
                    st = 'go'
                    xc = avg(yn)
                    vc = t.val[nx//2,:,nz//2]
                elif TEST == Z:
                    st = 'bo'
                    xc = avg(zn)
                    vc = t.val[nx//2,ny//2,:]
                
                if t.pos == TEST:
                    xc = avg(xc)
                
                plt.plot(xc, vc, st)
                plt.show()

if __name__ == '__main__':
    main()
