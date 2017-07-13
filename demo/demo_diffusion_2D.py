"""
Demonstrates periodicity for diffusion.
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

def main(show_plot=True, time_steps=1, plot_freq=1):

# =============================================================================
#
# Define problem
#
# =============================================================================

    NX = 32
    NY = NX
    NZ = 32

    xn = nodes(0, 1, NX)
    yn = nodes(0, 1, NY)
    zn = nodes(0, 1, NZ)

    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    grashof = 1.4105E+06
    prandtl = 0.7058
    rho   = zeros(rc)
    mu    = zeros(rc)
    kappa = zeros(rc)
    cap   = zeros(rc)
    rho  [:,:,:] = 1.0
    mu   [:,:,:] = 1.0 / sqrt(grashof)
    kappa[:,:,:] = 1.0 / (prandtl * sqrt(grashof))
    cap  [:,:,:] = 1.0

    # Time-stepping parameters
    dt  = 2.0e10        # time step
    ndt = time_steps  # number of time steps

    # Create unknowns; names, positions and sizes
    uf = Unknown('face-u-vel',  X, ru, DIRICHLET)
    vf = Unknown('face-v-vel',  Y, rv, DIRICHLET)
    wf = Unknown('face-w-vel',  Z, rw, DIRICHLET)
    t  = Unknown('temperature', C, rc, NEUMANN, per=(False, True, True))

    # This is a new test
    t.bnd[W].typ[:] = DIRICHLET
    t.bnd[W].val[:] = -0.0

    t.bnd[E].typ[:] = DIRICHLET
    t.bnd[E].val[:] = +0.0

    t_src = zeros(t.val.shape)
    t_src[  NX//8 :   NY//2,   NX//8 :   NY//2, : ] =  1.0
    t_src[5*NX//8 : 3*NY//4, 5*NX//8 : 3*NY//4, : ] = -4.0

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
        t.old[:]  = t.val[:]
        uf.old[:] = uf.val[:]
        vf.old[:] = vf.val[:]
        wf.old[:] = wf.val[:]

        # -----------------------
        # Temperature (enthalpy)
        # -----------------------
        calc_phi(t, (uf,vf,wf), (rho*cap), kappa, dt, (dx,dy,dz), None,
                 source = t_src)

# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:
            if ts % plot_freq == 0:
                plot.isolines(t.val, (uf, vf, wf), (xn, yn, zn), Z, levels=21)
                plot.gmv("tdc-staggered-%6.6d.gmv" % ts,
                         (xn, yn, zn), (uf, vf, wf, t))

if __name__ == '__main__':
    main()
