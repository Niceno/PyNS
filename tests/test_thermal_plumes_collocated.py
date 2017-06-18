"""
This script is to reproduce two-dimensional mixed convection case, with
the aim of testing the outflow boundary, particularly the "convective"
boundary condition which allows eddies to leave the domain.

There is also a staggered version of this script, called
"demo_plums_collocated.m".  It would be good to keep both version as
similar as possible to each other, to test the differences between
staggered and collocated arrangements always possible.
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

def main(show_plot=True, time_steps=1800, plot_freq=180):

# =============================================================================
#
# Define problem
#
# =============================================================================

    # Node coordinates
    xn = nodes(0, 10, 300)
    yn = nodes(0,  1,  40, 1/500, 1/500)
    zn = nodes(0,  3,   3)

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)

    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    rho   = zeros(rc)
    mu    = zeros(rc)
    kappa = zeros(rc)
    cap   = zeros(rc)
    rho  [:,:,:] = 1.
    mu   [:,:,:] = 0.1
    kappa[:,:,:] = 0.15
    cap  [:,:,:] = 1.0

    # Time-stepping parameters
    dt  = 0.003       # time step
    ndt = time_steps  # number of time steps

    # Create unknowns; names, positions and sizes
    uc = Unknown('cell-u-vel',  C, rc, DIRICHLET)
    vc = Unknown('cell-v-vel',  C, rc, DIRICHLET)
    wc = Unknown('cell-w-vel',  C, rc, DIRICHLET)
    uf = Unknown('face-u-vel',  X, ru, DIRICHLET)
    vf = Unknown('face-v-vel',  Y, rv, DIRICHLET)
    wf = Unknown('face-w-vel',  Z, rw, DIRICHLET)
    t  = Unknown('temperature', C, rc, NEUMANN)
    p  = Unknown('pressure',    C, rc, NEUMANN)

    # Specify boundary conditions
    uc.bnd[W].typ[:1,:,:] = DIRICHLET
    for k in range(0,nz):
        uc.bnd[W].val[:1,:,k] = par(1.0, yn);

    uc.bnd[E].typ[:1,:,:] = OUTLET
    uc.bnd[E].val[:1,:,:] = 1.0;

    for j in (B,T):
        uc.bnd[j].typ[:] = NEUMANN
        vc.bnd[j].typ[:] = NEUMANN
        wc.bnd[j].typ[:] = NEUMANN

    t.bnd[W].typ[:1,:,:] = DIRICHLET
    for k in range(0,nz):
        t.bnd[W].val[:1,:,k] = 1.0-yc;

    t.bnd[S].typ[:,:1,:] = DIRICHLET
    t.bnd[S].val[:,:1,:] = +1.0;
    t.bnd[N].typ[:,:1,:] = DIRICHLET
    t.bnd[N].val[:,:1,:] =  0.0

    # Specify initial conditions
    uc.val[:,:,:] = 1.0
    t.val[:,:,:] = 0

    # Copy the values to face velocities
    uf.val[:] = avg(X, uc.val)
    vf.val[:] = avg(Y, vc.val)
    wf.val[:] = avg(Z, wc.val)
    for j in (W,E):
        uf.bnd[j].val[:] = uc.bnd[j].val[:]
        vf.bnd[j].val[:] = avg(Y, vc.bnd[j].val[:])
        wf.bnd[j].val[:] = avg(Z, wc.bnd[j].val[:])
    for j in (S,N):
        uf.bnd[j].val[:] = avg(X, uc.bnd[j].val[:])
        vf.bnd[j].val[:] = vc.bnd[j].val[:]
        wf.bnd[j].val[:] = avg(Z, wc.bnd[j].val[:])
    for j in (B,T):
        uf.bnd[j].val[:] = avg(X, uc.bnd[j].val[:])
        vf.bnd[j].val[:] = avg(Y, vc.bnd[j].val[:])
        wf.bnd[j].val[:] = wc.bnd[j].val[:]

    obstacle = None

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
        uc.old[:] = uc.val[:]
        vc.old[:] = vc.val[:]
        wc.old[:] = wc.val[:]

        # -----------------------
        # Temperature (enthalpy)
        # -----------------------
        calc_t(t, (uf,vf,wf), (rho*cap), kappa, dt, (dx,dy,dz), obstacle)

        # ----------------------
        # Momentum conservation
        # ----------------------
        ef = zeros(rc), 150.0 * t.val, zeros(rc)

        calc_uvw((uc,vc,wc), (uf,vf,wf), rho, mu, dt, (dx,dy,dz), obstacle,
                 force = ef)

        # ---------
        # Pressure
        # ---------
        calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obstacle)

        # --------------------
        # Velocity correction
        # --------------------
        corr_uvw((uc,vc,wc), p, rho, dt, (dx,dy,dz), obstacle)
        corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), obstacle)

        # Check the CFL number too
        cfl = cfl_max((uc,vc,wc), dt, (dx,dy,dz))

# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:
            if ts % plot_freq == 0:
                plot.isolines(t.val, (uc, vc, wc), (xn, yn, zn), Z)
                plot.isolines(p.val, (uc, vc, wc), (xn, yn, zn), Z)
                plot.tecplot("tp-collocated-%6.6d" % ts, 
                             (xn, yn, zn), (uc, vc, wc, t, p))

if __name__ == '__main__':
    main()
