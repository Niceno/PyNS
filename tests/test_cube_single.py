#!/usr/bin/python
"""
Single cube placed in a square channel.
"""


# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.discretization import *
from pyns.display        import plot, write
from pyns.physical       import properties

def main(show_plot=True, time_steps=2000, plot_freq=20):

# =============================================================================
#
# Define problem
#
# =============================================================================

    # Node coordinates
    xn = nodes(0, 0.9, 90)
    yn = nodes(0, 0.6, 60)
    zn = nodes(0, 0.3, 30)

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)

    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    rho, mu, cap, kappa = properties.air(rc)

    # Time-stepping parameters
    dt  = 0.005       # time step
    ndt = time_steps  # number of time steps

    # Create unknowns; names, positions and sizes
    uf = Unknown("face-u-vel",  X, ru, DIRICHLET)
    vf = Unknown("face-v-vel",  Y, rv, DIRICHLET)
    wf = Unknown("face-w-vel",  Z, rw, DIRICHLET)
    p  = Unknown("pressure",    C, rc, NEUMANN)

    # Specify boundary conditions
    uf.bnd[W].typ[:1,:,:] = DIRICHLET
    for j in range(0, ny):
        uf.bnd[W].val[:1,j,:]  = par(0.4, zn)

    uf.bnd[E].typ[:1,:,:] = OUTLET

    for j in (S, N):
        uf.bnd[j].typ[:] = NEUMANN
        vf.bnd[j].typ[:] = NEUMANN
        wf.bnd[j].typ[:] = NEUMANN

    cube = zeros(rc)
    for j in range(22, 38):
        for i in range(22, 38):
            for k in range(0,16):
                cube[i,j,k] = 1

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
        uf.old[:] = uf.val[:]
        vf.old[:] = vf.val[:]
        wf.old[:] = wf.val[:]

        # ----------------------
        # Momentum conservation
        # ----------------------
        calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu,dt, (dx,dy,dz), cube)

        # ---------
        # Pressure
        # ---------
        calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), cube)

        # --------------------
        # Velocity correction
        # --------------------
        corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), cube)

        # Check the CFL number too
        cfl = cfl_max((uf,vf,wf), dt, (dx,dy,dz))

# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:
            if ts % plot_freq == 0:
                
                # Compute nodal velocities
                un, vn, wn = nodal_uvw((xn,yn,zn), (uf,vf,wf), cube) 
                
                # Plot everything
                plot.tecplot("cube-single-%6.6d" % ts, (xn,yn,zn), 
                             unknowns = (uf, vf, wf, p),
                             arrays   = (un, vn, wn) )

if __name__ == "__main__":
    main()
