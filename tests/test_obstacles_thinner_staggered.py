"""
#                                                       o ... scalars
#                          (n)                          - ... u velocities
#                                                       | ... v velocities
#       +-------+-------+-------+-------+-------+
#       |       |       |       |       |       |
#       |   o   -   o   -   o   -   o   -   o   | j=ny-2
#       |       |       |       |       |       |
#       +---|---+---|---+---|---+---|---+---|---+     j=ny-1
#       |       |       |       |       |       |
#       |   o   -   o   -   o   -   o   -   o   | ...
#       |       |       |       |       |       |
#  (w)  +---|---+---|---+---|---+---|---+---|---+    j=1        (e)
#       |       |       |       |       |       |
#       |   o   -   o   -   o   -   o   -   o   | j=1
#       |       |       |       |       |       |
#       +---|---+---|---+---|---+---|---+---|---+    j=0 (v-velocity)
#       |       |       |       |       |       |
#       |   o   -   o   -   o   -   o   -   o   | j=0   (scalar cell)
#       |       |       |       |       |       |
#       +-------+-------+-------+-------+-------+
#  y       i=0     i=1     ...     ...    i=nx-1      (scalar cells)
# ^            i=0      i=1    ...    i=nx-2      (u-velocity cells)
# |
# +---> x                  (s)
#
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
    xn = nodes(0, 1.25,  256)
    yn = nodes(0, 0.125,  32)
    zn = nodes(0, 0.125,  32)

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)

    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    rho, mu, cap, kappa = properties.air(rc)

    # Time-stepping parameters
    dt  = 0.005      # time step
    ndt = time_steps # number of time steps

    # Create unknowns; names, positions and sizes
    uf = Unknown("face-u-vel",  X, ru, DIRICHLET)
    vf = Unknown("face-v-vel",  Y, rv, DIRICHLET)
    wf = Unknown("face-w-vel",  Z, rw, DIRICHLET)
    p  = Unknown("pressure",    C, rc, NEUMANN)

    # Specify boundary conditions
    uf.bnd[W].typ[:1,:,:] = DIRICHLET
    uf.bnd[W].val[:1,:,:]  = 0.1 * outer( par(1.0, yn), par(1.0, zn) )

    uf.bnd[E].typ[:1,:,:] = OUTLET

    # Create obstacles
    plates = zeros(rc)

    class key:
        ip = -1
        im = -1
        jp = -1
        jm = -1
        kp = -1
        km = -1

    block = (key(), key(), key(), key());

    th = 5;
    block[0].im =  3*nx/16            # i minus
    block[0].ip =  block[0].im + th   # i plus
    block[0].jm =  0                  # j minus
    block[0].jp =  3*ny/4             # j plus
    block[0].km =  0                  # k minus
    block[0].kp =  3*ny/4             # k plus

    block[1].im =  5*nx/16            # i minus
    block[1].ip =  block[1].im + th   # i plus
    block[1].jm =  ny/4               # j minus
    block[1].jp =  ny                 # j plus
    block[1].km =  ny/4               # k minus
    block[1].kp =  ny                 # k plus

    block[2].im =  7*nx/16            # i minus
    block[2].ip =  block[2].im + th   # i plus
    block[2].jm =  0                  # j minus
    block[2].jp =  3*ny/4             # j plus
    block[2].km =  0                  # k minus
    block[2].kp =  3*ny/4             # k plus

    block[3].im =  9*nx/16            # i minus
    block[3].ip =  block[3].im + th   # i plus
    block[3].jm =  ny/4               # j minus
    block[3].jp =  ny                 # j plus
    block[3].km =  ny/4               # k minus
    block[3].kp =  ny                 # k plus

    for o in range(0, 4):
        for i in range(floor(block[o].im), floor(block[o].ip)):
            for j in range(floor(block[o].jm), floor(block[o].jp)):
                for k in range(floor(block[o].km), floor(block[o].kp)):
                    plates[i,j,k] = 1

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
        calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu, dt, (dx,dy,dz), 
                 obstacle = plates)

        # ---------
        # Pressure
        # ---------
        calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), 
               obstacle = plates)

        # --------------------
        # Velocity correction
        # --------------------
        corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), 
                 obstacle = plates)

        # Check the CFL number too
        cfl = cfl_max((uf,vf,wf), dt, (dx,dy,dz))

# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:
            if ts % plot_freq == 0:
                plot.isolines(p.val, (uf,vf,wf), (xn,yn,zn), Y)
                plot.isolines(p.val, (uf,vf,wf), (xn,yn,zn), Z)
                plot.gmv("obst-thinner-staggered-%6.6d" % ts, 
                         (xn,yn,zn), (uf,vf,wf,p))

if __name__ == "__main__":
    main()
