#!/usr/bin/env python3
"""
3D Impactor: Simplest case to illustrate how to implement Lagrangian Particle 
tracking techniques.

The fluid flow is first solved in the domain. 
As One-way coupling is being considered, this fluid velocity is then used to 
solve for the velocity of the particles. 
"""
from __future__ import division

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.discretization import *
from pyns.display        import plot, write 
from pyns.physical       import properties
from pyns.lagrangian     import *
from mpl_toolkits.mplot3d import Axes3D





import matplotlib.pyplot as plt 

def main(show_plot=True, time_steps=500, plot_freq=10):
    
    # Node coordinates
    xn = nodes(0, 1,    160)
    yn = nodes(0, 1,    160)
    zn = nodes(0, 0.25,   4)
    

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)
    
    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    rho, mu, cap, kappa = properties.air(rc)

    # Time-stepping parameters
    dt  = 0.15        # time step
    ndt = time_steps  # number of time steps

    # Create unknowns names, positions and sizes
    uf = Unknown("face-u-vel",  X, ru, DIRICHLET)
    vf = Unknown("face-v-vel",  Y, rv, DIRICHLET)
    wf = Unknown("face-w-vel",  Z, rw, DIRICHLET)
    p  = Unknown("pressure",    C, rc, NEUMANN)
    
    # Imposing the geometry by creating an obstacle
    #
    #
    #                             INLET
    #                       |                |
    #                       |                |
    #                       |                |
    #                       |                |
    #                       |                |
    #                       |                |
    #                       |                |
    #                       |                |
    #     __________________|                |___________________
    #   O                                                          O
    #   U                                                          U
    #   T                                                          T
    #   L                                                          L
    #   E                                                          E
    #   T ________________________________________________________ T
    
    obst = zeros(rc)
    
    for k in range(0,nz):
        for i in range(0,nx//3):
            for j in range(ny//4,ny):
                obst[i,j,k] = 1
                
    for k in range(0,nz):
        for i in range(2*nx//3,nx):
            for j in range(ny//4,ny):
                obst[i,j,k] = 1
    
    for k in range(0,nz):
        vf.bnd[N].val[nx//3:2*nx//3, 0, k] = -par(0.01, xn[nx//3:2*nx//3+1])
        uf.bnd[E].typ[0, :ny//4, k] = OUTLET
        uf.bnd[W].typ[0, :ny//4, k] = OUTLET
        
    for j in (B,T):
        uf.bnd[j].typ[:] = NEUMANN
        vf.bnd[j].typ[:] = NEUMANN
        wf.bnd[j].typ[:] = NEUMANN


    # Initialising the particles. 
    n = 10000  # number of particles 
    
    pt = initialiser(n, rho_p = 1000, d = 2.5e-6, verbose = False)

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
        ef = zeros(ru), zeros(rv), zeros(rw)

        calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu, dt, (dx,dy,dz), 
                 obstacle = obst)

        # ---------
        # Pressure
        # ---------
        calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), 
               obstacle = obst)

        # --------------------
        # Velocity correction
        # --------------------
        corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), 
                 obstacle = obst)

        # Check the CFL number too
        cfl = cfl_max((uf,vf,wf), dt, (dx,dy,dz))
        
        # Calculate the nodal velocities
        
        (un, vn, wn) = nodal_uvw((xn,yn,zn), (uf,vf,wf), 
                                  obstacle = obst)

        # ----------------------------
        # Lagrangian Particle Tracking
        # ----------------------------
    
        # Iterate through the number of particles to update their positions
        # and velocities. 

        calc_traj(pt, (un,vn,wn), rho, mu, (xn,yn,zn), (xc,yc,zc), dt, obst, n)
        
# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:
            if ts % plot_freq == 0:
                
               plot.gmv("impactor-%6.6d" % ts, (xn,yn,zn), 
                        unknowns = (uf, vf, wf, p),
                        arrays   = (un, vn, wn),
                        tracers  = (pt))
 
if __name__ == "__main__":
    main()      
  
