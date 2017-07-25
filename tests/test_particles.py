#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lagrangian Particle tracking test case.

Calculates the trajectories of particles for the impactor case. 
The fluid velocties are loaded in from the export Tecplot Results. 
 
"""


# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.discretization import *
from pyns.display        import plot, write 
from pyns.physical       import properties
from pyns.lagrangian     import *
import numpy as np 


def main(show_plot=True, time_steps=430000, plot_freq=1):
    
    # Loading in the Fluid Flow    
    f = np.loadtxt('fluent_fluid.txt')
    
    # Picking out the velocity components (1D arrays)
    u = f[:,3]
    v = f[:,4]
    w = f[:,5]
    
    # Node coordinates
    xn = nodes(0, 0.8,   640)
    yn = nodes(0, 0.1,   80)
    zn = nodes(0, 0.025,   4)
    
    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)
    
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)
    
    obst = zeros(rc)
    
    for k in range(0,nz):
        for i in range(0, 15*nx//32):
            for j in range(ny//4,ny):
                obst[i,j,k] = 1
                
    for k in range(0,nz):
        for i in range(17*nx//32, nx):
            for j in range(ny//4,ny):
                obst[i,j,k] = 1

    # Set physical properties
    rho, mu, cap, kappa = properties.air(rc)
    
    # Setting up the velocity arrays. 
    un = empty((nx+1,ny+1,nz+1))
    vn = empty((nx+1,ny+1,nz+1))
    wn = empty((nx+1,ny+1,nz+1))
    
    # Assigning values to the 3D velocity arrays (there 
    # might be a better way to do this) 

    d = -1 
    for k in range(0,nz+1):
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                d = d + 1
                un[i,j,k] = u[d]
            
    d = -1 
    for k in range(0,nz+1):
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                d = d + 1
                vn[i,j,k] = v[d]
        
    d = -1 
    for k in range(0,nz+1):
        for j in range(0,ny+1):
            for i in range(0,nx+1):
                d = d + 1
                wn[i,j,k] = w[d]
    
    # Time-stepping parameters
    dt  = 3e-06        # time step <= tau/2, where tau is the relaxtion time
    ndt = time_steps     # number of time steps

    # Initialising the particles. 
    n = 1000  # number of particles 

    pt = initialiser(n, rho_p=1000, d = 1.5e-6, verbose = False)
    
    # Output initial conditions to Tecplot
    plot.tecplot("particle-%6.6d" % 0,(xn,yn,zn), arrays = (un,vn,wn),
                 tracers = (pt))

    # ----------
    #
    # Time loop
    #
    # ----------
    for ts in range(1,ndt+1):

        # ----------------------------
        # Lagrangian Particle Tracking
        # ----------------------------
    
        # Iterate through the number of particles to update their positions
        # and velocities. 

        calc_traj(pt, (un,vn,wn), rho, mu, (xn,yn,zn), (xc,yc,zc), dt, obst, n)
        
        # =====================================================================
        #
        # Visualisation
        #
        # =====================================================================
        if show_plot:
            if ts % plot_freq == 0:
                write.time_step(ts)
                
                plot.tecplot("particle-%6.6d" % ts,
                             (xn,yn,zn),
                             arrays = (un,vn,wn),
                             tracers = (pt))

  
if __name__ == '__main__':
    main()      
  
