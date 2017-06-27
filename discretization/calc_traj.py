#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To calculate the trajectory of the particles, the coupled ODE's must be
solved. A Lagrange interpolation method is used to calculate the fluid velocity 
at the particles position and then a runge-kutta (order 4) is used to solve for 
the particles velocity.
"""
# Standard Python modules
from pyns.standard import *


#From PyNS modules
from pyns.constants import *
from pyns.discretization.adj_n_bnds   import adj_n_bnds
from pyns.discretization.adj_o_bnds   import adj_o_bnds
from pyns.discretization.obst_zero_val  import obst_zero_val
from pyns.operators.avg import avg_x
from pyns.operators.avg import avg_y
from pyns.operators.avg import avg_z
from pyns.operators.cat import cat_x
from pyns.operators.cat import cat_y
from pyns.operators.cat import cat_z
from pyns.lagrangian    import * 
import numpy as np 

def calc_traj(pt, uvwn, rho, mu, xyzn, xyzc, dt, obst, n):
    
    
    """
    Args:
      pt:   The created particles velocitie components and positions.
      uvwn:  Tuple holding the nodal velocities. 
      rho/mu: The fluids denisty and viscosity. 
      xyzn:  Tuple holding the nodal positions.
      xyzc:  Tuple holding the cell centres positions.
      dt:    Time step.
      obst: The obstacle.
      n: Number of particles
    

    Returns:
      Updated location and velocities of the particles
    """
    
    # Unpack tuples
    un,  vn,  wn  = uvwn
    xn,  yn,  zn  = xyzn
    xc,  yc,  zc  = xyzc
    
    # Extract the physical parameters of the Fluid. 
    rho = rho[0,0,0]
    mu = mu[0,0,0]
    
    obst = np.nonzero(obst)

    # Iterate through all the particles.
    for i in range(0,n):
        
        (iu, il) = closest_node(xn, pt[i].x)
        (ju, jl) = closest_node(yn, pt[i].y)
        (ku, kl) = closest_node(zn, pt[i].z)
        
        
        x = closest_cell(xc, pt[i].x)
        y = closest_cell(yc, pt[i].y)
        z = closest_cell(zc, pt[i].z)
        
        # Check if Particle is in an Obstacle cell, if True then set the particle's 
        # velocity to zero.
        
        if (x in obst[0]) and (y in obst[1]) and (z in obst[2]):
            (pt[i].u,pt[i].v,pt[i].w) = (0, 0, 0)
        
        # Check to see if the particle is still in the domain, if not, the particle
        # sticks to the surrounding boundary
        
        elif (min(xn) > pt[i].x or  pt[i].x > max(xn)) or \
             (min(yn) > pt[i].y or pt[i].y > max(yn)) or \
             (min(zn) > pt[i].z or pt[i].z > max(zn)):
                 
                 (pt[i].u, pt[i].v, pt[i].w) = (0, 0, 0)
        
        
        # If the particle hasn't ran into trouble and hit a boundary, we calculate 
        # its change in position and velocity. 
        
        else:
    
            # Find the Fluid velocities defined at the node.
            (u,v,w) = lagrange_interpol((un,vn,wn), (xn,yn,zn), 
                                        (pt[i].x, pt[i].y, pt[i].z), 
                                         iu, il, ju, jl, ku, kl)
        
     
            # Calculating the new particle velocities, using a fourth order Runge-kutta 
            # scheme. The y-component has the added gravity force. 
            
            
            pt[i].u = rk4(u, pt[i].u, rho, mu, dt)
            
            pt[i].v = rk4(v, pt[i].v, rho, mu, dt)
            
            pt[i].w = rk4(w, pt[i].w, rho, mu, dt)
            
        
            # Updating the positions of the particles. 
        
            pt[i].x = pt[i].x + pt[i].u * dt 
            
            pt[i].y = pt[i].y + pt[i].v * dt 
            
            pt[i].z = pt[i].z + pt[i].w * dt
        
        
    return pt[i].x, pt[i].y, pt[i].z, pt[i].u, pt[i].v, pt[i].w # end
    






