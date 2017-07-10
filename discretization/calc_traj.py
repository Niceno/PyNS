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
from pyns.lagrangian import * 

# =============================================================================
def calc_traj(pt, uvwn, rho, mu, xyzn, xyzc, dt, obst, n):
# -----------------------------------------------------------------------------        
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
    
    obst = lnot(obst)

    # ----------------------------------
    # Iterate through all the particles
    # ----------------------------------
    for p in range(0, n):
        
        (iu, il) = closest_node(xn, pt[p].x)
        (ju, jl) = closest_node(yn, pt[p].y)
        (ku, kl) = closest_node(zn, pt[p].z)
                
        i = closest_cell(xc, pt[p].x)
        j = closest_cell(yc, pt[p].y)
        k = closest_cell(zc, pt[p].z)
        
        # Check if particle is in an obstacle cell;
        # if true set the particle's velocity to zero.        
        if (i in obst[X]) and (j in obst[Y]) and (k in obst[Z]):
            (pt[p].u,pt[p].v,pt[p].w) = (0, 0, 0)
        
        # Check to see if the particle is still in the domain; 
        # if not, the particle sticks to the surrounding boundary        
        elif (min(xn) > pt[p].x or pt[p].x > max(xn)) or \
             (min(yn) > pt[p].y or pt[p].y > max(yn)) or \
             (min(zn) > pt[p].z or pt[p].z > max(zn)):                 
                 (pt[p].u, pt[p].v, pt[p].w) = (0, 0, 0)
                
        # If the particle hasn't ran into trouble and hit a boundary,
        # calculate its change in position and velocity.         
        else:
    
            # Find the Fluid velocities defined at the node.
            (u,v,w) = lagrange_interpol((un,vn,wn), (xn,yn,zn), 
                                        (pt[p].x, pt[p].y, pt[p].z), 
                                         iu, il, ju, jl, ku, kl)
             
            # Calculating the new particle velocities, 
            # using a fourth order Runge-kutta scheme. 
            # The y-component has the added gravity force. 
            pt[p].u = rk4(u, pt[p].u, rho[i,j,k], mu[i,j,k], dt)            
            pt[p].v = rk4(v, pt[p].v, rho[i,j,k], mu[i,j,k], dt)            
            pt[p].w = rk4(w, pt[p].w, rho[i,j,k], mu[i,j,k], dt)
                    
            # Updating the positions of the particles.         
            pt[p].x = pt[p].x + pt[p].u * dt            
            pt[p].y = pt[p].y + pt[p].v * dt             
            pt[p].z = pt[p].z + pt[p].w * dt
                
    return pt[p].x, pt[p].y, pt[p].z, pt[p].u, pt[p].v, pt[p].w  # end
    






