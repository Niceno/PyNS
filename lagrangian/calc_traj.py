#!/usr/bin/env python3
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
from pyns.physical.constants.gravitational import G

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
            
        (i0, i1, i2) = closest_node(xn, pt[p].x)
        (j0, j1, j2) = closest_node(yn, pt[p].y)
        (k0, k1, k2) = closest_node(zn, pt[p].z)
    
        
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
            
            # Check to see if the particle's closest node is a boundary node,
            # If so a first order interpolation scheme is implemented.
            if i2 == None or j2 == None or k2 == None:
                second_order = False
                 
            else:
                second_order = True

            # Calculate the interpolated fluid velocity at the particles center
            (u,v,w) = lagrange_interpol((un,vn,wn), (xn,yn,zn), 
                                        (pt[p].x, pt[p].y, pt[p].z), 
                                         i0, i1, i2, j0, j1, j2, k0, k1, k2, 
                                         second_order)
            
            # Calculating the new particle velocities, 
            # using a fourth order Runge-kutta scheme. 
            # The y-component has the added gravity term. 
            pt[p].u = rk4(u, pt[p].u, rho[i,j,k], mu[i,j,k], pt[p].d, pt[p].rho_p, dt)            
            pt[p].v = rk4(v, pt[p].v, rho[i,j,k], mu[i,j,k], pt[p].d, pt[p].rho_p, dt, G)            
            pt[p].w = rk4(w, pt[p].w, rho[i,j,k], mu[i,j,k], pt[p].d, pt[p].rho_p, dt)
                    
            # Updating the positions of the particles.         
            pt[p].x = pt[p].x + pt[p].u * dt            
            pt[p].y = pt[p].y + pt[p].v * dt             
            pt[p].z = pt[p].z + pt[p].w * dt
            
        
            
    return pt[p].x, pt[p].y, pt[p].z, pt[p].u, pt[p].v, pt[p].w  # end







