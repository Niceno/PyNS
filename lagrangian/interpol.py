#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple interpolator.
"""
from pyns.standard import *
 
from pyns.operators      import *

def interpol(vel, iu, il, ju, jl, ku, kl):
    """
    The fluid velocity at the particle's position is calculted as the 
    average of the surrounding defined velocities. 
    
    Args:
        vel: The component of the velocity which one wishes to calculate.
        iu, il: index of the upper and lower indices (x - direction) 
        ju, jl: index of the upper and lower indices (y - direction)
        ku, kl: index of the upper and lower indices (z - direction)
        
    Returns:
        The interpolated velocity. 
    """ 
    
    vel = (vel[iu, ju, ku] + vel[il, ju, ku] + vel[iu, jl, ku]
          + vel[il, jl, ku] + vel[iu, ju, kl] + vel[il, ju, ku]
          + vel[iu, jl, ku] + vel[il, jl, ku]) / 8 
                
    return vel
    

