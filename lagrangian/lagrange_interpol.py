#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lagrange interpolation.
"""


from pyns.standard import *
 
from pyns.operators      import *

def lagrange_interpol(uvwn, xyzn, xyzp, iu, il, ju, jl, ku, kl,):
    """
    The velocity at the particle's postition is calculted using a
    Lagrange interpolation method. 
    
    Args:
        vel: The component of the velocity which one wishes to calculate.
        xyzn: Tuple containing the position of the nodes. 
        xyzp: Tuple containing the particle's position. 
        iu, il: index of the upper and lower indices (x - direction) 
        ju, jl: index of the upper and lower indices (y - direction)
        ku, kl: index of the upper and lower indices (z - direction)
        
    Returns:
        The interpolated velocity at the particle's position. 
    """
    
    un,  vn,  wn  = uvwn
    xn,  yn,  zn  = xyzn
    xp,  yp,  zp  = xyzp
    
    # Need to find L_0 and L_1.
    # Check to see if the upper ndex is greater than the lower index. 
    if iu < il:
        iu, il = il, iu
    
    if ju < jl:
        ju, jl = jl, ju
    
    if ku < kl:
        ku, kl = kl, ku

        
    lx_l = (xp - xn[iu]) / (xn[il] - xn[iu])
    lx_u = (xp - xn[il]) / (xn[iu] - xn[il])
    
    
    ly_l = (yp - yn[ju]) / (yn[jl] - yn[ju])
    ly_u = (yp - yn[jl]) / (yn[ju] - yn[jl])
    
    lz_l = (zp - zn[ku]) / (zn[kl] - zn[ku])
    lz_u = (zp - zn[kl]) / (zn[ku] - zn[kl])
    
    
    u = ((lx_u * ly_u * lz_u) * (un[iu, ju, ku]) 
        + (lx_l * ly_u * lz_u) * (un[il, ju, ku])
        + (lx_u * ly_l * lz_u) * (un[iu, jl, ku])
        + (lx_l * ly_l * lz_u) * (un[il, jl, ku])
        + (lx_u * ly_u * lz_l) * (un[iu, ju, kl])
        + (lx_l * ly_u * lz_l) * (un[il, ju, kl])
        + (lx_u * ly_l * lz_l) * (un[iu, jl, kl])
        + (lx_l * ly_l * lz_l) * (un[il, jl, kl]))


    v = ((lx_u * ly_u * lz_u) * (vn[iu, ju, ku]) 
        + (lx_l * ly_u * lz_u) * (vn[il, ju, ku])
        + (lx_u * ly_l * lz_u) * (vn[iu, jl, ku])
        + (lx_l * ly_l * lz_u) * (vn[il, jl, ku])
        + (lx_u * ly_u * lz_l) * (vn[iu, ju, kl])
        + (lx_l * ly_u * lz_l) * (vn[il, ju, kl])
        + (lx_u * ly_l * lz_l) * (vn[iu, jl, kl])
        + (lx_l * ly_l * lz_l) * (vn[il, jl, kl]))
    
    w = ((lx_u * ly_u * lz_u) * (wn[iu, ju, ku]) 
        + (lx_l * ly_u * lz_u) * (wn[il, ju, ku])
        + (lx_u * ly_l * lz_u) * (wn[iu, jl, ku])
        + (lx_l * ly_l * lz_u) * (wn[il, jl, ku])
        + (lx_u * ly_u * lz_l) * (wn[iu, ju, kl])
        + (lx_l * ly_u * lz_l) * (wn[il, ju, kl])
        + (lx_u * ly_l * lz_l) * (wn[iu, jl, kl])
        + (lx_l * ly_l * lz_l) * (wn[il, jl, kl]))
 
    return u, v, w
    

