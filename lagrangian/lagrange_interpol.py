"""
1st order Lagrange interpolation.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *

# =============================================================================
def lagrange_interpol(uvwn, xyzn, xyzp, i0, i1, i2, j0, j1, j2, k0, k1, k2, second_order):
# -----------------------------------------------------------------------------  
    """
    The velocity at the particle's postition is calculted using a 1st order
    Lagrange interpolation method. 
    
    Args:
      vel:    The component of the velocity which one wishes to calculate.
      xyzn:   Tuple containing the position of the nodes. 
      xyzp:   Tuple containing the particle's position. 
      i0, i1, i2: Index of closest and neighbouring nodes (x - direction). 
      j0, j1, j2: Index of closest and neighbouring nodes (y - direction).
      k0, k1, k2: Index of closest and neighbouring nodes (z - direction).
      second_order: Truth variable which defines what order to interpolate the
                    velocity, depends on the value of i2, j2 and k2.  
        
    Returns:
      The interpolated velocity at the particle's position. 
    """
    
    # Unpack tuples
    un,  vn,  wn  = uvwn
    xn,  yn,  zn  = xyzn
    xp,  yp,  zp  = xyzp
    
    u = 0
    v = 0
    w = 0 

    if second_order is True:
          
        # Compute interpolation factors        
        lx_0 = ((xp - xn[i1]) / (xn[i0] - xn[i1]) 
              * (xp - xn[i2]) / (xn[i0] - xn[i2]))
        lx_1 = ((xp - xn[i0]) / (xn[i1] - xn[i0]) 
              * (xp - xn[i2]) / (xn[i1] - xn[i2]))
        lx_2 = ((xp - xn[i0]) / (xn[i2] - xn[i0]) 
             *  (xp - xn[i1]) / (xn[i2] - xn[i1]))
             
        # Place them in arrays for iterating over.
        lx = [lx_0, lx_1, lx_2]
        i = [i0, i1, i2]
        
        ly_0 = ((yp - yn[j1]) / (yn[j0] - yn[j1]) 
             *  (yp - yn[j2]) / (yn[j0] - yn[j2]))
        ly_1 = ((yp - yn[j0]) / (yn[j1] - yn[j0]) 
             *  (yp - yn[j2]) / (yn[j1] - yn[j2]))
        ly_2 = ((yp - yn[j0]) / (yn[j2] - yn[j0]) 
             *  (yp - yn[j1]) / (yn[j2] - yn[j1]))
             
        ly = [ly_0, ly_1, ly_2]
        j = [j0, j1, j2]
        
        lz_0 = ((zp - zn[k1]) / (zn[k0] - zn[k1]) 
             * (zp - zn[k2]) / (zn[k0] - zn[k2]))
        lz_1 = ((zp - zn[k0]) / (zn[k1] - zn[k0]) 
             *  (zp - zn[k2]) / (zn[k1] - zn[k2]))
        lz_2 = ((zp - zn[k0]) / (zn[k2] - zn[k0]) 
             *  (zp - zn[k1]) / (zn[k2] - zn[k1]))
             
        lz = [lz_0, lz_1, lz_2]
        k = [k0, k1, k2]

    
    else:

        # Compute interpolation factors        
        lx_0 = (xp - xn[i1]) / (xn[i0] - xn[i1])
        lx_1 = (xp - xn[i0]) / (xn[i1] - xn[i0])
        
        # Place them in arrays for iterating over.
        lx = [lx_0, lx_1]
        i = [i0, i1]
            
        ly_0 = (yp - yn[j1]) / (yn[j0] - yn[j1])
        ly_1 = (yp - yn[j0]) / (yn[j1] - yn[j0])
        
        ly = [ly_0, ly_1]
        j = [j0, j1]

        
        lz_0 = (zp - zn[k1]) / (zn[k0] - zn[k1])
        lz_1 = (zp - zn[k0]) / (zn[k1] - zn[k0])
        
        lz = [lz_0, lz_1]
        k = [k0, k1]
        
    # Using the interpolation factors and the values of the velocities at the
    # corresponding nodes, the interpolated velocities are computed. 
    for q in range(len(k)):
            for r in range(len(j)):
                for s in range(len(i)):
                    u += ((lx[s] * ly[r] * lz[q])) * (un[i[s], j[r], k[q]])
                    v += ((lx[s] * ly[r] * lz[q])) * (vn[i[s], j[r], k[q]])
                    w += ((lx[s] * ly[r] * lz[q])) * (wn[i[s], j[r], k[q]])
 
    return u, v, w  # end of function
    

