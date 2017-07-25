"""
A simple fourth order Runge-Kutta scheme to solve 
dv/dt = ((v-u)/ tau)
"""

# Standard Python modules
from pyns.standard import *

# ============================================================================= 
def rk4(vel, vel_p, rho, mu, d, rho_p, dt, gravity = None):
# -----------------------------------------------------------------------------  
    """
    Fourth-order Runge-Kutta method to solve for the particles velocity. 
    
    Args:
        vel: The interpolated velocity of the fluid at the particle's 
        position. 
        vel_p: The particle's velocity.
        rho/mu: The physical parameters of the fluid. 
        dt: Time step used for solving the fluid velocities.
        gravity: Is an optional parameter
        
    Returns:
        Updated velocity of the particle. 
    """
    
    # Defining the physical parameters of the particles.
     
    
    Re_p = (rho * (d) * (abs(vel - vel_p))) / mu
    
    f = 1 + 0.15 * Re_p**(0.687)
    
    # Relaxation Time
    tau = (rho_p) * ((d)**2) / (18 * mu)
    
    
    # Set up an array which will store the intermediate values.
    inter = []
    
    if gravity is None:
        
        # tau is the relaxation time
        k1 = f * (vel - vel_p) / tau
        inter = vel_p + k1 * (dt / 2)
            
        k2 = f * (vel - inter) / tau
        
        # Append to the intermediate values.
        inter = append(inter, vel_p + k2 * (dt / 2))
        
        k3 = f * (vel - inter[1]) / tau
        
        # Append to the intermediate values.
        inter = append(inter, vel_p + k3 * dt)
        
        k4 = f * (vel - inter[2]) / tau
        
        vel_p = vel_p + (( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0) * dt
        
        
    else:
        # tau is the relaxation time
        k1 = (f * (vel - vel_p) / tau) - gravity
        inter = vel_p + k1 * (dt / 2)
            
        k2 = (f * (vel - inter) / tau) - gravity
        
        # Append to the intermediate values.
        inter = append(inter, vel_p + k2 * (dt / 2))
        
        k3 = (f * (vel - inter[1]) / tau) - gravity
        
        # Append to the intermediate values.
        inter = append(inter, vel_p + k3 * dt)
        
        k4 = (f * (vel - inter[2]) / tau) - gravity
        
        vel_p = vel_p + (( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0) * dt
        

    return vel_p  
