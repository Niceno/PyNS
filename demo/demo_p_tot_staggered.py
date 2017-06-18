"""
Demonstrates the variation of projection algorythm which computes total
pressure, as a sum of all pressure corrections.

The total pressure being built up in this way counter-balances the
gravity term in momentum equations.

It seems that such an approach is important for bouyancy dominated flows.

Gravity term is under-relaxed here, but it works even without it.
"""

#!/usr/bin/python

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants          import *
from pyns.operators          import *
from pyns.discretization     import *
from pyns.display            import plot, write
from pyns.physical           import properties
from pyns.physical.constants import G

# =============================================================================
#
# Define problem
#
# =============================================================================

# Node coordinates
xn = nodes(0, 10,   80)
yn = nodes(0,  1,   20)
zn = nodes(0,  0.5,  5)

# Cell coordinates
xc = avg(xn)
yc = avg(yn)
zc = avg(zn)

# Cell dimensions
nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

# Set physical properties
rho   = zeros(rc)
mu    = zeros(rc)
kappa = zeros(rc)
cap   = zeros(rc)
rho   [:] =    1.15  # density              [kg/m^3]
mu    [:] =    0.1   # viscosity            [Pa s]

# Time-stepping parameters
dt  =   0.1  # time step
ndt = 200    # number of time steps

# Create unknowns; names, positions and sizes
uf = Unknown('face-u-vel',  X, ru, DIRICHLET)
vf = Unknown('face-v-vel',  Y, rv, DIRICHLET)
wf = Unknown('face-w-vel',  Z, rw, DIRICHLET)
p  = Unknown('pressure',    C, rc, NEUMANN)
p_tot = zeros(rc)

# Specify boundary conditions
uf.bnd[W].typ[:1,:,:] = DIRICHLET
for k in range(0,nz):
    uf.bnd[W].val[:1,:,k]  = par(0.1, yn)

uf.bnd[E].typ[:1,:,:] = OUTLET

for j in (B,T):
    uf.bnd[j].typ[:] = NEUMANN
    vf.bnd[j].typ[:] = NEUMANN
    wf.bnd[j].typ[:] = NEUMANN

adj_n_bnds(p)

obst = zeros(rc)

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
    g_v = -G * avg(Y, rho) * min(ts/100,1)

    ef = zeros(ru), g_v, zeros(rw)

    calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu,  \
             p_tot, dt, (dx,dy,dz), obst,
             force = ef)

    # ---------
    # Pressure
    # ---------
    calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obst)

    p_tot = p_tot + p.val

    # --------------------
    # Velocity correction
    # --------------------
    corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), obst)

    # Check the CFL number too
    cfl = cfl_max((uf,vf,wf), dt, (dx,dy,dz))

# =============================================================================
#
# Visualisation
#
# =============================================================================

    if ts % 20 == 0:
        plot.isolines(p_tot, (uf,vf,wf), (xn,yn,zn), Z)
        plot.isolines(p.val, (uf,vf,wf), (xn,yn,zn), Z)
