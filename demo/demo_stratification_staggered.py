"""
Solves flow in a channel with stable or unstable stratification.
The type of stratification is set by parameter "STRATIFICATION"

Uses non-Boussinesq model for buoyancy term; i.e. density depends on
temperature.

Gravity term is engaged gradually to avoid vortex at the outlet.
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

STRATIFICATION = 's'  # 's' or 'u'

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
cap   = zeros(rc)
kappa = zeros(rc)
rho   [:] =    1.25   # density              [kg/m^3]
mu    [:] =    0.1    # viscosity            [Pa s]
cap   [:] =    1.0    # capacity             [kg/m^3]
kappa [:] =    0.001  # thermal conductivity [W/mK]

# Time-stepping parameters
dt  =   0.2  # time step
ndt = 500    # number of time steps

# Create unknowns; names, positions and sizes
uf    = Unknown('face-u-vel',     X, ru, DIRICHLET)
vf    = Unknown('face-v-vel',     Y, rv, DIRICHLET)
wf    = Unknown('face-w-vel',     Z, rw, DIRICHLET)
p     = Unknown('pressure',       C, rc, NEUMANN)
t     = Unknown('temperature',    C, rc, NEUMANN)
p_tot = Unknown('total-pressure', C, rc, NEUMANN)

# Specify boundary conditions
uf.bnd[W].typ[:1,:,:] = DIRICHLET
for k in range(0,nz):
    uf.bnd[W].val[:1,:,k]  = par(0.1, yn)

uf.bnd[E].typ[:1,:,:] = OUTLET
uf.bnd[E].val[:1,:,:] = 0.1

for j in (B,T):
    uf.bnd[j].typ[:] = NEUMANN
    vf.bnd[j].typ[:] = NEUMANN
    wf.bnd[j].typ[:] = NEUMANN

t.val[:] = 70
t.bnd[W].typ[:1,:,:] = DIRICHLET
t.bnd[S].typ[:,:1,:] = DIRICHLET
t.bnd[N].typ[:,:1,:] = DIRICHLET
t.bnd[E].val[:1,:,:] = 70

if STRATIFICATION == 'u':
    dtemp = (60-80)/ny
    for k in range(0,nz):
        t.bnd[W].val[:1,:,k] = linspace(60-dtemp/2,80+dtemp/2,ny)
    t.bnd[S].val[:,:1,:] = 60
    t.bnd[N].val[:,:1,:] = 80

elif STRATIFICATION == 's':
    dtemp = (80-60)/ny
    for k in range(0,nz):
        t.bnd[W].val[:1,:,k] = linspace(80-dtemp/2,60+dtemp/2,ny)
    t.bnd[S].val[:,:1,:] = 80
    t.bnd[N].val[:,:1,:] = 60

adj_n_bnds(p)
adj_n_bnds(t)

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
    t.old[:]  = t.val[:]
    uf.old[:] = uf.val[:]
    vf.old[:] = vf.val[:]
    wf.old[:] = wf.val[:]

    # ---------------------------
    # Update physical properties
    # ---------------------------
    rho[:] = 1.06 + (t.val[:] - 60.0) * (0.99 - 1.06) / 20.0

    # -----------------------
    # Temperature (enthalpy)
    # -----------------------
    calc_t(t, (uf,vf,wf), (rho*cap), kappa, dt, (dx,dy,dz), obst)

    # ----------------------
    # Momentum conservation
    # ----------------------
    g_v = -G * avg(Y, rho)

    ext_f = zeros(ru), g_v, zeros(rw)

    calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu, dt, (dx,dy,dz), obst,
             pressure = p_tot,
             force    = ext_f)

    # ---------
    # Pressure
    # ---------
    calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obst)

    p_tot.val += p.val

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

    if ts % 50 == 0:
        plot.isolines(p_tot.val, (uf,vf,wf), (xn,yn,zn), Z)
        plot.isolines(t.val,     (uf,vf,wf), (xn,yn,zn), Z)
        plot.isolines(rho,       (uf,vf,wf), (xn,yn,zn), Z)
