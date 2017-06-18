"""
Demonstrates the membrane for Kerstin.  Two domains are computed
independently, but linked through boundary conditions in an inner loop
within each time step.  This seems the most practical approach of all
because the membrane model implementations are very obvious, at one
place, in the main function.  The convergence of the conditions at the
membrane seems to be rather fast too.
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

AIR = 0
H2O = 1

# Node coordinates for both domains
xn = (nodes(0, 0.1,  80), nodes(0,    0.1,  80))
yn = (nodes(0, 0.01, 20), nodes(0.01, 0.03, 30))
zn = (nodes(0, 0.1,  80), nodes(0,    0.1,  80))

# Cell coordinates
xc = (avg(xn[AIR]), avg(xn[H2O]))
yc = (avg(yn[AIR]), avg(yn[H2O]))
zc = (avg(zn[AIR]), avg(zn[H2O]))

# Cell dimensions
cell = [cartesian_grid(xn[AIR],yn[AIR],zn[AIR]),  \
        cartesian_grid(xn[H2O],yn[H2O],zn[H2O])]

nx,ny,nz, dx,dy,dz = (cell[AIR][0], cell[H2O][0]),  \
                     (cell[AIR][1], cell[H2O][1]),  \
                     (cell[AIR][2], cell[H2O][2]),  \
                     (cell[AIR][3], cell[H2O][3]),  \
                     (cell[AIR][4], cell[H2O][4]),  \
                     (cell[AIR][5], cell[H2O][5])
rc,ru,rv,rw =        (cell[AIR][6], cell[H2O][6]),  \
                     (cell[AIR][7], cell[H2O][7]),  \
                     (cell[AIR][8], cell[H2O][8]),  \
                     (cell[AIR][9], cell[H2O][9])

# Set physical properties
prop = [properties.air(rc[AIR]), properties.water(rc[H2O])]
rho, mu, cap, kappa = (prop[AIR][0], prop[H2O][0]),  \
                      (prop[AIR][1], prop[H2O][1]),  \
                      (prop[AIR][2], prop[H2O][2]),  \
                      (prop[AIR][3], prop[H2O][3])

# Time-stepping parameters
dt  =    0.004  # time step
ndt = 1000      # number of time steps

# Create unknowns; names, positions and sizes
uf    = [Unknown('face-u-vel',     X, ru[AIR], DIRICHLET),  \
         Unknown('face-u-vel',     X, ru[H2O], DIRICHLET)]
vf    = [Unknown('face-v-vel',     Y, rv[AIR], DIRICHLET),  \
         Unknown('face-v-vel',     Y, rv[H2O], DIRICHLET)]
wf    = [Unknown('face-w-vel',     Z, rw[AIR], DIRICHLET),  \
         Unknown('face-w-vel',     Z, rw[H2O], DIRICHLET)]
p     = [Unknown('pressure',       C, rc[AIR], NEUMANN),  \
         Unknown('pressure',       C, rc[H2O], NEUMANN)]
t     = [Unknown('temperature',    C, rc[AIR], NEUMANN),  \
         Unknown('temperature',    C, rc[H2O], NEUMANN)]
p_tot = [Unknown('total-pressure', C, rc[AIR], NEUMANN),  \
         Unknown('total-pressure', C, rc[H2O], NEUMANN)]

# Specify boundary conditions
for k in range(0,nz[AIR]):
    vf[AIR].bnd[N].val[:,:1,k]=-0.005
uf[AIR].bnd[E].typ[:1,:,:] = OUTLET
uf[AIR].bnd[E].val[:1,:,:] = 0.1

for k in range(0,nz[H2O]):
    uf[H2O].bnd[W].val[:1,:,k] = par(0.1, yn[H2O])
uf[H2O].bnd[E].typ[:1,:,:] = OUTLET
uf[H2O].bnd[E].val[:1,:,:] = 0.1

for c in range(AIR,H2O):
    for j in (B,T):
        uf[c].bnd[j].typ[:] = NEUMANN
        vf[c].bnd[j].typ[:] = NEUMANN
        wf[c].bnd[j].typ[:] = NEUMANN

t[AIR].bnd[S].typ[:,:1,:] = DIRICHLET
t[AIR].bnd[S].val[:,:1,:] = 60
t[AIR].bnd[N].typ[:,:1,:] = DIRICHLET
t[AIR].bnd[N].val[:,:1,:] = 70

t[H2O].bnd[W].typ[:1,:,:] = DIRICHLET
t[H2O].bnd[W].val[:1,:,:] = 80
t[H2O].bnd[S].typ[:,:1,:] = DIRICHLET
t[H2O].bnd[S].val[:,:1,:] = 70

t[AIR].val[:,:,:] = 50;
t[H2O].val[:,:,:] = 70;

for c in (AIR,H2O):
    adj_n_bnds(p[c])
    adj_n_bnds(t[c])

obst = [zeros(rc[AIR]), zeros(rc[H2O])]

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
    for c in (AIR,H2O):
        t[c].old[:]  = t[c].val
        uf[c].old[:] = uf[c].val
        vf[c].old[:] = vf[c].val
        wf[c].old[:] = wf[c].val

    # -----------------------
    # Temperature (enthalpy)
    # -----------------------

    # Compute old temperature of the membrane
    const = kappa[H2O][:, :1,:] * dy[AIR][:,-1:,:]  \
          / kappa[AIR][:,-1:,:] / dy[H2O][:, :1,:]
    t_mem_old = (const * t[H2O].val[:,:1,:] + t[AIR].val[:,-1:,:])  \
              / (1+const)

    while True:
        for c in (AIR,H2O):
            calc_t(t[c], (uf[c],vf[c],wf[c]), (rho[c]*cap[c]), kappa[c],  \
                   dt, (dx[c],dy[c],dz[c]), obst[c])

        # Compute new temperature of the membrane
        const = kappa[H2O][:, :1,:] * dy[AIR][:,-1:,:]  \
              / kappa[AIR][:,-1:,:] / dy[H2O][:, :1,:]
        t_mem = (const * t[H2O].val[:,:1,:] + t[AIR].val[:,-1:,:])  \
              / (1+const)

        # Update boundary conditions with membrane's temperature
        t[H2O].bnd[S].val[:,:1,:] = t_mem
        t[AIR].bnd[N].val[:,:1,:] = t_mem

        # Check if convergence has been reached
        eps = abs(t_mem-t_mem_old).max()
        print('Maximum temperature error in the membrane = %5.3e' % eps)
        if eps < 0.001:
            break

        t_mem_old = t_mem

    # ----------------------
    # Momentum conservation
    # ----------------------
    for c in (AIR,H2O):
        g_v = -G * avg(Y, rho[c])

        ext_f = zeros(ru[c]), g_v, zeros(rw[c])

        calc_uvw((uf[c],vf[c],wf[c]), (uf[c],vf[c],wf[c]), rho[c], mu[c],  \
                 dt, (dx[c],dy[c],dz[c]), obst[c],
                 pressure = p_tot[c],
                 force    = ext_f)

    # ---------
    # Pressure
    # ---------
    for c in (AIR,H2O):
        calc_p(p[c], (uf[c],vf[c],wf[c]), rho[c],  \
               dt, (dx[c],dy[c],dz[c]), obst[c])

        p_tot[c].val += p[c].val

    # --------------------
    # Velocity correction
    # --------------------
    for c in (AIR,H2O):
        corr_uvw((uf[c],vf[c],wf[c]), p[c], rho[c],  \
                 dt, (dx[c],dy[c],dz[c]), obst[c])

    # Check the CFL number too
    for c in (AIR,H2O):
        cfl = cfl_max((uf[c],vf[c],wf[c]), dt, (dx[c],dy[c],dz[c]))

# =============================================================================
#
# Visualisation
#
# =============================================================================

    if ts % 100 == 0:
        for c in (AIR,H2O):
            plot.isolines(t[c].val, (uf[c],vf[c],wf[c]), (xn[c],yn[c],zn[c]), Z)
            plot.isolines(p_tot[c], (uf[c],vf[c],wf[c]), (xn[c],yn[c],zn[c]), Z)
