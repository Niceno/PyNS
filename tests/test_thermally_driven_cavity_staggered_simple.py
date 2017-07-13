"""
This program solves thermally driven cavity at Ra = 1.0e6, in dimensional
and non-dimensional forms, for staggered variable arrangement.

Equations in dimensional form:

D(rho u)/Dt = nabla(mu (nabla u)^T) - nabla p + g
D(rho cp T)/Dt = nabla(lambda (nabla T)^T)

Equations in non-dimensional form, for natural convection problems

DU/Dt = nabla(1/sqrt(Gr) (nabla U)^T) - nabla P + theta
D theta/Dt = nabla(1/(Pr*sqrt(Gr)) (nabla theta)^T)

For thermally driven cavity, with properties of air at 60 deg:

nu   =  1.89035E-05;
beta =  0.003;
dT   = 17.126;
L    =  0.1;
Pr   = 0.709;

Characteristic non-dimensional numbers are:
Gr = 1.4105E+06
Ra = 1.0000E+06
"""

#!/usr/bin/python

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants      import *
from pyns.operators      import *
from pyns.discretization import *
from pyns.display        import plot, write
from pyns.physical       import properties
from pyns.solvers.norm   import norm
        
def main(show_plot=True, time_steps=12, plot_freq=1):

# =============================================================================
#
# Define problem
#
# =============================================================================

    xn = nodes(0, 1,  64, 1.0/256, 1.0/256)
    yn = nodes(0, 1,  64, 1.0/256, 1.0/256)
    zn = nodes(0, 0.1, 5)

    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    grashof = 1.4105E+06
    prandtl = 0.7058
    rho   = zeros(rc)
    mu    = zeros(rc)
    kappa = zeros(rc)
    cap   = zeros(rc)
    rho  [:,:,:] = 1.0
    mu   [:,:,:] = 1.0 / sqrt(grashof)
    kappa[:,:,:] = 1.0 / (prandtl * sqrt(grashof))
    cap  [:,:,:] = 1.0

    # Time-stepping parameters
    dt  = 2           # time step
    ndt = time_steps  # number of time steps

    # Create unknowns; names, positions and sizes
    uf    = Unknown('face-u-vel',     X, ru, DIRICHLET)
    vf    = Unknown('face-v-vel',     Y, rv, DIRICHLET)
    wf    = Unknown('face-w-vel',     Z, rw, DIRICHLET)
    t     = Unknown('temperature',    C, rc, NEUMANN)
    p     = Unknown('pressure',       C, rc, NEUMANN)
    p_tot = Unknown('total-pressure', C, rc, NEUMANN)

    # This is a new test
    t.bnd[W].typ[:] = DIRICHLET
    t.bnd[W].val[:] = -0.5

    t.bnd[E].typ[:] = DIRICHLET
    t.bnd[E].val[:] = +0.5

    for j in (B,T):
        uf.bnd[j].typ[:] = NEUMANN
        vf.bnd[j].typ[:] = NEUMANN
        wf.bnd[j].typ[:] = NEUMANN

    obstacle = None

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
        # Start inner iteration loop
        # ---------------------------
        
        # Allocate space for results from previous iteration
        t_prev = zeros(rc)
        u_prev = zeros(ru)
        v_prev = zeros(rv)
        w_prev = zeros(rw)

        for inner in range(1, 33):

            write.iteration(inner)
          
            # Store results from previous iteration
            t_prev[:] = t.val[:]
            u_prev[:] = uf.val[:]
            v_prev[:] = vf.val[:]
            w_prev[:] = wf.val[:]
          
            # Temperature (enthalpy)
            calc_t(t, (uf,vf,wf), (rho*cap), kappa, dt, (dx,dy,dz), obstacle,
                   advection_scheme = 'upwind',
                   under_relaxation = 0.75)
    
            # Momentum conservation
            ext_f = zeros(ru), avg(Y, t.val), zeros(rw)
    
            calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu, dt, (dx,dy,dz), obstacle,
                     pressure = p_tot,
                     force    = ext_f,
                     advection_scheme = 'upwind',
                     under_relaxation = 0.75)
    
            # Pressure
            calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obstacle,
                   verbatim = False)
    
            p_tot.val += 0.25 * p.val

            # Velocity correction
            corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), obstacle,
                     verbatim = False)

            # Print differences in results between two iterations
            t_diff = abs(t.val[:] - t_prev)
            u_diff = abs(uf.val[:] - u_prev)
            v_diff = abs(vf.val[:] - v_prev)
            w_diff = abs(wf.val[:] - w_prev)
            print("t_diff = ", norm(t_diff)/norm(rc))
            print("u_diff = ", norm(u_diff)/norm(ru))
            print("v_diff = ", norm(v_diff)/norm(rv))
            print("w_diff = ", norm(w_diff)/norm(rw))

        # --------------------------------------------------
        # Check the CFL number at the end of iteration loop
        # --------------------------------------------------
        cfl = cfl_max((uf,vf,wf), dt, (dx,dy,dz))

# =============================================================================
#
# Visualisation
#
# =============================================================================
        if show_plot:
            if ts % plot_freq == 0:
                plot.isolines(t.val, (uf, vf, wf), (xn, yn, zn), Z)
                plot.gmv("tdc-staggered-%6.6d" % ts, 
                         (xn, yn, zn), (uf, vf, wf, t))

if __name__ == '__main__':
    main()
