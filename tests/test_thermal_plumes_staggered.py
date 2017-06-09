#!/usr/bin/python

"""
This script is to reproduce two-dimensional mixed convection case, with
the aim of testing the outflow boundary, particularly the "convective"
boundary condition which allows eddies to leave the domain.

There is also a staggered version of this script, called
"demo_plums_collocated.m".  It would be good to keep both version as
similar as possible to each other, to test the differences between
staggered and collocated arrangements always possible.
"""

from standard import *

from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN, OUTLET
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.display.plot_isolines import plot_isolines
from scrins.discretization.adj_n_bnds import adj_n_bnds
from scrins.discretization.cartesian_grid import cartesian_grid
from scrins.discretization.nodes import nodes
from scrins.discretization.create_unknown import create_unknown
from scrins.discretization.cfl_max import cfl_max
from scrins.discretization.calc_p import calc_p
from scrins.discretization.calc_t import calc_t
from scrins.discretization.calc_uvw import calc_uvw
from scrins.discretization.corr_uvw import corr_uvw
from scrins.discretization.vol_balance import vol_balance
from scrins.display.print_time_step import print_time_step
from scrins.operators.avg import avg
from scrins.operators.par import par


def main(show_plot=True):
    """
    Docstring.
    """

    # =========================================================================
    #
    # Define problem
    #
    # =========================================================================

    # Node coordinates
    xn = nodes(0, 10, 300)
    yn = nodes(0, 1, 40, 1 / 500, 1 / 500)
    zn = nodes(0, 3, 3)

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)

    # Cell dimensions
    nx, ny, nz, dx, dy, dz, rc, ru, rv, rw = cartesian_grid(xn, yn, zn)

    # Set physical properties
    rho = zeros(rc)
    mu = zeros(rc)
    kappa = zeros(rc)
    cap = zeros(rc)
    rho[:, :, :] = 1.
    mu[:, :, :] = 0.1
    kappa[:, :, :] = 0.15
    cap[:, :, :] = 1.0

    # Time-stepping parameters
    dt = 0.003      # time step
    ndt = 1500         # number of time steps

    # Create unknowns; names, positions and sizes
    uf = create_unknown('face-u-vel', X, ru, DIRICHLET)
    vf = create_unknown('face-v-vel', Y, rv, DIRICHLET)
    wf = create_unknown('face-w-vel', Z, rw, DIRICHLET)
    t = create_unknown('temperature', C, rc, NEUMANN)
    p = create_unknown('pressure', C, rc, NEUMANN)

    # Specify boundary conditions
    uf.bnd[W].typ[:1, :, :] = DIRICHLET
    for k in range(0, nz):
        uf.bnd[W].val[:1, :, k] = par(1.0, yn);

    uf.bnd[E].typ[:1, :, :] = OUTLET
    uf.bnd[E].val[:1, :, :] = 1.0;

    for j in (B, T):
        uf.bnd[j].typ[:] = NEUMANN
        vf.bnd[j].typ[:] = NEUMANN
        wf.bnd[j].typ[:] = NEUMANN

    t.bnd[W].typ[:1, :, :] = DIRICHLET
    for k in range(0, nz):
        t.bnd[W].val[:1, :, k] = 1.0 - yc;

    t.bnd[S].typ[:, :1, :] = DIRICHLET
    t.bnd[S].val[:, :1, :] = +1.0;
    t.bnd[N].typ[:, :1, :] = DIRICHLET
    t.bnd[N].val[:, :1, :] = 0.0

    adj_n_bnds(t)
    adj_n_bnds(p)

    # Specify initial conditions
    uf.val[:, :, :] = 1.0
    t.val[:, :, :] = 0

    obst = zeros(rc)

    # =========================================================================
    #
    # Solution algorithm
    #
    # =========================================================================

    # ----------
    #
    # Time loop
    #
    # ----------
    for ts in range(1, ndt + 1):

        print_time_step(ts)

        # -----------------
        # Store old values
        # -----------------
        t.old[:] = t.val[:]
        uf.old[:] = uf.val[:]
        vf.old[:] = vf.val[:]
        wf.old[:] = wf.val[:]

        # -----------------------
        # Temperature (enthalpy)
        # -----------------------
        calc_t(t, (uf, vf, wf), (rho * cap), kappa, dt, (dx, dy, dz), obst)

        # ----------------------
        # Momentum conservation
        # ----------------------
        ef = zeros(ru), 150.0 * avg(Y, t.val), zeros(rw)

        calc_uvw((uf, vf, wf), (uf, vf, wf), rho, mu,
                 zeros(rc), ef, dt, (dx, dy, dz), obst)

        # ---------
        # Pressure
        # ---------
        calc_p(p, (uf, vf, wf), rho, dt, (dx, dy, dz), obst)

        # --------------------
        # Velocity correction
        # --------------------
        corr_uvw((uf, vf, wf), p, rho, dt, (dx, dy, dz), obst)

        # Compute volume balance for checking
        err = vol_balance((uf, vf, wf), (dx, dy, dz), obst)
        print('Maximum volume error after correction: %12.5e' % abs(err).max())

        # Check the CFL number too
        cfl = cfl_max((uf, vf, wf), dt, (dx, dy, dz))
        print('Maximum CFL number: %12.5e' % cfl)

    # =========================================================================
    #
    # Visualisation
    #
    # =========================================================================

        if ts % 150 == 0:
            plot_isolines(t.val, (uf, vf, wf), (xn, yn, zn), Z)
            plot_isolines(p.val, (uf, vf, wf), (xn, yn, zn), Z)


if __name__ == '__main__':
    main()
