"""
Function to plot a two-dimensional slice through solution.  It always plots
velocity fields on top of isolines.
"""

# Standard Python modules
from standard import *

# ScriNS modules
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.operators.avg import avg
from scrins.operators.avg_x import avg_x 
from scrins.operators.avg_y import avg_y
from scrins.operators.avg_z import avg_z
from scrins.operators.cat import cat
from scrins.operators.cat_x import cat_x
from scrins.operators.cat_y import cat_y
from scrins.operators.cat_z import cat_z


# =============================================================================


def plot_isolines(phi, uvw, xyzn, dir):
    # -----------------------------------------------------------------------------
    """
    Args:
      phi:  Unknown to be plotted.
      uvw:  Tuple containing three velocity components.  They can be either
            collocated or staggered.
      xyzn: Tuple containing one-dimensional arrays with "x", "y" and "z"
            coordinates.
      dir:    Direction for cutting, can be X, Y or Z.

    Returns:
      none!
    """

    # Unpack tuples
    u, v, w = uvw
    xn, yn, zn = xyzn

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)

    # Collocated velocity components
    if u.pos == C:
        uc = u.val
        vc = v.val
        wc = w.val
    else:
        uc = avg(
            X, cat_x((u.bnd[W].val[:1, :, :], u.val, u.bnd[E].val[:1, :, :])))
        vc = avg(
            Y, cat_y((v.bnd[S].val[:, :1, :], v.val, v.bnd[N].val[:, :1, :])))
        wc = avg(
            Z, cat_z((w.bnd[B].val[:, :, :1], w.val, w.bnd[T].val[:, :, :1])))

    # Pick coordinates for plotting (xp, yp) and values for plotting
    if dir == Y:
        jp = floor(yc.size / 2)
        xp, yp = meshgrid(xc, zc)
        zp = transpose(phi[:, jp, :], (1, 0))
        up = transpose(uc[:, jp, :], (1, 0))
        vp = transpose(wc[:, jp, :], (1, 0))
    if dir == Z:
        kp = floor(zc.size / 2)
        xp, yp = meshgrid(xc, yc)
        zp = transpose(phi[:, :, kp], (1, 0))
        up = transpose(uc[:, :, kp], (1, 0))
        vp = transpose(vc[:, :, kp], (1, 0))

    # Set levels and normalize the colors
    levels = linspace(zp.min(), zp.max(), 11)
    norm = cm.colors.Normalize(vmax=zp.max(), vmin=zp.min())

    plt.figure()
    plt.gca(aspect='equal')
    plt.contour(xp, yp, zp, levels, cmap=plt.cm.rainbow, norm=norm)
    plt.quiver(xp, yp, up, vp)
    plt.axis([min(xn), max(xn), min(yn), max(yn)])
    plt.show()

    return  # end of function
