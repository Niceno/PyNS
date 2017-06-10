"""
This constant sets iterative linear solver (such as "bicgstab") tolerance.

If very tight values are needed (often the case for pressure solution), one
can use their product (TOL*TOL).  To loosen the tolerance up, sqrt(TOL) is
possible.
"""

# Solver tolerance
TOL = 1.0e-8
