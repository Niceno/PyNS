"""
Defines class of the type "Matrix".  "Matrix" is the object which holds the 
non-zero entries in a system matrix.  Essentially it stores a bundle of 
non-zero diagonals named after compass directions.
"""

# Standard Python modules
from pyns.standard import *

# PyNS modules
from pyns.constants import *
from pyns.operators import *
from pyns.display   import write

class Matrix:
    
    # =========================================================================
    def __init__(self, res):
    # -------------------------------------------------------------------------

        # -----------------------------------------------------------------
        # Number of diagonals.  
        #
        # Now set to seven, like default, which is: W, E, S, N, B, T and C
        #
        # In the future, we might think of:
        # 13: W, E, S, N, B, T, C, WW, EE, SS, NN, BB, TT 
        # 15: W, E, S, N, B, T, C, BSW, BSE, BNW, BNE, TSW, TSE, TNW, TNE 
        # 21: all of the above
        # -----------------------------------------------------------------
        self.ndiag = 7
        
        self.W = zeros(res)
        self.E = zeros(res)
        self.S = zeros(res)
        self.N = zeros(res)
        self.B = zeros(res)
        self.T = zeros(res)
        self.C = zeros(res)
    
        return  # end of function
