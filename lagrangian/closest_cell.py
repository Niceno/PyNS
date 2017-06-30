#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the closest cell centre to the particles position. 
It is important when checking whether the particle is inside the obstacle or 
not.
"""

# Importing the numpy module 
import numpy as np

def closest_cell(cells, particle_position):
    """
    Args:
        cells: These are the cells coordinates in the x, y or z direction.
        particle: This is the particle location in the x, y or z direction.
        
    Returns:
        The index of the cell which the particle is closest to.
    """
    idxc = (np.abs(cells-particle_position)).argmin()
        
    return idxc

