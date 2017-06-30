#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates the closest node to the particles position. 
It is important when interpolating the particles velocity. 
"""

# Importing the numpy module
import numpy as np 

def closest_node(nodes, particle):
    """
    Args:
        nodes: These are the nodes coordinates in the x, y or z direction.
        particle: This is the particle location in the x, y or z direction.
        
    Returns:
        The index of the cell which the particle is closest to, in the x, y or 
        z direction.
    """
    
    # Returns the index of the closest node.
    idx = (np.abs(nodes-particle)).argmin()
    
    # For the interpolater we need to find the index of the node on the 
    # other side to the closest node.
    #
    # o = particle
    # | = nodes
    #
    #
    # ---|------|------|------|------|----o-|------|------|------|------|
    #                               idx1   idx

    
    if nodes[idx] > particle:
        idx1 = idx - 1 
    else:
        idx1 = idx + 1 
        
    return idx, idx1

