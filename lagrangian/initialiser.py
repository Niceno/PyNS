#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates the particles with certain initial position and velocities. 
"""
from pyns.standard import *

from pyns.lagrangian import *
from pyns.discretization import *

import random

def initialiser(x, show_initial= []):
    """
    Args:
        x:   The number particle's that one wishes to generate. 
        show_initial: Boolen, when True the particle's initial properties are 
                      printed. 
    
    Returns:
      Particles.
    """
    pt = [Particles(random.uniform(0.4,0.6), random.uniform(0.1,0.9), 
                    random.uniform(0.01,0.24), 0, 0, 0) for i in range(0,x)]
    
    if show_initial:
        for i in range(0,x):
            print(pt[i].x, pt[i].y, pt[i].z, pt[i].u, pt[i].v, pt[i].w)
        
    return pt
