#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particles: This Class object of a Particle.
TO DO: Add constraints to how close particles can get to each other. 
       Add definition of the tau variable
"""
# Standard Python modules
from pyns.standard import *

from pyns.constants import *
from pyns.discretization import *



class Particles(object):
    
    number = []
    
    def __init__(self, x, y, z, u, v, w, rho_p, d):
        """
        Initialising a particle.
        """
        self.x = x
        self.y = y
        self.z = z
        self.u = u
        self.v = v
        self.w = w
        self.rho_p = rho_p
        self.d = d 
      
        self.number.append(self)

   