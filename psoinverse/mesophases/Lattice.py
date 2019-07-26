# Imports
import numpy as np

class Lattice(object):
    """ Object representing a crystallographic basis vector lattice"""
    
    def __init__(self, basis, **kwargs):
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = basis
        super().__init__(**kwargs)
    
    
