# Lattice Module. 
# Partially derived from crystals package available at:
#   https://github.com/LaurentRDC/crystals.git
# Generallized to handle 1-, 2-, or 3-dimensional systems

# Imports
from abc import ABC, abstractmethod
import numpy as np

class Lattice(ABC):
    """ Abstract base class of n-dimensional lattices """
    
    def __init__(self, **kwargs):
        super().__init__(kwargs)
    
    @abstractmethod
    def __repr__(self):
        return NotImplemented
    
    @abstractmethod
    def __hash__(self):
        return super().hash()
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            reutrn np.allclose(self.lattice_vectors, other.lattice_vectors, atol=1e-4)
        return NotImplemented
    
    @abstractmethod
    def __array__(self, *args, **kwargs):
        pass
    
class Lattice3D(Lattice):
    """ Object representing a crystallographic basis vector lattice"""
    
    def __init__(self, basis, **kwargs):
        a1, a2, a3 = basis
        self.a1 = np.asarray(a1, dtype=np.float)
        self.a2 = np.asarray(a2, dtype=np.float)
        self.a3 = np.asarray(a3, dtype=np.float)
        super().__init__(**kwargs)
    
    def __repr__(self):
        s = "< Lattice object with parameters {:.3f}A, {:.3f}, {:.3f}, {:.3f}, {:.3f} >"
        return s.format(*self.lattice_parameters)
    
    def __hash__(self):
        return hash(self.lattice_parameters) | super().hash()
        
    def __array__(self, *args, **kwargs):
        """ returns a 3x3 float array. Each row is a lattice vector. """
        return np.array(self.lattice_vectors, *args, *kwargs)
    
    @classmethod
    def from_parameters(cls, a, b, c, alpha, beta, gamma):
        
    
    
    
