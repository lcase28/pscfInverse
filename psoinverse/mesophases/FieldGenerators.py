# Imports
from abc import ABC, abstractmethod
from crystals import Lattice, LatticeSystem
import numpy as np

#class FieldGenerator(ABC):
#    """ Abstract base class for field generators """
#    
#    def __init__(**kwargs):
#        lat = kwargs.get("lattice", None)
#        if lat is not None:
#            self.lattice = 
#        super().__init(**kwargs)
#    
#    @abstractmethod
#    def WriteFile(**kwargs):
#        pass
#    
#    def DensityFields(self, **kwargs):
    

class FieldGenerator(object):
    """ Generator class for 3D k-grid density fields of diblock systems """
    
    def __init__(self, lattice, particlePositions, **kwargs):
        self.lattice = lattice
        super().__init__(**kwargs)
    
    def to_kgrid(self, ngrid):
        
        for (i,x) in enumerate(ngrid):
            if i == 0:
                ngrid[i] = x/2
            else:
                ngrid[i] = x - 1
            
        nwaves = np.prod(ngrid)
        field
        
        # primary loop for n-dimensional generation
        for G in enumerate(itertools.product(*[range(x) for x in ngrid])):
            # fill in loop operations - G is miller indices in n-dimensions.
            
    
