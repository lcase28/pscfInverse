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
    
def sphere_form_factor(qR):
    out = 3 * (np.sin(qR) - qR * np.cos(qR)) / qR**3
    return out

class FieldGenerator(object):
    """ Generator class for 3D k-grid density fields of diblock systems """
    
    def __init__(self, lattice, **kwargs):
        self.lattice = lattice
        self.particles = kwargs.get("particlePositions",None)
        self.nparticles = kwargs.get("N_particles")
        self.nspecies = kwargs.get("N_monomer")
        self.smear = kwargs.get("sigma_smear")
        self.dim = kwargs.get("dim")
        self.crystal_system = kwargs.get("crystal_system")
        self.group_name = keargs.get("group_name")
        self.ngrid = kwargs.get("ngrid")
        self.formfactor = kwargs.get("formfactor", sphere_form_factor)
        super().__init__(**kwargs)
    
    @classmethod
    def from_file(cls, fname):
        with open(fname) as f:
            for line in f.readlines():
                line = line.strip()
                splitline = line.split(maxsplit=1)
                if len(splitline) > 1:
                    key = splitline[0].strip()
                    data = splitline[1].strip().split()
                    #interpret key for data parsing
                    #   most keys should be directly passable to __init__
                    #   some require additional treatment (lattice, particle positions)
                    
                
    
    # TODO: Figure out how to generate 2D, 1D initial guesses
    def to_kgrid(self, frac):
        ngrid = self.ngrid
        # Shift grid for k-grid
        kgrid = np.zeros_like(ngrid)
        for (i,x) in enumerate(ngrid):
            if i == 0:
                kgrid[i] = x/2
            else:
                kgrid[i] = x - 1
            
        nwaves = np.prod(kgrid + np.ones_like(kgrid))
        rho = np.zeros(nwave, self.nspecies)
        vol = self.lattice.volume
        Rsph = ((3 * frac[-1]) / (4 * np.pi * self.nparticles))**(1./3)
        a, b, c, alpha, beta, gamma = self.lattice.lattice_parameters
        
        # primary loop for n-dimensional generation
        t = 0
        for G in itertools.product(*[range(x) for x in kgrid]):
            # fill in loop operations - G is miller indices in n-dimensions.
            G = np.array(G) #convert tuple to array
            brillouin = self.miller_to_brillouin(G, ngrid)
            t = t + 1
            if np.array_equiv(brillouin, np.zeros_like(brillouin)):
                rho[t,:] = frac[:] # not sure if this will work correctly
            else:
                # Actual wave-form calculations needed here.
                R, I = self.sum_ff(brillouin)
                GdotR = np.multiply(brillouin, [1/a, 1/b, 1/c])
                qR = 2 * np.pi * Rsph * np.dot(GdotR, GdotR)**0.5
                rho[t, -1] = const * R * self.formfactor(qR) * np.exp(-self.smear * qR**2 / 2)
                rhoTemp = -rho[t, -1] / np.sum(frac[:-1])
                for i in range(self.nspecies-1):
                    rho[t, i] = rhoTemp * frac[i]
                
                
    def sum_ff(G):
        R = 0
        I = 0
        for i in range(nparticles):
            qR = 2 * np.pi * np.dot(G, self.particles[i,:])
            R = R + np.cos(qR)
            I = I + np.sin(qR)
        
        return R, I
    
    def miller_to_brillouin(G, grid):
        """
        Convert miller indices to brillouin zone (Aliasing)
        
        """
        out = np.zeros_like(G, dtype=int)
        dim = length(grid)
        for i in [2,3]:
            if dim >= i:
                if G[i] > grid[i]/2:
                    out[i] = G[i] - grid[i]
                else:
                    out[i] = G[i]
        return out
    
    
