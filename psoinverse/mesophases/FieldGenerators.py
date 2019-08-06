# Imports
from abc import ABC, abstractmethod
from crystals import Lattice, LatticeSystem
import numpy as np
from psoinverse.util.stringTools import str_to_num
import itertools

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
    
    def __init__(self, **kwargs):
        self.lattice = kwargs.get("lattice", Lattice.from_parameters(1,1,1,90,90,90))
        self.particles = kwargs.get("particlePositions",None)
        self.nparticles = kwargs.get("N_particles")
        self.nspecies = kwargs.get("N_monomer")
        self.smear = kwargs.get("sigma_smear")
        self.dim = kwargs.get("dim")
        self.crystal_system = kwargs.get("crystal_system")
        self.nCellParam = kwargs.get("N_cell_param")
        self.cellParam = kwargs.get("cell_param")
        self.group_name = kwargs.get("group_name")
        self.ngrid = kwargs.get("ngrid")
        self.outfile = kwargs.get("output_filename", "rho_kgrid")
        self.formfactor = kwargs.get("formfactor", sphere_form_factor)
        super().__init__()
    
    @classmethod
    def from_file(cls, fname):
        with open(fname) as f:
            kwargs = {}
            for line in f.readlines():
                line = line.strip()
                splitline = line.split(maxsplit=1)
                if len(splitline) > 1:
                    key = splitline[0].strip()
                    if key != "group_name":
                        data = splitline[1].strip().split()
                        data = data[0] if len(data) == 1 else data
                    else:
                        data = splitline[1].strip()
                    #interpret key for data parsing
                    print(key)
                    print(data)
                    print("\n")
                    if key == "lattice_const":
                        kwargs.update(lattice=Lattice.from_parameters(*[float(d) for d in data]))
                    elif key == "particlePositions":
                        pos = np.array([float(d) for d in data])
                        pos = np.reshape(pos, (-1, 3))
                        kwargs.update([(key,pos)])
                    else:
                        if len(data) > 1 and not isinstance(data, str):
                            try:
                                data = np.array([str_to_num(d) for d in data])
                            except(ValueError, TypeError):
                                pass
                        else:
                            try:
                                data = str_to_num(data)
                            except(ValueError,TypeError):
                                pass
                        kwargs.update([(key, data)])
        return cls(**kwargs)
    
    # TODO: Figure out how to generate 2D, 1D initial guesses
    def to_kgrid(self, frac):
        coreindex = 0
        ngrid = self.ngrid
        print("ngrid = ", ngrid)
        # Shift grid for k-grid
        kgrid = np.zeros_like(ngrid)
        for (i,x) in enumerate(ngrid):
            if i == 0:
                kgrid[i] = x/2
            else:
                kgrid[i] = x - 1
        print("kgrid = ",kgrid)
            
        nwaves = str_to_num(np.prod(kgrid + np.ones_like(kgrid)))
        print("nwaves = ", nwaves)
        rho = np.zeros((nwaves, self.nspecies))
        print("rho init size: ", rho.size)
        vol = self.lattice.volume
        print("volume = ", vol)
        const = frac[coreindex] / self.nparticles
        print("frac / nparticles = ", const)
        Rsph = ((3 * frac[coreindex] * vol) / (4 * np.pi * self.nparticles))**(1./3)
        print("Rsph = ", Rsph)
        a, b, c, alpha, beta, gamma = self.lattice.lattice_parameters
        print("Lattice params: ", a, b, c, alpha, beta, gamma)
        
        # primary loop for n-dimensional generation
        t = -1
        for G in itertools.product(*[range(x+1) for x in kgrid]):
            # fill in loop operations - G is miller indices in n-dimensions.
            G = np.array(G) #convert tuple to array
            brillouin = self.miller_to_brillouin(G, ngrid)
            t = t + 1
            print("\nt = ", t)
            print("miller: ", G)
            print("brillouin: ", brillouin)
            if np.array_equiv(brillouin, np.zeros_like(brillouin)):
                rho[t,:] = frac[:] # not sure if this will work correctly
            else:
                # Actual wave-form calculations needed here.
                R, I = self.sum_ff(brillouin)
                print("R = ", R)
                GdotR = np.multiply(brillouin, [1/a, 1/b, 1/c])
                print("GdotR = ", GdotR)
                qR = 2 * np.pi * Rsph * np.dot(GdotR, GdotR)**0.5
                print("qR = ", qR)
                rho[t, coreindex] = const * R * self.formfactor(qR) * np.exp(-(self.smear**2 * qR**2 / 2))
                rhoTemp = -rho[t, coreindex] / np.sum(frac[1:])
                for i in range(self.nspecies-1):
                    rho[t, i+1] = rhoTemp * frac[i+1]
                for (i,r) in enumerate(rho[t,:]):
                    if r == -0.0:
                        rho[t,i] = 0.0
            print("rho[{}] = {}".format(t, rho[t,:]))
        
        print("rho final size: ", rho.size)
        return rho
                
    def to_file(self, frac, fileRoot=None):
        if fileRoot is not None:
            fileRoot = fileRoot.strip().rstrip("/")
            fname = fileRoot + "/" + self.outfile.strip().strip("'")
        else:
            fname = self.outfile.strip().strip("'")
        with open(fname, "w") as f:
            f.write("    format 1 0")
            f.write("\ndim\n\t{}".format(self.dim))
            f.write("\ncrystal_system\n\t{}".format(self.crystal_system))
            f.write("\nN_cell_param\n\t{}".format(self.nCellParam))
            
            cellParamString = "\ncell_param\n"
            for i in range(self.nCellParam):
                cellParamString += "\t{}"
            try:
                f.write(cellParamString.format(*self.cellParam))
            except(TypeError):
                f.write(cellParamString.format(self.cellParam))
            
            f.write("\ngroup_name\n\t{}".format(self.group_name))
            f.write("\nN_monomer\n\t{}".format(self.nspecies))
            
            ngridString = "\nngrid\n"
            for i in range(self.dim):
                ngridString += "\t{}"
            try:
                f.write(ngridString.format(*self.ngrid))
            except(TypeError):
                f.write(ngridString.format(self.ngrid))
            
            rho = self.to_kgrid(frac)
            baseString = "({:.4E},0.0000E+00)"
            rowString = "\n" + baseString
            for i in range(self.nspecies-1):
                rowString += "  " + baseString
            for i in range(len(rho[:,0])):
                f.write(rowString.format(*rho[i,:]))
        return 0
        
    def sum_ff(self, G):
        R = 0
        I = 0
        for i in range(self.nparticles):
            qR = 2 * np.pi * np.dot(G, self.particles[i,:])
            R = R + np.cos(qR)
            I = I + np.sin(qR)
        
        return R, I
    
    def miller_to_brillouin(self, G, grid):
        """
        Convert miller indices to brillouin zone (Aliasing)
        
        """
        out = np.zeros_like(G, dtype=int)
        out[0] = G[0]
        dim = self.dim-1
        for i in [1,2]:
            if dim >= i:
                if G[i] > grid[i]/2:
                    out[i] = G[i] - grid[i]
                else:
                    out[i] = G[i]
        return out
    
    
