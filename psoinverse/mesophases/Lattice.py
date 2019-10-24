# Lattice Module. 
# Partially derived from crystals package available at:
#   https://github.com/LaurentRDC/crystals.git
# Generallized to handle 1-, 2-, or 3-dimensional systems

# Imports
from abc import ABC, abstractmethod
from copy import deepcopy
from enum import Enum, unique
import numpy as np

@unique
class LatticeSystem3D(Enum):
    """
    Enumeration of basic 3D lattice systems (not Bravais)
    """
    cubic = 1
    hexagonal = 2
    rhombohedral = 3
    tetragonal = 4
    orthorrhombic = 5
    monoclinic = 6
    triclinic = 7

    
class Lattice(Lattice):
    """ Object representing a crystallographic basis vector lattice"""
    
    ## TODO: Scrub inputs for value/type constraints
    def __init__(self, dim, basis, **kwargs):
        """
        Generate a lattice object.
        
        Params:
        -------
        dim: int, dim >= 2
            Number of dimensions in the lattice.
        basis: numpy.ndarray
            Matrix (size: dim-by-dim) representation of basis vectors.
            Coordinates of each vector are defined relative to the
            'dim'-dimensional standard basis
        """
        self.dim = dim;
        self.basis = basis;
        self.reciprocal
    
    @classmethod
    def latticeFromParameters(cls, dim, **kwargs):
        """
        Generate a lattice object for 2-D or 3-D lattice.
        
        Params:
        -------
        dim: int, in set {2, 3}
            The dimensionality of the resulting lattice
        
        Keyword Params:
        ---------------
        a:  float, required
            Magnitude of first basis vector.
        b:  float, required
            Magnitude of second basis vector.
        c:  float, (only if dim == 3)
            Magnitude of third basis vector.
        alpha:  float, only if dim == 3
            Angle (in degrees, range (0, 180) ) between vector b and c.
        beta:   float, only if dim == 3
            Angle (in degrees, range (0, 180) ) between vector a and c.
        gamma:  float, required
            Angle (in degrees, range (0, 180) ) between vectors a and b.
        
        Returns:
        --------
        lat: Lattice
            Lattice defined by given parameters.
            
        Raises:
        -------
        TypeError: If input arguments are non-numeric.
        ValueError: If dim not one of {2, 3}, or invalid lattice parameters.
        """
        basis = cls.basisFromParameters(dim, **kwargs)
        return cls(dim, basis)
            
    @classmethod
    def basisFromParameters(cls, dim, **kwargs):
        """
        Return a set of basis vectors.
        
        Params:
        -------
        dim: int, in set {2, 3}
            The dimensionality of the resulting lattice
        
        Keyword Params:
        ---------------
        a:  float, required
            Magnitude of first basis vector.
        b:  float, required
            Magnitude of second basis vector.
        c:  float, (only if dim == 3)
            Magnitude of third basis vector.
        alpha:  float, only if dim == 3
            Angle (in degrees, range (0, 180) ) between vector b and c.
        beta:   float, only if dim == 3
            Angle (in degrees, range (0, 180) ) between vector a and c.
        gamma:  float, required
            Angle (in degrees, range (0, 180) ) between vectors a and b.
        
        Returns:
        --------
        basis: numpy.ndarray, 'dim'-by-'dim'
            Lattice defined by given parameters.
            
        Raises:
        -------
        TypeError: If input arguments are non-numeric.
        ValueError: If dim not one of {2, 3}, or invalid lattice parameters.
        """
        # extract all parameters
        a = kwargs.get("a",None)
        b = kwargs.get("b",None)
        c = kwargs.get("c",None)
        alpha = kwargs.get("alpha",None)
        beta = kwargs.get("beta",None)
        gamma = kwargs.get("gamma",None)
        
        # check inputs
        if (a is None) or (b is None) or (gamma is None):
            raise TypeError("Required lattice parameter is missing.")
            
        if dim == 3 and ((c is None) or (alpha is None) or (beta is None)):
            raise TypeError("Missing lattice parameter for 3D lattice")
        
        # initialize basis
        basis = np.zeros((dim,dim))
        
        # Complete common calculations
        gammaRad = deg2rad(gamma)
        basis[0,0] = a
        basis[1,0] = b*np.cos(gammaRad)
        basis[1,1] = b*np.sin(gammaRad)
        
        # Additional 3D calculations
        if dim == 3:
            alphaRad = np.deg2rad(alpha)
            betaRad = np.deg2rad(beta)
            basis[2,0] = c*np.cos(betaRad)
            basis[2,1] = c*np.cos(alphaRad)*sin(gammaRad)
            basis[2,2] = np.sqrt( c**2 - basis[2,0]**2 - basis[2,1]**2)
        return basis
    
    ## Properties
    @property
    def latticeParameters(self):
        """
            Lattice parameters as list of lengths and angles.
            3D: [a, b, c, alpha, beta, gamma]
            2D: [a, b, gamma]
        """
        a = self.basis[0,:]
        b = self.basis[1,:]
        aMag = np.linalg.norm(a)
        bMag = np.linalg.norm(b)
        gamma = np.rad2deg( np.arccos( np.dot(a,b) / (aMag*bMag) ) )
        if self.dim == 2:
            return [aMag, bMag, gamma]
        elif self.dim == 3:
            c = self.basis[2,:]
            cMag = np.linalg.norm(c)
            alpha = np.rad2deg( np.arccos( np.dot(b,c) / (bMag*cMag) ) )
            beta = np.rad2deg( np.arccos( np.dot(a,c) / (aMag*cMag) ) )
            return [aMag, bMag, cMag, alpha, beta, gamma]
        else
            return NotImplemented
    
    @latticeParameters.setter
    def latticeParameters(self, **newParams):
        """ Update the lattice to have the given parameters """
        self.basis = self.__class__.basisFromParameters(self.dim, **newParams)
    
    @property
    def latticeVectors(self):
        """ Return lattice vectors as numpy.ndarray """
        return deepcopy(self.basis)
    
    @latticeVectors.setter
    def latticeVectors(self, newBasis):
        """
            Update the basis vectors to match newBasis
            
            Params:
            -------
            newBasis: numpy.ndarray, self.dim square
            
            Raises:
            -------
            ValueError: If self.basis.size != newBasis.size
        """
        if self.basis.size != newBasis.size:
            raise ValueError("New Basis does not have same dimension as lattice")
        self.basis = newBasis
    
    @property
    def volume(self):
        """ Area (2D) or Volume (3D) enclosed by basis vectors """
        return np.linalg.det(self.basis)
    
    @property
    def metricTensor(self):
        """ The real-space metric tensor """
        return np.matmul(self.basis, self.basis.T)
    
    @property
    def reciprocal(self):
        """ Lattice object for reciprocal lattice. """
        return self.__class__.(self.dim, self.reciprocalBasis)
    
    @property
    def reciprocalVectors(self):
        """ Basis vectors of reciprocal lattice as numpy.ndarray,
            with rows corresponding to each vector on standard basis """
        reciprocalMetricTensor = np.linalg.inv(self.metricTensor)
        return np.matmul(reciprocalMetricTensor, self.basis)
    
    ## Instance Methods
    def changeFromBasis(self, vect, reference):
        """
            Return vector coordinates relative to THIS lattice basis.
            
            Parameters
            ----------
            vect : float, 1-D list-like
                The vector, in fractional coordinates relative to
                old basis.
            reference : Lattice
                The lattice whose vectors represent current basis
                of the given vector coordinates.
            
            Returns
            -------
            newVect : float, numpy.ndarray
                Coordinates of the vector relative to basis vectors
                of THIS lattice.
        """
        v = np.reshape( np.asarray(vect, dtype=np.float64), (self.dim, 1) )
        COB = np.matmul(self.basis, reference.basis.T)
        return np.matmul(COB, v)
    
    def changeToBasis(self, vect, target):
        """
            Return vector coordinates relative to target lattice basis.
            
            Parameters
            ----------
            vect : float, 1-D list-like
                The vector, in fractional coordinates relative to
                THIS basis.
            target : Lattice
                The lattice whose vectors represent the basis in which
                to cast the vector.
            
            Returns
            -------
            newVect : float, numpy.ndarray
                Coordinates of the vector relative to basis vectors
                of target lattice.
        """
        return target.changeFromBasis(vect, self)
        
    
    ## "Private" internal methods
    def __repr__(self):
        s = "< Lattice object with parameters {:.3f}A, {:.3f}, {:.3f}, {:.3f}, {:.3f} >"
        return s.format(*self.lattice_parameters)
    
    def __hash__(self):
        return hash(self.lattice_parameters) | super().hash()
        
    def __array__(self, *args, **kwargs):
        """ returns a 3x3 float array. Each row is a lattice vector. """
        return np.array(self.lattice_vectors, *args, *kwargs)
    
    
#class Lattice(ABC):
#    """ Abstract base class of (2+)-dimensional lattices """
#    
#    def __init__(self, **kwargs):
#        super().__init__(kwargs)
#    
#    ## Properties
#    @property
#    @abstractmethod
#    def latticeParameters(self):
#        """
#        Implementing class should return FULL set of lattice parameters.
#        (i.e. for any 3D system, 6 parameters should be returned, even if
#        fewer are needed to describe the lattice based on system defintion)
#        
#        Returns:
#        --------
#        LatticeParameters: list-like
#            The parameters of the lattice, in conventional order.
#        """
#        return NotImplemented
#    
#    @latticeParameters.setter
#    @abstractmethod
#    def latticeParameters(self, params):
#        
#    
#    ## Class Methods
#    @classmethod
#    @abstractmethod
#    def latticeFromParameters(cls, *args, **kwargs):
#        return NotImplemented
#    
#    ## "Private" Methods
#    @abstractmethod
#    def __repr__(self):
#        return NotImplemented
#    
#    @abstractmethod
#    def __hash__(self):
#        return super().hash()
#    
#    def __eq__(self, other):
#        if isinstance(other, self.__class__):
#            reutrn np.allclose(self.lattice_vectors, other.lattice_vectors, atol=1e-4)
#        return NotImplemented
#    
#    @abstractmethod
#    def __array__(self, *args, **kwargs):
#        pass
    
