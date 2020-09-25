"""
Module containing basic wrapper classes of point and velocity.

Each class is a simple extension of numpy.ndarray, modified to ensure
uniform treatment throughout this program.
"""
from enum import Enum, unique
from functools import total_ordering
import numbers
import numpy as np

__all__ = [Velocity, Point, OptimizationType,FitnessComparator,FITNESS_SELECTOR,set_optimization_type]

class Velocity(object):
    """
    Velocity object to enforce velocity maxima.
    """
    def __init__(self, template, maximum=None):
        self.__dim = np.asarray(template).size
        self.__coord = None
        if maximum is not None:
            self.__cap = np.absolute(np.asarray(maximum).flatten())
        else:
            self.__cap = np.inf * np.ones(self.__dim)
        if not self.__cap.size == self.__dim:
            raise(ValueError("template and maximum must have same size"))
        self.components = template
        
    def reverseComponent(self, index):
        if index >= self.dimensions:
            raise(IndexError("Index {} outside bounds for velocity with dimension {}".format(index,self.dimensions)))
        cd = self.components
        cd[index] *= -1.
        self.components = cd
    
    @property
    def dimensions(self):
        return self.__dim
        
    @property
    def components(self):
        return np.array(self.__coord)
        
    @components.setter
    def components(self, newVal):
        newVal = np.asarray(newVal).flatten()
        if not self.dimensions == newVal.size
            raise(ValueError("Velocity must maintain dimensionality"))
        negFlags = np.sign(newVal)
        newAbs = np.absolute(newVal)
        newCap = np.minimum(newAbs, self.__cap)
        self.__coord = negFlags * newCap
        
@total_ordering
class Point(object):
    """
    Point object containing both a position and a fitness.
    """
    
    def __init__(self, coords, fitness):
        self.__coord = np.asarray(coords).flatten()
        self.__dim = self.__coord.size
        self.__fit = fitness
        
    @property
    def dimensions(self):
        return self.__dim
    
    @property
    def components(self):
        return np.array(self.__coord)
    
    @components.setter
    def components(self, newVals):
        newVal = np.asarray(newVals).flatten()
        if not self.dimensions == newVal.size():
            raise(ValueError("Points must maintain dimensionality"))
        self.__coord = newVal
        
    @property
    def fitness(self):
        return self.__fit
        
    @fitness.setter
    def fitness(self, f):
        self.__fit = f
        
    def __gt__(self,other):
        if not isinstance(other, Point):
            raise(TypeError("Points must be compared to other points"))
        return self.fitness > other.fitness
        
    def __lt__(self,other):
        if not isinstance(other, Point):
            raise(TypeError("Points must be compared to other points"))
        return self.fitness < other.fitness
        
    def __eq__(self,other):
        if not isinstance(other, Point):
            raise(TypeError("Points must be compared to other points"))
        return self.fitness == other.fitness
    
    def __repr__(self):
        return "{}({},{})".format(self.__class__.__name__, self.__coord, self.__fit)
        
    def __str__(self):
        return "{}({},{})".format(self.__class__.__name__, self.__coord, self.__fit)
        
    def same_point(self,other):
        return self == other and np.allclose(self.components, other.components)

@unique
class OptimizationType(Enum):
    minimize = 1
    maximize = 2

class FitnessComparator(object):
    """ Class to enable uniform comparison between points according to selected optimization target """
    
    def __init__(self, optimType=OptimizationType.maximize):
        """
        Parameters
        ----------
        optimType : psovectors.OptimizationType
            The optimization type desired.
        """
        self.__optimType = optimType
    
    @property
    def badFitness(self):
        """ For the given optimization type, return the worst possible fitness value. 
        Any point with this fitness will never be favored as the best point (except in a tie) """
        if self.optimizationType == OptimizationType.maximize:
            return -np.inf
        return np.inf
    
    @property
    def optimizationType(self):
        """ The OptimizationType enforced by the FitnessComparator instance. """
        return self.__opimType
    
    def betterPoint(self, point1, point2):
        """
        Return the point deemed "better" according to the set opimization type.
        Thus, return the point with lower fitness if optimizationType is OptimizationType.minimize.
        
        If the two points are "equal" point1 is returned. No explicit accommodation is made for
        floating point precision. For points with very close fitnesses, the returned point may 
        depend on roundoff errors.
        
        Parameters
        ----------
        point1, point2 : Point
            The points to compare
            
        Returns
        -------
        The better point according to optimization type.
        """
        if self.__optimType is OptimizationType.minimize:
            if point1.fitness <= point2.fitness:
                return point1
            else:
                return point2
        elif self.__optimType is OptimizationType.maximize:
            if point1.fitness >= point2.fitness:
                return point1
            else:
                return point2
        else:
            raise(RuntimeError("Unrecognized optimization type set"))
    
    def bestPoint(self, pointSet):
        """ 
        Return the point from pointSet with the best fitness.
        
        If two or more points tie for "best" point, the first one encountered is returned.
        
        Parameters
        ----------
        pointSet : iterable of Point objects
            The set of points being mutually compared.
        
        Returns
        -------
        The point from pointSet with the best fitness.
        """
        tmp = None
        for (i,p) in enumerate(pointSet):
            if i == 0:
                tmp = p
            else:
                tmp = self.betterPoint(tmp, p)
        return tmp
            
## Shared Variables for Library

FITNESS_SELECTOR = FitnessComparator()

## Methods to Modify Shared States

def set_optimization_type(optim_type):
    """
    Changes core module variables to reflect new optimization type.
    
    Parameters
    ----------
    optim_type : OptimizationType
        The new optimization targeting type.
    
    Raises
    ------
    ValueError
        If optim_type is not a member of OptimizationType
    """
    if not isinstance(optim_type, OptimizationType):
        raise(ValueError("{} is not a member of {}".format(optim_type,OptimizationType)))
    FITNESS_SELECTION = FitnessComparator(optim_type)
