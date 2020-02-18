"""
Module containing basic wrapper classes of point and velocity.

Each class is a simple extension of numpy.ndarray, modified to ensure
uniform treatment throughout this program.
"""
from functools import total_ordering
import numbers
import numpy as np

class Velocity(PsoVector):
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
class Point(PsoVector):
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
    
