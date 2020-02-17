"""
Module defining pso variable sets and points
"""

# Library imports
from psoinverse.PSO.psovariables import PsoVariable
from psoinverse.PSO.psovectors import Point

# Third-party imports
import numpy as np


class PsoVariableSet(object):
    """
    Class defining basic variable set management
    """
    
    def __init__(self, variables):
        """
        Parameters
        ----------
        variables : iterable of PsoVariable (or sub-class) instances
            The variables being grouped in the set.
        """
        self.__nvar = 0
        self.__variables = []
        for v in variables:
            self.__nvar += 1
            self.__variables.append(v)
        
    def checkBounds(self, testposition):
        val = np.asarray(testposition).flatten()
        if not val.size == self.__nvar:
            raise(ValueError("Checking bounds on position of invalid size"))
        res = []
        for i in range(self.__nvar):
            res.append( self.__variables[i].checkBounds(val[i]) )
        return res
    
    def checkAcceptance(self, testposition):
        val = np.asarray(testposition).flatten()
        if not val.size == self.__nvar:
            raise(ValueError("Checking bounds on position of invalid size"))
        res = []
        for i in range(self.__nvar):
            res.append( self.__variables[i].checkAcceptance(val[i]) )
        return res
        
    ## TODO: Write TryUpdate and SetUpdate methods (possibly separate for pso vs true values)
        
    @property
    def psoValues(self):
        return np.asarray( [v.psoValue for v in self.__variables] )
    
    @psoValues.setter
    def psoValues(self, newVals):
        newVals = np.asarray(newVals).flatten()
        if not newVals.size == self.__nvar:
            raise(ValueError("New value set for psoValues is improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.psoValue = newVals[i]
    
    @property
    def trueValues(self):
        return np.asarray( [v.trueValue for v in self.__variables])
    
    @trueValues.setter
    def trueValues(self, newVals):
        newVals = np.asarray(newVals).flatten()
        if not newVals.size == self.__nvar:
            raise(ValueError("New value set for trueValues is improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.trueValue = newVals[i]
    
    @property
    def psoLowerBounds(self):
        return np.asarray( [v.psoLowerBound for v in self.__variables])
    
    @psoLowerBounds.setter
    def psoLowerBounds(self, newVals):
        newVals = np.asarray(newVals).flatten()
        if not newVals.size == self.__nvar:
            raise(ValueError("New value set for psoLowerBounds is improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.psoLowerBound = newVals[i]
    
    @property
    def psoUpperBounds(self):
        return np.asarray( [v.psoUpperBound for v in self.__variables])
    
    @psoUpperBounds.setter
    def psoUpperBounds(self, newVals):
        newVals = np.asarray(newVals).flatten()
        if not newVals.size == self.__nvar:
            raise(ValueError("New value set for psoUpperBounds is improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.psoUpperBound = newVals[i]
    
    @property
    def psoBounds(self):
        return [self.psoLowerBounds, self.psoUpperBounds]
    
    @psoBounds.setter
    def psoBounds(self, newbnds):
        lowbnds = np.asarray(newbnds[0]).flatten()
        upbnds = np.asarray(newbnds[1]).flatten()
        if not lowbnds.size == upbnds.size:
            raise(ValueError("Dimensionality of bounds must agree"))
        if not lowbnds.size == self.__nvar:
            raise(ValueError("New value sets for psoBounds are improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.psoBounds = [lowbnds[i], upbnds[i]]
    
    @property
    def trueLowerBounds(self):
        return np.asarray( [v.trueLowerBound for v in self.__variables])
    
    @trueLowerBounds.setter
    def trueLowerBounds(self, newVals):
        newVals = np.asarray(newVals).flatten()
        if not newVals.size == self.__nvar:
            raise(ValueError("New value set for trueLowerBounds is improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.trueLowerBound = newVals[i]
            
    @property
    def trueUpperBounds(self):
        return np.asarray( [v.trueUpperBound for v in self.__variables])
    
    @trueUpperBounds.setter
    def trueUpperBounds(self, newVals):
        newVals = np.asarray(newVals).flatten()
        if not newVals.size == self.__nvar:
            raise(ValueError("New value set for trueUpperBounds is improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.trueUpperBound = newVals[i]
            
    @property
    def trueBounds(self):
        return [self.trueLowerBounds, self.trueUpperBounds]
        
    @trueBounds.setter
    def trueBounds(self, newbnds):
        lowbnds = np.asarray(newbnds[0]).flatten()
        upbnds = np.asarray(newbnds[1]).flatten()
        if not lowbnds.size == upbnds.size:
            raise(ValueError("Dimensionality of bounds must agree"))
        if not lowbnds.size == self.__nvar:
            raise(ValueError("New value sets for trueBounds are improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.trueBounds = [lowbnds[i], upbnds[i]]
            
    
        
