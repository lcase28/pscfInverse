""" 
Classes to manage Polymer-specific PSO variables in a uniform,
simulator-independent way.
"""

from psoinverse.PSO.SearchSpace import DictPoint, SearchBounds

from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from enum import Enum, unique
import numpy as np
import re

@unique
class VariableTypes(Enum):
    BlockFraction = 1
    BlockLength = 2
    ChainLength = 3
    Chi = 4
    
class MesophaseVariable(ABC):
    """ abstract base for a mesophase variable """
    
    @property
    @abstractmethod
    def psoValue(self):
        pass
    
    @psoValue.setter
    @abstractmethod
    def psoValue(self, val):
        pass
    
    @property
    @abstractmethod
    def scftValue(self):
        pass
    
    @scftValue.setter
    @abstractmethod
    def scftValue(self,val):
        pass
    
    @property
    @abstractmethod
    def psoBounds(self):
        pass
    
    @psoBounds.setter
    @abstractmethod
    def psoBounds(self, val):
        pass
    
    @property
    @abstractmethod
    def scftBounds(self):
        pass
    
    @scftBounds.setter
    @abstractmethod
    def scftBounds(self, val):
        pass
        
    @property
    @abstractmethod
    def keyword(self):
        pass
    
    @abstractmethod
    def keywordMatch(self,key):
        pass
    
    @property
    @abstractmethod
    def flag(self):
        pass
    
class BlockFractionVariable(MesophaseVariable):
    """ A polymer block fraction """
    
    __keywordForm = "BlockFraction[(](?P<PNum>[0-9]+)[)][(](?P<BNum>[0-9]+)[)]"
    
    def __init__(self, polyNum, blockNum, val = 0.5, **kwargs):
        """ 
            Constructor for BlockFractionVariable
            
            Parameters
            ----------
            polyNum : int
                The polymer ID. IDs assumed to start at 0.
            blockNum : int
                The block ID within the indicated polymer.
                IDs are assumed to start at 0.
            val : optional real scalar. Default 0.5.
                The SCFT value of the block fraction.
            lower : keyword, optional, real scalar
                The lower bound of SCFT values. Default 0.
                Must be in range [0,1).
            upper : keyword, optional, real scalar
                The upper bound of SCFT values. Default 1.
                Must be in range (0,1]
        """
        self.PolymerNumber = polyNum
        self.BlockNumber = blockNum
        self.value = val
        
        lower = kwargs.get("lower", 0.0)
        upper = kwargs.get("upper", 1.0)
        bnds = kwargs.get("bounds", None)
        if bnds is None:
            self.lowBnd = lower
            self.hiBnd = upper
        else:
            try:
                lower, upper = bnds
            except(ValueError):
                raise(ValueError("Bounds must contain a lower and upper value"))
            else:
                self.lowBnd = lower
                self.hiBnd = upper
    
    @property
    def polymer(self):
        return self.PolymerNumber
    
    @property
    def block(self):
        return self.BlockNumber
    
    @property
    def psoValue(self):
        return self.value
    
    @psoValue.setter
    def psoValue(self,val):
        self.scftValue = val
    
    @property
    def scftValue(self):
        return self.value
    
    @scftValue.setter
    def scftValue(self,val):
        bnds = self.scftBounds
        lb = bnds[0]
        ub = bnds[1]
        newVal = min(ub, max(lb, val) )
        self.value = newVal
    
    @property
    def psoBounds(self):
        return np.array([self.lowBnd, self.hiBnd])
    
    @psoBounds.setter
    def psoBounds(self, val):
        self.scftBounds = val
    
    @property
    def scftBounds(self):
        return np.array([self.lowBnd, self.hiBnd])
    
    @scftBounds.setter
    def scftBounds(self, val):
        try:
            lower, upper = val
        except ValueError:
            raise(ValueError("Pass an iterable with a lower and upper bound."))
        else:
            self.lowBnd = max(lower, 0.0)
            self.hiBnd = min(upper, 1.0)
        
    @property
    def keyword(self):
        return "BlockFraction(" + str(self.PolymerNumber) + ")(" + str(self.BlockNumber) + ")"
    
    def keywordFormMatch(self, key):
        return re.fullmatch(self.__keywordForm,key) is not None
    
    def keywordMatch(self, key):
        return self.keywordFormMatch(key) and self.keyword == key
    
    @property
    def flag(self):
        return VariableTypes.BlockFraction
        
class ChiVariable(MesophaseVariable):
    """ A Chi*N interaction value """
    
    __keywordForm = "Chi[(](?P<Mon1>[0-9]+)[)][(](?P<Mon2>[0-9]+)[)]"
    
    def __init__(self, Monomer1, Monomer2, val = 0.0, **kwargs):
        """ 
            Constructor for ChiVariable.
            
            Parameters
            ----------
            Monomer1 : int
                The first monomer ID.
                Monomer indexing should start at 1.
            Monomer2 : int
                The second monomer ID.
                Not equal to Monomer1.
                Monomer indexing should start at 1.
            val : optional real scalar. Default 0.0.
                The SCFT value of the block fraction.
            lower : keyword, optional, real scalar
                The lower bound of SCFT values. Default 0.
            upper : keyword, optional, real scalar
                The upper bound of SCFT values. Default 100.
            
            Raises
            ------
            ValueError
                When MonomerIDs < 1, or when Monomer2 == Monomer1
        """
        self.Monomer1 = min(Monomer1,Monomer2)
        self.Monomer2 = max(Monomer1,Monomer2)
        
        # Ensure monomer id numbering starts at 1
        diff = self.Monomer2 - self.Monomer1
        if self.Monomer1 < 1:
            raise(ValueError("Monomer IDs must be >= 1"))
        if diff < 1:
            raise(ValueError("Monomer IDs must differ by at least 1"))
        
        self.value = val
        
        self.lowBnd = kwargs.get("lower", 0.0)
        self.hiBnd = kwargs.get("upper", 100.0)
    
    @property
    def monomerIDs(self):
        """ 
            Monomer ID numbers, with indexing starting at 1 
            Monomer1 < Monomer2
        """
        return self.Monomer1, self.Monomer2
    
    @property
    def psoValue(self):
        return self.value
    
    @psoValue.setter
    def psoValue(self,val):
        self.scftValue = val
    
    @property
    def scftValue(self):
        return self.value
    
    @scftValue.setter
    def scftValue(self,val):
        bnd = self.scftBounds
        lb = bnd[0]
        ub = bnd[1]
        val = min( ub, max( lb, val ) )
        self.value = val
    
    @property
    def psoBounds(self):
        return np.array([self.lowBnd, self.hiBnd])
    
    @psoBounds.setter
    def psoBounds(self, val):
        self.scftBounds = val
    
    @property
    def scftBounds(self):
        return np.array([self.lowBnd, self.hiBnd])
    
    @scftBounds.setter
    def scftBounds(self, val):
        try:
            lower, upper = val
        except ValueError:
            raise(ValueError("Pass an iterable with a lower and upper bound."))
        else:
            self.lowBnd = lower
            self.hiBnd = upper
        
    
    @property
    def keyword(self):
        return "Chi(" + str(self.Monomer1) + ")(" + str(self.Monomer2) + ")"
    
    def keywordFormMatch(self, key):
        return re.fullmatch(self.__keywordForm,key) is not None
    
    def keywordMatch(self,key):
        return self.keywordFormMatch and key == self.keyword
    
    @property
    def flag(self):
        return VariableTypes.Chi

class VariableSet(object):
    """ A collection for Mesophase variables """
    
    def __init__(self, Vars, **kwargs):
        """
            Parameters
            ----------
            Vars : iterable of MesophaseVariable Objects
                The variables in the set
        """
        self.Variables = OrderedDict( [ (v.keyword, v) for v in Vars ] )
    
    @property
    def psoPoint(self):
        """
            Return a PSO dictPoint which would allow PSO
            to be run on the variable set. Returned Point
            will contain current values of all variables.
            
            NOTE: in the returned dictPoint, the fitness
            is defaulted to a value of zero which must be
            updated by the caller.
        """
        kws = []
        coords = []
        for (k, v) in self.Variables.items():
            kws.append(k)
            coords.append(v.psoValue)
        pt = DictPoint( Fitness = 0.0, \
                        Coords = coords, \
                        Scale = 1.0, \
                        keys = kws )
        return pt
        
    @psoPoint.setter
    def psoPoint(self, newpt):
        """ 
            Update the variable values to those stored in newpt
            
            If newpt contains keys not found in VariableSet, 
            these are silently ignored.
            
            If any keys in VariableSet are not in newpt, 
            these are not updated.
            
            If any point values are outside accepted range
            of corresponding variable, they are silently
            rounded to closes bound.
            
            Parameters
            ----------
            newpt : psoinverse.PSO.SearchSpace.DictPoint
                A DictPoint for which the keys match the
                keys of variables stored here.
        """
        pointDict = newpt.get_dict()
        for (k,v) in pointDict.items():
            myVar = self.Variables.get(k,None)
            if myVar is not None:
                myVar.psoValue = v
        
    @property
    def psoBounds(self):
        low = []
        hi = []
        for (k,v) in self.Variables.items():
            bnds = v.psoBounds
            low.append(bnds[0])
            hi.append(bnds[1])
        return SearchBounds(lower=low, upper=hi)
        
    def items(self):
        for (k,c) in deepcopy(self.Variables.items()):
            yield c
        
