""" 
Classes to manage Polymer-specific PSO variables in a uniform,
simulator-independent way.
"""
# Standard Library Imports
from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from enum import Enum, unique
import re

# Third Party Imports
import numpy as np

# Project Imports
import psoinverse.polymer.parameters as params
from psoinverse.pso.variables import VariableBase
from psoinverse.pso.containers import PsoVariableSet

@unique
class VariableTypes(Enum):
    BlockFraction = 1
    BlockLength = 2
    ChainLength = 3
    Chi = 4
    
class MesophaseVariable(VariableBase):
    """ abstract base for a mesophase variable """
    
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    
    def _pso_to_true(self,val):
        return val
    
    def _true_to_pso(self,val):
        return val
    
    @property
    @abstractmethod
    def parameter(self):
        """
        The polymer parameter (or list of parameters) set by the variable.
        """
        pass
    
class BlockFractionVariable(MesophaseVariable):
    """ A polymer block fraction """
    
    def __init__(self, polyNum, blockNum, value = 0.5, lower = 0.0, upper = 1.0):
        """ 
        Constructor for BlockFractionVariable
        
        Parameters
        ----------
        polyNum : int
            The polymer ID. IDs assumed to start at 0.
        blockNum : int
            The block ID within the indicated polymer.
            IDs are assumed to start at 0.
        val : real, optional,  Default 0.5.
            An initial value of the block fraction.
        lower : keyword, optional, real scalar
            The lower bound of SCFT values. Default 0.
            Must be in range [0,1).
        upper : keyword, optional, real scalar
            The upper bound of SCFT values. Default 1.
            Must be in range (0,1]
        """
        self.__polymer_number = polyNum
        self.__block_number = blockNum
        
        if value < 0.0 or value > 1.0:
            raise(ValueError("Block Fraction must be on [0.0,1.0]; gave {}".format(value)))
        if lower < 0.0 or lower >= 1.0:
            raise(ValueError("Fraction lower bound must be on [0.0, 1.0); gave {}".format(lower)))
        if upper <= 0.0 or upper > 1.0:
            raise(ValueError("Fraction upper bound must be on (0.0, 1.0]; gave {}".format(upper)))
        
        lbl = "bFrac({},{})".format(polyNum,blockNum)
        
        super().__init__(value, lower, upper, label = lbl)
    
    @property
    def polymer(self):
        return self.PolymerNumber
    
    @property
    def block(self):
        return self.BlockNumber
    
    @property
    def parameter(self):
        return params.BlockLength(self.polymer,self.block,self.trueValue)
    
class DiblockFractionVariable(MesophaseVariable):
    """ A diblock copolymer composition """
    
    def __init__(self, polyNum, value = 0.5, lower = 0.0, upper = 1.0):
        """ 
        Constructor for BlockFractionVariable
        
        Parameters
        ----------
        polyNum : int
            The polymer ID. IDs assumed to start at 0.
        val : real, optional,  Default 0.5.
            An initial value of the block fraction.
        lower : keyword, optional, real scalar
            The lower bound of SCFT values. Default 0.
            Must be in range [0,1).
        upper : keyword, optional, real scalar
            The upper bound of SCFT values. Default 1.
            Must be in range (0,1]
        """
        self.__polymer_number = polyNum
        
        if value < 0.0 or value > 1.0:
            raise(ValueError("Block Fraction must be on [0.0,1.0]; gave {}".format(value)))
        if lower < 0.0 or lower >= 1.0:
            raise(ValueError("Fraction lower bound must be on [0.0, 1.0); gave {}".format(lower)))
        if upper <= 0.0 or upper > 1.0:
            raise(ValueError("Fraction upper bound must be on (0.0, 1.0]; gave {}".format(upper)))
        
        lbl = "dbFrac({})".format(polyNum)
        
        super().__init__(value, lower, upper, label = lbl)
    
    @property
    def polymer(self):
        return self.__polymer_number
    
    @property
    def parameter(self):
        out = []
        out.append(params.BlockLength(self.polymer,0,self.trueValue))
        out.append(params.BlockLength(self.polymer,1,1.0-self.trueValue))
        return out
    
class ChiVariable(MesophaseVariable):
    """ A Chi*N interaction value """
    
    def __init__(self, Monomer1, Monomer2, val = 0.0, lower=0.0, upper=100.0):
        """ 
        Constructor for ChiVariable.
        
        Parameters
        ----------
        Monomer1 : int
            The first monomer ID.
            Monomer indexing should start at 0.
        Monomer2 : int
            The second monomer ID.
            Not equal to Monomer1.
            Monomer indexing should start at 0.
        val : optional real scalar. Default 0.0.
            The SCFT value of the block fraction.
        lower : keyword, optional, real scalar
            The lower bound of SCFT values. Default 0.
        upper : keyword, optional, real scalar
            The upper bound of SCFT values. Default 100.
        
        Raises
        ------
        ValueError
            When MonomerIDs < 0, or when Monomer2 == Monomer1
        """
        self.__monomer1 = min(Monomer1,Monomer2)
        self.__monomer2 = max(Monomer1,Monomer2)
        
        # Ensure monomer id numbering starts at 1
        diff = self.monomer2 - self.monomer1
        if self.monomer1 < 0:
            raise(ValueError("Monomer IDs must be >= 0"))
        if diff < 1:
            raise(ValueError("Monomer IDs must differ by at least 1"))
        
        lbl = "ChiN({},{})".format(self.monomer1, self.monomer2)
        
        super().__init__(val,lower,upper,label = lbl)
    
    @property
    def monomer1(self):
        return self.__monomer1
    
    @property
    def monomer2(self):
        return self.__monomer2
    
    @property
    def monomerIDs(self):
        """ 
        Monomer ID numbers, with indexing starting at 0
        Monomer1 < Monomer2
        """
        return self.monomer1, self.monomer2
    
    @property
    def parameter(self):
        return params.ChiN(self.monomer1, self.monomer2, self.trueValue)

class PolymerVariableSet(PsoVariableSet):
    """ A collection for polymeric variables """
    
    def __init__(self, Vars):
        """
        Parameters
        ----------
        Vars : iterable of MesophaseVariable Objects
            The variables in the set
        """
        for v in Vars:
            if not isinstance(v,MesophaseVariable):
                msg = "Variable {} of type {} is not allowed in PolymerVariableSet."
                raise(TypeError(msg.format(v,type(v))))
        super().__init__(Vars)
    
    @property
    def parameters(self):
        out = []
        for v in self.variables:
            p = v.parameter
            if type(p) == type([]):
                for t in p:
                    out.append(t)
            else:
                out.append(p)
        return out
    
    
