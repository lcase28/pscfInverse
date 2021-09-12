""" 
Classes to manage Polymer-specific PSO variables in a uniform,
simulator-independent way.
"""
# Standard Library Imports
from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from enum import Enum, unique
import itertools
import re

# Third Party Imports
import numpy as np

# Project Imports
import pscfinverse.polymer.parameters as params
from pscfinverse.pso.variables import VariableBase
from pscfinverse.pso.containers import PsoVariableSet

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
    
    @abstractmethod
    def getEquation(self):
        """
        Return a description of the linear equation defined by the variable.
        
        Each variable should be linearizable such that it can be included
        as part of a system of linear equations of the form Ax=b.
        
        Returns
        -------
        param_coeffs : list of PolymerParameter Objects
            All polymer parameter coefficients in the linear equation.
            Each PolymerParameter object will represent a parameter
            in the equation, and the "value" attribute of the parameter
            should be its coefficient in the linear equation, rather than
            a value of the parameter itself (as most variables are likely
            under-specified on their own).
        value : real
            The constant value in the linear equation. The corresponding
            element in b of the linear system.
        """
        pass
    
class BlockLengthVariable(MesophaseVariable):
    """ 
    A variable directly controlling the total length of one or more blocks.
    """
    
    def __init__(self, blocks, value = 0.5, lower = 0.0, upper = 1.0):
        """ 
        Constructor for BlockFractionVariable
        
        If multiple blocks are specified, the variable controls the total
        length of the blocks as a sum. If only one is specified, the
        variable effectively defines its length directly.
        
        Default values of optional parameters are chosen with the assumption
        that the chain the blocks are within is scaled to length 1.
        
        Labels are formatted as: 
                'Length[[_p#b#]...]'
        where '#' stands in for the polymer ('p#') or block ('b#') ID numbers,
        and where the structure '_p#b#' is repeated for each block in the set.
        
        Parameters
        ----------
        blocks : tuple or list of tuple
            The polymer block(s) considered in this variable.
            Each tuple contains two elements: the ID number of
            the polymer, and the ID number of the block in that
            polymer, taking the form (polymerID, blockID).
            For both polymerID and blockID, indexing starts at 0.
        val : real, optional,  Default 0.5.
            An initial value of the block fraction.
        lower : keyword, optional, real scalar
            The lower bound of SCFT values. Default 0.
        upper : keyword, optional, real scalar
            The upper bound of SCFT values. Default 1.
        """
        if type(blocks) == type((1,1)):
            self.__blocks = [blocks]
        elif type(blocks) == type([]):
            self.__blocks = blocks
        else:
            raise(TypeError("Parameter blocks ({}) must be a tuple or list.".format(blocks)))
        
        if value < 0.0:
            raise(ValueError("Block Length must be positive; gave {}".format(value)))
        if lower < 0.0 or upper <= 0.0:
            raise(ValueError("Length bounds must be positive; gave {}, {}".format(lower,upper)))
        
        lbl = "Length"
        for b in self.__blocks:
            lbl = lbl + "_p{}b{}".format(b[0],b[1])
        
        super().__init__(value, lower, upper, label = lbl)
    
    @property
    def polymer(self):
        return self.PolymerNumber
    
    @property
    def block(self):
        return self.BlockNumber
    
    def getEquation(self):
        paramset = []
        for b in self.__blocks:
            polyID, blockID = b
            paramset.append( params.BlockLength( polyID, blockID, 1 ) )
        return paramset, self.trueValue
    
class BlockRatioVariable(MesophaseVariable):
    """
    Describes the ratio of the total lengths of block sets.
    
    Scaling of the variable is logarithmic, such that the pso
    value is ln(Len1/Len2)
    """
    
    def __init__(self, blocks1, blocks2, value = 0.0, lower = -4.5, upper = 4.5):
        """
        Initialize a BlockRatioVariable instance.
        
        BlockRatioVariables Represent ratios on a logarithmic scale.
        If L1 is the total length of blocks1 and Len2 is the total length
        of blocks2, then the "value" of this variable would be
        
            self.value = LN( Len1 / Len2 )
        
        where LN is the natural logarithm. Lower and upper bounds follow
        this same scaling.
        
        Parameters
        ----------
        blocks1 : tuple or list of tuples
            The first (numerator) set of block lengths to sum. Each tuple 
            contains the polymer ID in index 0 and the block ID
            in index 1. Indexing starts at 0.
        blocks2 : tuple or list of tuples
            The second (denominator) set of block lengths to sum.
            Tuples structurd as described for blocks1
        value : real
            An initial value for the 
        """
        if type(blocks1) == type((1,1)):
            self._blocks1 = [blocks1]
        elif type(blocks1) == type([]):
            self._blocks1 = blocks1
        else:
            raise(TypeError("Block sets must be lists or tuples."))
        
        if type(blocks2) == type((1,1)):
            self._blocks2 = [blocks2]
        elif type(blocks2) == type([]):
            self._blocks2 = blocks2
        else:
            raise(TypeError("Block sets must be lists or tuples."))
        
        lbl = "LengthRatio"
        for b in self._blocks1:
            lbl = lbl + "_p{}b{}".format(b[0],b[1])
        lbl += "_r"
        for b in self._blocks2:
            lbl = lbl + "_p{}b{}".format(b[0],b[1])
        
        super().__init__(value, lower, upper, label = lbl)
    
    def getEquation(self):
        paramset = []
        for b in self._blocks1:
            polyID, blockID = b
            paramset.append( params.BlockLength( polyID, blockID, 1 ) )
        for b in self._blocks2:
            polyID, blockID = b
            paramset.append( params.BlockLength( polyID, blockID, -np.exp( self.trueValue ) ) )
        return paramset, 0.0

class ChiVariable(MesophaseVariable):
    """ A Chi*N interaction value """
    
    def __init__(self, Monomer1, Monomer2, value = 0.0, lower=0.0, upper=100.0):
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
        value : optional real scalar. Default 0.0.
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
        
        lbl = "ChiN_m{}_m{}".format(self.monomer1, self.monomer2)
        
        super().__init__(value,lower,upper,label = lbl)
    
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
    
    def getEquation(self):
        return [params.ChiN(self.monomer1, self.monomer2, 1)], self.trueValue

class KuhnVariable(MesophaseVariable):
    """
    A Kuhn length variable.
    
    This variable can act as a stand-in for any flexibility-measuring value
    such as statistical segment lengths or persistence lengths which may be 
    required by a particular SCFT software.
    """
    
    def __init__(self, monomer, value = 1, lower = 0.2, upper = 5.0):
        self._mon = monomer
        lbl = "Kuhn_m{}".format(self._mon)
        super().__init__(value, lower, upper, label = lbl)
    
    def getEquation(self):
        return [params.KuhnLength(self._mon, 1.0)], self.trueValue

class PolymerVariableSet(PsoVariableSet):
    """ A collection for polymeric variables """
    
    def __init__(self, Vars, constraints):
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
        self._next_param = itertools.count(start=1)
        self._next_equation = itertools.count(start=1)
        self._param_count = 0
        self._equation_count = 0
        self._param_map = {}
        self._param_set = {}
        self._A_matr = np.zeros((0,0))
        self._b_vect = np.zeros(0)
        self._x_vect = np.zeros(0)
        # Evaluate System will allocate all data required for updates later.
        self.constraints = constraints
        super().__init__(Vars)
        self._evaluate_system()
    
    def _add_parameter(self,key,param):
        self._param_count = next(self._next_param)
        newIndex = self._param_count - 1
        self._param_map.update({key:newIndex})
        self._param_set.update({key:deepcopy(param)})
        newcol = np.zeros((self._equation_count,1))
        self._A_matr = np.concatenate((self._A_matr, newcol), axis=1)
        self._x_vect.resize(self._param_count)
    
    def _add_equation(self):
        self._equation_count = next(self._next_equation)
        newrow = np.zeros((1,self._param_count))
        self._A_matr = np.concatenate((self._A_matr, newrow), axis=0)
        self._b_vect = np.append(self._b_vect, 0.0)
    
    def _build_equation_system(self,*args):
        """
        Build the current system of linear equation.
        
        Each input value should be an iterable set of
        PolymerVariable objects.
        """
        equationindex = 0
        for varlist in args:
            for var in varlist:
                # First time through, this will build up the equation list
                if (equationindex+1) > self._equation_count:
                    self._add_equation()
                param_coefs, value = var.getEquation() 
                for p in param_coefs:
                    # first time through, this will build the parameter map
                    key = p.key
                    if key not in self._param_map:
                        self._add_parameter(key,p)
                    # look up the column index for the parameter
                    paramindex = self._param_map.get(key)
                    self._A_matr[equationindex][paramindex] = p.value
                self._b_vect[equationindex] = value
                equationindex += 1
        
    def _evaluate_system(self):
        """
        Use the current values of the polymer variables to evaluate parameter values.
        
        All variables defined for Polymer systems can be taken as part of a system
        of linear equations of the form Ax=b. In this method, coefficients in A and
        values in b are collected from the variables in the set to determine
        """
        # Build system of equations
        self._build_equation_system(self.variables, self.constraints)
        # Solve Linear equations
        Ainv = np.linalg.inv(self._A_matr)
        self._x_vect = np.matmul( Ainv, self._b_vect, out=self._x_vect )
        
        # Collect Parameter values
        for (key,ind) in self._param_map.items():
            param = self._param_set.get(key)
            param.value = self._x_vect[ind]
    
    def startUpdate(self, *args, **kwargs):
        out = super().startUpdate(*args,**kwargs)
        self._evaluate_system()
        return out
    
    def cancelUpdate(self, *args, **kwargs):
        super().cancelUpdate(*args, **kwargs)
        self._evaluate_system()
    
    @property
    def parameters(self):
        out = []
        for (key,param) in self._param_set.items():
            out.append(param)
        return out
    
    
