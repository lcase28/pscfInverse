""" 
Classes to manage Polymer-specific PSO variables in a uniform,
simulator-independent way.
"""

from abc import ABC, abstractmethod
import enum import Enum, unique
import re

@unique
class VariableTypes(Enum):
    
    BlockFraction = 1
    BlockLength = 2
    ChainLength = 3
    Chi = 4
    

class MesophaseVariable(ABC):
    """ abstract base for a mesophase variable """
    
    @abstractmethod
    @property
    def psoValue(self):
        """ Value, applying any selected scaling/projection for PSO calculation """
        pass
    
    @abstractmethod
    @psoValue.setter
    def psoValue(self, val):
        pass
    
    @abstractmethod
    @property
    def scftValue(self):
        pass
    
    @abstractmethod
    @scftValue.setter
    def scftValue(self,val):
        pass
        
    @abstractmethod
    @property
    def keyword(self):
        pass
    
    @abstractmethod
    def keywordMatch(self,key):
        pass
    
    @abstractmethod
    @property
    def flag(self):
        pass
    
class BlockFractionVariable(MesophaseVariable):
    """ A polymer block fraction """
    
    __keywordForm = "BlockFraction[(](?P<PNum>[0-9]+)[)][(](?P<BNum>[0-9]+)[)]"
    
    def __init__(self, polyNum, blockNum, val = 0.0):
        self.PolymerNumber = polyNum
        self.BlockNumber = blockNum
        self.value = val
    
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
        self.value = val
    
    @property
    def scftValue(self):
        return self.value
    
    @scftValue.setter
    def scftValue(self):
        self.value = val
    
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
    
    def __init__(self, Monomer1, Monomer2, val = 0.0):
        self.Monomer1 = Monomer1
        self.Monomer2 = Monomer2
        self.value = val
    
    @property
    def monomerIDs(self):
        return self.Monomer1, self.Monomer2
    
    @property
    def psoValue(self):
        return self.value
    
    @psoValue.setter
    def psoValue(self,val):
        self.value = val
    
    @property
    def scftValue(self):
        return self.value
    
    @scftValue.setter
    def scftValue(self):
        self.value = val
    
    @property
    def keyword(self):
        return "Chi(" + str(self.Monomer1) + ")(" + str(self.Monomer2) + ")"
    
    def keywordFormMatch(self, key):
        return re.fullmatch(self.__keywordForm,key) is not None
    
    def keywordMatch(self,key):
        return self.keywordFormMatch and key == self.keyword
    
    @property
    def flag(self):
        return VariableTypes.BlockFraction
        


