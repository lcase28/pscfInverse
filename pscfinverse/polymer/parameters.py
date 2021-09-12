"""
This module contains simple data wrappers to be used by simulator
mesophase objects. These classes represent the basic, direct 
parameters that are being manipulated, and allow creation of 
complex "variable" objects, without adding complexity to 
interface with calculation managers.

The idea is that each class here can map 1:1 with
parameters that would be used in definitions of SCFT problems.

The classes defined here are meant to serve two purposes:

    1.  The parameters provide a uniform interface with which
        SCFT interfaces can work.
    2.  They are used in the construction of linear systems of
        equations to translate PSO variable values and constraints
        to usable polymer values.
"""
# Standard Library imports
from abc import ABC, abstractmethod
## TODO: Create set of common error-checks to be used locally (i.e. check for valid ID numbers)

class PolymerParameter(ABC):
    @abstractmethod
    def __init__(self, value):
        self._value = value
    
    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self,val):
        self._value = val
    
    @property
    @abstractmethod
    def key(self):
        """
        Return a hashable key to uniquely identify the parameter.
        
        Generally, a good key would be a tuple containing a string
        identifying the parameter type (shared by all instances of
        the same parameter class), followed by integer identifiers.
        These keys should be constructed such that two parameter
        objects specifying the same part of the polymer system would
        have the same key, even if holding different values.
        
        These keys are typically used for resolution of systems of
        equations with complex variables.
        """
        pass
    
    @property
    @abstractmethod
    def label(self):
        """
        Return a string label for the parameter.
        
        Like the key property, the form of the label should be unique
        to the parameter type, with necessary identifiers (polymer ID
        numbers, monomer ID numbers, etc) included such that two
        instances specifying the same parameter will have the same
        label.
        """
        pass
        
class BlockLength(PolymerParameter):
    def __init__(self, polymerid, blockid, value):
        self.__pid = polymerid
        self.__blockid = blockid
        super().__init__(value)
        
    @property
    def polymerID(self):
        return self.__pid
    
    @property
    def blockID(self):
        return self.__blockid
    
    @property
    def key(self):
        return ("BlockLength", self.polymerID, self.blockID)
    
    @property
    def label(self):
        return "Block_Len_p{}_b{}".format(self.polymerID,self.blockID)
    
    def __str__(self):
        return "<BlockLength Parameter, polymer {}, block {}, value {}>".format(self.polymerID, self.blockID, self.value)
    
class ChiN(PolymerParameter):
    def __init__(self, mon1, mon2, value):
        self.__m1 = min(mon1,mon2)
        self.__m2 = max(mon1,mon2)
        
        # Ensure monomer id numbering starts at 0
        diff = self.monomer2 - self.monomer1
        if self.monomer1 < 0:
            raise(ValueError("Monomer IDs must be >= 0"))
        if diff < 1:
            raise(ValueError("Monomer IDs must differ by at least 1"))
        super().__init__(value)
    
    @property
    def monomer1(self):
        return self.__m1
    
    @property
    def monomer2(self):
        return self.__m2
    
    @property
    def monomerIDs(self):
        """ 
        Monomer ID numbers, with indexing starting at 1 
        Monomer1 < Monomer2
        """
        return self.monomer1, self.monomer2
    
    @property
    def key(self):
        return ( "ChiN", self.monomer1, self.monomer2 )
    
    @property
    def label(self):
        return "ChiN_m{}_m{}".format(*self.monomerIDs)
    
    def __str__(self):
        return "<ChiN Parameter, monomers ({},{}), value {}>".format(self.monomer1, self.monomer2, self.value)

class KuhnLength(PolymerParameter):
    def __init__(self, monomerid, value):
        self.__mid = monomerid
        super().__init__(value)
        
    @property
    def monomerID(self):
        return self.__mid
    
    @property
    def key(self):
        return ( "KuhnLength", self.monomerID )
    
    @property
    def label(self):
        return "Kuhn_m{}".format(self.monomerID)
    
    def __str__(self):
        return "<KuhnLength Parameter, monomer {}, value {}>".format(self.monomerID, self.value)

class PolymerBlendFraction(PolymerParameter):
    def __init__(self, polymerid, value):
        self.__pid = polymerid
        super().__init__(value)
    
    @property
    def polymerID(self):
        return self.__pid
    
    @property
    def key(self):
        return ( "PolymerBlendFraction", self.polymerID )
    
    @property
    def label(self):
        return "Blend_Frac_p{}".format(self.polymerID)
    
    def __str__(self):
        return "<PolymerBlendFraction Parameter, polymer {}, value {}>".format(self.polymerID, self.value)

