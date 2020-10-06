"""
This module contains simple data wrappers to be used by simulator
mesophase objects. These classes represent the basic, direct 
parameters that are being manipulated, and allow creation of 
complex "variable" objects, without adding complexity to 
interface with calculation managers.

The idea is that each class here can map 1:1 with
parameters that would be used in definitions of SCFT problems.
"""

## TODO: Create set of common error-checks to be used locally (i.e. check for valid ID numbers)

class BlockLength(object):
    
    def __init__(self, polymerid, blockid, value):
        self.__pid = polymerid
        self.__blockid = blockid
        self.__value = value
        
    @property
    def polymerID(self):
        return self.__pid
    
    @property
    def blockID(self):
        return self.__blockid
    
    @property
    def value(self):
        return self.__value
    
class ChiN(object):
    
    def __init__(self, mon1, mon2, value):
        self.__m1 = min(mon1,mon2)
        self.__m2 = max(mon1,mon2)
        self.__value = value
        
        # Ensure monomer id numbering starts at 0
        diff = self.monomer2 - self.monomer1
        if self.monomer1 < 0:
            raise(ValueError("Monomer IDs must be >= 0"))
        if diff < 1:
            raise(ValueError("Monomer IDs must differ by at least 1"))
    
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
    def value(self):
        return self.__value

class KuhnLength(object):
    
    def __init__(self, monomerid, value):
        self.__mid = monomerid
        self.__value = value
        
    @property
    def monomerID(self):
        return self.__mid
    
    @property
    def value(self):
        return self.__value

class PolymerBlendFraction(object):
    
    def __init__(self, polymerid, value):
        self.__pid = polymerid
        self.__value = value
    
    @property
    def polymerID(self):
        return self.__pid
    
    @property
    def value(self):
        return self.__value

