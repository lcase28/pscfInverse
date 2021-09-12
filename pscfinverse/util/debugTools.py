""" Module containing universal class to manage debug mode """
from enum import Enum, unique
import pathlib

@unique
class ReportLevel(Enum):
    none = 0
    basic = 1
    avg = 2
    detail = 3
    

class Debug(object):
    """
    Stores debugging state for full pscfinverse module.
    Called thoughout project stack with various printouts.
    Only prints those at set detail.
    """
    
    _threshold = 0
    
    _logLocation = None
    
    @classmethod
    @property
    def threshold(cls):
        return cls._threshold
    
    @classmethod
    @threshold.setter
    def threshold(cls, val):
        if val < 0:
            cls._threshold = ReportLevel.none
        elif val > 3:
            cls._threshold = ReportLevel.detail
        else:
            cls._threshold = val
    
    @classmethod
    @property
    def logLocation(cls):
        return cls._logLocation
    
    @classmethod
    @logLocation.setter
    def logLocation(cls, val):
        self._logLocation = val
    
    @classmethod
    def write(cls, val, string):
        if val >= cls.threshold:
            cls.__write(string)
    
    @classmethod
    def clearLog(cls):
        if cls._logLocation is not None:
            with open(cls._logLocation,'w') as f:
                f.write('Debug Report\n')
            
    @classmethod
    def __write(cls, s):
        if cls._logLocation is None:
            print(s)
        else:
            with open(cls.logLocation,'a') as f:
                f.write(s)
