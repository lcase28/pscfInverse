"""
Module defining methodology of PSO variable management.
"""

# Library imports

# Third Party Imports
from abc import ABC, abstractmethod
import copy
from functools import total_ordering
from enum import Enum, unique
import numpy as np

@unique
class BoundTypes(Enum):
    reflective = 1
    periodic = 2
    fly = 3

class PsoVariable(ABC):
    """
    Base class for all PSO variables. Defines the basic
    operations expected of all variables going to be used
    in PSO optimizations.
    
    Abstract in that it does not define a scaling behavior
    """
    __next_ID = 0
    
    def __init__(self, value, lower, upper, \
                psoValues=True, \
                boundsType=BoundTypes.reflective, \
                label=None, \
                initLower = None, \
                initUpper = None, \
                initIsPsoValue = True)
        """
        Initialize a basic PsoVariable object.
        
        Inheriting classes should perform actions required to 
        enable translations between pso and true values, then
        call this initializer to set underlying values.
        
        Parameters
        ----------
        value : real
            The current value of the variable.
        lower : real
            The lower bound of the variable
        upper : real > lower
            The upper bound of the variable
        psoValues : bool, optional (default=True)
            If True (default), the previous 3 parameters are taken
            to refer to the 'psoValue' and 'psoBounds' of 
            the variable.
            If False, the previous 3 parameters are taken to 
            refer to the 'trueValue' and 'trueBounds' of
            the variable.
        boundsType : psovariables.BoundTypes
            Defines the bounding behavior of the variable.
            See psovariables.BoundTypes for details.
        label : string (Optional)
            A label to use when displaying values
        """
        # Define ID
        self.__id = self.__class__.__next_ID
        self.__class__.__next_ID += 1
        
        # create underlying memebers
        self.__psoLowerBound = None
        self.__psoUpperBound = None
        self.__psoValue = None
        
        # set bounds style
        self.boundStyle = boundsType
        
        # underlying data is the psoValue
        if psoValues:
            self.psoBounds = [lower, upper]
            self.psoValue = value
        else:
            self.trueBounds = [lower, upper]
            self.trueValue = value
        
        # Set label
        if label is not None:
            self.__label = label
        else:
            self.__label = "Variable{}".format(self.__id)
        
    # Practical methods
    def checkBounds(self, val, psoValue = True):
        if psoValue:
            lb = self.psoLowerBound
            ub = self.psoUpperBound
        else:
            lb = self.trueLowerBound
            ub = self.trueUpperBound
        return val >= lb and val <= ub
        
    def checkAcceptance(self, val, psoValue = True):
        cb = self.checkBounds(val, psoValue)
        if self.boundStyle == BoundTypes.reflective:
            return cb
        else:
            return True
    
    # abstract Methods
    @abstractmethod
    def __pso_to_true(self, val):
        """
        Convert the scaled pso value to unscaled
        true value. This must perform the inverse
        operation as __true_to_pso(...) such that 
        val == self.__true_to_pso(self.__pso_to_true(val))
        is tautological (**within numerical precision**)
        
        Default behavior is to apply no scaling. Inheriting
        classes can select this behavior with a super() call.
        
        This operation must be defined on an absolute basis,
        and should not be dependent on the bounds of the variable.
        If a scaling based on bounds is desired in a derived class,
        be sure that all property implementations are adjusted accordingly.
        """
        return val
        
    @abstractmethod
    def __true_to_pso(self, val):
        return val
    
    # Properties
    @property
    def psoValue(self):
        return self.__psoValue
    
    @psoValue.setter
    def psoValue(self, newValue):
        if self.boundStyle == BoundTypes.reflective:
            if self.checkBounds(newValue):
                self.__psoValue = newValue
            else:
                warnMsg = "{} rejected position update to {}"
                warnMsg = warnMsg.format(self.__label, newValue)
                raise(RuntimeWarning(warnMsg))
        elif self.boundStyle == BoundTypes.periodic:
            shiftVal = newValue - self.psoLowerBound
            modVal = shiftVal % self.psoRange
            self.__psoValue = modVal + self.psoLowerBound
        elif self.boundStyle == BoundTypes.fly:
            self.__psoValue = newValue
        
    @property
    def psoBounds(self):
        return [self.__psoLowerBound, self.__psoUpperBound]
    
    @psoBounds.setter
    def psoBounds(self,newbnds):
        lowbnd = newbnds[0]
        upbnd = newbnds[1]
        if lowbnd < upbnd:
            self.__psoLowerBound = lowbnd
            self.__psoUpperBound = upbnd
        else:
            raise(ValueError("Attempted to set invalid bounds {}".format(newbnds)))
    
    @property
    def psoLowerBound(self):
        return self.__psoLowerBound
    
    @psoLowerBound.setter
    def psoLowerBound(self, newval):
        if newval < self.psoUpperBound:
            self.__psoLowerBound = newval
        else:
            errmsg = "Attempted to set invalid lower bound {} >= {}"
            errmsg = errmsg.format(newval, self.psoUpperBound)
            raise(ValueError(errmsg))
        
    @property
    def psoUpperBound(self):
        return self.__psoUpperBound
        
    @psoUpperBound.setter
    def psoUpperBound(self, newval):
        if newVal > self.psoLowerBound:
            self.__psoUpperBound = newval
        else:
            errmsg = "Attempted to set invalid upper bound {} <= {}"
            errmsg = errmsg.format(newval, self.psoLowerBound)
            raise(ValueError(errmsg))
    
    @property
    def trueValue(self):
        return self.__pso_to_true(self.psoValue)
    
    @trueValue.setter
    def trueValue(self,newVal):
        self.psoValue = self.__true_to_pso(newVal)
    
    @property
    def trueLowerBound(self):
        return self.__pso_to_true(self.psoLowerBound)
    
    @trueLowerBound.setter
    def trueLowerBound(self, newval):
        self.psoLowerBound = self.__true_to_pso(newval)
    
    @property
    def trueUpperBound(self):
        return self.__pso_to_true(self.psoUpperBound)
    
    @trueUpperBound.setter
    def trueUpperBound(self, newval):
        self.psoUpperBound = self.__true_to_pso(newval)
    
    @property
    def trueBounds(self):
        return [self.trueLowerBound, self.trueUpperBound]
    
    @trueBound.setter
    def trueBounds(self, newvals):
        psoLowBnd = self.__true_to_pso(newvals[0])
        psoUpBnd = self.__true_to_pso(newvals[1])
        self.psoBounds = [psoLowBnd, psoUpBnd]
            
class UnscaledVariable(PsoVariable):
    """
    Derived variable class accepting PsoVariable's default,
    unscaled behavior. psoValue and trueValue are treated as interchangeable
    """
    
    def __init__(self,*args,**kwargs):
        super().__init__(*args, **kwargs)
    
    def __pso_to_true(self,val):
        return val
    
    def __true_to_pso(self,val):
        return val
    
class LinearScaledVariable(PsoVariable):
    """
    A constant scaling factor is used between true and pso value.
    
    psoValue = A * trueValue
    trueValue = psoValue / A
    
    for a scalar constant A
    """
    
    def __init__(self, scale, *args, **kwargs):
        self.__scalefactor = scale
        super().__init__(*args, **kwargs)
        
    def __pso_to_true(self, val):
        return val / self.__scalefactor
        
    def __true_to_pso(self, val):
        return val * self.__scalefactor
    
class LinearOffsetVariable(LinearScaledVariable):
    """
    Conversion factors in scaling are given as
    
    p = A * (t - r)
    t = (p / A) + r
    
    """
    
    def __init__(self, offset, *args, **kwargs):
        self.__scaleoffset = offset
        super().__init__(*args, **kwargs)
    
    def __pso_to_true(self, val):
        return super().__pso_to_true(val) + self.__scaleoffset
    
    def __true_to_pso(self, val):
        return super().__true_to_pso(val - self.__scaleoffset)
    
class LogScaledVariable(PsoVariable):
    """
    A the true value is scaled logarithmically for the pso search.
    
    psoValue = ln(trueValue)
    trueValue = exp(psoValue)
    """
    
    def __init__(self, scale, *args, **kwargs):
        self.__scalefactor = scale
        super().__init__(*args, **kwargs)
        
    def __pso_to_true(self, val):
        return np.exp(val)
        
    def __true_to_pso(self, val):
        return np.log(val)
    



