"""
Module defining pso variable sets and points
"""

# Library imports
from psoinverse.PSO.psovariables import PsoVariable
from psoinverse.PSO.psovectors import Point

# Third-party imports
import numpy as np
from enum import Enum

class ContainerStatusTypes(Enum):
    stable = 1
    uninitialized = 2
    unstable = 3
    error = 4
    unknown = 5

class ContainerStatusError(Exception):
    """
    Raised when an action on a container which is not allowed
    in its current state
    
    Attributes
    ----------
        caller_class
            The Class of the container raising the error
        caller_status : psoinverse.PSO.psocontainers.ContainerStatusTypes
            The status of the caller at the time of the error
        trigger_action
            The attempted action (as a method object) which triggered
            the error to be thrown
        explanation : string
            Explanation of why the action is prohibited in this instance.
    """
    def __init__(self, src, status, action, message):
        self.caller_class = src
        self.caller_status = status
        self.trigger_action = action
        self.explanation = message
    
class PsoVariableSet(object):
    """
    Class defining basic variable set management
    """
    
    def __init__(self, variables):
        """
        Once initialized, variable set is fixed.
        
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
        
        # Attributes for use during reversible value updates
        self.__soft_update_flag = False   # flag indicating soft update not yet confirmed.
        self.__old_pso_data = None        # old variable values (psoValue)
        self.__old_true_data = None       # old variable values (trueValue)
        self.__update_acceptance = None   # position acceptance flags for update
        self.__update_bounds_check = None # position in-bounds flags for update
        
    # Methods for monitored updates of variable values
    
    def checkBounds(self, testposition, psoVals=True):
        val = np.asarray(testposition).flatten()
        if not val.size == self.__nvar:
            raise(ValueError("Checking bounds on position of invalid size"))
        res = []
        for i in range(self.__nvar):
            res.append( self.__variables[i].checkBounds(val[i], psoVals) )
        return res
    
    def checkAcceptance(self, testposition, psoVals=True):
        val = np.asarray(testposition).flatten()
        if not val.size == self.__nvar:
            raise(ValueError("Checking bounds on position of invalid size"))
        res = []
        for i in range(self.__nvar):
            res.append( self.__variables[i].checkAcceptance(val[i], psoVals) )
        return res
        
    def tryUpdate(self, source, psoVals=True):
        """
            Perform a soft update on variable values.
            This will update the values in the variables, but internally
            store the pre-update values to allow for reversion to pre-update
            state.
            
            See PsoVariableSet.confirmUpdate for "fixing" the update.
            
            See PsoVariableSet.cancelUpdate for "reversion" of the update.
            
            Parameters
            ----------
            source : array-like of real
                The new values to set the variables to. Order assumed to match
                that returned from PsoVariableSet.[pso/true]Values property.
            psoVals : bool (optional)
                If True (default), assumes that the values in source correspond
                to the psoValue of the variables. If False, assumes the trueValue.
                
            Returns
            -------
            success : numpy.ndarray of bool
                Array  of flags indicating if the updated value was accepted (True)
                or rejected (False) for each variable.
                
            Raises
            ------
            ContainerStatusError
                If tryUpdate has already been called without confirming or cancelling
                the update, and tryUpdate is called again.
            ValueError
                If the given data is not the same size as the variable set.
        """
        src = np.array(source).flatten()
        if not src.size == self.__nvar:
            msg = "Cannot update {}. Size mismatch with source, {}"
            raise(ValueError(msg.format(self.__class__.__name__, source)))
        if self.__soft_update_flag:
            raise(ContainerStatusError(self.__class__, self.status, self.tryUpdate, \
                        "Must confirm or cancel update before attempting new update"))
        
        self.__old_pso_data = self.psoValues
        self.__old_true_data = self.trueValues
        self.__update_bounds_check = self.checkBounds(src)
        success = self.checkAcceptance(src)
        self.__update_acceptance = self.checkAcceptance(src)
        
        if psoVals:
            self.psoValues = src
        else:
            self.trueValues = src
        
        self.__soft_update_flag = True
        
        return success
        
    def confirmUpdate(self):
        """
        When a soft update has been started, this method fixes the update as the stable
        state of the object.
        
        Returns
        -------
        success : numpy.ndarray of bool
            True if the updated value was accepted. False if it was rejected.
            This value is the same as that returned by tryUpdate.
            
        Raises
        ------
        ContainerStatusError
            If no soft update has been initiated with tryUpdate prior to the call
        """
        if not self.__soft_update_flag:
            raise(ContainerStatusError(self.__class__, self.status, self.confirmUpdate, \
                        "Must initialize a soft update before confirming one."))
        success = self.__update_acceptance
        self.__clear_soft_backup()
        return success
        
    def cancelUpdate(self):
        """
        When a soft update has been started, this method aborts the soft update and
        restores the object to its status before the update.
        
        If no soft update has been started, This method silently performs no action.
        """
        if not self.__soft_update_flag:
            return
        src = self.__old_data
        self.__clear_soft_backup()
        self.psoValues = src
        
    def forceUpdate(self, source, psoVals=True):
        """
            Applies the update and confirms it. Has the same effect as consecutive calls to
            tryUpdate and confirmUpdate.
            
            Parameters
            ----------
            source : array-like of real
                The new values to set the variables to. Order assumed to match
                that returned from PsoVariableSet.[pso/true]Values property.
            psoVals : bool (optional)
                If True (default), assumes that the values in source correspond
                to the psoValue of the variables. If False, assumes the trueValue.
                
            Returns
            -------
            success : numpy.ndarray of bool
                Array  of flags indicating if the updated value was accepted (True)
                or rejected (False) for each variable.
                
            Raises
            ------
            ContainerStatusError
                If a soft update has been started but not confirmed or canceled.
            ValueError
                If the given data is not the same size as the variable set.
        """
        if self.__soft_update_flag:
            raise(ContainerStatusError(self.__class__, self.status, self.forceUpdate, \
                        "Must confirm or cancel update before attempting new update"))
        success = self.tryUpdate(source,psoVals)
        suc = self.confirmUpdate()
        return success
    
    def __clear_soft_backup(self):
        self.__update_acceptance = None
        self.__update_bounds_check = None
        self.__old_pso_data = None
        self.__old_true_data = None
        self.__soft_update_flag = False
        
    # Properties, and unmonitored modifiers
    
    @property
    def dimensions(self):
        tmp = self.stablePsoValues
        return tmp.size
    
    @property
    def status(self):
        if self.__soft_update_flag:
            return ContainerStatusTypes.unstable
        return ContainerStatusTypes.stable
    
    @property
    def inBounds(self):
        """ True if all stable pso values are in bounds for the variables. """
        tmp = True
        for b in self.checkBounds(self.stablePsoValues):
            tmp = tmp and b
        return tmp
    
    @property
    def stablePsoValues(self):
        """ The most recent stable psoValue data. Constant return value between confirmUpdate() calls. """
        if self.__soft_update_flag:
            return np.array(self.__old_pso_data)
        else:
            return self.psoValues
    
    @property
    def psoValues(self):
        return np.asarray( [v.psoValue for v in self.__variables] )
    
    @psoValues.setter
    def psoValues(self, newVals):
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "psoValues.setter", "Unable to update unstable values"))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "trueValues.setter", "Unable to update unstable values."))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "psoLowerBounds.setter", "Unable to update bounds of unstable set"))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "psoUpperBounds.setter", "Unable to update bounds of unstable set"))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "psoBounds.setter", "Unable to update bounds of unstable set"))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "trueLowerBounds.setter", "Unable to update bounds of unstable set"))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "trueUpperBounds.setter", "Unable to update bounds of unstable set."))
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
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, "trueBounds.setter", "Unable to update bounds of unstable set."))
        lowbnds = np.asarray(newbnds[0]).flatten()
        upbnds = np.asarray(newbnds[1]).flatten()
        if not lowbnds.size == upbnds.size:
            raise(ValueError("Dimensionality of bounds must agree"))
        if not lowbnds.size == self.__nvar:
            raise(ValueError("New value sets for trueBounds are improper size"))
        for i in range(self.__nvar):
            v = self.__variables[i]
            v.trueBounds = [lowbnds[i], upbnds[i]]
            
class PsoPositionData(object):
    """
    Wrapper object for the PsoVariableSet which also stores
    current and best positions and fitnesses.
    """
    
    ## TODO: Implement queue-based position history storage with max capacity.
    
    def __init__(self, varSet):
        """
        Partial initialization of object
        """
        if not isinstance(varSet, PsoVariableSet):
            raise(TypeError("{} is not a valid type for a variable set".format(varSet.__class__.__name__)))
        self.__variableSet = varSet
        self.__status = ContainerStatusTypes.uninitialized
        self.__fitness = None
        self.__best_point = None
        self.__old_fitness = None
        
    @classmethod
    def initializedInstance(cls, varSet, fitness):
        instance = cls(varSet)
        pt = instance.__variableSet.psoValues
        instance.forceUpdate(pt, fitness)
        instance.setCurrentAsBest()
        return instance
        
    # Updating methods
    
    def tryUpdate(self, src, psoValues=True):
        """
        Reversibly update the current position.
        
        Parameters
        """
        if self.status == ContainerStatusTypes.unstable:
            raise(ContainerStatusError(self.__class__, self.status, self.tryUpdate, \
                        "Cannot try new update until prior try is confirmed or cancelled."))
        success = self.__variableSet.tryUpdate(src, psoValues)
        self.__old_fitness = self.__fitness
        self.__fitness = None
        self.__status = ContainerStatusTypes.unstable
        return success
        
    def confirmUpdate(self, fitness):
        """
        Accept currently attempted update and set fitness to given value
        """
        if not self.status == ContainerStatusTypes.unstable:
            raise(ContainerStatusError(self.__class__, self.status, self.confirmUpdate, \
                    "Must initiate a soft update before confirming the soft update."))
        success = self.__variableSet.confirmUpdate()
        self.__fitness = fitness
        self.__old_fitness = None
        self.__status = ContainerStatusTypes.stable
        return success
        
    def cancelUpdate(self):
        """ Reject active update attempt """
        if not self.status == ContainerStatusTypes.unstable:
            raise(ContainerStatusError(self.__class__, self.status, self.cancelUpdate, \
                    "Must initiate a soft update before cancelling the soft update."))
        self.__variableSet.cancelUpdate()
        self.__fitness = self.__old_fitness
        self.__status = ContainerStatusTypes.stable
        
    def forceUpdate(self, src, fitness, psoVals=True):
        self.tryUpdate(src, psoVals)
        success = self.confirmUpdate(fitness, psoVals)
        return success
    
    def setCurrentAsBest(self):
        """ Set the current position as the best position seen """
        if not self.status == ContainerStatusTypes.stable:
            raise(ContainerStatusError(self.__class__, self.status, self.setCurrentAsBest, \
                    "The container must be stable before setting the best position."))
        self.__best_fitness = self.currentFitness
        self.__best_pso_coords = self.currentPsoValues
        self.__best_true_coords = self.currentTrueValues
        
    # Properties
    
    @property
    def dimensions(self):
        return self.__variableSet.dimensions
    
    @property
    def stablePsoValues(self):
        if self.__status == ContainerStatusTypes.uninitialized:
            msg = "No stable state has yet been set. No stable data Available."
            raise(ContainerStatusError(self.__class__, self.__status, "stablePsoValues", msg))
        return self.__variableSet.stablePsoValues
    
    @property
    def stableTrueValues(self):
        if self.__status == ContainerStatusTypes.uninitialized:
            msg = "No stable state has yet been set. No stable data Available."
            raise(ContainerStatusError(self.__class__, self.__status, "stableTrueValues", msg))
        return self.__variableSet.stableTrueValues
    
    @property
    def stableFitness(self):
        if self.__status == ContainerStatusTypes.uninitialized:
            msg = "No stable state has yet been set. No stable data Available."
            raise(ContainerStatusError(self.__class__, self.__status, "stableFitness", msg))
        if self.__status == ContainerStatusTypes.unstable:
            return self.__old_fitness
    
    @property
    def currentPsoValues(self):
        return self.__variableSet.psoValues
        
    @property
    def currentTrueValues(self):
        return self.__variableSet.trueValues
    
    @property
    def currentFitness(self):
        return self.__fitness
    
    @property
    def bestPsoValues(self):
        return np.array(self.__best_pso_coords)
        
    @property
    def bestTrueValues(self):
        return np.array(self.__best_true_coords)
    
    @property
    def bestFitness(self):
        return self.__best_fitness
    
    @property
    def stablePoint(self):
        """ A point containing the most recent stable position and fitness. Constant return between confirmUpdate calls. """
        if self.__status == ContainerStatusTypes.uninitialized:
            msg = "No stable state has yet been set. No stablePoint Available."
            raise(ContainerStatusError(self.__class__, self.__status, "stablePoint", msg))
        return Point(self.stablePsoValues, self.stableFitness)
    
    @property
    def currentPoint(self):
        """ A Point object with current fitness and current pso values. """
        return Point(self.currentPsoValues, self.currentFitness)
    
    @property
    def bestPoint(self):
        """ A Point object with current fitness and current pso values. """
        return Point(self.bestPsoValues, self.bestFitness)
    
    @property
    def variableSet(self):
        """ 
            The underlying variable set used for position management.
            The returned value should be treated as read-only; modification
            may lead to undefined behavior.
        """
        return self.__variableSet
        
# TODO: Implement container class to wrap velocities and manage updates.

