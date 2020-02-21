"""Module contains core classes for Agent implementation."""

# Third-party imports
from abc import ABC, abstractmethod
from enum import Enum, Flag, unique
from copy import deepcopy
import numpy as np

# Project imports
from .psocontainers import PsoPositionData, PsoVariableSet
from .psovectors import Point, Velocity, OptimizationType, FitnessComparator

@unique
class AGENT_STATUS_TYPES(Flag):
    STABLE = 0
    STARTUP = auto()
    UPDATE = auto()
    TRY = auto()
    EVALUATE = auto()
    CONFIRM = auto()
    ERROR = auto()
    INITIALIZED = ~STARTUP
    UPDATE_ALL = UPDATE | TRY | EVALUATE | CONFIRM
    
class Agent(ABC):
    """Abstract base class for PSO Agents."""
    
    __NextID = 0  # Static counter for generating unique agent IDs

    def __init__(self, varSet, velocitySource, fitnessCompare=None, randGen = None, **kwargs):
        """
            Constructor for Agent base class.
            
            When generating a set of agents for a swarm, the same inputs can be used
            to all constructors because any input not treated as read-only is cloned 
            with deepcopy to ensure logical isolation of the different instances.
            
            Parameters
            ----------
            varSet : PsoVariableSet (or subclass)
                The variables being considered in the PSO.
            velocitySource : Velocity instance
                A template velocity for the agent.
            randGen : numpy.random.RandomState, optional
                A RandomState object to use to initialize positions.
                If not included, values in ptSrc and velSrc are used.
            lowerSpawnRange : array-like of real, optional
                The lower bound of the range in which to initialize the agent.
                Ignored if randGen is excluded.
            upperSpawnRange : array-like of real, optional
                The upper bound of the range in which to initialize the agent.
                Required if lowerSeedRange is included. Ignored if randGen is excluded.
            velocitySpawnBound : array-like of real, optional
                The maximum magnitude of velocity in each dimension on initialization.
        """
        # Deep copy critical inputs
        varSet = deepcopy(varSet)
        velocitySource = deepcopy(velocitySource)
        
        # Initialize status flag
        self.__status = AGENT_STATUS_TYPES.STARTUP
        
        # Fitness comparison object
        if fitnessCompare is None:
            fitnessCompare = FitnessComparator()
        self.__fitCompare = fitnessCompare
        
        # Initialize to zero steps taken
        self.steps = -1

        # Grab a unique ID for this agent
        self.id = Agent.__NextID
        Agent.__NextID += 1 # Update the static counter to deliver unique IDs
        
        # Partial Initialization of position and velocity data.
        self.__positions = PsoVariableSet(varSet)
        self.__velocity = velocitySource
        self.__old_velocity = deepcopy(self.__velocity)
        self.__old_point = None
        
        # Select whether or not to randomize initial position and velocity
        if randGen is None:
            # Do not randomize... use given values
            initPosition = self.__positions.currentPsoValues
            initVelocity = self.__velocity.components
        else:
            # Randomize position and velocity
            initPosition, initVelocity = self.__randomize(randGen, *args, **kwargs)
        
        # Finish Initializing
        self.tryUpdate(initPosition, initVelocity)
        
    def __randomize(self, randGen, \
                    lowerSpawnRange=None, \
                    upperSpawnRange=None, \
                    velocitySpawnBound = None):
        """
            Return a new instance at randomized initial position and velocity.
            
            Parameters
            ----------
            randGen : numpy.random.Random
                Pre-seeded random number generator to use for initial values
            lowerSpawnRange : array-like of real, optional
                The lower bound of the range in which to initialize agent positions.
                If None (Default), the lower bound of the varSet will be used.
                If included, values should correspond to "psoValue" quantities.
            upperSpawnRange : array-like of real, optional.
                The upper bound of the range in which to initialize agent positions.
                If None (Default), the lower bound of the varSet will be used.
                If included, values should correspond to "psoValue" quantities.
            velocitySpawnBound : array-like of real, optional.
                The maximum magnitude of velocity components to use, such that initial
                velocity values will be in range [-bound, bound]. If None (default),
                a value equal to 1/10 of the chosen spawn range width ( = upper - lower ).
                If a negative value is given, the magnitude is used.
                If included, values should be "psoValue" scaled quantities.
            
            Returns
            -------
            randomPosition : numpy.ndarray
                A randomly generated position array with size equal to dimensionality of the search
            randomVelocity : numpy.ndarray
                A randomly generated velocity array with size equal to dimensionality of the search
        """
        # Select lower and upper bounds for position selection
        if lowerSpawnRange is None:
            lowerSpawnRange = self.__positions.variableSet.psoLowerBounds
        else:
            lowerSpawnRange = np.array(lowerSpawnRange).flatten()
        if upperSpawnRange is None:
            upperSpawnRange = self.__positions.variableSet.psoUpperBounds
        else:
            upperSpawnRange = np.array(upperSpawnRange).flatten()
        if not upperSpawnRange.size == lowerSpawnRange.size:
            msg = "Spawn range dimensions disagree.\n\tLower: {}\n\tUpper: {}"
            raise(ValueError(msg.format(lowerSpawnRange,upperSpawnRange)))
        if not upperSpawnRange.size == self.__positions.dimensions:
            msg = "Spawn range dimensions don't match problem dimensionality ({}).\n\tLower: {}\n\tUpper: {}"
            reaise(ValueError(msg.format(self.__positions.dimensions, lowerSpawnRange, upperSpawnRange)))
        
        # Randomize position
        randomPosition = randGen.uniform(lowerSpawnRange, upperSpawnRange)
        
        # Determine velocity constraints
        if velocitySpawnBound is None:
            velocitySpawnBound = upperSpawnRange - lowerSpawnRange
            velocitySpawnBound *= 0.1
        else:
            velocitySpawnBound = np.array(velocitySpawnBound).flatten()
        if not velocitySpawnBound.size == upperSpawnRange.size:
            msg = "Velocity and Position spawn range dimensions disagree.\n\tVelocity: {}\n\tPosition: {}"
            raise(ValueError(msg.format(velocitySpawnRange,upperSpawnRange)))
        velocitySpawnBound = np.absolute(velocitySpawnBound)
        randomVelocity = randGen.uniform(-velocitySpawnBound, velocitySpawnBound)
        
        return randomPosition, randomVelocity
    
    @classmethod
    def createSeveral(cls, nAgent, parallelManager, *args, **kwargs):
        out = []
        for i in range(nAgent):
            out.append(cls(*args, **kwargs))
        #TODO: enable parallel evaluation
        for a in out:
            a.evaluateUpdate()
        for a in out:
            a.confirmUpdate()
    
    def __enter_update_mode(self):
        self._endErrorState()
        # remove record of error-retained update flags
        self.__endState(AGENT_STATUS_TYPES.UPDATE_ALL) 
        self.__startState(AGENT_STATUS_TYPES.UPDATE)
        # Store old values
        self.__old_velocity.components = self.__velocity.components
        self.__old_point = self.__positions.currentPoint
        # Additional flag
        self.__has_evaluated = False
        
    def __exit_update_mode(self):
        self.__endState(AGENT_STATUS_TYPES.UPDATE)
        # as soon as first update succeeds, remove startup flag. No effect if error occurs
        self.__endState(AGENT_STATUS_TYPES.STARTUP)
        self.__old_point = None
        self.__new_fitness = None
        self.__steps += 1
    
    def tryUpdate(self, newPos, newVel):
        """
        Set the agent to update mode, and attempt to update the position and velocity.
        
        For any dimension where the new position is rejected, the velocity is negated.
        
        If the accepted new position is out of the search bounds, an error state is triggered.
        
        Parameters
        ----------
        newPos : array-like of real
            The coordinates of the new PsoPosition
        newVel : array-like of real
            The components of the velocity vector.
            
        Returns
        -------
        True if values successfully updated and in bounds.
        False if an error occurred or the position is out of bounds.
        """
        self.__enter_update_mode()
        self.__startState(AGENT_STATUS_TYPES.TRY)
        # Update velocity
        self.__velocity.components = newVel
        # Attempt to update position
        accepted = self.__position.tryUpdate(newPos)
        for (i,b) in enumerate(accepted):
            if not b:
                self.__velocity.reverseComponent(i)
        if not self.inBounds:
            self._startErrorState()
            return False
        # If no error thrown, remove UPDATE_VECTORS flag to indicate successful update.
        self.__endState(AGENT_STATUS_TYPES.TRY)
        return True
        
    def evaluateUpdate(self):
        """ 
        Calculate the fitness for the position set during tryUpdate().
        
        If the agent is in an error state (which was cleared at start of the update), no evaluation is done.
        
        Returns
        -------
        fitness : numerical
            The fitness the agent will have if the update is confirmed.
        success : bool
            True if evaulation succeeded. False if evaluation failed or did not occur at all.
        
        Raises
        ------
        RuntimeError if an update was not started before calling this method.
        """
        if not self.inUpdate():
            raise(RuntimeError("No active update pending."))
        self.__has_evaluated = True
        self.success = False
        # If error state previously triggered and not ended, do not evaluate.
        if not self.inErrorState:
            self.__startState(AGENT_STATUS_TYPES.EVALUATE)
            self.__new_fitness, success = self.__evaluate
            if not success:
                self._startErrorState()
            self.__endState(AGENT_STATUS_TYPES.EVALUATE)
        return self.testFitness, success
    
    def confirmUpdate(self):
        """ 
        Finalize and accept the update. If fitness not yet calculated, calculate it then confirm. 
        """
        if not self.inUpdate():
            raise(RuntimeError("No active update pending."))
        if not self.__has_evaluated:
            self.evaluateUpdate()
        self.__positions.confirmUpdate(self.testFitness)
        self.__exit_update_mode()
    
    def cancelUpdate(self):
        if not self.inUpdate():
            raise(RuntimeError("No active update pending."))
        self.__positions.cancelUpdate(self.__new_fitness)
        self.__velocity.components = self.__old_velocity.components
        self.__exit_update_mode()
    
    @property
    def inErrorState(self):
        return self.__status & AGENT_STATUS_TYPES.ERROR
    
    @property
    def inUpdate(self):
        return self.__status & AGENT_STATUS_TYPES.UPDATE
        
    def _startError(self):
        self.__status_on_error
        
    def _startState(self, state):
        self.__status = self.__status | state
        
    def _endState(self, state):
        if not self.inErrorSate:
            self.__status = self.__status & ~state
    
    @property
    def stablePoint(self):
        return self.__positions.stablePoint
    
    @property
    def bestPoint(self):
        return self.__positions.bestPoint
    
    @property
    def stableVelocity(self):
        if not self.inUpdate:
            return deepcopy(self.__velocity)
        else:
            return deepcopy(self.__old_velocity)
    
    @property
    def point(self):
        return Point(self.__positions.currentPsoValues, self.testFitness)
    
    @property
    def velocity(self):
        return deepcopy(self.__velocity)
    
    @property
    def testFitness(self):
        if self.inErrorState:
            return self.__fitCompare.badFitness
        else:
            return self.__new_fitness
    
    # Derived classes must override, update fitness, and call base using super calls to generate the MRO.
    @abstractmethod
    def evaluate(self):
        if self.PBest is None:
            # Initialization condition
            self.PBest = deepcopy(self.Location)
        # Now that fitness has been updated, compare to PBest
        if self.Location > self.PBest and self.seekMax:
            self.PBest.fill_from(self.Location)
        elif self.Location < self.PBest and not self.seekMax:
            self.PBest.fill_from(self.Location)
        return True
    
    def __str__(self):
        s = "Agent {}:\n\tScaled: {}\n\tCoords: {}\n\tVelcty: {}\n\tFitnes: {}"
        s = s.format(self.id, self.get_coords(), self.Location.Coords, self.Velocity, self.Location.Fitness)
        return s

class FunctionAgent(Agent):
    def __init__(self, fn, **kwargs):
        self.Function = fn
        
        super().__init__(**kwargs)
        
        self.evaluate()

    def evaluate(self):
        self.Location.Fitness = self.Function(self.get_coords())

        return super().evaluate()


