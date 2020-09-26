"""Module contains core classes for Agent implementation."""
# Standard Library Imports
from abc import ABC, abstractmethod
from copy import deepcopy
from enum import Enum, Flag, unique
import pathlib

# Third-party imports
import numpy as np

# Project imports
from psoinverse.pso.containers import PsoPositionData, PsoVariableSet
from psoinverse.pso.core import Point, Velocity, FitnessComparator, FITNESS_SELECTOR
from psoinverse.util.iotools import writeCsvLine, checkPath

class Agent(ABC):
    """Abstract base class for PSO Agents."""
    
    __NextID = 0  # Static counter for generating unique agent IDs

    def __init__(self, varSet, velocitySource, calcManager, parentRoot=None):
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
        calcManager : CalculationManager
            Manager (shared among agents) which manages the parallel launching of calculations.
        parentRoot : pathlib.Path or string (optional)
            The root directory of the Swarm's run.
            Each agent will create their own sub-directory to this root 
            in which to write log files.
        """
        # Grab a unique ID for this agent
        self.__id = Agent.__NextID
        Agent.__NextID += 1 # Update the static counter to deliver unique IDs
        
        # Deep copy critical inputs
        varSet = deepcopy(varSet)
        velocitySource = deepcopy(velocitySource)
        
        # Initialize status flags
        self.__updating = False
        self.__next_step = 0
        
        # Partial Initialization of position and velocity data.
        self.__positions = PsoPositionData(varSet)
        self.__velocity = velocitySource
        self.__old_velocity = deepcopy(self.__velocity)
        
        self.__calc_manager = calcManager
        
        # Set Agent Root
        if parentRoot is None:
            parentRoot = pathlib.Path.cwd()
        if not isinstance(parentRoot, pathlib.Path):
            parentRoot = pathlib.Path(parentRoot)
        root = parentRoot / "agent{}".format(self.id)
        root, success = checkPath(root)
        if not success:
            raise(ValueError("Given root path failed to resolve."))
        self.root = root
        
        self._start_logs()
    
    @classmethod
    def createSeveral(cls, nAgent, *args, **kwargs):
        """
        Generate a full set of uninitialized Agents.
        
        A shortcut for generating a swarm of agents from the same input set.
        
        Parameters
        ----------
        nAgent : int > 0
            The number of agents to generate.
        *args : ordered parameters
            Any ordered parameters required for the class's __init__ method.
            See __init__ of the Agent-derived class for details.
        **kwargs : keyworded parameters
            Any keyword parameters required for the class's __init__ method.
            see __init__ of the Agent-derived class for details.
        
        Returns
        -------
        agentSet : list
            A list containing nAgent uninitialized Agents.
        
        Raises
        ------
        ValueError :
            If nAgent <= 0.
        """
        if nAgent <= 0:
            raise(ValueError("Must generate a positive number of agents."))
        out = []
        for i in range(nAgent):
            out.append(cls(*args, **kwargs))
        return out
    
    @staticmethod
    def randomizeSeveral(agents, randGen, lowerBound=None, upperBound=None, velocityBound=None):
        """
        Start an update with randomized values for each in a set of agents.
        
        Caller can optionally specify bounds (psoValue) in which to spawn the position.
        If bounds are excluded, the full search bounds are used.
        
        If velocityBound is not specified, 1/10 of the spawn range is used.
        
        Sequential calls to <CLASS>.createSeveral() and <CLASS>.randomizeSeveral will
        produce a randomly initialized set of Agents with initial calculations started.
        This set will be ready to pass directly to a swarm.
        
        Parameters
        ----------
        agents : iterable set of Agent-derived instances
            The set of agents to randomize.
        randGen : numpy.random.RandomState
            Seeded random number generator from which to draw values
        lowerBound : array-like of real, Optional
            The lower bound of the positional spawn range (psoValue).
            Must be within the bounds permitted by the variable set.
        upperBound : array-like of real, Optional
            The upper bound of the positional spawn range (psoValue).
            Must be within the bounds permitted by the variable set.
        velocityBound : array-like of real, Optional
            The upper bound of the velocity magnitude. 
            Velocity is initialized such that -velocityBound < velocity < velocityBound.
        
        Returns
        -------
        agents : list of Agent-derived instances
            A list containing the same set of Agent objects passed to the method.
        
        Raises
        ------
        ValueError : 
            If spawn bounds conflict with each other (lower > upper),
            If search space bounds are violated (lower or upper is outside of search space).
            If given bounds do not match dimensionality of the search space.
        """
        out = []
        for a in agents:
            a.randomize(randGen, lowerBound, upperBound, velocityBound)
            out.append(a)
        return out
    
    @property
    def id(self):
        return self.__id
    
    @property
    def lastStep(self):
        """ 
        The most recent step for which a stable position is available.
        
        On creation, before initial point is calculated, returns -1. 
        """
        out = self.nextStep - 1
        return out
    
    @property
    def nextStep(self):
        """ The next step the agent will complete. """
        return self.__next_step
    
    def hasStep(self, step):
        """ Determine if the Agent has a stable position for the specified step.
        
        A return value of True indicates that the agent has completed the specified
        step, and that a stable and best point will be known for that step.
        
        Parameters
        ----------
        step : int >= 0
            The step number of interest. Step 0 is considered initialization.
        
        Returns
        -------
        hasStep : boolean
            True if the step has been completed. False otherwise.
            If step < 0 (an invalid step number) silently returns False.
        """
        if step < 0:
            return False
        if self.nextStep > step:
            return True
        else:
            return False
    
    def stablePointAtStep(self, step):
        """ Returns a Point with agent's position and fitness at step. """
        return self.__positions.stablePointAtStep(step)
    
    def bestPointAtStep(self, step):
        """ Returns agent's best-seen point at step. """
        return self.__positions.bestPointAtStep(step)
    
    @property
    @abstractmethod
    def calculationFinished(self):
        """ Determine if calculations required for fitness determination are complete.
        
        ABSTRACT: Derived classes must override. 
        
        The method should check that all calculations required for fitness determination
        are completed. If methods were sent to the parallel computation manager,
        the async_result object returned on submission should be checked for
        calculation completion.
        
        Returns
        -------
        calcDone : boolean
            True if calculations are finished and ready for fitness determination.
            False otherwise.
        """
        pass
    
    @property
    def readyForStartUpdate(self):
        """ Determine if the agent's previous step has completed and new update can start. 
        
        Overriding methods should include a call to super and incorporate its return, 
        A general approach for this could be ( <child_class_ready> and <parent_class_ready> )
        
        Returns
        -------
        isReady : boolean
            True if the agent is ready to begin a new update.
        """
        if not self.updating:
            return True
        else:
            return False
    
    @property
    def readyForFinishUpdate(self):
        """ Determine if a call to finishUpdate(...) is appropriate for the agent.
        
        When true is returned, all preconditions have been met to safely be able to 
        complete an update and end up in a stable condition.
        
        Generally, derived classes should not override this method. Generally, derived
        classes will vary their calculation approaches, which would properly be reflected
        in the self.calculationFinished property, which is checked here.
        
        Overriding methods should include a call to super and incorporate its result.
        
        Returns
        -------
        isReady : boolean
            True if agent is ready to finish the current update.
        """
        if self.updating and self.calculationFinished:
            return True
        else:
            return False
    
    @property
    def stable(self):
        flag = self.readyForStartUpdate and self.nextStep > 0
        return flag
    
    @property
    def updating(self):
        return self.__updating
    
    def randomize(self, randGen, lowerBound=None, upperBound=None, velocityBound=None):
        """
        Start an update with randomized values.
        
        Caller can optionally specify bounds (psoValue) in which to spawn the position.
        If bounds are excluded, the full search bounds are used.
        
        If velocityBound is not specified, 1/10 of the spawn range is used.
        
        Parameters
        ----------
        randGen : numpy.random.RandomState
            Seeded random number generator from which to draw values
        lowerBound : array-like of real, Optional
            The lower bound of the positional spawn range (psoValue).
            Must be within the bounds permitted by the variable set.
        upperBound : array-like of real, Optional
            The upper bound of the positional spawn range (psoValue).
            Must be within the bounds permitted by the variable set.
        velocityBound : array-like of real, Optional
            The upper bound of the velocity magnitude. 
            Velocity is initialized such that -velocityBound < velocity < velocityBound.
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError : 
            If spawn bounds conflict with each other (lower > upper),
            If search space bounds are violated (lower or upper is outside of search space).
            If given bounds do not match dimensionality of the search space.
        """
        # Verify Agent is able to update
        if not self.readyForStartUpdate:
            raise(RuntimeError("Agent {} updated but not ready for update.".format(self.id)))
        
        # Check and choose bounds
        if lowerBound is not None:
            lowerBound = np.array(lowerBound)
            if not len(lowerBound) == self.__positions.dimensions:
                msg = "Dimension mismatch on lower bound: given {}, needs {}"
                raise(ValueError(msg.format(len(lowerBound),self.__positions.dimensions)))
            for b in self.__positions.variableSet.checkBounds(lowerBound):
                if not b:
                    raise(ValueError("Search bound violation on lower bound"))
        else:
            lowerBound = self.__positions.variableSet.psoLowerBounds
        
        if upperBound is not None:
            upperBound = np.array(upperBound)
            if not len(upperBound) == self.__positions.dimensions:
                msg = "Dimension mismatch on upper bound: given {}, needs {}"
                raise(ValueError(msg.format(len(upperBound),self.__positions.dimensions)))
            for b in self.__positions.variableSet.checkBounds(upperBound):
                if not b:
                    raise(ValueError("Search bound violation on upper bound"))
        else:
            upperBound = self.__positions.variableSet.psoUpperBounds
        
        boundRange = upperBound - lowerBound
        
        for i in range(len(lowerBound)):
            if lowerBound[i] >= upperBound[i]:
                msg = "Bounds conflict in dimension {}: lower({}) > upper({})"
                raise(ValueError(msg.format(i,lowerBound[i],upperBound[i])))
        
        if velocityBound is not None:
            velocityBound = np.array(velocityBound)
            if not len(velocityBound) == self.__positions.dimensions:
                msg = "Dimension mismatch on velocity bound: given {}, needs {}"
                raise(ValueError(msg.format(len(velocityBound),self.__positions.dimensions)))
        else:
            velocityBound = (1./10.) * (upperBound - lowerBound)
        
        # Calculate random position
        dim = self.__positions.dimensions
        tempVals = randGen.rand(dim)
        newPos = tempVals * boundRange + lowerBound
        
        # Calculate random velocity
        tempVals = 2. * randGen.rand(dim) - np.ones(dim)
        newVel = tempVals * velocityBound
        
        # Start Update with Random Values
        self.startUpdate(newPos, newVel)
    
    def startUpdate(self, newPos, newVel):
        """
        Start process to update the position, velocity and fitness.
        
        For any dimension where the new position is rejected, the velocity is negated.
        
        Parameters
        ----------
        newPos : array-like of real
            The coordinates of the new PsoPosition.
        newVel : array-like of real
            The components of the velocity vector.
            
        Returns
        -------
        True if values successfully updated and in bounds.
        False if an error occurred or the position is out of bounds.
        """
        # Update velocity
        self.__old_velocity.components = self.__velocity.components
        self.__velocity.components = newVel
        # Attempt to update position
        accepted = self.__positions.startUpdate(newPos)
        for (i,b) in enumerate(accepted):
            if not b:
                self.__velocity.reverseComponent(i)
        if not self.__positions.inBounds:
            raise(NotImplementedError("Case of open search bounds not implemented."))
        self._log_update_attempt(newPos,newVel,accepted)
        self._setup_calculations(self.__calc_manager)
        self.__updating = True
        return True
    
    @abstractmethod
    def _setup_calculations(self, calcManager):
        """
        Prepare for calculations, and pass calculation tasks to calcManager.
        
        ABSTRACT: must override.
        
        Overriding method should adjust any internal states required for calculations
        to be performed. Calculations should then be passed (along with ALL required
        data as arguments) to the calcManager to be run in parallel. Be sure to store
        async_result objects from calcManager to check for completion of calculations.
        
        If overriding method's computations are low-weight and do not require
        parallelization, overhead of the calcManager can be bypassed by performing
        the calculations here, and storing results. In this case, be sure to override
        self.calculationsFinished() to always return True.
        
        Parameters
        ----------
        calcManager : CalculationManager
            The calculation manager responsible for launching calculations in parallel.
        """
        pass
        
    def finishUpdate(self):
        """ 
        Finalize and accept the update, and set new fitness.
        """
        if not self.readyForFinishUpdate:
            if not self.updating:
                raise(RuntimeError("No active update to finish for agent {}.".format(self.id)))
            else:
                raise(RuntimeError("Calculations incomplete for agent {} update.".format(self.id)))
        fitness = self._cleanup_calculations(self.__calc_manager)
        if self.hasStep(0):
            oldBest = self.__positions.bestFitness
            # oldBest listed first below to favor maintaining current best position during tie.
            betterFit = FITNESS_SELECTOR.betterFitness(oldBest, fitness)
            if betterFit == oldBest:
                success = self.__positions.finishUpdate(fitness)
            else:
                success = self.__positions.finishUpdate(fitness, newBest=True)
        else:
            success = self.__positions.finishUpdate(fitness, newBest=True)
        self.__updating = False
        self.__next_step += 1
        self._log_step()
        return success
    
    @abstractmethod
    def _cleanup_calculations(self, calcManager):
        """
        Perform any state cleanup following calculation completion and return fitness.
        
        ABSTRACT: Derived types must override.
        
        Overriding method should undo any state changes completed for calculations, 
        and return to a stable state. In addition to restoring a stable state,
        the fitness should be determined from the results of calculations passed
        to the calculation manager.
        
        If override of self._setup_calculations() performed calculations without
        calculation manager, and handled state cleanup, this method can simply
        return the fitness value determined in earlier call.
        
        Parameters
        ----------
        calcManager : CalculationManager
            The manager to which earlier calculations were passed.
        
        Returns
        -------
        fitness : real, numeric
            The new fitness value.
        """
        pass
    
    def cancelUpdate(self):
        """
        Abort the current update and restore to previous state.
        
        Both position and velocity are restored. Step is not incremented.
        
        Calculations submitted to the calculation manager will still be
        completed. To avoid backup in calculation manager, use this utility
        sparingly.
        """
        if not self.updating:
            raise(RuntimeError("No active update pending."))
        self.__positions.cancelUpdate()
        self.__velocity.components = self.__old_velocity.components
        self._unset_calculations()
    
    @abstractmethod
    def _unset_calculations(self):
        """
        Restore state to pre-update conditions.
        
        ABSTRACT: Derived types must override.
        
        Overriding method should undo any state changes completed for calculations, 
        and return to the last stable state.
        
        Presently, it is impossible to cancel calculations submitted to calculation
        manager.
        """
        pass
    
    @property
    def variableSet(self):
        return self.__positions.variableSet
    
    @property
    def stablePoint(self):
        return self.__positions.stablePoint
    
    @property
    def bestPoint(self):
        return self.__positions.bestPoint
    
    @property
    def stableVelocity(self):
        if not self.updating:
            return deepcopy(self.__velocity)
        else:
            return deepcopy(self.__old_velocity)
    
    @property
    def velocity(self):
        return deepcopy(self.__velocity)
    
    def _start_logs(self):
        self.__pso_vals_fname = self.root/"psoValues.csv"
        self.__tru_vals_fname = self.root/"trueValues.csv"
        self.__velocity_fname = self.root/"velocities.csv"
        self.__updates_fname = self.root/"updateAttempts.csv"
        
        # PsoValue data
        lbl = ["step"]
        lbl.append("bestStep")
        lbl.append("fitness")
        for v in self.variableSet.labels:
            lbl.append(v)
        writeCsvLine(self.__pso_vals_fname,lbl,'w')
        
        # TrueValue Data - use same labels as PsoValue data
        writeCsvLine(self.__tru_vals_fname,lbl,'w')
        
        # Velocity data
        lbl = ["step"]
        for v in self.variableSet.labels:
            lbl.append(v)
        writeCsvLine(self.__velocity_fname,lbl,'w')
        
        # Update attempts
        lbl = ["step"]
        lbl.append("datatype")
        for v in self.variableSet.labels:
            lbl.append(v)
        writeCsvLine(self.__updates_fname,lbl,'w')
    
    def _log_step(self):
        # PsoValue data
        dat = [self.lastStep]
        dat.append(self.__positions.bestStep)
        dat.append(self.__positions.stableFitness)
        for p in self.__positions.stablePsoValues:
            dat.append(p)
        writeCsvLine(self.__pso_vals_fname,dat,'a')
        #TrueValue Data
        dat = [self.lastStep]
        dat.append(self.__positions.bestStep)
        dat.append(self.__positions.stableFitness)
        for p in self.__positions.trueValues:
            dat.append(p)
        writeCsvLine(self.__tru_vals_fname,dat,'a')
        # Velocity data
        dat = [self.lastStep]
        for v in self.__velocity.components:
            dat.append(v)
        writeCsvLine(self.__velocity_fname,dat,'a')
    
    def _log_update_attempt(self, position, velocity, accepted):
        step = self.nextStep
        datatypes = ['POSITION','VELOCITY','ACCEPTED']
        datasources= [position, velocity, accepted]
        for (lbl,src) in zip(datatypes,datasources):
            dat = [step,lbl]
            for p in src:
                dat.append(p)
            writeCsvLine(self.__updates_fname,dat,'a')
    
    def __str__(self):
        s = "Agent {}:\n\tStable: {}\n\tBest: {}\n\tVelocity: {}"
        s = s.format(self.id, self.stablePoint, self.bestPoint, self.stableVelocity)
        return s

class FunctionAgent(Agent):
    def __init__(self, fn, *args, **kwargs):
        self.__function = fn
        self.__last_launch = None
        
        super().__init__(*args, **kwargs)
    
    @property
    def calculationFinished(self):
        """
        Returns True if the update calculations are complete.
        """
        if self.__last_launch is None:
            return False
        return self.__last_launch.ready()
        
    def _setup_calculations(self, calcManager):
        """
        Prepare for calculations, and pass to the calculation manager.
        """
        vals = self.variableSet.trueValues
        self.__last_launch = calcManager.addTask(self.__function, [vals])
        
    def _cleanup_calculations(self, calcManager):
        """
        Read fitness from calculation result and return. 
        """
        fitness = self.__last_launch.get()
        self.__last_launch = None
        return fitness
    
    def _unset_calculations(self):
        self.__last_launch = None

