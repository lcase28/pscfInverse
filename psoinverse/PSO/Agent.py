"""Module contains core classes for Agent implementation."""

# Third-party imports
from abc import ABC, abstractmethod
from copy import deepcopy
import numpy as np

# Project imports
from .SearchSpace import Point #SimulationPoint as Point
from .SearchSpace import SearchBounds

class Agent(ABC):
    """Abstract base class for PSO Agents."""
    
    __NextID = 0  # Static counter for generating unique agent IDs

    def __init__(self, boundaries, ptSrc, velSrc, \
                    randGen = None, batchRunner=None, **kwargs):
        """
        Constructor for Agent base class.
        
        Parameters
        ----------
        boundaries : SearchSpace.SearchBounds or sub-class
            The desired search space.
        ptSrc : psoinverse.PSO.SearchSpace.Point or derivative
            A point to use as template for agent's point
        velSrc : 1-D list-like
            Template for velocities variable.
        randGen : numpy.random.RandomState, optional
            A RandomState object to use to initialize positions.
            If not included, values in ptSrc and velSrc are used.
        """
        # Initialize to zero steps taken
        self.steps = 0

        # Grab a unique ID for this agent
        self.id = Agent.__NextID
        Agent.__NextID += 1 # Update the static counter to deliver unique IDs

        self.boundaries = boundaries
        self.seekMax = kwargs.get("seekMax",True)
        self.Location = deepcopy(ptSrc)
        self.Velocity = np.array(velSrc)
        self._runner = batchRunner
        
        if randGen is None:
            randomizeValues = False
        else:
            randomizeValues = True
        
        # Initialize Position
        if randomizeValues:
            # general case
            spawnRange = None
            rng = self.__choose_spawn_range(boundaries, spawnRange)
            dims = len(self.boundaries.upper)
            LB = self.boundaries.lower
            initCoords = (randGen.rand(dims) * rng + LB)
            self.Location.Coords = initCoords
        
        # Initialize Velocity
        if randomizeValues:
            self.Velocity = 2. * (randGen.rand(dims) - 0.5)
            self.Velocity *= 0.1 * rng

        # Store PBest - best Point() location discovered in all history
        self.PBest = None
        self._inerror = False
        self._validPbest = False # Pbest not yet set.
        # Prepare for calculations to be run for initial fitness.
        # Swarm will trigger the actual calculation and evaluation.
        self._setup_calculations()
        
    @property
    def validPbest(self):
        return self._validPbest
    
    @property
    def inErrorState(self):
        return self._inerror
        
    def _startErrorState(self):
        self._inerror = True
        
    def _endErrorState(self):
        self._inerror = False
        
    def __choose_scaling_source(self, boundaries, spawnRange, useScale):
        # Select scaling source
        if useScale is None:
            if spawnRange is not None and isinstance(spawnRange, SearchBounds):
                scale = spawnRange.getScale()
                if not len(scale) == len(self.boundaries.upper):
                    raise(ValueError("Dimension mismatch between spawnRange and boundaries"))
            else:
                scale = self.boundaries.getScale()
        else:
            scale = np.array(useScale).astype(float)
            if not len(scale) == len(self.boundaries.upper):
                raise(ValueError("Dimension mismatch between useScale and boundaries"))
        return scale
    
    def __choose_spawn_range(self, boundaries, spawnRange):
        # Select Particle spawn range source
        if spawnRange is not None and isinstance(spawnRange, SearchBounds):
            rng = spawnRange.getRange()
            if not len(rng) == len(self.boundaries.upper):
                raise(ValueError("Dimension mismatch between spawnRange and boundaries"))
        else:
            rng = self.boundaries.getRange()
        return rng
        
    def get_coords(self):
        return self.Location.get_scaled_coords()
    
    def startUpdate(self, neighbors, integrator, acceleration=None):
        self._integrate_motion(neighbors, integrator, acceleration)
        self._setup_calculations()
    
    def _integrate_motion(self, neighbors, integrator, acceleration=None):
        """ Update the agent's Position (and Velocity) based on the current Position and neighbor information.
        """
        ## Update the position according to PSO dynamics. Retain in tmp array for boundary checks / constraints
        self.steps += 1
        attempt = Point()
        attempt.fill_from(self.Location)
        attempt.Coords, new_velocity = integrator.update(self, neighbors)

        # Reflective boundaries
        if self.boundaries is not None:
            scaled_attempt = attempt.get_scaled_coords()
            for i, in_bounds in enumerate(self.boundaries.inBounds(scaled_attempt)):
                if not in_bounds:
                    print("REFLECT AGENT {}, V COMPONENT {}".format(self.id,i))
                    print("OLD POSITION = {}".format(self.Location.get_scaled_coords()))
                    print("TRIAL POSITION = {}".format(scaled_attempt))
                    new_velocity[i] *= -1.0
                    attempt.Coords[i] = self.Location.Coords[i]
                    print("NEW POSITION = {}".format(attempt.get_scaled_coords()))

        # Replace internal Location and Velocity with new values
        vtemp = deepcopy(self.Velocity)
        ltemp = deepcopy(self.Location.Coords)

        self.Velocity = new_velocity
        self.Location.Coords = attempt.Coords

        # Fitness is incremented during Evaluate()
        self.Location.Fitness = 0
        self._endErrorState() # set flag to normal condition.
                                # self.evaluate will switch back
                                # if required, but base class
                                # assumes normal state will be
                                # reached unless otherwise
                                # specified.
        ## If choosing to revert on failure, store new and old position and velocity
        ## presently don't do this, so values not stored.
        return True
    
    # Derived classes must override, perform any steps to prepare for fitness evaluation
    #   then send functions and argument to self._runner for parallel evaluation if using.
    # Function sent to self._runner must not require modifications to the Agent's data.
    #   If Data must be returned, task id returned from self._runner should be stored and
    #   used in self.finishUpdate to retrieve the return value.
    @abstractmethod
    def _setup_calculations(self):
        pass
    
    # Derived classes must override, update fitness, and call base using super calls to generate the MRO.
    def finishUpdate(self):
        """
        Basic wrapper method for process in finishing an update procedure
        """
        self._evaluate_fitness()
        self._evaluate_Pbest()
        
        
    ## Derived classes must override.
    ## Overriding method must retrieve any data produced by self._runner execution of
    ##  methods in order to update the fitness of the Agent.
    @abstractmethod
    def _evaluate_fitness(self):
        """
            Derived classes must override. Overriding method must update
            the Agent's fitness.
            
            If the overriding method is unable to calculate a valid fitness value,
            an arbitrarily high (numpy.inf, if Agent.seekMax=False) 
            or low (numpy.ninf, if Agent.seekMax=True) value should be set
            for location fitness, which will leave the location ignored.
            
            In addition to arbitrarily extreme fitness, overriding method should 
            call private method self._startErrorState() to set an Agent-level
            flag indicating that the current location is invalid.
        """
        pass
        
    def _evaluate_Pbest(self):
        """
            Pbest Behavior:
            When the Agent's current location is valid (no errors, has fitness):
                During initialization: 
                    Pbest is set to current location.
                Otherwise: 
                    Pbest is updated to current location if current
                    location has better fitness value.
            When Agent's current location is invalid (error, no valid fitness):
                During initialization: 
                    Pbest is set to current location,
                    flag is set to return False from Agent.validPbest property.
                Other times, while Agent.validPbest==False: 
                    Pbest is set to current location. 
                    Agent.validPbest==False flag is retained.
                Other times, when Agent.validPbest==True:
                    A valid location has previously been found.
                    This valid Pbest is retained, and current location ignored.
            
            Returns
            -------
            success : bool
                True if evaluation was successful.
                False if an error occurred.
        """
        if not self._validPbest:
            # Edge condition: No valid Pbest yet found.
            if self.PBest is None:
                # Initialization condition
                self.PBest = deepcopy(self.Location)
            else:
                # Searching for initial Pbest
                self.PBest.fill_from(self.Location)
                
            if self.inErrorState:
                self._validPbest = False
            else:
                self._validPbest = True
                
        if not self.inErrorState:
            # Now that fitness has been updated, compare to PBest
            if self.Location > self.PBest and self.seekMax:
                self.PBest.fill_from(self.Location)
            elif self.Location < self.PBest and not self.seekMax:
                self.PBest.fill_from(self.Location)
            successFlag = True
        else:
            successFlag = False
        return successFlag
    
    def __str__(self):
        s = "Agent {}:\n\tScaled: {}\n\tCoords: {}\n\tVelcty: {}\n\tFitnes: {}"
        s = s.format(self.id, self.get_coords(), self.Location.Coords, self.Velocity, self.Location.Fitness)
        return s

class FunctionAgent(Agent):
    def __init__(self, fn, **kwargs):
        self.Function = fn
        
        super().__init__(**kwargs)
        
        self.evaluate()
    
    def _setup_calculations(self):
        """
        Current implementation does not make use of parallelization.
        """
        return True
        
    def _evaluate_fitness(self):
        """
        Batch parallelization is bypassed. Instead, fitness is evaluated here.
        """
        self.Location.Fitness = self.Function(self.get_coords())


