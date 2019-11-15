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

    def __init__(self, boundaries, ptSrc, velSrc, randGen = None, **kwargs):
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
        self.evaluate()
        self._inerror = False
        
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
    
    def update(self, neighbors, integrator, acceleration=None):
        """ Update the agent's Position (and Velocity) based on the current Position and neighbor information.

            ::Returns:: the agent's personal best fitness value (not Position)
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

        # Replace internal Location and Velocity with new values
        move = True
        if move:
            vtemp = deepcopy(self.Velocity)
            ltemp = deepcopy(self.Location.Coords)

            self.Velocity = new_velocity
            self.Location.Coords = attempt.Coords

            # Fitness is incremented during Evaluate()
            self.Location.Fitness = 0
            # Use self.Evaluate as a validator
            if not self.evaluate():
                debug("Agent {} Failed!".format(self.id))
                self.Velocity = vtemp
                self.Location.Coords = ltemp
                return False

        return True

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


