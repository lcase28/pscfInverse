"""Module contains core classes for Agent implementation."""

# Third-party imports
import numpy as np
import traceback
import os
import shutil
from abc import ABC, abstractmethod

# Project imports
from .SearchSpace import Point #SimulationPoint as Point
from .SearchSpace import SearchBounds

class Agent(ABC):
    """Abstract base class for PSO Agents."""
    
    __NextID = 0  # Static counter for generating unique agent IDs

    def __init__(self, boundaries, **kwargs):
        """
        Constructor for Agent base class.
        
        Parameters
        ----------
        boundaries : SearchSpace.SearchBounds or sub-class
            The desired search space.
        
        Keyword Arguments
        -----------------
        v0 : 1D list-like, numerical values, optional
            Initial velocity scale; all velocities will be
            randomized [-v0, v0].
            If v0 is None, the initial velocity will be set
            to 1/10 of the boundary range in each direction.
        spawnRange : SearchSpace.SearchBounds or sub-class, optional
            The range in which to spawn (initialize) the agent.
        useScale : 1D list-like, numerical values, optional
            The scaling to use on the Agent's positions.
        initPosition : SearchSpace.Point or subclass, or 1D list-like, optional
            A forced initial position for Agent. Intended for use
            in test and validation cases. spawnRange will be ignored.
        initVelocity : 1D list-like, optional
            A forced initial velocity for Agent. Intended for use
            in test and validation cases. v0 will be ignored.
        """
        # Initialize to zero steps taken
        self.steps = 0

        # Grab a unique ID for this agent
        self.id = Agent.__NextID
        Agent.__NextID += 1 # Update the static counter to deliver unique IDs

        self.boundaries = boundaries
        
        # Parse kwargs
        useScale = kwargs.get("useScale",None)
        spawnRange = kwargs.get("spawnRange",None)
        v0 = kwargs.get("v0", None)
        self.seekMax = kwargs.get("seekMax",True)
        # optional arguments to enable fixed-starting position initialization for debugging
        # and unit tests -- Not intended for general use.
        initPosition = kwargs.get("position", None)
        initVelocity = kwargs.get("velocity",None)
        
        # determine scaling of point
        scale = self.__choose_scaling_source(boundaries, spawnRange, useScale)
        
        # Initialize Position
        if initPosition is None:
            # general case
            rng = self.__choose_spawn_range(boundaries, spawnRange)
            dims = len(self.boundaries.upper)
            LB = self.boundaries.lower
            initCoords = (np.random.rand(dims) * rng + LB) / scale
            self.Location = Point(Coords=initCoords, Fitness= -1e8, Scale=scale)
        else:
            # debug / testing case
            if isinstance(initPosition, Point):
                self.location = Point()
                self.location.fill_from(initPosition)
            else:
                try:
                    initCoords = np.array(initPosition)
                except:
                    raise
                if not len(initCoords) == len(self.boundaries.upper):
                    raise(ValueError("initPosition does not match dimension of boundaries"))
                
                self.Location = Point(Coords=initCoords, Fitness=-1e8, Scale=scale)
        
        # Initialize Velocity
        if initVelocity is None:
            # general case
            if v0 is None:
                self.Velocity = 2. * (np.random.rand(dims) - 0.5)
                self.Velocity *= 0.1 * rng / scale
            else:
                self.Velocity = 2 * v0 * (np.random.rand(dims) - 0.5) / scale
        else:
            # debug / testing case
            try:
                self.Velocity = np.array(initVelocity)
            except:
                raise
            
            if not len(self.Velocity) == len(self.boundaries.upper):
                raise(ValueError("initVelocity does not match dimension of boundaries"))

        # Store PBest - best Point() location discovered in all history
        self.PBest = Point()
        self.PBest.fill_from(self.Location)
        
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
            from copy import deepcopy
            vtemp = deepcopy(self.Velocity)
            ltemp = deepcopy(self.Location.Coords)

            self.Velocity = new_velocity
            #output("Agent {}: V = {}".format(self.id, self.Velocity))
            self.Location.Coords = attempt.Coords

            # Fitness is incremented during Evaluate() since we use super() calls to generate the MRO
            self.Location.Fitness = 0
            # Use self.Evaluate as a validator
            if not self.evaluate():
                debug("Agent {} Resetting!".format(self.id))
                self.Velocity = vtemp
                self.Location.Coords = ltemp
                return False

        return True

    # Derived classes must override, update fitness, and call base using super calls to generate the MRO.
    @abstractmethod
    def evaluate(self):
        # Now that fitness has been updated, compare to PBest
        if self.Location > self.PBest and self.seekMax:
            self.PBest.fill_from(self.Location)
        elif self.Location < self.PBest and not self.seekMax:
            self.PBest.fill_from(self.Location)
        else:
            pass
            # When we store the PBest point, it gets a reference to self.Location.sim
            # If sim_tmp = None (PBest never set), then a new record will be generated automatically in Agent.Evaluate()
            #self.Location.simulations = sim_tmp
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


