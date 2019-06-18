## @package Integrators
# A set of classes to perform position and velocity updates on PSO Agents.
#
# The Integrators do not assign new positions and velocities, but are meant to be called by
# the agents to simply perform the update calculations.
#
# Experimental component to ease the implementation and testing of various update schemes.

# Third Party imports
import numpy as np
from abc import ABC, abstractmethod

# Local Library imports
from .SearchSpace import Point
from .Swarm import Agent

## Abstract Base class for PSO Integrators
#
# Inheriting classes must override method Update
class Integrator(ABC):
    
    ## Constructor
    def __init__(self, **kwargs):
        super().__init__()
        
        seekMax = kwargs.get("seekMax",False)
        try:
            self.seekMax = bool(seekMax)
        except(TypeError, ValueError):
            self.seekMax = False
            print("DEBUG: seekMax set to default value of False")
    
    ## Return updated Position and Velocity for agent with Neighbors
    #
    # @param Position The Current Position of the Agent (SearchSpace.Point)
    # @param Velocity The velocity of the Agent (numpy.array or array-like)
    # @param Neighbors An iterable set of neighboring agents
    # @return new_position The new position of the agent (if accepted)
    # @return new_velocity The updated velocity of the agent
    #
    # @throws ValueError if Position, Velocity, or neighbor Positions do not
    #                    agree in dimension
    @abstractmethod
    def update(self, Position, Velocity, Neighbors, **kwargs):
        return Position, Velocity
    
    
## "Standard" PSO Update Scheme
class StandardIntegrator(Integrator):
    ## Constructor
    def __init__(self, chi=None, c1=None, c2=None, **kwargs):
        super().__init__(**kwargs)
        try:
            self.chi = float(chi)
        except(TypeError, ValueError):
            self.chi = 0.729
        try:
            self.c1 = float(c1)
        except(TypeError, ValueError):
            self.c1 = 2.05
        try:
            self.c2 = float(c2)
        except(TypeError, ValueError):
            self.c2 = 2.05
    
    ## Return updated position and velocity for agent with Neighbors
    #
    # @param Position The Current Position of the Agent (SearchSpace.Point)
    # @param Velocity The velocity of the Agent (numpy.array or array-like)
    # @param Neighbors An iterable set of neighboring agents
    # @return new_position The new position of the agent (if accepted)
    # @return new_velocity The updated velocity of the agent
    #
    # @throws ValueError if Position, Velocity, or neighbor Positions do not
    #                    agree in dimension
    def update(self, target, Neighbors, **kwargs):
        nbest = self.get_nbest(Neighbors) # best among neighbors
        
        # Get data from target
        Position = target.Location
        Velocity = target.Velocity
        pbest = target.PBest
        
        # Random forcing terms
        e1 = np.random.rand(len(Position.Coords))
        e2 = np.random.rand(len(Position.Coords))
        
        # Inertia??
        res = kwargs.get("acceleration",None) # Temporary result storage
        if res is not None:
            res = np.array(res)
            acc = acc.astype(float)
        else:
            acc = np.zeros_like(e1)
        
        # new Velocity
        val1 = self.c1 * e1 * (pbest.Coords - Position.Coords)
        val2 = self.c2 * e2 * (nbest.Coords - Position.Coords)
        new_velocity = self.chi * (Velocity + val1 + val2)
        
        # new Position
        new_position = Position.Coords + new_velocity
        
        # return results
        return new_position, new_velocity
    
    ## Return best position seen by any neighbor.
    #
    # @param Neighbors An iterable set of Agent objects representing neighbors
    def get_nbest(self, Neighbors):
        nbest = None
        for n in Neighbors:
            if not isinstance(n,Agent):
                raise(TypeException("Neighbors List must contain only subclasses of Agent"))
            if nbest is None:
                nbest = n.Location
            else:
                if self.seekMax:
                    if n.Location > nbest:
                        nbest = n.Location
                else:
                    if n.Location < nbest:
                        nbest = n.Location
        
        return nbest
    
