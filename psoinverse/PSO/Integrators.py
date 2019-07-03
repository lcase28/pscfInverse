"""
A set of classes to perform position and velocity updates on PSO Agents.

The Integrators do not assign new positions and velocities, but are meant to be called by
the agents to simply perform the update calculations and return the new positions.

Experimental component to ease the implementation and testing of various update schemes.
"""

# Third Party imports
from abc import ABC, abstractmethod
import numpy as np

# Local Library imports
from .SearchSpace import Point
#from .Swarm import Agent

# Inheriting classes must override method Update
class Integrator(ABC):
    """
    Abstract base class for Integrator objects.
    
    Defines the basic functionality of an Integrator.
    
    Abstract Methods
    ----------------
    update
        Calculates the new position and velocity of a target Agent
    """
    
    def __init__(self, **kwargs):
        """
        Constructor for abstract Integrator Class
        
        Keyword Arguments
        -----------------
        seekMax : Boolean, optional, Default False
            If True, the Integrator considers higher Fitness better
            If False, lower fitness is considered better
            If given value can not be cast as bool, seekMax set to
            default and warning is printed.
        """
        super().__init__()
        
        seekMax = kwargs.get("seekMax",False)
        try:
            self.seekMax = bool(seekMax)
        except(TypeError, ValueError):
            self.seekMax = False
            print("DEBUG: seekMax set to default value of False")
    
    @abstractmethod
    def update(self, Target, Neighbors, **kwargs):
        """
        [Abstract] Calculate and return the updated position and velocity for Target
        
        Parameters
        ----------
        Target : Swarm.Agent or subclass
            The Agent being updated
        Neighbors : 1D List-like containing Agents
            The neighbors of Target
        
        Returns
        -------
        new_position : 1D numpy array
            The new position of Target after update
        new_velocity : 1D numpy array
            The new velocity of Target after update
        """
        return Position, Velocity
    
## "Standard" PSO Update Scheme
class StandardIntegrator(Integrator):
    ## Constructor
    #
    # @param chi Constriction factor value. (Default: 0.729)
    # @param c1 Personal Best weighting factor. (Default: 2.05)
    # @param c2 Neighbor Best weighting factor. (Default: 2.05)
    # @param seekMax (Optional, key-word) Boolean. Default: False
    #                If True, Neighbor best is the MAXIMUM neighbor value.
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
        
        seed = kwargs.get('RandSeed',None)
        self.__randGen = np.random.RandomState(seed)
    
    ## Return updated position and velocity for agent with Neighbors
    #
    # Operates on and returns the scaled (O(1)) coordinates and velocities.
    #
    # Assumes that Agent.Location.Coords and Agent.Velocity are equivallently scaled.
    #
    # @param target The Agent being updated in the call (will not be modified in call)
    # @param Neighbors An iterable set of neighboring agents
    # @param acceleration Additional velocity increment (Optional, key-word)
    # @return new_position The new position of the agent (if accepted) (np.array)
    # @return new_velocity The updated velocity of the agent (np.array)
    #
    def update(self, target, Neighbors, **kwargs):
        nbest = self.get_nbest(Neighbors) # best among neighbors
        # Get data from target
        Position = target.Location
        Velocity = target.Velocity
        pbest = target.PBest
        
        # Random forcing terms
        e1 = self.__randGen.rand(len(Position.Coords))
        e2 = self.__randGen.rand(len(Position.Coords))
        
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
        new_velocity = self.chi * (Velocity + val1 + val2) + acc
        
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
            #if not isinstance(n,Agent):
            #    raise(TypeException("Neighbors List must contain only subclasses of Agent"))
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
    
