
# Third Party imports
from abc import ABC, abstractmethod
import numpy as np

# Local Library imports
from psoinverse.pso.core import Point, FITNESS_SELECTOR
from psoinverse.pso.agent import Agent

class Integrator(ABC):
    """
    Abstract base class for Integrator objects.
    
    Defines the basic functionality of an Integrator.
    
    Abstract Methods
    ----------------
    integrate
        Calculates the new position and velocity of a target Agent
    """
    
    def __init__(self, randGen):
        """
        Constructor for abstract Integrator Class
        """
        self.__randGen = randGen
        self.__fit_compare = FITNESS_SELECTOR
    
    @abstractmethod
    def integrate(self, target, neighbors):
        """
        [Abstract] Calculate and return the updated position and velocity for Target.
        
        Derived classes should override this definition and implement an update scheme.
        
        Derived classes should use self.chooseBest when selecting best positions among
        neighbors to ensure consistent comparison schemes.
        
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
        pass
    
    @property
    def randomGenerator(self):
        return self.__randGen
    
    @property
    def fitCompare(self):
        return self.__fit_compare
    
## "Standard" PSO Update Scheme
class StandardIntegrator(Integrator):
    """
    Integrator implementing "standard" pso updates.
    """
    
    def __init__(self, randGen, chi=0.729, c1=2.05, c2=2.05):
        """
        Constructor for StandardIntegrator.
        
        Parameters
        ----------
        randGen : numpy.random.RandomState
            The random number generator.
        chi : float, 0 < chi < 1
            Constriction factor value. (Default: 0.729)
        c1 : float 0 < c1
            Personal Best weighting factor. (Default: 2.05)
        c2 : float, 0 < c2
            Neighbor Best weighting factor. (Default: 2.05)
        """
        super().__init__(randGen)
        self.chi = float(chi)
        self.c1 = float(c1)
        self.c2 = float(c2)
    
    @property
    def chi(self):
        return self.__chi
    
    @chi.setter
    def chi(self, val):
        val = float(val)
        if val <= 0. or val > 1.:
            raise(ValueError("Chi must be in the range (0,1]: gave {}".format(val)))
        self.__chi = val
    
    @property
    def c1(self):
        return self.__c1
    
    @c1.setter
    def c1(self,val):
        val = float(val)
        if val < 0:
            raise(ValueError("C1 must be non-negative: gave {}".format(val)))
        self.__c1 = val
        
    @property
    def c2(self):
        return self.__c2
    
    @c2.setter
    def c2(self,val):
        val = float(val)
        if val < 0:
            raise(ValueError("C2 must be non-negative: gave {}".format(val)))
        self.__c2 = val
    
    def integrate(self, target, neighbors):
        """
        Return updated position and velocity for agent with neighbors.
        
        Assumes that position and velocity components are equivallently scaled.
        
        Parameters
        ----------
        target : Agent
            The Agent being updated in the call (will not be modified in call).
        Neighbors : iterable of Agent
            An iterable set of neighboring agents.
        
        Returns
        -------
        new_position : numpy.array
            The new position of the agent (if accepted).
        new_velocity : numpy.array
            The updated velocity of the agent.
        """
        # get data from target
        step = target.lastStep
        position = target.stablePoint
        velocity = target.stableVelocity
        pbest = target.bestPoint
        # get data from neighbors
        nlist = [ n.bestPointAtStep(step) for n in neighbors ]
        nbest = self.fitCompare.bestPoint(nlist)
        
        # Random forcing terms
        e1 = self.randomGenerator.rand(len(position.components))
        e2 = self.randomGenerator.rand(len(Position.components))
        
        # new Velocity
        val1 = self.c1 * e1 * (pbest.components - position.components)
        val2 = self.c2 * e2 * (nbest.components - position.components)
        new_velocity = self.chi * (velocity.components + val1 + val2)
        
        # new Position
        new_position = position.components + new_velocity
        
        # return results
        return new_position, new_velocity

