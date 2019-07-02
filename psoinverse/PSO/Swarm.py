import numpy as np
import networkx as nx
#from functools import total_ordering
import time
import traceback
import os
import shutil

# moving Point class to alternate module, splitting functionality
# This should allow identical functionality throughout this module
#  because functionality was fully recreated in SearchSpace Module
from .SearchSpace import Point #SimulationPoint as Point
from .SearchSpace import SearchBounds
from .Integrators import Integrator

# Helper function for debugging - just wraps the process of spitting out a string to a file
def debug(line):
    #with open("debug.out", 'a') as f:
    #    f.write("{}\n".format(line))
    print("DEBUG: {}".format(line))

def output(line):
    from datetime import datetime
    with open("runtime.out", 'a') as f:
        f.write("{} :: {}\n".format(datetime.now().strftime("%d/%m/%y %H:%M:%S"), line))

def log_exception():
    debug("Exception raised: {}".format(traceback.format_exc()))

# =========================================================
#
# Swarm Class
#
# =========================================================
class Swarm(object):
    def __init__(self, graph, agents, integrator):
        assert len(graph) == len(agents), "Number of nodes in graph doesn't match number of agents"
        
        self.Agents = agents
        self.Graph = graph
        self.integrator = integrator
        self.stepsTaken = 0
        
        # record-keeping
        self.History, self.Best = [], []
        
        # switch
        self.SeekMax = integrator.seekMax
    
    def step(self):
        for i in range(len(self.Agents)):
            neighbors = [self.Agents[a] for a in self.Graph.neighbors(i)]
            self.Agents[i].update(neighbors, self.integrator)
        self.stepsTaken += 1
    
    def printState(self):
        for a in self.Agents:
            print("\n{}".format(a))
    
    # Return global best of whole swarm
    def get_gbest(self):
        gbest_pt = self.Agents[0].PBest
        best_agent = self.Agents[0]

        for a in self.Agents:
            if self.SeekMax and a.PBest > gbest_pt:
                gbest_pt = a.PBest
                best_agent = a
            elif not self.SeekMax and a.PBest < gbest_pt:
                gbest_pt = a.PBest
                best_agent = a
            else:
                pass

        return best_agent

    def print_neighbor_graph(self, fname="SwarmNetwork.png"):
        #import networkx as nx
        import pylab as py

        g = nx.Graph()
        [g.add_node(a.id) for a in self.Agents]
        [[g.add_edge(a.id, b.id) for b in a.neighbors] for a in self.Agents]

        # draw_graphviz doesn't work in networkx 1.11. Replace with 3 lines below.
        #nx.draw_graphviz(g)
        from networkx.drawing.nx_agraph import graphviz_layout
        pos = graphviz_layout(g)
        nx.draw(g, pos)

        py.savefig(fname)

    def write_output(self, outputbasedir):
        from copy import deepcopy

        # Maintain a history of the entire set of coordinates and fitnesses for all agents
        #self.History.append([list(agent.get_coords()) + [agent.Location.Fitness] for agent in self.Agents])

        # Maintain a list of the best agent at each step.
        # First, find the current best in the whole swarm.
        gbest_agent = self.get_gbest()
        # Only add to the list of Best coordinates if the gbest has been updated
        if gbest_agent.Location == gbest_agent.PBest or len(self.Best)<1:
            self.Best.append(deepcopy(gbest_agent.PBest))
            # Copy the simulation files to a gbest location
            gbest_agent.PBest.copy_simulations_to(outputbasedir+"GBest_step{}/".format(gbest_agent.steps))

        # Output all agent data
        for a in self.Agents:
            filename=outputbasedir+"History_Agent{}.dat".format(a.id)
            if a.steps <= 1:
                f=open(filename,'w')
                f.write("# ncoords = {}\n".format(len(a.Location.keys)))
                f.write("# step, fitness, {}, {}\n".format(", ".join(s for s in a.Location.keys), ", ".join("V_{}".format(s) for s in a.Location.keys)))
            else:
                f=open(filename,'a')
            f.write("{} {} {} {}\n".format(
                                    a.steps, a.Location.Fitness,
                                    " ".join(repr(s) for s in a.get_coords()),
                                    " ".join(repr(s) for s in a.Velocity)
                                    ))
            f.close()

        filename=outputbasedir+"History_GBest.dat"
        if gbest_agent.steps <= 1:
            f=open(filename,'w')
            f.write("# ncoords = {}\n".format(len(a.Location.keys)))
            f.write("# step, fitness, {}, AgentID\n".format(", ".join(s for s in a.Location.keys)))
        else:
            f=open(filename,'a')
        f.write("{} {} {} {}\n".format(
                                gbest_agent.steps, gbest_agent.PBest.Fitness,
                                " ".join(repr(s) for s in gbest_agent.PBest.get_scaled_coords()),
                                gbest_agent.id
                                ))
        f.close()

# =========================================================
#
# Agent Class (Abstract)
#
# =========================================================
class Agent(object):
    # STATIC MEMBERS
    NextID = 0  # Static counter for generating unique agent IDs

    def __init__(self, boundaries, **kwargs): #v0=None, spawnRange=None, useScale=None):
        """
        Constructor for Agent base class.
        
        Parameters
        ----------
        boundaries : SearchSpace.SearchBounds or sub-class
            The desired search space.
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
        self.id = Agent.NextID
        Agent.NextID += 1 # Update the static counter to deliver unique IDs

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

# =========================================================
#
# FunctionAgent Class (Derives From: Agent)
#
# =========================================================
class FunctionAgent(Agent):
    def __init__(self, fn, **kwargs):
        self.Function = fn
        
        super().__init__(**kwargs)
        
        self.evaluate()

    def evaluate(self):
        self.Location.Fitness = self.Function(self.get_coords())

        return super().evaluate()



