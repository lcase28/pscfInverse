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
from .SearchSpace import SimulationPoint as Point
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

    def __init__(self, boundaries, v0=None):
        """
            v0 is the initial velocity scale; all velocities will be randomized [-v0, v0]
              If v0 is None, the initial velocity will be set to 1/10 of the boundary range in each direction
            boundaries are the bounds of the search domain and should be a numpy array with shape (dims, 2).
        """
        # Initialize to zero steps taken
        self.steps = 0

        # Grab a unique ID for this agent
        self.id = Agent.NextID
        Agent.NextID += 1 # Update the static counter to deliver unique IDs

        #debug("Agent class init, ID = {}".format(self.id))

        # Store the seach boundaries and generate scale
        self.boundaries = boundaries
        scale = np.array([b[1] - b[0] for b in boundaries])

        #print "AGENT {}, boundaries = {}".format(self.id, scale)

        # The coords used in the update are scaled by the range of the search domain
        # Generate a Point object to store and manipulate the coordinate of the agent in phase space
        dims = len(boundaries)
        self.Location = Point()
        self.Location.Coords = (np.random.rand(dims) * (boundaries[:, 1] - boundaries[:, 0]) + boundaries[:, 0]) / scale
        self.Location.Scale = scale
        #print "AGENT {}, initial coordinates = {}, scaled coords = {} ".format(self.id, self.Location.Coords, self.Location.get_scaled_coords())
        self.Location.Fitness = -1e8

        # Randomize velocity uniformly on [-v0, v0]. Velocity is added to scaled coordinates.
        if v0 is None:
            self.Velocity = 2. * (np.random.rand(dims) - 0.5)
            self.Velocity *= 0.1 * np.array([b[1] - b[0] for b in boundaries]) / scale
        else:
            self.Velocity = 2 * v0 * (np.random.rand(dims) - 0.5) / scale

        #print "AGENT {}, initial velocity = {} ".format(self.id, self.Velocity)

        # Store PBest - best Point() location discovered in all history
        self.PBest = Point()
        self.PBest.fill_from(self.Location)
        #assert self.PBest.Fitness <  -1e7, "pbest fitness is too positive."

        # Initialize empty neighbor list and set fitness to zero. Will be maximized positive.
        self.neighbors = []
        #self.Location.Fitness = 0  # Do not do this - in some cases the fitness metric can be any sign and we don't want to pin the initial GBest/PBest to zero

        #debug("Agent {}: coords   = {}".format(self.id, self.get_coords()))
        #debug("Agent {}: velocity = {}".format(self.id, self.Velocity))
        #debug("Agent {}: fitness  = {}".format(self.id, self.Location.Fitness))

    def get_coords(self):
        return self.Location.get_scaled_coords()

    def get_nbest(self, neighbors):
        # Find the best among the neighbors (including self)
        best = self.PBest
        for n in neighbors:
            if n.PBest > best:
                best = n.PBest
        return best

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
            for i, out_of_bounds in enumerate(
                    np.logical_or(scaled_attempt < self.boundaries[:, 0], scaled_attempt > self.boundaries[:, 1])):
                if out_of_bounds:
                    print("REFLECT AGENT {}, V COMPONENT {}".format(self.id,i))
                    print("OLD POSITION = {}".format(self.Location.get_scaled_coords()))
                    print("TRIAL POSITION = {}".format(scaled_attempt))
                    #print
                    #print
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

    ## Disconnect self and a neighbor 'a' from each other
    #def disconnect(self, a):
    #    self.neighbors.remove(a)
    #    a.neighbors.remove(self)

    ## Connect self and a neighbor 'a' to each other
    #def connect(self, a):
    #    if a not in self.neighbors:
    #        self.neighbors.append(a)
    #    if self not in a.neighbors:
    #        a.neighbors.append(self)

    # Derived classes must override, update fitness, and call base using super calls to generate the MRO.
    def evaluate(self):
        # Now that fitness has been updated, compare to PBest
        if self.Location > self.PBest:
            # Re-use the old PBest simulation record
            #sim_tmp = {k: v for k, v in self.PBest.simulations.items()}

            self.PBest.fill_from(self.Location)

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

        super(FunctionAgent, self).__init__(**kwargs)

        self.PBest.Fitness = -1

        self.evaluate()

    def evaluate(self):
        self.Location.Fitness = self.Function(self.get_coords())

        return super(FunctionAgent, self).evaluate()



