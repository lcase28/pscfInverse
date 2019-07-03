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
from .Agent import Agent

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



