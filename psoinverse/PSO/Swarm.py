from copy import deepcopy
import numpy as np
import networkx as nx
import pathlib
#from functools import total_ordering
import time
import traceback
import os
import shutil

# moving Point class to alternate module, splitting functionality
# This should allow identical functionality throughout this module
#  because functionality was fully recreated in SearchSpace Module
from .SearchSpace import Point #SimulationPoint as Point
#from .SearchSpace import SearchBounds
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

def checkPath(root):
    root = root.resolve()
    if root.exists() and root.is_dir():
        return root.resolve(), True
    try:
        if not root.exists():
            root.mkdir(parents=True)
        if not root.is_dir():
            root = root.parent
        return root.resolve(), True
    except (FileNotFoundError, FileExistsError, RuntimeError):
        return None, False

# =========================================================
#
# Swarm Class
#
# =========================================================
class Swarm(object):
    def __init__(self, graph, agents, integrator, root = None, logName = "swarm_log", autoLog = True):
        assert len(graph) == len(agents), "Number of nodes in graph doesn't match number of agents"
        self.Agents = agents
        self.Graph = graph
        self.integrator = integrator
        self.stepsTaken = 0
        # record-keeping
        self.autoLog = autoLog
        if root is None:
            self.root = pathlib.Path.cwd()
            flag = True
        else:
            self.root, flag  = checkPath(root)
        if flag:
            self.logFile = self.root / logName
        else:
            raise(ValueError("Invalid root passed to swarm"))
        #self.History, self.Best = [], []
        # switch
        self.SeekMax = integrator.seekMax
        self._inerror = False
        self.bestAgent = self.get_gbest()
        if self.autoLog:
            self.logStatus(True)
        
    
    @property
    def inErrorState(self):
        return self._inerror
    
    def step(self):
        self.stepsTaken += 1
        for i in range(len(self.Agents)):
            neighbors = [self.Agents[a] for a in self.Graph.neighbors(i)]
            self.Agents[i].update(neighbors, self.integrator)
        if self.autoLog:
            self.logStatus(False)
    
    def printState(self):
        for a in self.Agents:
            print("\n{}".format(a))
    
    # Return global best of whole swarm
    def get_gbest(self):
        """ Returns the agent with the best historical best. None if all in error. """
        gbest_pt = None #self.Agents[0].PBest
        best_agent = None #self.Agents[0]

        for a in self.Agents:
            if not a.inErrorState:
                if best_agent is not None and gbest_pt is not None:
                    if self.SeekMax and a.PBest > gbest_pt:
                        gbest_pt = a.PBest
                        best_agent = a
                    elif not self.SeekMax and a.PBest < gbest_pt:
                        gbest_pt = a.PBest
                        best_agent = a
                else:
                    gbest_pt = a.PBest
                    best_agent = a
                    
        if best_agent is None:
            self._inerror = True

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

    def statusString(self, includeDefinitions = False):
        s = ""
        if includeDefinitions:
            s += "Swarm Characteristics:\n"
            s += "\tNumber Agents : {}\n".format(len(self.Agents))
            s += "\tGraph Type : {}\n".format(self.Graph)
            s += "\tIntegrator : {}\n".format(self.integrator)
            s += "\tSeek Max : {}\n".format(self.SeekMax)
            s += "\tRoot Path : {}\n".format(self.root)
            s += "\tLog File : {}\n".format(self.logFile)
        # Output status
        s += "Swarm Status:\n"
        s += "\tStep Number : {}\n".format(self.stepsTaken)
        s += "\tIn Error State : {}\n".format(self.inErrorState)
        if self.inErrorState:
            s += "\tBest Agent ID : {}\n".format(self.bestAgent)
        else:
            s += "\tBest Agent ID : {}\n".format(self.bestAgent.id)
        s += "\tAgent Data (Current):\n"
        s += "\tID_Num\t\tFitness\t\tPosition\n"
        formstring = "\t{}\t\t{}\t\t{}\n"
        for a in self.Agents:
            posStr = "".join(["{:.4f}, ".format(i) for i in a.Location.Coords])
            s += formstring.format(a.id, a.Location.Fitness, posStr)
        s += "\tAgent Data (P_Best):\n"
        s += "\tID_Num\t\tFitness\t\tPosition\n"
        formstring = "\t{}\t\t{}\t\t{}\n"
        for a in self.Agents:
            posStr = "".join(["{:.4f}, ".format(i) for i in a.PBest.Coords])
            s += formstring.format(a.id, a.PBest.Fitness, posStr)
        return s
    
    def logStatus(self, includeDefinitions = False):
        with self.logFile.open(mode='a') as f:
            f.write(self.statusString(includeDefinitions))
        

