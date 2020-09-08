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
    def __init__(self, graph, agents, integrator, \
                root = None, logName = "swarm_log", autoLog = True):
        """ 
            Constructor for swarm class.
            
            Parameters
            ----------
            graph : networkx.Graph (undirected)
                A graph describing the communication links in the swarm.
                Nodes should be integers corresponding to the Agent IDs.
            agents : list of psoinverse.PSO.Agent.Agent
                A list containing all agents in the swarm. List indexes
                should correspond with the Agent IDs, such that 
                Agent_n can be accessed through agents[n].
            integrator : psoinverse.PSO.Integrators.Integrator child class object
                An object derived from the Integrator abstract base class
                representing the update procedure for the swarming behavior.
            root : pathlib.path, optional.
                A path object pointing to the Swarm's root directory.
                All output from Swarm will be placed in this directory
                (or in sub-directories, if so implemented). If directory
                does not exist, it will be created.
                Default Value is the current working directory
            logName : string, optional
                A valid filename (on user's operating-system) relative to root.
                Whenever told to write to log, swarm will append
                to the file of this name.
            autoLog : bool, optional
                If True, Swarm will log status on creation and after each step.
                If False, log will only be written to at user command.
                Default is True.
            
            Raises
            ------
            ValueError :
                If number of nodes in graph does not match number of agents.
                If root is unable to be resolved or created.
        """
        if not len(graph) == len(agents): 
            raise(ValueError("Number of nodes in graph doesn't match number of agents"))
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
            raise(ValueError("Unable to resolve root in Swarm"))
        self.SeekMax = integrator.seekMax
        self._inerror = False
        self.bestAgent = self.get_gbest()
        if self.autoLog:
            self.logStatus(True)
        
    
    @property
    def inErrorState(self):
        """ 
            Returns true if the swarm is in an error state. 
            
            Most likely error state cause is that all agents
            entered an error state.
        """
        return self._inerror
    
    def step(self):
        self.stepsTaken += 1
        for i in range(len(self.Agents)):
            neighbors = [self.Agents[a] for a in self.Graph.neighbors(i)]
            self.Agents[i].update(neighbors, self.integrator)
        if self.autoLog:
            self.logStatus(False)
    
    def printState(self):
        """ Depreciated. logStatus and statusString should be used. """
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

    def write_output(self, outputbasedir):
        """ Legacy code. Method not currently supported. See logStatus and statusString. """
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
            s += "\tGraph Type : {!s}\n".format(type(self.Graph))
            s += "\tIntegrator : {}\n".format(self.integrator)
            s += "\tSeek Max : {}\n".format(self.SeekMax)
            s += "\tRoot Path : {}\n".format(self.root)
            s += "\tLog File : {}\n".format(self.logFile)
        # Output status
        agentTableHeader = "\tID_Num\tFitness\t\t\tPosition\n"
        agentTableLine = "\t{}\t{:E}\t{}\n"
        s += "Swarm Status:\n"
        s += "\tStep Number : {}\n".format(self.stepsTaken)
        s += "\tIn Error State : {}\n".format(self.inErrorState)
        if self.inErrorState:
            s += "\tBest Agent ID : {}\n".format(self.bestAgent)
        else:
            s += "\tBest Agent ID : {}\n".format(self.bestAgent.id)
        s += "\tAgent Data (Current):\n"
        s += agentTableHeader
        for a in self.Agents:
            posStr = "".join(["{:E}, ".format(i) for i in a.Location.Coords])
            s += agentTableLine.format(a.id, a.Location.Fitness, posStr)
        s += "\tAgent Data (P_Best):\n"
        s += agentTableHeader
        for a in self.Agents:
            posStr = "".join(["{:E}, ".format(i) for i in a.PBest.Coords])
            s += agentTableLine.format(a.id, a.PBest.Fitness, posStr)
        return s
    
    def logStatus(self, includeDefinitions = False):
        with self.logFile.open(mode='a') as f:
            f.write(self.statusString(includeDefinitions))
        

