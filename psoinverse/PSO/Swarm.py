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
from psoinverse.util.stringTools import str_to_num, wordsGenerator
from psoinverse.util.parallelBatches import LocalBatchRunner

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
                root = None, batchRunner = None, \
                logName = "swarm_log", autoLog = True):
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
        tempAgent = self.Agents[0]
        self.dim = len(tempAgent.Location.Coords)
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
        if batchRunner is None:
            self._runner = LocalBatchRunner(1)
        else:
            self._runner = batchRunner
        self._finish_step(is_startup=True)
        
    
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
        self._start_step()
        self._finish_step()
    
    def _start_step(self):
        for i in range(len(self.Agents)):
            neighbors = [self.Agents[a] for a in self.Graph.neighbors(i)]
            self.Agents[i].startUpdate(neighbors, self.integrator)
        
    def _finish_step(self, is_startup=False):
        self._runner.runBatch()
        for i in range(len(self.Agents)):
            self.Agents[i].finishUpdate()
        self.bestAgent = self.get_gbest()
        if self.autoLog:
            self.logStatus(is_startup)
    
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
            s += "Swarm_Characteristics\n"
            s += "\tNumber_Agents\t\t{}\n".format(len(self.Agents))
            s += "\tGraph_Type\t\t{!s}\n".format(type(self.Graph))
            s += "\tIntegrator\t\t{}\n".format(self.integrator)
            s += "\tDimensions\t\t{}\n".format(self.dim)
            s += "\tSeek_Max\t\t{}\n".format(self.SeekMax)
            s += "\tRoot_Path\t\t{}\n".format(self.root)
            s += "\tLog_File\t\t{}\n".format(self.logFile)
            s += "END_Swarm_Characteristics\n"
        # Start New Step update
        s += "\nStart_Step_Number\t{}\n".format(self.stepsTaken)
        # Output status
        s += "\tSwarm_Status\n"
        s += "\t\tIn_Error_State\t\t{}\n".format(self.inErrorState)
        if self.inErrorState:
            s += "\t\tBest_Agent_ID\t\t{}\n".format(self.bestAgent)
        else:
            s += "\t\tBest_Agent_ID\t\t{}\n".format(self.bestAgent.id)
        s += "\tEND_Swarm_Status\n"
        # Current Positions
        positTableHeader = "\t\tID_Num\tFitness\t\tPosition\n"
        positTableLine = "\t\t{}\t{: E}\t{}\n"
        listForm = "{: E}\t"
        s += "\tCurrent_Positions\n"
        s += positTableHeader
        for a in self.Agents:
            posStr = "".join([listForm.format(i) for i in a.Location.Coords])
            s += positTableLine.format(a.id, a.Location.Fitness, posStr)
        s += "\tEND_Current_Positions\n"
        # Best Positions
        s += "\tBest_Positions\n"
        s += positTableHeader
        for a in self.Agents:
            posStr = "".join([listForm.format(i) for i in a.PBest.Coords])
            s += positTableLine.format(a.id, a.PBest.Fitness, posStr)
        s += "\tEND_Best_Positions\n"
        # Velocities
        agentTableHeader = "\t\tID_Num\tVelocity\n"
        agentTableLine = "\t\t{}\t{}\n"
        s += "\tVelocities\n"
        s += agentTableHeader
        for a in self.Agents:
            posStr = "".join([listForm.format(i) for i in a.Velocity])
            s += agentTableLine.format(a.id, posStr)
        s += "\tEND_Velocities\n"
        s += "END_Step\n"
        
        return s
    
    def logStatus(self, includeDefinitions = False):
        with self.logFile.open(mode='a') as f:
            f.write(self.statusString(includeDefinitions))
        
class LogReadingException(ValueError):
    """ Special Class for errors while parsing a SwarmLog file """
    
    def __init__(self, expected_key, got_key):
        self.expectedKey = expected_key
        self.gotKey = got_key
        ErrorMsg = "Error reading Log File: Expected key {}; Got {}"
        super().__init__(ErrorMsg.format(self.expectedKey,self.gotKey))

class SwarmLog(object):
    """
    Class to read a swarm's log file for analysis.
    """
    
    def __init__(self, fname):
        self.fileSource = fname
        initFlag = False
        self.steps = []
        self.nsteps = 0
        with fname.open(mode='r') as f:
            words = wordsGenerator(f)
            for word in words:
                key = word
                if not initFlag:
                    if key == 'Swarm_Characteristics':
                        self._readCharacteristics(words)
                        initFlag = True
                    else:
                        raise(LogReadingException("Swarm_Characteristics",key))
                else:
                    if key == 'Start_Step_Number':
                        self._readStep(words)
                    else:
                        print(self.nsteps)
                        raise(LogReadingException("Start_Step_Number",key))
    
    def _readCharacteristics(self, words):
        key = next(words)
        while not key == 'END_Swarm_Characteristics':
            if key == 'Number_Agents':
                self.nagent = str_to_num(next(words))
            elif key == 'Graph_Type':
                data = next(words)
                data = data.join(next(words))
                self.graphtype = data.join(next(words))
            elif key == 'Integrator':
                data = next(words)
                self.integrator = data.join([next(words) for i in range(3)])
            elif key == 'Dimensions':
                self.dim = str_to_num(next(words))
            elif key == 'Seek_Max':
                self.seekmax = bool(next(words))
            elif key == 'Root_Path':
                self.rootpath = pathlib.Path(next(words))
            elif key == 'Log_File':
                self.logfile = pathlib.Path(next(words))
            else:
                raise(ValueError("Unexpected keyword in Characteristics section: {}".format(key)))
            key = next(words)
    
    def _readStep(self, words):
        # read step number
        num = str_to_num(next(words))
        if not num == self.nsteps:
            raise(LogReadingException(self.nsteps, num))
        
        step = self.StepData(num, self.dim, self.nagent, words)
        self.steps.append(step)
        
        self.nsteps = self.nsteps + 1
    
    def tocsv(self, swarmFile, agentFile, root=None):
        if root is None:
            root = pathlib.Path.cwd()
        
        self.writeSwarmCSV(root/swarmFile)
        
        self.writeAgentCSV(root/agentFile)
        
    def writeAgentCSV(self,fname):
        header = "step,agent,cfit,bfit"
        lineForm = "\n{},{},{:.6E},{:.6E}"
        dat = [0, 0, 0.0, 0.0]
        for i in range(self.dim):
            lineForm += ",{:.6E},{:.6E},{:.6E}"
            header += ",cpos_{},bpos_{},vel_{}".format(i,i,i)
            dat.append(0.0)
            dat.append(0.0)
            dat.append(0.0)
        with fname.open(mode='w') as f:
            f.write(header)
            for i in range(self.nsteps):
                step = self.steps[i]
                dat[0] = i
                for j in range(self.nagent):
                    dat[1] = j
                    cfit, cpos = step.getAgentCurrent(j)
                    bfit, bpos = step.getAgentBest(j)
                    vel = step.getAgentVelocity(j)
                    dat[2] = cfit
                    dat[3] = bfit
                    dat[4:3:] = cpos
                    dat[5:3:] = bpos
                    dat[6:3:] = vel
                    f.write(lineForm.format(*dat))
    
    def writeSwarmCSV(self,fname):
        header = "step,inError,bagent,bfit"
        lineform = "\n{},{},{},{:.6E}"
        dat = [0, 0, 0, 0.0]
        for i in range(self.dim):
            lineform += ",{:.6E}"
            header += ",bpos_{}".format(i)
            dat.append(0.0)
        with fname.open(mode='w') as f:
            f.write(header)
            for i in range(self.nsteps):
                step = self.steps[i]
                dat[0] = i
                if step.inError:
                    dat[1] = 1
                else:
                    dat[1] = 0
                dat[2] = step.bestAgent
                bfit, bpos = step.getAgentBest(step.bestAgent)
                dat[3] = bfit
                dat[4:] = bpos
                f.write(lineform.format(*dat))
                
                    
    
    class StepData(object):
        """ Helper class which reads and stores the data for a single step """
        
        def __init__(self, idNum, dim, nagent, words): #curPos, bestPos, curVel, errState, bestAgent):
            self.id = idNum
            self.dim = dim
            self.nagent = nagent
            
            # First section is General Swarm info
            key = next(words)
            if key == "Swarm_Status":
                self._readSwarmStatus(words)
            else:
                raise(LogReadingException("Swarm_Status",key))
            
            # Then read all agent data
            key = next(words)
            if key == "Current_Positions":
                self._readCurrentPositions(words)
            else:
                raise(LogReadingException("Current_Positions",key))
            
            key = next(words)
            if key == "Best_Positions":
                self._readBestPositions(words)
            else:
                raise(LogReadingException("Best_Positions",key))
            
            key = next(words)
            if key == "Velocities":
                self._readVelocities(words)
            else:
                raise(LogReadingException("Velocities",key))
            
            # check for end of step data
            key = next(words)
            if not key == "END_Step":
                raise(LogReadingException("END_Step",key))
        
        def getAgentCurrent(self, agentID):
            if agentID < self.nagent:
                dat = self.agentData[agentID]
            else:
                raise(ValueError("Invalid agentID, {}. Expected below {}.".format(agentID, self.nagent)))
            
            return dat.get("cfit"), dat.get("cpos")
        
        def getAgentBest(self, agentID):
            if agentID < self.nagent:
                dat = self.agentData[agentID]
            else:
                raise(ValueError("Invalid agentID, {}. Expected below {}.".format(agentID, self.nagent)))
            
            return dat.get("bfit"), dat.get("bpos")
        
        def getAgentVelocity(self, agentID):
            if agentID < self.nagent:
                dat = self.agentData[agentID]
            else:
                raise(ValueError("Invalid agentID, {}. Expected below {}.".format(agentID, self.nagent)))
            
            return dat.get("vel")
            
            
        
        # Private Functions for file parsing
        
        def _readSwarmStatus(self, words):
            # First read error state value
            key = next(words)
            if key == "In_Error_State":
                self.inError = bool(next(words))
            else:
                raise(LogReadingException("In_Error_State",key))
            
            # Next read best agent
            key = next(words)
            if key == "Best_Agent_ID":
                self.bestAgent = str_to_num(next(words))
            else:
                raise(LogReadingException("Best_Agent_ID",key))
            
            # check for end flag
            key = next(words)
            if not key == "END_Swarm_Status":
                raise(LogReadingException("END_SwarmStatus",key))
        
        def _readPositionTableHeaders(self, words):
            ExpHeads = ["ID_Num","Fitness","Position"]
            for expKey in ExpHeads:
                key = next(words)
                if not key == expKey:
                    raise(LogReadingException(expKey,key))
        
        def _readPositionLine(self, words):
            # read ID number
            num = str_to_num(next(words))
            
            # read Fitness
            fitness = str_to_num(next(words))
            
            # read Position
            pos = np.zeros(self.dim)
            for i in range(self.dim):
                pos[i] = str_to_num(next(words))
            
            return num, fitness, pos
            
        def _readVelocityTableHeaders(self,words):
            ExpHeads = ["ID_Num","Velocity"]
            for expKey in ExpHeads:
                key = next(words)
                if not key == expKey:
                    raise(LogReadingException(expKey,key))
            
        def _readVelocityLine(self, words):
            # read ID number
            num = str_to_num(next(words))
            
            # read Velocity
            pos = np.zeros(self.dim)
            for i in range(self.dim):
                pos[i] = str_to_num(next(words))
            
            return num, pos
            
                
        def _readCurrentPositions(self,words):
            # Check for table headings
            self._readPositionTableHeaders(words)
            
            self.agentData = []
            # Read Agents' Data
            for i in range(self.nagent):
                num, ft, pos = self._readPositionLine(words)
                if num == i:
                    dat = {"ID": num, "cfit":ft, "cpos": pos}
                    self.agentData.append(dat)
                else:
                    raise(LogReadingException(i,num))
            
            #Check for section end flag
            key = next(words)
            if not key == "END_Current_Positions":
                raise(LogReadingException("END_Current_Position",key))
        
        def _readBestPositions(self, words):
            # assumed that self.agentData has been initialized
            self._readPositionTableHeaders(words)
            
            for i in range(self.nagent):
                num, ft, pos = self._readPositionLine(words)
                if num == i:
                    dat = {"ID": num, "bfit": ft, "bpos": pos}
                    self.agentData[i].update(dat)
                else:
                    raise(LogReadingException(i,num))
            
            key = next(words)
            if not key == "END_Best_Positions":
                raise(LogReadingException("End_Best_Position",key))
        
        def _readVelocities(self, words):
            self._readVelocityTableHeaders(words)
            
            for i in range(self.nagent):
                num, vel = self._readVelocityLine(words)
                if i == num:
                    dat = {"ID": num, "vel": vel}
                    self.agentData[i].update(dat)
                else:
                    raise(LogReadingException(i, num))
            
            key = next(words)
            if not key == "END_Velocities":
                raise(LogReadingException("END_Velocities",key))
        
        
