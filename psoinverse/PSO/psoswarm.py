
# Standard Imports
from copy import deepcopy
import pathlib
import time

# Third-party Imports
import numpy as np
import networkx as nx

# Library imports
from psoinverse.pso.core import Point, FITNESS_COMPARATOR
from psoinverse.pso.integrators import Integrator
from psoinverse.pso.agent import Agent

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

def allHaveStep(neighbors, step):
    """ True if all agents have completed step """
    for n in neighbors:
        if not n.hasStep(step):
            return False
    return True

# =========================================================
#
# Swarm Class
#
# =========================================================
class Swarm(object):
    def __init__(self, graph, agents, integrator, root = None):
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
            All agents should have had their initialization started
            (Calculations need not have finished, but should be 
            awaiting step 0 calculation results)
        integrator : psoinverse.PSO.Integrators.Integrator child class object
            An object derived from the Integrator abstract base class
            representing the update procedure for the swarming behavior.
        root : pathlib.path, optional.
            A path object pointing to the Swarm's root directory.
            All output from Swarm will be placed in this directory
            (or in sub-directories, if so implemented). If directory
            does not exist, it will be created.
            Default Value is the current working directory
        
        Raises
        ------
        ValueError :
            If number of nodes in graph does not match number of agents.
            If root is unable to be resolved or created.
        """
        if not len(graph) == len(agents): 
            raise(ValueError("Number of nodes in graph doesn't match number of agents"))
        self.agents = agents
        self.graph = graph
        self.integrator = integrator
        self.__fit_compare = FITNESS_COMPARATOR
        self.__next_step = 0
        # record-keeping
        if root is None:
            self.root = pathlib.Path.cwd()
            flag = True
        else:
            self.root, flag  = checkPath(root)
        self.__start_log_file()
    
    @property
    def nextStep(self):
        """ The next step requiring processing and logging. """
        return self.__next_step
        
    @property
    def lastStep(self):
        """
        Return the number of steps fully completed within the swarm.
        """
        return self.nextStep - 1
        
    def hasStep(self, step):
        """
        Check if all agents have completed the given step number.
        """
        for a in self.agents:
            if not a.hasStep(step):
                return False
        return True
    
    @property
    def maxStep(self):
        """
        Return the highest step achieved by any agent.
        """
        out = -1
        for a in self.agents:
            if out < a.lastStep:
                out = a.lastStep
        return out
    
    @property
    def minStep(self):
        """ Return the lowest step completed by any agent. """
        out = -1
        for (i,a) in enumerate(self.agents):
            if out > a.lastStep or i == 0:
                out = a.lastStep
        return out
    
    def bestHistoricalAgentId(self, step):
        """
        Return the ID number of the agent with the best historical position at step.
        """
        bpos = [ a.bestPointAtStep(step) for a in self.agents ]
        best, index = self.__fit_compare.bestPointIndex(bpos)
        return index
    
    def bestCurrentAgentId(self, step):
        cpos = [ a.stablePointAtStep(step) for a in self.agents ]
        best, index = self.__fit_compare.bestPointIndex(bpos)
        return index
    
    def agentNeighbors(self, index):
        return [ self.agents[i] for i in self.graph.neighbors(index) ]
    
    def run(self, nstep, tryCooldown=0.1):
        """
        Run the swarm out to nstep steps.
        
        The swarm will run the agents out asynchronously
        until all have completed nstep steps.
        
        Because agents may take time to calculate,
        pauses will be taken between attempts to advance
        each agent. This delay time is tunable with the
        variable tryCooldown. If agent calculations are expected
        to take substantial time, longer cooldown periods may be
        used. For agents expected to calculate rapidly, smaller
        times will speed the overall program.
        
        Parameters
        ----------
        nstep : int
            The maximum number of steps to take.
            This number must include steps taken in previous
            calls to swarm.run(), as those are counted in step
            totals.
        tryCooldown : positive, real number (optional)
            The amount of time (in seconds) to pause in between 
            attempts to advance the agents.
            Default = 0.1
        
        Returns
        -------
        runtime : float
            Total seconds elapsed between the start of the method
            and its termination.
        proctime : float
            Total system and user CPU time in the main process
            (excludes cooldown time and cpu time in parallel
            calculations).
        """
        starttime = time.time()
        startproc = time.process_time()
        while not self.hasStep(nstep):
            self.finishAgents()
            self.startAgents(nstep)
            self.logData(nstep)
            time.sleep(tryCooldown)
        endtime = time.time()
        endproc = time.process_time()
        runtime = endtime - starttime
        proctime = endproc - startproc
        return runtime, proctime
    
    def startAgents(self, nstep):
        """
        Attempt to start an update for each agent in the swarm.
        
        In order for an update to be started for an agent, the 
        following conditions are checked:
            1.  Agent.readyForStartUpdate must return True.
            2.  The Agent's neighbors must all have completed the the agent's last step
            3.  The Agent must not have completed nstep steps.
        
        Parameters
        ----------
        nstep : int
            The maximum number of steps to allow each Agent.
        
        Returns
        -------
        startCount : int
            The number of agents starting updates during call.
        """
        startCount = 0
        for (i,a) in enumerate(self.agents):
            if not a.readyForStartUpdate:
                continue # not ready for update
            step = a.lastStep
            if step == nstep:
                continue # completed all steps
            neighbors = self.agentNeighbors(i)
            if not allHaveStep(neighbors, step):
                continue # neighbors are not caught up
            newPos, newVel = self.integrator.integrate(a,neighbors)
            a.startUpdate(newPos,newVel)
            startCount += 1
        return startCount
    
    def finishAgents(self):
        """
        Attempt to finish an update for each agent in the swarm.
        
        Returns
        -------
        finishCount : int
            The number of agents finishing an update during the call.
        """
        finishCount = 0
        for (i,a) in enumerate(self.agents):
            if a.readyForFinishUpdate:
                a.finishUpdate()
                finishCount += 1
        return finishCount
    
    def logData(self, nstep):
        """
        Parse and record swarm-level data.
        
        Return
        ------
        didUpdate : boolean
            True if the next step was available and recorded.
            False otherwise
        """
        if not self.hasStep(self.nextStep):
            return False
        logdat = []
        step = self.nextStep
        logdat.append( step )
        logdat.append( self.minStep )
        logdat.append( self.maxStep )
        logdat.append( self.bestHistoricalAgentId(step) )
        logdat.append( self.bestCurrentAgentId(step) )
        for (i,a) in enumerate(self.agents):
            logdat.append( a.lastStep )
        self.__write_list_to_log(logdat)
        self.__next_step += 1
        return True
    
    def __start_log_file(self):
        dat = []
        dat.append("step")
        dat.append("minStep")
        dat.append("maxStep")
        dat.append("histBestAgent")
        dat.append("currBestAgent")
        for (i,a) in enumerate(self.agents):
            dat.append("agent{}step".format(i))
        self.__write_list_to_log(dat)
    
    def __write_list_to_log(self,data):
        logFile = self.root / "swarmLog.csv"
        strline = ','.join(map(str,data))
        strline = strline + '\n'
        with open(logFile,'a') as f:
            f.write(strline)

