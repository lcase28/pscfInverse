# Library imports
from psoinverse.PSO.Agent import Agent
from psoinverse.PSO.SearchSpace import Point
from psoinverse.mesophases.phaseManagement import MesophaseManager

# temporary imports (will be changing form later)
# TODO: change these
from psoinverse.PSO.Swarm import checkPath

# Other imports
from copy import deepcopy
import numpy as np
import pathlib

class ScftAgent(Agent):
    """ Derived Agent Class for SCFT phase search. """
    
    def __init__(self, phaseManager, randGen, root = None, \
                logFileName = 'AgentLog', autoLog = True, \
                batchRunner, *args, **kwargs):
        """
            Constructor for ScftAgent
            
            Parameters
            ----------
            phaseManager : MesophaseManager object
                Fully initialized phase manager class
            randGen : np.random.RandomState
                Seeded random number generator used to initialize position
            root : pathlib.Path
                Path to directory where agent can place any files generated.
                Defaults to pathlib.Path.cwd()
            logFileName : string, valid filename
                The file to which the Agent will write its log.
            autoLog : bool
                Whether or not the agent should automatically write out to its
                log, or wait for caller to do so. Default: True
        """
        self.phaseManager = deepcopy(phaseManager)
        if root is None:
            self.root = pathlib.Path.cwd()
        else:
            self.root, success = checkPath(root)
            if not success:
                raise(ValueError("Given root path failed to resolve."))
        self.logFile = self.root / logFileName
        self.autoLog = autoLog
        location = self.phaseManager.psoPoint
        bounds = self.phaseManager.psoBounds
        velSrc = np.zeros_like(location.Coords)
        # TODO: Figure out better way to force clean file on startup.
        self._startup = True
        super().__init__(bounds, location, velSrc, randGen, batchRunner, seekMax = True)
        self._startup = False
        
    def _setup_calculations(self):
        stepRoot = self.root / "step{}".format(self.steps)
        success = self.phaseManager.startUpdate(stepRoot, self._runner, self.Location)
        if not success:
            print("scftAgent phasemanager setup Fail")
            self._startErrorState()
    
    def finishUpdate(self):
        super().finishUpdate()
        self._tryLog()
    
    def _evaluate_fitness(self):
        stepRoot = self.root / "step{}".format(self.steps)
        self.phaseManager.finishUpdate(stepRoot, self._runner)
        if not self.phaseManager.consistent:
            print("scftAgent phaseManager inconsistent.")
            self._startErrorState()
        self.Location.Fitness = self.phaseManager.fitness
    
    @property
    def statusString(self):
        s = "Agent {}:\n".format(self.id)
        s += "Step {}:\n".format(self.steps)
        s += "Location:\n{}\n".format(self.Location)
        s += "Personal Best:\n{}\n".format(self.PBest)
        s += "Phase Manager Status:\n{}\n".format(self.phaseManager.statusString)
        return s
    
    def writeLog(self, refresh=False):
        s = self.statusString
        if refresh:
            with self.logFile.open(mode='w') as f:
                f.write(s)
        else:
            with self.logFile.open(mode='a') as f:
                f.write(s)
    
    def _tryLog(self):
        if self.autoLog:
            self.writeLog(self._startup)
        
