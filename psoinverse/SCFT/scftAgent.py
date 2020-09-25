# Library imports
from psoinverse.PSO.Agent import Agent
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
    
    def __init__(
        self, 
        varSet, 
        phaseManage, 
        velocitySource, 
        calcManager, 
        parentRoot = None):
        """
        Constructor for ScftAgent
        
        Parameters
        ----------
        varSet : PolymerVariableSet
            The variables being searched in this optimization.
        phaseManage : MesophaseManager
            Fully initialized phase manager class.
        velocitySource : Velocity
            A template Velocity object.
            Values will be randomized on initial update.
        runner : CalculationManager
            A parallel calculation manager for running calculations.
        parentRoot : pathlib.Path
            The parent directory of the Agent's output root.
            Generally, this is the Swarm's root directory. Each Agent
            will place their output in a subdirectory of this root
            named 'parentRoot/agent{id}/'.
        """
        if not isinstance(varSet, PolymerVariableSet):
            raise(TypeError("ScftAgent requires a PolymerVariableSet, not {}".format(type(varSet))))
        if not isinstance(phaseManage, MesophaseManager):
            raise(TypeError("ScftAgent requires MesophaseManager, not {}".format(type(phaseManage))))
        self.__phaseManager = deepcopy(phaseManage)
        super().__init__(varSet, velocitySource, calcManager, parentRoot)
    
    @property
    def calculationFinished(self):
        """
        Determine if all phase calculations are complete.
        """
        return self.__phaseManager.calculationFinished()
        
    @property
    def readyForStartUpdate(self):
        out1 = self.__phaseManager.readyForStartUpdate
        out2 = super().readyForStartUpdate
        return out1 and out2
    
    @property
    def readyForFinishUpdate(self):
        out1 = self.__phaseManager.readyForFinishUpdate()
        out2 = super().readyForFinishUpdate()
        return out1 and out2
        
    def _setup_calculations(self, calcManager):
        stepRoot = self.root / "step{}".format(self.nextStep)
        success = self.phaseManager.startUpdate(stepRoot, calcManager)
    
    def _cleanup_calculations(self, calcManager):
        stepRoot = self.root / "step{}".format(self.nextStep)
        self.phaseManager.finishUpdate(stepRoot)
        fitness = self.phaseManager.fitness
        return fitness
    
    def _unset_calculations(self):
        self.phaseManager.cancelUpdate()
    
        
