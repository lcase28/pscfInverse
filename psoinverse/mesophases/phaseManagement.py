"""
Module contains base class definitions for candidate mesophase management.

Declarations here allow PSO algorithm objects to uniformly access mesophase
data and simulations from any given SCFT simulator by abstracting specifics
of the simulator formatting from the interface through class hierarchies.
"""
from psoinverse.mesophases.mesophaseVariables import VariableSet

from abc import ABC, abstractmethod
from copy import deepcopy
import numpy as np
import os
from pathlib import Path
import io

def checkPath(root):
    """
        Checks if path exists, creates it if not.
        Checks if path is a directory, if not takes the
        parent directory.
        Fully resolves path.
        
        Parameters
        ----------
        root : pathlib.Path
            The directory path to be resolved.
        
        Returns
        -------
        resolvedRoot : pathlib.Path
            The resolved path to the specified directory.
        flag : bool
            True if the path resolved without error.
            False if path resolution threw one of:
            FileNotFoundError, FileExistsError, RuntimeError,
            which capture the expected throws from the Path
            object during resolution.
    """
    root = root.resolve()
    if root.exists() and root.is_dir():
        return root.resolve(), True
    # path either does not exist or is not directory
    try:
        if not root.exists():
            root.mkdir(parents=True)
        if not root.is_dir():
            root = root.parent
        return root.resolve(), True
    except (FileNotFoundError, FileExistsError, RuntimeError):
        return None, False
    
class MesophaseBase(ABC):
    """ 
    Abstract base class for mesophase wrapper classes.
    Derived classes manage a single mesophase for a single SCFT simulator.
    """
    
    @abstractmethod
    def __init__(self, ID):
        """
        All inheriting classes should define their own __init__
        First call in __init__ should be  call to the base class
        constructor to set the ID/name of the phase and set initial
        error state.
        """
        self.name = ID
        self._energy = np.inf
        self._updating = False
        self._lastLaunch = None
        self._updateCount = -1
    
    @property
    def stable(self):
        """
        True if not updating, not waiting on calculations, and has valid fitness.
        
        Return value is analogous to self.readyForStartUpdate, but will return
        False on initialization, before the first calculation has been run.
        """
        out = self._updateCount >= 0 and self.readyForStartUpdate
        return out
    
    @property
    def updating(self):
        return self._updating
    
    @property
    def calculationFinished(self):
        """
        True if the latest calculation is complete and results available.
        
        If not in the process of an update, returns False.
        """
        if self._lastLaunch is None:
            return False
        else:
            return self._lastLaunch.ready()
    
    @property
    def readyForStartUpdate(self):
        if self._lastLaunch is None and not self.updating:
            return True
        else:
            return False
    
    @property
    def readyForFinishUpdate(self):
        flag = self.updating and self.calculationFinished
        return flag
        
    def startUpdate(self, VarSet, root, runner):
        """ 
        Update phase variables, and launch a simulation.
        
        If deriving classes override, they should end
        with a call to super().update, and allow super
        to run calls to setParams and _evaluate.
        
        Parameters
        ----------
        VarSet : psoinverse.mesophases.mesophaseVariables.VariableSet
            The set of all variables to be updated, with their 
            current values
        root : pathlib.Path (OS-dependent type)
            The root directory of the run. All files from the
            simulation are to be placed in this directory.
            If the directory does not exist, it will be created.
        
        Returns
        -------
        flag : bool
            True if updated without error.
            False otherwise.
        """
        if not self.readyForStartUpdate:
            return False
        
        success = self.setParams(VarSet)
        if not success:
            return False
            
        root, success = self._checkPath(root)
        if not success:
            raise(ValueError("Invalid root directory provided"))
            
        args, success = self._setup_calculations(root, runner)
        if not success:
            return False
        self._updating = True
        self._lastLaunch = runner.addTask(self._launchSim, args)
        return True
    
    @abstractmethod
    def _setup_calculations(self, root):
        """
        Perform any setup required before launching the calculation.
        
        Make any required state changes to the object prior to launching
        the calculation. Put together a list of arguments for to pass
        to self._launchSim and return this list (will be passed to 
        self._launchSim using *args).
        
        Parameters
        ----------
        root : pathlib.Path
            The root directory of the run. All files from the
            simulation are to be placed in this directory.
            update() method should have already ensured this
            path exists prior to call to _evaluate().
        
        Returns
        -------
        args : list
            List of arguments for self._launchSim
        flag : boolean
            True if setup occurred without issue.
        """
        pass
    
    @abstractmethod
    def _launchSim(self):
        """
        Finalize setup and run the simulation.
        
        This method will be passed to the parallel calculation manager.
        During this method, any additional calculation setup can be done, 
        and the calculation should be run. This call should not make any
        modifications to the object state.
        
        Any arguments should be returned from self._setup_calculations(...).
        
        Return values are discouraged, but can be accessed through
        the multiprocessing.AsyncResult attribute self._lastLaunch.
        """
        pass
    
    def finishUpdate(self, root):
        if not self.readyForFinishUpdate:
            return False
        flag = True
        root, success = self._checkPath(root)
        if not success:
            raise(ValueError("Invalid root directory provided"))
        ener, success = self._evaluate_energy(root)
        if not success:
            flag = False
        # If reaches this point, state valid
        self._updateCount += 1
        self._updating = False
        self._lastLaunch = None
        self._energy = ener
        return flag
    
    @abstractmethod
    def _evaluate_energy(root):
        """
        Determine energy of the phase from calculation results.
        
        Derived classes must override. Overriding method must
        read results from the SCFT calculations to determine
        the free energy of the phase.
        
        Parameters
        ----------
        root : pathlib.Path
            The root directory containing calculation results.
        
        Returns
        -------
        energy : numeric or numpy.inf
            The energy of the phase based on SCFT.
            If an error occurs and an energy is unable to be
            determined, numpy.inf should be returned as an arbitrarily
            high free energy.
        """
        pass
        
    def _checkPath(self, root):
        return checkPath(root)
    
    @abstractmethod
    def setParams(self, VarSet, **kwargs):
        """
        Update the specified set of variables.
        
        Parameters
        ----------
        VarSet : psoinverse.mesophases.mesophaseVariables.VariableSet
            The set of MesophaseVariable objects to be updated.
        
        Returns
        -------
        flag : bool
            True if parameters updated without issue.
            False otherwise
        
        NOTE:
        Built-in Variable Types can be found in the mesophaseVariables
        module.
        
        NOTE: 
        A given variable type may or may not be implemented
        by a derived Mesophase. See derived classes for details.
        """
        pass
    
    @property
    def energy(self):
        """
        The energy of the mesophase as of the most recent simulation.
        
        Returns
        -------
        E : real or numpy.inf
            If the simulation failed to converge, or another
            error occurred, returns numpy.inf.
            Else returns the simulated energy.
        """
        return self._energy
    
    @property
    def phaseName(self):
        return self.name

class MesophaseManager(object):
    """
    Manages a set of Mesophase objects.
    
    The manager contains a target phase and a set of candidates.
    At each update, the manager distributes the new values to the
    competitor phases to initiate updates. It ensures that all phases
    have updated before new updates can start. It then uses the energy
    of the target and candidate phases to determine the "fitness" based
    on the relative energetic stability of the target relative to
    the competitors.
    """
    
    def __init__(self, candidates, target):
        """
        Initialize the MesophaseManager.
        
        Parameters
        ----------
        candidates : dict
            Keys are name of phase.
            Values are MesophaseBase-like objects.
            The set of initialized candidate phases.
            All should be capable of launching simulations
            as-is.
        target : MesophaseBase-like
            The phase object to be treated as the "target" phase.
            As with candidates, should be initialized and ready to run.
        """
        self.candidates = candidates
        self.target = target
        self.target_name = self.target.phaseName
        ## Ensure target is not in candidate phases
        self.candidates.pop( self.target_name, None )
        self._fitness = np.NINF
    
    @property
    def stable(self):
        """
        Specifies if all competing phases are stable.
        """
        flag = self.target.stable
        for (k,c) in self.candidates.items():
            flag = flag and c.stable
        return flag
    
    @property
    def updating(self):
        """
        True if any candidate is still updating.
        """
        flag = self.target.updating
        for (k,c) in self.candidates.items():
            flag = flag or c.updating
        return flag
    
    @property
    def calculationFinished(self):
        """
        True if all candidates are finished calculating
        """
        flag = self.target.calculationFinished
        for (k,c) in self.candidates.items():
            flag = flag and c.calculationFinished
        return flag
    
    @property
    def readyForStartUpdate(self):
        """
        True if all candidates are ready for an update.
        """
        flag = self.target.readyForStartUpdate
        for (k,c) in self.candidates.items():
            flag = flag and c.readyForStartUpdate
        return flag
    
    @property
    def readyForFinishUpdate(self):
        """
        True if all candidates are ready for finishUpdate call.
        """
        flag = self.target.readyForFinishUpdate
        for (k,c) in self.candidates.items():
            flag = flag and c.readyForFinishUpdate
        return flag
    
    def startUpdate(self, varSet, root, runner):
        """
        Update the Mesophase variable values, setup calculations,
        and pass calculation commands to runner.
        
        If the object is not in a state to begin an update, the 
        request is rejected, and the method returns False.
        
        Parameters
        ----------
        varSet : PolymerVariableSet
            The set of variables (with updated values) to modify
            for this set of calculations.
        root : pathlib.Path
            Path to the root directory of the run. All files
            from this update will be written to this path, or
            or a child path.
        runner : util.parallelBatches.LocalBatchRunner
            The job manager to which all calculations should be
            submitted for parallel processing.
        
        Returns
        -------
        flag : bool
            Returns True if the update was successful. 
            Returns False otherwise.
        """
        if not self.readyForStartUpdate:
            return False
        
        root, success = checkPath(root) # resolve root path
        if not success:
            raise(ValueError("PhaseManager root path {} is invalid.".format(root)))
        
        flag = True
        phaseRoot = root/self.target.phaseName
        success = self.target.startUpdate(varSet, phaseRoot, runner)
        if not success:
            flag = False  # if target fails, error state
        for (k, c) in self.candidates.items():
            phaseRoot = root/c.phaseName
            success = c.startUpdate(varSet, phaseRoot, runner)
            flag = flag or success
        
        return flag
        
    def finishUpdate(self, root):
        root, success = checkPath(root) # resolve root path
        if not success:
            raise(ValueError("PhaseManager root path {} is invalid.".format(root)))
        
        flag = True
        phaseRoot = root/self.target.phaseName
        success = self.target.finishUpdate(phaseRoot)
        flag = flag or success
        for (k, c) in self.candidates.items():
            phaseRoot = root/c.phaseName
            success = c.finishUpdate(phaseRoot)
            flag = flag or success
        
        self._fitness = self._evaluate_fitness
        
        return flag
    
    def _evaluate_fitness()
        tgtE = self.target.energy
        # Mesophase energies are defined relative to homogeneous
        # disorder. Initializing 'fit' to -tgtE effectively acts
        # to consider the disordered phase as a candidate.
        fit = -tgtE  
        for (k,c) in self.candidates.items():
            test = c.energy - tgtE
            if fit > test:
                fit = test
        return fit
        
    @property
    def fitness(self):
        """ 
        The relative energetic stability of the target phase.
        
        This result should be considered the optimization fitness
        to be maximized.
        """
        return self._fitness
    
    @property
    def statusString(self):
        s = self.__class__.__name__ + ", "
        s += "Fitness = {}, ".format(self.fitness)
        formstring = "{}({}) = {:E}"
        c = self.target
        s += formstring.format("Tgt", c.phaseName, c.energy)
        for (k, c) in self.candidates.items():
            s += formstring.format("Comp", c.phaseName, c.energy)
        s += "\n"
        return s
            
    
