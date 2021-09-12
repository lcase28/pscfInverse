"""
Module contains base class definitions for candidate mesophase management.

Declarations here allow PSO algorithm objects to uniformly access mesophase
data and simulations from any given SCFT simulator by abstracting specifics
of the simulator formatting from the interface through class hierarchies.
"""
# Standard Library Imports
from abc import ABC, abstractmethod
from copy import deepcopy
import os
from pathlib import Path
import io

# Third Party Library Imports
import numpy as np

# Project Imports
from pscfinverse.polymer.variables import PolymerVariableSet
from pscfinverse.util.iotools import checkPath

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
            if self._lastLaunch.ready():
                if self._lastLaunch.successful():
                    return True
                else:
                    self._lastLaunch.get()
            else:
                return False
    
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
        VarSet : pscfinverse.mesophases.mesophaseVariables.VariableSet
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
            
        args, success = self._setup_calculations(root)
        if not success:
            return False
        self._updating = True
        self._lastLaunch = runner.addTask(self._launchSim, args)
        #self._launchSim(*args)
        return True
    
    @abstractmethod
    def _setup_calculations(self, root):
        """
        Perform any setup required before launching the calculation.
        
        Make any required state changes to the object prior to launching
        the calculation. Put together a list of arguments for to pass
        to self._launchSim and return this list (will be passed to 
        self._launchSim using *args). This list should include any instance
        data, as self._launchSim is a static method and does not receive a
        copy of instance data implicitly.
        
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
    
    @staticmethod
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
        VarSet : pscfinverse.mesophases.mesophaseVariables.VariableSet
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
    
    The manager contains a target phase and a set of competitors.
    At each update, the manager distributes the new values to the
    competitor phases to initiate updates. It ensures that all phases
    have updated before new updates can start. It then uses the energy
    of the target and candidate phases to determine the "fitness" based
    on the relative energetic stability of the target relative to
    the competitors.
    """
    
    def __init__(self, competitors, target):
        """
        Initialize the MesophaseManager.
        
        All considered phases (target and competitors) must have unique 
        values for their phaseName property to ensure integrity of the 
        data output methodology. Duplicate phase names will raise a
        ValueError.
        
        Parameters
        ----------
        competitors : list of MesophaseBase sub classes
            The set of initialized candidate phases.
            All should be capable of launching simulations
            as-is.
        target : MesophaseBase sub-class or list of MesophaseBase sub class
            The phase objects to be treated as the "target" phases.
            As with competitors, should be initialized and ready to run.
        """
        ## Ensure target is not in candidate phases
        names = []
        duplicationMsg = "Duplication of phase name {} not allowed"
        for c in competitors:
            if c.phaseName in names:
                raise(ValueError(duplicationMsg.format(c.phaseName)))
            else:
                names.append(c.phaseName)
        if isinstance(target,MesophaseBase):
            targetlist = [target]
        elif isinstance(target,list):
            targetlist = target
        else:
            raise(TypeError("Invalid type for input target: {}".format(type(target))))
        for t in targetlist:
            if t.phaseName in names:
                raise(ValueError(duplicationMsg.format(t.phaseName)))
            else:
                names.append(t.phaseName)
                
        self._competitors = competitors
        self._target = targetlist
        self._fitness = np.NINF
    
    @property
    def competitors(self):
        """ A List of competing candidate phases. """
        return self._competitors
    
    @property
    def target(self):
        """ The list of target phases. """
        return self._target
    
    @property
    def stable(self):
        """
        Specifies if all competing phases are stable.
        """
        flag = True
        for t in self.target:
            flag = flag and t.stable
        for c in self.competitors:
            flag = flag and c.stable
        return flag
    
    @property
    def updating(self):
        """
        True if any candidate is still updating.
        """
        flag = False
        for t in self.target:
            flag = flag or t.updating
        for c in self.competitors:
            flag = flag or c.updating
        return flag
    
    @property
    def calculationFinished(self):
        """
        True if all competitors are finished calculating
        """
        flag = True
        for t in self.target:
            flag = flag and t.calculationFinished
        for c in self.competitors:
            flag = flag and c.calculationFinished
        return flag
    
    @property
    def readyForStartUpdate(self):
        """
        True if all candidates are ready for an update.
        """
        flag = True
        for t in self.target:
            flag = flag and t.readyForStartUpdate
        for c in self.competitors:
            flag = flag and c.readyForStartUpdate
        return flag
    
    @property
    def readyForFinishUpdate(self):
        """
        True if all candidates are ready for finishUpdate call.
        """
        flag = True
        for t in self.target:
            flag = flag and t.readyForFinishUpdate
        for c in self.competitors:
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
        for t in self.target:
            phaseRoot = root/t.phaseName
            success = t.startUpdate(varSet, phaseRoot, runner)
            if not success:
                flag = False  # if target fails, error state
        for c in self.competitors:
            phaseRoot = root/c.phaseName
            success = c.startUpdate(varSet, phaseRoot, runner)
            flag = flag or success
        
        return flag
        
    def finishUpdate(self, root):
        root, success = checkPath(root) # resolve root path
        if not success:
            raise(ValueError("PhaseManager root path {} is invalid.".format(root)))
        
        flag = True
        for t in self.target:
            phaseRoot = root/t.phaseName
            success = t.finishUpdate(phaseRoot)
            flag = flag or success
        for c in self.competitors:
            phaseRoot = root/c.phaseName
            success = c.finishUpdate(phaseRoot)
            flag = flag or success
        
        self._fitness = self._evaluate_fitness()
        
        return flag
    
    def _evaluate_fitness(self):
        startphase = self.target[0]
        tgtE = startphase.energy
        # Consider most stable target for fitness
        for t in self.target:
            test = t.energy
            if test < tgtE:
                tgtE = test
        # Mesophase energies are defined relative to homogeneous
        # disorder. Initializing 'fit' to -tgtE effectively acts
        # to consider the disordered phase as a candidate.
        fit = -tgtE  
        for c in self.competitors:
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
    def labels(self):
        """ 
        A list of names of competing phases. 
        
        Each candidate label is marked to indicate status as a target (tgt)
        or competitor (cmp) such that each label takes the form 
        "<lbl>_<phaseName>" where <lbl> represents its status flag (tgt or cmp),
        and <phaseName> is the unique phase name assigned to the candidate.
        
        The target phases are always placed at the start of the list.
        """
        lbl = []
        formstr = "tgt_{}"
        for v in self.target:
            lbl.append(formstr.format(v.phaseName))
        formstr = "cmp_{}"
        for v in self.competitors:
            lbl.append(formstr.format(v.phaseName))
        return lbl
    
    @property
    def energies(self):
        """ 
        A list of phase energies at the most recent update. 
        
        The energy at each list position corresponds to the candidate phase
        in the same position of the MesophaseManager.labels property.
        """
        out = []
        for v in self.target:
            out.append(v.energy)
        for v in self.competitors:
            out.append(v.energy)
        return out
    
    @property
    def statusString(self):
        s = self.__class__.__name__ + ", "
        s += "Fitness = {}, ".format(self.fitness)
        formstring = "{}({}) = {:E}"
        for t in self.target:
            s += formstring.format("Tgt", t.phaseName, t.energy)
        for c in self.competitors:
            s += formstring.format("Cmp", c.phaseName, c.energy)
        s += "\n"
        return s
            
    
