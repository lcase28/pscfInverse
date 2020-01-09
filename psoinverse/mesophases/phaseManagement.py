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
    def __init__(self, ID, **kwargs):
        """
            All inheriting classes should define their own __init__
            First call in __init__ should be  call to the base class
            constructor to set the ID/name of the phase and set initial
            error state.
        """
        self.name = ID
        self._seterror()
        
    def update(self, VarSet, root, **kwargs):
        """ 
            Update phase variables, launch a simulation,
            and update phase energy from result.
            
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
        success = self.setParams(VarSet)
        if not success:
            self._seterror()
            return False
            
        root, success = self._checkPath(root)
        if not success:
            raise(ValueError("Invalid root directory provided"))
            
        ener, success = self._evaluate(root)
        if not success:
            self._seterror()
            return False
            
        # only reaches here if all else valid
        self._validState = True
        self._energy = ener
        return True
    
    def _seterror(self):
        self._validState = False
        self._energy = np.inf
    
    def _checkPath(self, root):
        return checkPath(root)
    
    @property
    def validState(self):
        """
            True if the Mesophase is in a stable (successfully 
            resolved and converged) state.
            False if the Mesophase is not up-to-date, or has 
            otherwise encountered an error.
        """
        return self._validState
    
    @abstractmethod
    def _evaluate(self, root, **kwargs):
        """
            Launch a simulation of the mesophase and parse results.
            
            Parameters
            ----------
            root : pathlib.Path
                The root directory of the run. All files from the
                simulation are to be placed in this directory.
                update() method should have already ensured this
                path exists prior to call to _evaluate().
            
            Returns
            -------
            energy : real or numpy.inf
                The energy returned by the simulation, relative to
                a hypothetical homogeneous disordered state.
                If an error occurred, and no energy available, 
                numpy.inf returned (arbitrarily high energy)
            flag : bool
                True if simulation converged without issue. 
                False if an error occurred and no energy could
        """
        pass
    
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
            E : real or np.NaN
                If the simulation failed to converge, returns np.NaN.
                Else returns the simulated energy.
        """
        if self.validState:
            return self._energy
        else:
            return np.inf
    
    @property
    def phaseName(self):
        return self.name

class MesophaseManager(object):
    """
        Manages a set of Mesophase objects based on universal set of 
        variable parameters. Responsible for coordinating launch, parsing,
        and comparison of multiple mesophase simulations
        
        *** should be SCFT-Simulator-Independent ***
    """
    
    def __init__(self, candidates, target, variables, **kwargs):
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
            variables : psoinverse.mesophases.mesophaseVariables.VariableSet
                The set of variables being handled in this PSO run.
        """
        self.candidates = candidates
        self.target = target
        self.target_name = self.target.phaseName
        self.variables = variables
        ## Ensure target is not in candidate phases
        self.candidates.pop( self.target_name, None )
        ## Flag to indicate whether successful simulations have 
        ##  been run on the current variable state.
        self._consistent = False
    
    @property
    def consistent(self):
        """ 
        Indicates whether or not the manager is in a consistent state.
        
        Returns True if successful simulations have been run on the 
        current set of variable values. Returns False otherwise.
        """
        return self._consistent
        
    @property
    def psoPoint(self):
        return self.variables.psoPoint
        
    @psoPoint.setter
    def psoPoint(self, val):
        self.variables.psoPoint = val
        self._consistent = False
    
    @property
    def psoBounds(self):
        return self.variables.psoBounds
    
    def update(self, root, newPoint = None, **kwargs):
        """
            Update the Mesophase variable values, run simulations,
            and recalculate fitness.
            
            Parameters
            ----------
            root : pathlib.Path
                Path to the root directory of the run. All files
                from this update will be written to this path, or
                or a child path.
                NOTE: Not implemented to handle concurrency.
            newPoint : psoinverse.PSO.SearchSpace.DictPoint
                A PSO point object (such as that returned from 
                MesophaseManager.psoPoint) containing the updated
                values for this update.
                If newPoint is None, simulations will be run at 
                current variable values
            
            Returns
            -------
            flag : bool
                Returns True if the update was successful. 
                Returns False otherwise.
        """
        if newPoint is not None:
            self.psoPoint = newPoint
        
        flag = self._evaluate(root)
        if not flag:
            print("PhaseManager eval Fail")
            print(self.statusString)
            self._errstate()
        else:
            print("PhaseManager update success")
            self._errstate(False)
        
        return flag
        
    # TODO: Revise to allow parallelization of phases??
    def _evaluate(self, root):
        root, success = checkPath(root) # resolve root path
        if not success:
            raise(ValueError("PhaseManager root path {} is invalid.".format(root)))
        phaseRoot = root/self.target.phaseName
        success = self.target.update( VarSet = self.variables, \
                                 root = phaseRoot )
        if not success:
            return False  # if target fails, error state
        ovrSuccess = False
        for (k, c) in self.candidates.items():
            phaseRoot = root/c.phaseName
            success = c.update( VarSet = self.variables, \
                                root = phaseRoot )
            ovrSuccess = ovrSuccess or success
        if not ovrSuccess:
            print("phaseManager Candidates Fail")
            #return np.nan, False # error state if all candidates fail
        #self._errstate(False) # simplify by leaving this process to update()
        return True
    
    def _errstate(self, flag=True):
        if flag:
            self._consistent = False
        else:
            self._consistent = True
    
    @property
    def fitness(self):
        """ 
            Calculate and return the fitness of the mesophase set 
            
            If object is in an inconsistent state, returns numpy.inf
        """
        WarnText = "Mesophase error state encountered " + \
                    "at an unexpected time. %s entered " + \
                    "invalid state outside of update. " + \
                    "Inconsistent behavior may result."
        unexpTarg = "Target Phase"
        unexpCand = "Candidate Phase"
        if not self.consistent:
            return np.NINF
        if not self.target.validState:
            self._errstate()
            raise(RuntimeWarning(WarnText % unexpTarg))
            return np.NINF
        tgtE = self.target.energy
        # Mesophase energies are defined relative to homogeneous
        # disorder. Initializing 'fit' to -tgtE effectively acts
        # to consider the disordered phase as a candidate.
        fit = -tgtE  
        for (k,c) in self.candidates.items():
            if c.validState:
                test = c.energy - tgtE
                if fit < test:
                    fit = test
        if fit is None:
            self._errstate()
            raise(RuntimeWarning(WarnText % unexpCand))
            return np.NINF
        return fit
    
    @property
    def statusString(self):
        s = self.__class__.__name__ + "\n"
        s += "Consistent = {!s}\n".format(self.consistent)
        s += "Fitness = {}\n".format(self.fitness)
        s += "Competing Phases Summary:\n"
        s += "Phase\tTgt\tValid\tEnergy\n"
        formstring = "{}\t{}\t{!s}\t{:E}\n"
        c = self.target
        s += formstring.format(c.phaseName, "Y", c.validState, c.energy)
        for (k, c) in self.candidates.items():
            s += formstring.format(c.phaseName, "N", c.validState, c.energy)
        return s
            
    
