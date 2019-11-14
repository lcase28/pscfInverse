# LIBRARY IMPORTS
from psoinverse.mesophases.phaseManagement import MesophaseBase
from psoinverse.mesophases.mesophaseVariables import MesophaseVariable
from psoinverse.mesophases.mesophaseVariables import VariableTypes as varType
import psoinverse.SCFT.PSCF.FileManagers as filemanagers
from filemanagers.fieldfile import CoordFieldFile, WaveVectFieldFile
from filemanagers.paramfile import ParamFile
from filemanagers.outfile import OutFile

# EXTERNAL IMPORTS
import io
import numpy as np
import os

class PSCFMesophase(MesophaseBase):
    """ A Mesophase manager for simulating in PSCF. """
    
    def __init__(self, ID, kgrid, param, **kwargs):
        """
            Initialize the PSCF Mesophase using pre-formatted files.
            
            IMPLEMENTATION NOTE: Presently, no checks are done to ensure 
            that all template files are compatible (that system definitions match).
            These checks are assumed to have been done by the user.
            
            In parameter definitions, 'fileManager' refers to psoinverse.SCFT.PSCF.FileManagers
            
            Parameters
            ----------
            ID : string
                The name of the phase (should be unique within the run, 
                but this is not enforced)
            kgrid : filemanager.fieldfile.WaveVectFieldFile
                A pre-populated Wave-vector format field file.
            param : filemanager.paramfile.ParamFile
                A pre-populated param file. Should minimally contain definitional
                fields (format, MONOMERS, CHAINS, SOLVENTS, ..., BASIS) as well as
                the desired ITERATE fields. It should contain no other fields.
            
            
        """
        self.kgrid = kgrid
        self.param = param
        super().__init__(ID)
    
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
        super().update(VarSet,root)
    
    def setParams(self, VarSet):
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
        self._psoVars = VarSet
        # TODO: Add aditional variable options.
        for v in VarSet.items():
            f = v.flag
            if f == varType.BlockFraction:
                # TODO: Add capability to have any block act as floating fraction
                temp = deepcopy(self.param.block_length[v.polymer])
                N = np.sum(np.array(temp))
                newLen = v.scftValue * N
                Nshift = temp[v.block] - newLen
                temp[v.block] = newLen
                temp[-1] = temp[-1] + Nshift
                self.param.block_length[v.polymer] = temp
            elif f == varType.Chi:
                # TODO: Add capability to handle T-dependent chi
                m1, m2 = v.monomerIDs
                self.param.chi[m2][m1] = v.scftValue
            else:
                # TODO: Adjust this to return False rather than raising error for unknowns
                raise(NotImplementedError("Variable of Type " + str(f) + " not implemented"))
        return True
        
    def _evaluate(self, root):
        
        
    @property
    def energy(self):
        """
            The energy of the mesophase as of the most recent simulation.
            
            THIS IS PRESENTLY UNDER CONSTRUCTION. TEMPORARILY USING QUADRATIC
            FITNESS FUNCTION.
            
            Returns
            -------
            E : real or np.NaN
                If the simulation failed to converge, returns np.NaN.
                Else returns the simulated energy.
        """
        # TODO: Implement actual energy value.
        tot = 0.0
        for v in self._psoVars:
            tot = tot + v.scftVal**2
        return np.sqrt(max(tot, 0.0))
    
    
