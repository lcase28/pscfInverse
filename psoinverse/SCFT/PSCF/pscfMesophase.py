# LIBRARY IMPORTS
from psoinverse.SCFT.PSCF.FileManagers.fieldfile import CoordFieldFile, WaveVectFieldFile
from psoinverse.SCFT.PSCF.FileManagers.paramfile import ParamFile
from psoinverse.SCFT.PSCF.FileManagers.outfile import OutFile
from psoinverse.mesophases.FieldGenerators import FieldCalculator
from psoinverse.mesophases.phaseManagement import MesophaseBase
from psoinverse.mesophases.mesophaseVariables import MesophaseVariable
from psoinverse.mesophases.mesophaseVariables import VariableTypes as varType

# EXTERNAL IMPORTS
from copy import deepcopy
import io
import numpy as np
import os
import subprocess as sub    # Requires Python 3.5 --> System dependence

class PSCFMesophase(MesophaseBase):
    """ A Mesophase manager for simulating in PSCF. """
    
    def __init__(self, ID, kgrid, param, fieldGen, **kwargs):
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
                the desired ITERATE fields.
            fieldGen : psoinverse.mesophases.FieldGenerators.FieldGenerator
                A pre-initialized FieldGenerator object, set to generate fields
                for this phase.
        """
        self.kgrid = kgrid
        self.param = param
        self.fieldGen = fieldGen
        # lastLaunch will be used to hold a subprocess.CompletedProcess instance
        #   Will store system status data for the most recent run
        #   Previous runs currently not deemed useful.
        self.lastLaunch = None
        self.out = None
        super().__init__(ID)
    
    @classmethod
    def instanceFromFiles(cls, ID, kgridFile, paramFile, fieldGenFile, **kwargs):
        """
            Return an instance of PSCFMesophase initialized from given files.
            
            NOTE: No checks are done to ensure the specified files exist.
            If file-related exceptions occur, these will cascade.
            
            IMPLEMENTATION NOTE: Presently, no checks are done to ensure 
            that all template files are compatible (that system definitions match).
            These checks are assumed to have been done by the user.
            
            In parameter definitions, 'fileManager' refers to psoinverse.SCFT.PSCF.FileManagers
            
            Parameters
            ----------
            ID : string
                The name of the phase (should be unique within the run, 
                but this is not enforced)
            kgrid : pathlib.Path
                An already-resolved path object to a kgrid field file.
            param : pathlib.Path
                An already-resolved path object to a param file.
            fieldGen : pathlib.path
                An already-resolved path object to a fieldGen source file.
        """
        kgrid = WaveVectFieldFile(kgridFile.resolve())
        param = ParamFile(paramFile.resolve())
        fieldGen = FieldCalculator.from_file(fieldGenFile)
        return cls(ID, kgrid, param, fieldGen)
    
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
        return super().update(VarSet,root)
    
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
                #print(self.param.chi,m1,m2)
                self.param.chi[m2-1][m1-1] = v.scftValue
            else:
                # TODO: Adjust this to return False rather than raising error for unknowns
                raise(NotImplementedError("Variable of Type " + str(f) + " not implemented"))
        return True
        
    def _evaluate(self, root):
        """
            Launch a pscf simulation of the phase using current parameters.
            Parse the results and return the energy.
            
            ** Presently will not handle multiple blocks of same monomer
            type. i.e. will not handle ABA-type structures.
            
            Parameters
            ----------
            root : pathlib.Path
                The root directory of the run. All simulation files are
                placed in this location. It is assumed that path has 
                already been resolved.
            
            Returns
            -------
            energy : real
                The energy returned by the simulation.
            flag : bool
                True if simulation converged without issue.
                False if an error occurred and no energy found.
        """
        monFrac = self._getMonomerFractions()
        ngrid = self.param.ngrid
        newField = self.fieldGen.to_kgrid(monFrac,ngrid)
        self.kgrid.fields = newField
        kgridFile = root / 'rho_kgrid_in'
        self.kgrid.write(kgridFile.open(mode='x'))
        paramFile = root / 'param'
        self.param.write(paramFile.open(mode='x'))
        # TODO: Implement actual pscf launch and output parsing
        self._launchSim(root)
        self._readOutput(root)
        return 1.0, True
    
    def _launchSim(self, root):
        """
        Handle Launching of the simulation in given root directory
        
        Store completed process object in self.lastLaunch
        """
        launchTemplate = "pscf {} > simulationLog"
        launchCmd = launchTemplate.format("param")
        self.lastLaunch = sub.run(launchCmd, cwd=root)
        
        
    def _readOutput(self, root):
        
        
    
    def _getMonomerFractions(self):
        """ Return the volume fraction of all monomers """
        nmonomer = self.param.N_monomer
        nchain = self.param.N_chain
        frac = np.zeros(nmonomer)
        Ntot = 0.0
        for i in range(nchain):
            nblock = self.param.N_block[i]
            for j in range(nblock):
                mon = self.param.block_monomer[i][j] - 1
                Nb = self.param.block_length[i][j]
                Ntot += Nb
                frac[mon] += Nb
        for i in range(nmonomer):
            frac[i] = frac[i] / Ntot
        return frac
    
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
        for v in self._psoVars.items():
            tot = tot + v.scftValue**2
        return np.sqrt(max(tot, 0.0))
    
    
