# LIBRARY IMPORTS
from psoinverse.polymer.phases import MesophaseBase
from psoinverse.polymer.variables import MesophaseVariable
import psoinverse.polymer.parameters as parameters
import psoinverse.util.contexttools as contexttools

# Related Libraries
from pscfFieldGen.generation import (
    read_input_file, 
    generate_field_file,
    seed_calculator
)
from pscfFieldGen.filemanagers import PscfParam
from pscfFieldGen.filemanagers.pscf import (
    WaveVectFieldFile, 
    ParamFile, 
    OutFile, 
    SymFieldFile, 
    CoordFieldFile 
)

# EXTERNAL IMPORTS
from copy import deepcopy
import io
import numpy as np
import os
import pathlib
import subprocess as sub

class PscfMesophase(MesophaseBase):
    """ A Mesophase manager for simulating in PSCF. """
    
    def __init__(self, ID, paramWrap, fieldGen, coreOptions=None):
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
        param : PscfParam or pscf.ParamFile from pscfFieldGen.filemanagers
            A pre-populated param file. Should minimally contain definitional
            fields (format, MONOMERS, CHAINS, SOLVENTS, ..., BASIS) as well as
            the desired ITERATE fields.
        fieldGen : pscfFieldGen.fieldGenerators.FieldCalculator
            A pre-initialized FieldGenerator object, set to generate fields
            for this phase.
        coreOptions : list-like of int
            A list of monomer ID's which can be chosen as the "core" monomer.
            At each evaluation, the monomer in this list with the lowest overall
            volume fraction will be set as the core in initial guesses.
        """
        errmsg = "PscfMesophase requires paramWrap of type {}, or {}; gave {}"
        if isinstance(paramWrap,PscfParam):
            self.paramWrap = paramWrap
            self.param = paramWrap.file
        elif isinstance(paramWrap,ParamFile):
            self.param = paramWrap
            self.paramWrap = PscfParam(paramWrap)
        else:
            raise(TypeError(errmsg.format(PscfParam, ParamFile, type(paramWrap))))
        self.fieldGen = fieldGen
        if coreOptions is None:
            self.coreOptions = [ i for i in range(self.param.N_monomer) ]
        else:
            self.coreOptions = coreOptions
        super().__init__(ID)
    
    @classmethod
    def fromFieldGenFile(cls, ID, infile, coreOptions=None):
        """
        Use a pscfFieldGen input file to initialize a PscfMesophase object.
        
        The rules for the input file are identical to that of running pscfFieldGen.
        Operations are performed within the directory of the input file,
        creating no problems with distributed input files for multiple PscfMesophases.
        
        Parameters
        ----------
        ID : string
            The name of the phase (should be unique within the run, 
            but this is not enforced)
        infile : pathlib.Path
            The path to the file 
        
        Returns
        -------
        newPhase : PscfMesophase
            An instance of PscfMesophase based on the input file.
        """
        infile = infile.resolve()
        newdir = str(infile.parent)
        with contexttools.cd(newdir):
            param, calc, outfile, cmon = read_input_file(infile)
        out = cls(ID, param, calc, coreOptions)
        return out
    
    def startUpdate(self, VarSet, root, runner):
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
        return super().startUpdate(VarSet,root,runner)
    
    def _setup_calculations(self, root):
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
                numpy.nan if flag = False
            flag : bool
                True if simulation converged without issue.
                False if an error occurred and no energy found.
        """
        seed_calculator(self.fieldGen, self.paramWrap)
        return [root, self.paramWrap, self.fieldGen, self.coreOptions], True
    
    def finishUpdate(self, root):
        return super().finishUpdate(root)
    
    def _evaluate_energy(self, root):
        """
        Return energy of the phase.
        
        Parameters
        ----------
        root : pathlib.Path
            The root directory of the run. 
            All simulation files are
            placed in this location. 
            It is assumed that path has 
            already been resolved.
        
        Returns
        -------
        energy : real
            The energy of the mesophase.
            numpy.nan if flag = False
        flag : bool
            True if simulation converged without issue.
            False if an error occurred and no energy found.
        """
        # Check output file
        outfilename = root/"{}out".format(self.param.output_prefix)
        outfile = OutFile(outfilename.resolve())
        # check for finite values before any calculations
        critVals = [outfile.final_error, outfile.f_Helmholtz, outfile.f_homo]
        if not np.all(np.isfinite(critVals)):   
            # infinite energy if error or energy not finite value (+/-inf or NaN)
            return np.inf, False
        # Check convergence
        metError = outfile.final_error <= outfile.error_max
        metIter = outfile.iterations <= outfile.max_itr
        converged = metError and metIter
        if not converged:
            return np.inf, False
        # Check field similarity
        try:
            infield = "in.rho" # This specific file is currently a requirement.
            startFile = root/infield
            startField = SymFieldFile(startFile.resolve())
            endFile = root/"{}rho".format(self.param.output_prefix)
            endField = SymFieldFile(endFile.resolve())
            fieldSim = startField.fieldSimilarity(endField)
            if np.min(fieldSim) <= 0.7:
                return np.inf, False
        except ValueError as Err:
            if self.param.dim > 1:
                print("Error checking field similarity for in {}.".format(root))
                print(Err)
                return np.inf, False
        # Calculate Energy
        ener = outfile.f_Helmholtz - outfile.f_homo
        return ener, True
        
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
        """
        for p in VarSet.parameters:
            if isinstance(p,parameters.BlockLength):
                #temp = deepcopy(self.param.block_length[v.polymer])
                #N = np.sum(np.array(temp))
                #newLen = v.scftValue * N
                #Nshift = temp[v.block] - newLen
                #temp[v.block] = newLen
                #temp[-1] = temp[-1] + Nshift
                #self.param.block_length[v.polymer] = temp
                temp = self.param.block_length[p.polymerID]
                temp[p.blockID] = p.value
                self.param.block_length[p.polymerID] = temp
            elif isinstance(p,parameters.ChiN):
                m1, m2 = p.monomerIDs
                self.param.chi[m2][m1] = p.value
            elif isinstance(p, parameters.KuhnLength):
                mon = p.monomerID
                self.kuhn[mon] = p.value
            elif isinstance(p, parameters.PolymerBlendFraction):
                self.phi_chain[p.polymerID] = p.value
            else:
                raise(NotImplementedError("{} is not an implemented parameter".format(type(p))))
        return True
    
    @staticmethod
    def _launchSim(root, paramWrap, fieldGen, coreOptions):
        """
        Handle Launching of the simulation in given root 
        directory
        
        Store completed process object in self.lastLaunch
        
        Parameters
        ----------
        root : pathlib.Path
            The root directory of the run. 
            All simulation files are
            placed in this location. 
            It is assumed that path has 
            already been resolved.
        
        Returns
        -------
        flag : bool
            True if simulation converged without issue.
            False if an error occurred and no energy found.
        """
        param = paramWrap.file
        # Create Input Files
        monFrac = paramWrap.getMonomerFractions()
        #ngrid = param.ngrid
        #lattice = paramWrap.getLattice()
        core_mon = coreOptions[0]
        for i in coreOptions:
            if monFrac[i] < monFrac[core_mon]:
                core_mon = i
        #w = paramWrap.getInterfaceWidth(core_mon)
        #newField = fieldGen.to_kgrid(monFrac,ngrid,interfaceWidth=w, coreindex=core_mon, lattice=lattice)
        #kgrid = paramWrap.cleanFieldFile()
        #kgrid.fields = newField
        kgridName = param.fieldTransforms[0][1] # Pull init guess filename from param file
        kgridFile = root / kgridName
        generate_field_file(paramWrap,fieldGen,kgridFile,core=core_mon)
        #with kgridFile.open(mode='w') as f:
        #    kgrid.write(f)
        #kgrid.write(kgridFile.open(mode='w'))
        paramFile = root / 'param'
        with paramFile.open(mode='w') as f:
            param.write(f)
        #self.param.write(paramFile.open(mode='w'))
        # Launch Calculations
        infile = paramFile #root / "param"
        outfile = root / "simLog"
        try:
            with open(infile) as fin:
                with open(outfile,'w') as fout:
                    lastLaunch = sub.run("pscf",stdin=fin,stdout=fout, cwd=root)
            lastLaunch.check_returncode()
        except sub.CalledProcessError as err:
            with open(root/"errorLog") as ferr:
                ferr.write("Simuation error in: {}\n{}".format(root,err))
            return False
        return True
        
