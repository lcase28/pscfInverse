# LIBRARY IMPORTS
#from psoinverse.SCFT.PSCF.FileManagers.fieldfile import CoordFieldFile, WaveVectFieldFile
#from psoinverse.SCFT.PSCF.FileManagers.paramfile import ParamFile
#from psoinverse.SCFT.PSCF.FileManagers.outfile import OutFile
#from psoinverse.mesophases.FieldGenerators import FieldCalculator
from psoinverse.polymer.phases import MesophaseBase
from psoinverse.polymer.variables import MesophaseVariable
import psoinverse.polymer.parameters as parameters
import psoinverse.util.contexttools as contexttools

# Related Libraries
from pscfFieldGen.generation import FieldCalculatorBase, read_input_file, generate_field_file
from pscfFieldGen.filemanagers import PscfParam
from pscfFieldGen.filemanagers.pscf import WaveVectFieldFile, ParamFile, OutFile

# EXTERNAL IMPORTS
from copy import deepcopy
import io
import numpy as np
import os
import pathlib
import subprocess as sub    # Requires Python 3.5 --> System dependence

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
        self.coreOptions = coreOptions
        super().__init__(ID)
    
    @classmethod
    def fromFieldGenFile(cls, ID, infile):
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
        out = cls(ID, param, calc)
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
        ngrid = self.param.ngrid
        lattice = self.paramWrap.getLattice()
        self.fieldGen.seedCalculator(ngrid,lattice)
        return root, True
    
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
        outfilename = root/"out"
        outfile = OutFile(outfilename.resolve())
        metError = outfile.final_error <= outfile.error_max
        metIter = outfile.iterations <= outfile.max_itr
        converged = metError and metIter
        if not converged:
            return np.inf, False
        ener = outfile.f_Helmholtz - outfile.f_homo
        flag = True
        ## TODO: Check Field Similarity with symmetrized files for quicker I/O
        # Check field similarity
        #startFile = root/"rho_kgrid_in"
        #startField = WaveVectFieldFile(startFile.resolve())
        #endFile = root/"rho_kgrid"
        #endField = WaveVectFieldFile(endFile.resolve())
        #self.fieldSim = startField.fieldSimilarity(endField)
        #if np.min(self.fieldSim) <= 0.8:
        #    ener = np.inf
        #    flag = False
        return ener, flag
        
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
                self.param.chi[m2-1][m1-1] = p.value
            elif isinstance(p, parameters.KuhnLength):
                mon = p.monomerID
                self.kuhn[mon] = p.value
            elif isinstance(p, parameters.PolymerBlendFraction):
                self.phi_chain[p.polymerID] = p.value
            else:
                raise(NotImplementedError("{} is not an implemented parameter".format(type(p))))
        return True
    
    def _launchSim(self, root):
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
        # Create Input Files
        monFrac = self.paramWrap.getMonomerFractions()
        ngrid = self.param.ngrid
        lattice = self.paramWrap.getLattice()
        core_mon = self.coreOptions[0]
        for i in self.coreOptions:
            if monFrac[i] < monFrac[core_mon]:
                core_mon = i
        w = self.paramWrap.getInterfaceWidth(core_mon)
        newField = self.fieldGen.to_kgrid(monFrac,ngrid,interfaceWidth=w, coreindex=core_mon, lattice=lattice)
        kgrid = self.paramWrap.cleanFieldFile()
        kgrid.fields = newField
        kgridName = self.param.fieldTransforms[0][1] # Pull init guess filename from param file
        kgridFile = root / kgridName
        kgrid.write(kgridFile.open(mode='w'))
        paramFile = root / 'param'
        self.param.write(paramFile.open(mode='w'))
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
        
