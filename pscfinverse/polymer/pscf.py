# LIBRARY IMPORTS
from pscfinverse.polymer.phases import MesophaseBase
from pscfinverse.polymer.variables import MesophaseVariable
import pscfinverse.polymer.parameters as parameters
import pscfinverse.util.contexttools as contexttools

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

class PscfIterationData:
    """
    Helper class to store data on an individual iteration step.
    """
    def __init__(self, itr, tot_res, res, stress, cell_param, f_helm):
        self.itr = itr
        self.error = tot_res
        self.res = res
        if len(stress) == len(cell_param):
            self.N_cell_param = len(cell_param)
            self.stress = stress
            self.cell_param = cell_param
        else:
            raise(ValueError("stress and cell_param must have the same number of elements."))
        self.f_helm = f_helm
    
    @classmethod
    def from_file(cls, f):
        """
        Initialize by reading directly from open, readable file f.
        
        If no next iteration is found, return None.
        """
        out = None
        foundEntry = False
        while not foundEntry:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            dat = line.split(maxsplit=1)
            if len(dat) > 0:
                if dat[0].strip() == "Iteration":
                    foundEntry = True
                    itr = int(dat[1].strip())
        def parseline(src):
            src = src.strip()
            dat = src.split("=",maxsplit=1)
            src = dat[1]
            src = src.strip()
            src = src.split()
            pout = []
            for s in src:
                pout.append(float(s.strip()))
            return pout
        if foundEntry:
            # SCF + stress residuals
            line = f.readline()
            if not line:
                return None
            tot_res = parseline(line)
            # SCF residuals
            line = f.readline()
            if not line:
                return None
            res = parseline(line)
            # stress
            line = f.readline()
            if not line:
                return None
            stress = parseline(line)
            # cell param
            line = f.readline()
            if not line:
                return None
            cell_param = parseline(line)
            # f_helmholtz
            line = f.readline()
            if not line:
                return None
            helm = parseline(line)
            out = cls(itr,tot_res[0],res[0],stress,cell_param,helm[0])
        return out
    
    @property
    def csv_headers(self):
        formstr = "{},{}"
        s = formstr.format("itr","error")
        s = formstr.format(s,"residuals")
        for i in range(self.N_cell_param):
            s = formstr.format(s,"stress_{}".format(i))
        for i in range(self.N_cell_param):
            s = formstr.format(s,"cell_param_{}".format(i))
        s = formstr.format(s,"f_Helmholtz")
        return "{}\n".format(s)
    
    @property
    def csv_data(self):
        formstr = "{},{}"
        s = formstr.format(self.itr, self.error)
        s = formstr.format(s,self.res)
        for i in range(self.N_cell_param):
            s = formstr.format(s,self.stress[i])
        for i in range(self.N_cell_param):
            s = formstr.format(s,self.cell_param[i])
        s = formstr.format(s,self.f_helm)
        return "{}\n".format(s)
    
    def write_csv(self, target, headers=False):
        """
        Write self to a CSV file.
        
        if headers is True, also write a header line.
        """
        if headers:
            target.write(self.csv_headers)
        target.write(self.csv_data)


class PscfIterateLog:
    """
    Class to read the log from a PSCF ITERATE command.
    
    NOTE: Log must have been saved to file.
    
    NOTE: Presently does not handle SWEEP operations.
    """
    
    def __init__(self, fname):
        """
        Initialize object by reading from file named fname.
        
        Argument expected_iterations is used to pre-allocate
        arrays for efficiency. Should be an integer if included.
        """
        self._iter = []
        with open(fname) as f:
            while True:
                newiter = PscfIterationData.from_file(f)
                if newiter is None:
                    break
                self._iter.append(newiter)
        self._itr = np.array([i.itr for i in self._iter])
        self._error = np.array([i.error for i in self._iter])
        self._res = np.array([i.res for i in self._iter])
        self._stress = np.array([i.stress for i in self._iter])
        self._cell = np.array([i.cell_param for i in self._iter])
        self._helm = np.array([i.f_helm for i in self._iter])
        
        stackvect = np.column_stack( (self._itr, self._error, self._res, self._helm) )
        left_stack = stackvect[:,0:-1]
        right_stack = stackvect[:,-1:]
        self._full = np.concatenate( (left_stack, self._stress, self._cell, right_stack), axis=1)
    
    @property
    def niter(self):
        return len(self._iter)
    
    @property
    def start_iteration(self):
        """ Return the iteration number that was given at the start of the log (0 or 1) """
        return self._itr[0]
    
    @property
    def end_iteration(self):
        """ Return the iteration number that was given for the last found iteration. """
        return self._itr[-1]
    
    @property
    def itr(self):
        return np.array(self._itr)
    
    @property
    def error(self):
        return np.array(self._error)
    
    @property
    def residual(self):
        return np.array(self._res)
    
    @property
    def stress(self):
        return np.array(self._stress)
    
    @property
    def N_cell_param(self):
        return self._iter[0].N_cell_param
    
    @property
    def cell_param(self):
        return np.array(self._cell)
    
    @property
    def f_Helmholtz(self):
        return np.array(self._helm)
    
    def csv_headers(self, prepend=""):
        prepend = prepend.strip()
        prepend = prepend.strip(",")
        return "{},{}\n".format(prepend,self._iter[0].csv_headers)
    
    def csv_data(self, prepend=""):
        s = ""
        prepend = prepend.strip()
        prepend = prepend.strip(",")
        for itr in self._iter:
            s = s + "{},{}\n".format(prepend,itr.csv_data)
        return s
    
    def csv_string(self, prepend_header="", prepend_data=""):
        head = self.csv_headers(prepend_header)
        data = self.csv_data(prepend_data)
        return "{}{}".format(head,data)
    
    def write_csv(self, target, prepend_data="", write_header=False, prepend_header=""):
        if write_header:
            target.write(self.csv_string(prepend_header, prepend_data))
        else:
            target.write(self.csv_data(prepend_data))
    
    def __len__(self):
        return self.niter
    
    def __array__(self):
        return np.array(self._full)
    

class PscfMesophase(MesophaseBase):
    """ A Mesophase manager for simulating in PSCF. """
    
    def __init__(self, ID, paramWrap, fieldGen, coreOptions=None, disorderCheckFilename=None, disorderCheckTolerance=0.001):
        """
        Initialize the PSCF Mesophase using pre-formatted files.
        
        IMPLEMENTATION NOTE: Presently, no checks are done to ensure 
        that all template files are compatible (that system definitions match).
        These checks are assumed to have been done by the user.
        
        In parameter definitions, 'fileManager' refers to pscfinverse.SCFT.PSCF.FileManagers
        
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
        disorderCheckFilename : string (optional)
            The name, within the SCFT run directory,
            of the coordinate grid file to check total
            variation range in final disorder check.
            Final check skipped if omitted.
        disorderCheckTolerance : float (optional)
            The minimum value of field variation that the field
            must exceed to be accepted. If the total variation in
            all components' fields are below this threshold,
            the result is rejected. Default is 0.001
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
        self._disorder_check_filename = disorderCheckFilename
        self._disorder_check_tolerance = disorderCheckTolerance
        super().__init__(ID)
    
    @classmethod
    def fromFieldGenFile(cls, ID, infile, coreOptions=None, disorderCheckFilename=None, disorderCheckTolerance=0.001):
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
        disorderCheckFilename : string (optional)
            The name, within the SCFT run directory,
            of the coordinate grid file to check total
            variation range in final disorder check.
            Final check skipped if omitted.
        disorderCheckTolerance : float (optional)
            The minimum value of field variation that the field
            must exceed to be accepted. If the total variation in
            all components' fields are below this threshold,
            the result is rejected.
        
        Returns
        -------
        newPhase : PscfMesophase
            An instance of PscfMesophase based on the input file.
        """
        infile = infile.resolve()
        newdir = str(infile.parent)
        with contexttools.cd(newdir):
            param, calc, outfile, cmon = read_input_file(infile)
        out = cls(ID, param, calc, coreOptions,disorderCheckFilename)
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
        # Check for disorder by field variation range
        if self._disorder_check_filename is not None:
            dcfname = root / self._disorder_check_filename
            gridfile = CoordFieldFile(dcfname.resolve())
            if np.amax(gridfile.valueRanges()) < self._disorder_check_tolerance:
                return np.inf, False
        # Calculate Energy
        ener = outfile.f_Helmholtz - outfile.f_homo
        return ener, True
        
    def setParams(self, VarSet):
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
                self.param.kuhn[mon] = p.value
            elif isinstance(p, parameters.PolymerBlendFraction):
                self.param.phi_chain[p.polymerID] = p.value
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
        
