import sys,os
import numpy as np
from ParamFactories import Section
import os


# ==============================================================
# Candidate Class
# ==============================================================
# This class allows a pfac to have mesophase-specific parts to be overwritten while keeping H and cell persistent from the last run
class Candidate(object):
    def __init__(self, key, seed_field, cell, **params):
        """
        :param key: string for phase key
        :param seed_field: string for location of seed field
        :param cell: ParameterFactories.Section object containing cell parameters
        :param params: Dictionary of parameters for overriding in ParameterFactory for this candidate phase
        """
        self.cell = cell
        self.key = key
        # cell is aliased to modelsBK_model1BK_cell
        self.parameters = dict(cell=cell, SeedField=seed_field)
        self.parameters.update(params)
        self.hamiltonian = 10000.

    def get_cell_lengths(self, simulationpath):
        # If we're not fetching from a finished simulation path then
        # just return the celllengths contained in the parameters.
        if simulationpath is None or self.parameters['VariableCell'] == 'false':
            return self.cell.get("CellLengths")

        # Get the number of spatial dimensions
        dims = int(self.cell.get("Dim"))
        cell_lengths = np.zeros(dims)
        # Get column names
        with open(simulationpath+"/operators.dat") as f:
            line = f.readline()
            columns = line.split()[1:]
            colmap = dict(zip(columns, range(len(columns))))

        # Load operators.dat
        opdata = np.loadtxt(simulationpath+"/operators.dat")

        # Load the final box tensor
        if dims > 1:
            k=1
            h=np.zeros(shape=(dims,dims))
            for i in range(dims):
                for j in range(dims):
                    h[i][j] = opdata[-1][colmap['CellTensor{}'.format(k)]]
                    k=k+1

            # Compute the length of each cell vector
            for i in range(dims):
                cell_lengths[i] = np.linalg.norm(h[i])
        else:
            cell_lengths[0] = opdata[-1][colmap['CellTensor']]

        # Update this object's record of the latest cell lengths
        # TODO:
        # Do not update the class/parameters record of celllengths if they are nonsense
        self.cell.set(CellLengths=cell_lengths)

        # TODO: extend to compute angles also. Currently low priority because all candidates have angles set by symmetry.
        return cell_lengths


    def get_H(self, parent=None, simulationpath=None):
        """
        Run an SCFT calculation with this candidate at the given point in parameter space and return its free energy

        :param point: :class:`Point` representing a point in parameter space
        :param parent: An `SCFTAgent` that can run SCFT simulations (implements RunSCFT_AutoDT)
        :return: Intensive free energy of the candidate morphology at `point`
        """

        if parent == None and simulationpath != None:
            # Find H from a simulation already ran
            self.hamiltonian = np.loadtxt(simulationpath+"/operators.dat")[-1][1]
            return self.hamiltonian


        candidate_sim = parent.Location.simulations[self.key]
        try:
            lastCell = self.get_cell_lengths(candidate_sim)
        except:
            print 'cant retrieve cellLengths for sim {} for candidate with key: {}'.format( candidate_sim,
                                                                                  self.key)
            dims = int(candidate_sim.get('Dim').values()[0])
            lastCell = 7.0*np.ones(dims)

        # if one of the lastCell elements is zero, change it to 5.
        try:
            # in try block because the first time it is executed dims is empty. 
            dims = int(candidate_sim.get('Dim').values()[0])
            for i in range(dims):
                if lastCell[i]==0:
                    lastCell[i]=5.0
                    print 'changed one of the lastCell elements in sim {} from zero to 5.'.format(candidate_sim)
        except:
            pass
        # sanitize the last cell input. If above 20, reduce it to 7.0
        if np.isnan(lastCell[0]):
            lastCell = 7.0*np.ones(dims)
        elif lastCell[0]>20.0:
            lastCell = 7.0/lastCell[0]*lastCell
            print 'sanitized cell input for {} phase in {} sim from large value'.format(self.key,candidate_sim)
        elif lastCell[0]< 1.0:
            print 'sanitized cell input for {} phase in {} sim from small value'.format(self.key,candidate_sim)
            lastCell = 3.0/lastCell[0]*lastCell

        self.cell.get("cell1").set(CellLengths=lastCell)

        parent.RunSCFT(parent.Location, self.key, **self.parameters)

        # If candidate_sim was None before, then it may have been set during the RunSCFT call.
        candidate_sim = parent.Location.simulations[self.key]

        try:
            hamiltonian_value = FTS.Load(candidate_sim, "operators.dat")[1][-1]

            print "Agent {}: H(chain = [{},{}])_{} = {}".format(parent.id,alpha0,alpha1,
                                                                  candidate_sim.get("mesophase_tag"),
                                                                      hamiltonian_value)
        except:
            print "couldn't get the hamiltonian value to load in ", candidate_sim
            sys.exit()
            pass

        return hamiltonian_value


# ==============================================================
# MesophaseFactory Class
# ==============================================================
class MesophaseFactory(object):
    @staticmethod
    def generate_candidate(mesophase_key):
        """
        Generate a Candidate mesophase

        :param mesophase_key:
        :return: Candidate object with initial guesses that should be suitable for converging each phase
        """

        seedfield = "phases/{}/fields_k.in".format(mesophase_key)
        extra_params = dict(ReadInputFields="true")

        if mesophase_key == "DIS":
            num_dimensions = 1
            spacegroup = dict(SpaceGroupIndex=2)
            cell_side_lengths = [1.0]
            angles = [90]
            threads = 1
            spacegroup = dict(SpaceGroupIndex=2)
            npw_array = [32]
            extra_params["ReadInputFields"]="false"
            seedfield=None
            for i in range(1, 2):
                extra_params["modelsBK_model1BK_initfieldsBK_initfield{}BK_InitType".format(i)] = "Zero"
            extra_params["VariableCell"]  = "false"

        if mesophase_key == "LAM":
            num_dimensions = 1
            spacegroup = dict(SpaceGroupIndex=2)
            cell_side_lengths = [4.0]
            angles = [90]
            npw_array = [64]
            threads = 1
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "HEX":
            num_dimensions = 2
            spacegroup = dict(SpaceGroupIndex=17)
            cell_side_lengths = [3.87] * 2
            angles = [120]
            npw_array = [32,32]
            threads = 2
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "SDM":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=227)
            cell_side_lengths = [6.7] * 3
            angles = [90,90,90]
            npw_array = [32,32,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "DDM":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=224)
            cell_side_lengths = [6.5] * 3
            angles = [90,90,90]
            npw_array = [32,32,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "A15":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=223)
            cell_side_lengths = [9.15] * 3
            angles = [90, 90, 90]
            npw_array = [48,48,48]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "GYR":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=230)
            cell_side_lengths = [8.4669] * 3
            angles = [90, 90, 90]
            npw_array = [32,32,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "O70":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=70)
            cell_side_lengths = [3.763, 7.526, 13.03]
            angles = [90, 90, 90]
            npw_array = [16,32,64]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "BCC":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=229)
            cell_side_lengths = [6.5] * 3
            angles = [90, 90, 90]
            npw_array = [32,32,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "FCC":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=225)
            cell_side_lengths = [7.0] * 3
            angles = [90, 90, 90]
            npw_array = [32,32,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "SGM":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=136)
            cell_side_lengths = [15.6,15.6,8.6]
            angles = [90, 90, 90]
            npw_array = [64,64,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "HCP":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=194)
            cell_side_lengths = [5.0,5.0,8.25]
            angles = [90, 90, 120]
            npw_array = [32,32,32]
            threads = 4
            extra_params["VariableCell"]  = "true"

        if mesophase_key == "LVS":
            num_dimensions = 3
            spacegroup = dict(SpaceGroupIndex=227)
            cell_side_lengths = [13.5]*3
            angles = [90, 90, 90]
            npw_array = [48,48,48]
            threads = 4
            extra_params["VariableCell"]  = "true"

        return Candidate(mesophase_key, seedfield,
                         cell=Section(Dim=num_dimensions, CellLengths=cell_side_lengths, CellAngles=angles, SpaceGroupIndex=spacegroup['SpaceGroupIndex'], NPW=npw_array, Symmetrize="on"),
                         OMPThreads=threads, **extra_params)
