import numpy as np
import threading
from PSO.Swarm import LocalSCFTAgent
from PSO.mesophase_factory import MesophaseFactory
import os
import sys
from SimulationLauncher import getSimulationStatus

from pprint import pprint

# Import PolyFTSIO module for field IO
sys.path.insert(0, "/home/kdelaney/Work/research/code/Polymers/FullFTS/PolyFTS_develop/tools/") # Add search path for PolyFTSIO module
#sys.path.insert(0, "/Users/kdelaney/Work/research/code/Polymers/FullFTS/PolyFTS_develop/tools/") # Add search path for PolyFTSIO module
from PolyFTSIO import *
#import plot.FieldPlot as ftpl

# Adapted by Kris Delaney from earlier code by Sean Paradiso and Mihir Khadilkar

# ==============================================================
# MorphologyAgent Class
#  Take a list of phases and compute a single-agent fitness from free energies of all phases
# ==============================================================
class MorphologyAgent(LocalSCFTAgent):
    def __init__(self, target, alternative_candidates, **kwargs):
        """
        Create a MorphologyAgent.

        :param target: Candidate for the target mesophase.
        :param alternative_candidates: list of alternative `Candidate` phases.
        :param kwargs: PSO parameters. Check out `SCFTAgent`.__init__ for details/options
        """

        self.target = MesophaseFactory.generate_candidate(target)
        self.candidates = [MesophaseFactory.generate_candidate(c) for c in alternative_candidates]

        #SeedFields = {c.key : c.parameters['SeedField'] for c in self.candidates + [self.target]}

        # Here we have created self.target and self.candidates, which are objects of type class Candidate.
        # But to the baseclass LocalSCFTAgent we pass only strings of the simulation names to be used as keys.
        #
        # We need to override any function that can call RunSCFT and pass the settings from the Candidate object to the launcher.

        super(MorphologyAgent, self).__init__(simulation_keys=[target] + alternative_candidates, **kwargs)


    def ComputeFitness(self):
        """
        Computes the fitness as the difference in free energy between the target
        mesophase and the lowest free energy mesophase among the alternate candidates.
        :return: H_target - min(H_candidates)
        """

        # Run all SCFT simulations
        # In some cases, we need to check for failures and assign a bad fitness, and only failures
        # of the target phase may matter. In that case we can avoid running all of the other candidates.
        #
        # Run target first
        self.RunSCFT(self.Location, self.target.key, **self.target.parameters)

        # Checks to be made
        # * target -> DIS
        # * target -> supergroup
        # * target sim status == FAILED
        # * target cell exploded to very large / very small values

        # Check whether the target morphology melted to DIS and assign a bad fitness if it did
        if checkDISCollapse(self.Location.simulations[self.target.key]):
            print " [FAIL]: Target melted to DIS for agent {}\n".format(self.id)
            return -1E6

        # Read field data in k space
        fielddata = PolyFTSFieldReader()
        simpath = self.Location.simulations[self.target.key]
        if os.path.isfile(simpath+"/fields_k.bin"):
            fielddata.readFields(simpath+"/fields_k.bin")
        elif os.path.isfile(simpath+"/fields_k.dat"):
            fielddata.readFields(simpath+"/fields_k.dat")
        else:
            print " [FAIL]: Error reading target simulation density file\n"
            return -1E6

        # Check that SCFT converged to the same phase that was seeded
        if checkStructureSimilarity(fielddata, self.target.parameters['SeedField']):
            print " [FAIL]: Target converged to a different phase for agent {}".format(self.id)
            print "         Path = {}\n".format(self.Location.simulations[self.target.key])
            return -1E6

        # Check for status == failed
        if getSimulationStatus(self.Location.simulations[self.target.key]) == 'FAILED':
            print " [FAIL]: Target simulation failed for agent {}\n".format(self.id)
            return -1E6

        cell_last = getCellMaxVecLength(fielddata)
        if cell_last > 30.0:
            print " [FAIL]: The target phase simulation cell grew beyond 30 Rg\n"
            return -1E6

        #
        #
        # Target did not have any obvious failures. Run the comparison candidates.
        for c in self.candidates:
            self.RunSCFT(self.Location, c.key, **c.parameters)

        # We will now compute the fitness based on free energies. Fetch for all phases.
        # We don't need to check for failures - any structure that emerges with a competitive free energy
        # should be used to assess the fitness of the target.
        H_values = {}
        for c in [self.target] + self.candidates:
            H_values[c.key] = self.GetH(c.key)

        # Prepare for next runs:
        #   Update guess cell (within limits) if a given structure didn't fail and it didn't have a supergroup collapse



        # PSO *maximizes* the objective function
        # define dH to be positive when target has lower free energy than most stable alternative candidate
        dH = sorted([H_values[c.key] for c in self.candidates])[0] - H_values[self.target.key]
        print " All free energies = {};  Fitness = {}\n".format(H_values, dH)
        return dH


###################################

def checkDISCollapse(simpath):
    fielddata = PolyFTSFieldReader()
    if os.path.isfile(simpath+"/density.bin"):
        fielddata.readFields(simpath+"/density.bin")
    elif os.path.isfile(simpath+"/density.dat"):
        fielddata.readFields(simpath+"/density.dat")
    else:
        print "Error reading density file"
        return True

    total_variance = 0.
    for i in range(fielddata.nfields):
        total_variance += np.var(fielddata.AllFields[i])

    return total_variance < 1E-2


def checkStructureSimilarity(fielddata,ref_field=None,threshold=0.75):
    if ref_field==None:
        ref_field = "phases/"+sim.get('mesophase_tag')+"/fields_k.dat"
    fielddata2 = PolyFTSFieldReader()
    try:
        fielddata2.readFields(ref_field)
    except:
        print "Reference field = {}".format(ref_field)
        sys.exit("Error reading reference field")

    # Normalize field 0 from phase 1
    vec1 = fielddata.AllFields[0]
    vec1 = vec1 / np.linalg.norm(vec1)

    # Normalize field 0 from phase 2
    vec2 = fielddata2.AllFields[0]
    vec2 = vec2 / np.linalg.norm(vec2)

    return abs(np.dot(vec1,vec2))<threshold


def getCellMaxVecLength(fielddata):
    # Process the cell tensor.
    # Take the length of each cell vector and return the max.
    maxlen = 0.
    for i in range(fielddata.Dim):
        #print "Cell vector {} = {}".format(i, fielddata.hcell[i])
        #print "  Length = {}".format(np.linalg.norm(fielddata.hcell[i]))
        maxlen = max(maxlen, np.linalg.norm(fielddata.hcell[i]))

    return maxlen
