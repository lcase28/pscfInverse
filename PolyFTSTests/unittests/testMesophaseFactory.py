#!/usr/bin/env python
from ParamFactories import *
import SimulationLauncher as SL
from PSO.mesophase_factory import MesophaseFactory

PolyFTSbin = "~/bin/PolyFTSGPU.x"
#PolyFTSbin = "~/bin/PolyFTS.x"

# Create a parameter factory from the template
pfac = PolyFTSFactory()
pfac.ParseParameterFile('params_template.in')

# Set the chiN and the block fractions
pfac.set(chiN12=24.0)
#pfac.set(BlockFractions1=[0.2])
#pfac.set(SimulationBK_NumBlocks=1) # Reduce the number of blocks to speed up the tests

# Generate all mesophases for a diblock sweep
phasekeys = ['BCC','HEX','GYR','O70','LAM','DIS']
phaseobjects = dict()
# Generate candidate mesophase objects
for phase in phasekeys:
    phaseobjects[phase] = MesophaseFactory.generate_candidate(phase)

resultfile = open("diblocksweep.dat","w")
resultfile.write("# fA ")
for phase in phasekeys:
    resultfile.write("{} ".format(phase))
resultfile.write("\n")

# Block fraction sweep
for fA in np.linspace(0.1, 0.5, 21):
    rootdir="DiblockPhaseSweep_chiN{}/fA_{}/".format(pfac.get("chiN12"),fA)
    resultfile.write("{} ".format(fA))
    pfac.set(BlockFractions1=[fA])
    print "Running phase sweep point = ",rootdir
    for phase in phasekeys:
        if os.path.isdir(rootdir+phase):
            print "Directory exists. Not running."
        else:
            #print "Phase key = ",phase
            print "- Input cell lengths = ",phaseobjects[phase].get_cell_lengths(None)
            #print "- Parameters = ",phaseobjects[phase].parameters
            #print "- Calling the launcher with phase's arguments: "
            SL.LaunchSimulation(pfac, rootdir+phase, PolyFTSbin, **phaseobjects[phase].parameters)
        # Fetch outputs from this run
        resultfile.write("{} ".format(phaseobjects[phase].get_H(simulationpath=rootdir+phase)))
        resultfile.flush()
        print "- Updated cell lengths for next run = ",phaseobjects[phase].get_cell_lengths(rootdir+phase)
        print "\n"
    print "\n\n"
    resultfile.write("\n")

resultfile.close()
