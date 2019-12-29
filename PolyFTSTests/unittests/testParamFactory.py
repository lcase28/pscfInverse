#!/usr/bin/env python

from ParamFactories import *
import SimulationLauncher as SL

#PolyFTSbin = "~/bin/PolyFTSGPU.x"
PolyFTSbin = "~/bin/PolyFTS.x"

pfac = PolyFTSFactory()
pfac.ParseParameterFile('params_template.in')

############# SINGLE PARAMETER SET #################
# Example of setting a key entry (using a short-hand Alias from ParameterFactories.py)
pfac.set(chiN12=30.0)
## Example of getting a key entry
#chiABN = pfac.get("chiN12")
#print "New chiABN = ",chiABN
#print "Writing param file to --> ",pfac.PrintCurrent("testlocalparam.in")
# Can also get the entire inputfile as a string using pfac.__str__()

# Test running this simulation (single)
#SL.LaunchSimulation(pfac, "./agent0_testrun", PolyFTSbin)

############ MULTPLE PARAMETERS IN A SWEEP ###########
# Test adding a 2D sweep
pfac.Sweep(chiN12=np.linspace(10.,20.,6))
pfac.Sweep(BlockFractions1=np.linspace(2,5,4)*0.1)
print("Number of simulations = ",pfac.GetTotalSimulations())
#pfac.PrintGraph() # Writes a tree graph representation of the sweep to Graph.png

#########################
# Test launching. This function will be called by the simulationlauncher class, but we can conceivably simplify and launch directly.
# This function is for launching SWEEPS
pfac.LaunchSimulations(SL.LaunchSimulation,PolyFTSbin,"diblocksweep/")

print("Done testing")
