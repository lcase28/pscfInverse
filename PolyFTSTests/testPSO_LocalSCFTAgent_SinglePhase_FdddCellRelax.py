#!/usr/bin/env python

# Test the PSO algorithm by optimizing both the Rosenbrock function and the Greiwank function
import glob, os
from ParamFactories import *
import SimulationLauncher as SL
from PSO.Swarm import *
import networkx as nx
import numpy as np
import PSO.SwarmPlot as spl
import time
import shutil

# Base directory for storing all simulations and outputs
basedir="PSOTest_Fddd_CellRelax/"
# Start with a clean set of directories if run already exists
if os.path.isdir(basedir+"PBests"):
    shutil.rmtree(basedir+"PBests")
for dirname in glob.glob(basedir+"step*"):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
for dirname in glob.glob(basedir+"GBest*"):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)

# Prepare for PolyFTS runs
PolyFTSbin = "~/bin/PolyFTSGPU.x"

pfac = PolyFTSFactory()
pfac.ParseParameterFile('params_template.in')
pfac.set(VariableCell=False) # Turn off automatic variable cell
pfac.set(modelsBK_model1BK_operatorsBK_CalcStressTensor=False) # Speed up operators computation
pfac.set(SpaceGroupIndex=70) # Enforce O70 symmetry group
pfac.set(NPW=[16, 32, 64]) # Tweak resolution to be useful for the orthorhombic cell
# Run where Fddd is very stable
pfac.set(BlockFractions1=0.44)
pfac.set(chiN12=11.5)
# Optimal target at these parameters:
# [3.836 7.6416 13.445]

# Set up the simulation launcher
launcher = SL.SimulationLauncher(pfac, PolyFTSbin, basedir)

# Parameters for PSO search
PSOparams = OrderedDict()
boundaries = np.array([[3,8],[3,8],[8,15]])
# Note the trailing _idx for setting elements of an array in
# the parameter factory
PSOparams["CellLengths_0"] = boundaries[0]
PSOparams["CellLengths_1"] = boundaries[1]
PSOparams["CellLengths_2"] = boundaries[2]

# TODO: Next step is to do MorphologyAgent but to optimize over variables that are not constrained (e.g., fA in diblock and chiN)
#       Then we don't need the constraints setting.

# TODO: REMOVE ME - accelerating the SCFT just to test plotting and output
#pfac.set(simulationBK_NumBlocks=5)

# Set up the agents
nagents = 10
agent_list=[]
for i in range(nagents):
    a = LocalSCFTAgent(["O70"], launcher=launcher, PSOparameters=PSOparams, SeedFields={"O70": "phases/O70/fields_k.in"})
    agent_list.append(a)
    print "Added agent {}".format(a.id)
# Initialize the swarm
cycle_graph = nx.cycle_graph(nagents)
swarm = Swarm(cycle_graph, agent_list)
# Output initial swarm configuration
swarm.write_output(basedir)

# Run PSO
nsteps = 50
plotter = spl.SwarmPlotter("Fddd_CellOpt_PSO", 3, nsteps, labels=["Lx (Rg)", "Ly (Rg)", "Lz (Rg)"], axislimits=boundaries, OptimalCoord=[3.836,7.6416,13.445], makemovie=True, OutputDir=basedir, fps=8)
for i in range(nsteps):
    plotter.updatePlot(agent_list, swarm.get_gbest().PBest, i)
    # Update all the agents
    starttime = time.time()
    for agent in swarm.Agents:
        agent.update()
    swarm.write_output(basedir)
    print "Step {} Finished in {} hours".format(i, (time.time() - starttime) / 3600)
plotter.finishMovie()
