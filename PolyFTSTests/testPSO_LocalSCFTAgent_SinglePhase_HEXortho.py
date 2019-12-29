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
basedir="PSOTest_HEXortho_CellOpt/"
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
PolyFTSbin = "~/bin/PolyFTS.x"
pfac = PolyFTSFactory()
pfac.ParseParameterFile('params_HEXortho_template.in')
pfac.set(VariableCell=False) # Turn off automatic variable cell
pfac.set(modelsBK_model1BK_operatorsBK_CalcStressTensor=False) # Speed up operators computation
# Run where Fddd is very stable
pfac.set(BlockFractions1=0.36)
pfac.set(chiN12=14.5)

# Set up the simulation launcher
launcher = SL.SimulationLauncher(pfac, PolyFTSbin, basedir)

# Parameters for PSO search
PSOparams = OrderedDict()
boundaries = np.array([[2.5,8.5],[2.5,8.5]])
PSOparams["CellLengths_0"] = boundaries[0]
PSOparams["CellLengths_1"] = boundaries[1]


# Set up the agents
nagents = 10
agent_list=[]
for i in range(nagents):
    a = LocalSCFTAgent(["HEXortho"], launcher=launcher, PSOparameters=PSOparams, SeedFields={"HEXortho": "phases/HEXortho/fields_k.in"})
    agent_list.append(a)
    print "Added agent {}".format(a.id)
# Initialize the swarm
cycle_graph = nx.cycle_graph(nagents)
swarm = Swarm(cycle_graph, agent_list)
# Output initial swarm configuration
swarm.write_output(basedir)

# Run PSO
nsteps = 50
# Note that we have two degenerate OptimalCoord solutions presented below
plotter = spl.SwarmPlotter("HEXortho_CellOpt_PSO", 2, nsteps, labels=["Lx (Rg)", "Ly (Rg)"], axislimits=boundaries, OptimalCoord=[[3.9802e+00, 6.8953e+00],[6.8953e+00, 3.9802e+00]], makemovie=True, OutputDir=basedir, fps=8)
for i in range(nsteps):
    plotter.updatePlot(agent_list, swarm.get_gbest().PBest, i)
    # Update all the agents
    starttime = time.time()
    for agent in swarm.Agents:
        agent.update()
    swarm.write_output(basedir)
    print "Step {} Finished in {} hours".format(i, (time.time() - starttime) / 3600)
plotter.finishMovie()
