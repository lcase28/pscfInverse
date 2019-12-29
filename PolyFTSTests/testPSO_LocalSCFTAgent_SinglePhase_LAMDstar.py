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
basedir="PSOTest_LAM_Dstar/"
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
pfac.ParseParameterFile('params_LAM_template.in')
# Set up the simulation launcher
launcher = SL.SimulationLauncher(pfac, PolyFTSbin, basedir)

# Parameters for PSO search
PSOparams = OrderedDict()
PSOparams["CellLengths"] = [2.5,5.5]

# Set up the agents
nagents = 10
agent_list=[]
for i in range(nagents):
    a = LocalSCFTAgent(["LAM"], launcher=launcher, PSOparameters=PSOparams)
    agent_list.append(a)
    print "Added agent {} to swarm".format(a.id)
# Initialize the swarm
cycle_graph = nx.cycle_graph(nagents)
swarm = Swarm(cycle_graph, agent_list)
# Output initial swarm configuration
swarm.write_output(basedir)

# Run PSO
nsteps = 50
plotter = spl.SwarmPlotter("LAMDstar", 1, nsteps, labels=["Lx (Rg)"], axislimits=[2.5,5.5], OptimalCoord=4.0431e+00, makemovie=True, OutputDir=basedir, fps=8)
for i in range(nsteps):
    plotter.updatePlot(agent_list, swarm.get_gbest().PBest, i)
    # Update all the agents
    starttime = time.time()
    for agent in swarm.Agents:
        agent.update()
    swarm.write_output(basedir)
    print "Step {} Finished in {} hours".format(i, (time.time() - starttime) / 3600)
plotter.finishMovie()
