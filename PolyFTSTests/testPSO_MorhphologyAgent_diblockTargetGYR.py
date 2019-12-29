#!/usr/bin/env python

import glob, os
from ParamFactories import *
import SimulationLauncher as SL
from PSO.Swarm import *
from PSO.morphology_agent import MorphologyAgent
import networkx as nx
import numpy as np
import PSO.SwarmPlot as spl
import time
import shutil

# Base directory for storing all simulations and outputs
basedir="PSOTest_FindTargetMorphology_ABdiblockTargetGYR/"
# Start with a clean set of directories if run already exists
if os.path.isdir(basedir+"PBests"):
    shutil.rmtree(basedir+"PBests")
for dirname in glob.glob(basedir+"step*"):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
for dirname in glob.glob(basedir+"GBest*"):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)

# Function to render the AB phase boundaries in the chiN, fA domain.
def PlotABPhaseDiagram(mplaxis):
  # Precomputed phase boundary lines
  GL = np.array([[20.00,0.3745], [19.00,0.3789], [18.00,0.3837], [17.00,0.3893], [16.00,0.3955], [15.00,0.4030], [14.00,0.4110], [13.63,0.4144]])
  OL = np.array([[13.63,0.4144], [13.50,0.4160], [13.00,0.4222], [12.00,0.4369], [11.50,0.4468], [11.00,0.4598], [10.75,0.4702], [10.495,0.5]])
  GO = np.array([[13.63,0.4144], [13.50,0.4151], [13.00,0.4183], [12.00,0.4270], [11.75,0.4299], [11.58,0.4325]])
  HG = np.array([[20.00,0.3378], [19.00,0.3417], [18.00,0.3465], [17.00,0.3522], [16.00,0.3593], [15.00,0.3683], [14.00,0.3800], [13.00,0.3958], [12.00,0.4187], [11.75,0.4266], [11.58,0.4325]])
  HO = np.array([[11.58,0.4325], [11.50,0.4349], [11.00,0.4532], [10.75,0.4667], [10.495,0.5]])
  BH = np.array([[20.00,0.2436], [19.00,0.2518], [18.00,0.2614], [17.00,0.2722], [16.00,0.2853], [15.00,0.2999], [14.00,0.3183], [13.00,0.3410], [12.00,0.3721], [11.50,0.3950], [11.00,0.4245], [10.75,0.4470], [10.495,0.50 ]])
  DB = np.array([[20.00,0.2117], [19.00,0.2218], [18.00,0.2328], [17.00,0.2456], [16.00,0.2610], [15.00,0.2776], [14.00,0.2986], [13.00,0.3256], [12.00,0.3628], [11.50,0.3856], [11.00,0.4182], [10.75,0.4435], [10.495,0.5]])
  # Plot boundaries
  for b in [GL, OL, GO, HG, HO, BH, DB]:
      # fA < 0.5
      mplaxis.plot(b.T[1],b.T[0],c='#444444',alpha=0.6,ls='--')
      # fA > 0.5
      mplaxis.plot(1.-b.T[1],b.T[0],c='#444444',alpha=0.6,ls='--')


# Prepare for PolyFTS runs
PolyFTSbin = "~/bin/PolyFTSGPU.x"
pfac = PolyFTSFactory()
pfac.ParseParameterFile('params_template.in')

# Set up the simulation launcher
launcher = SL.SimulationLauncher(pfac, PolyFTSbin, basedir)

# Parameters for PSO search
PSOparams = OrderedDict()
boundaries = np.array([[0.2,0.5],[10.0,18.0]])
# Note the trailing _idx for setting elements of an array in
# the parameter factory
PSOparams["BlockFractions1"] = boundaries[0]
PSOparams["chiN12"]          = boundaries[1]

###########################
# TODO: ENSURE THAT INITIAL CELL VECTOR IS CHAINED AND CONSIDER CHAINING FIELDS TOO.
# TODO: consider moving PSO updates from Agent to Swarm. Swarm keeps a neighbor list.
###########################

# Set up the agents
# Morphology target
target = "GYR"
# Alternative candidate library
candidate_list = ["LAM", "DIS", "HEX", "O70"]
nagents = 10
#target = "LAM"
#candidate_list = ["DIS"]
#nagents = 10
agent_list=[]
for i in range(nagents):
    a = MorphologyAgent(target, candidate_list, launcher=launcher, PSOparameters=PSOparams)
    agent_list.append(a)
    print "Added agent {} to swarm".format(a.id)
# Initialize the swarm
cycle_graph = nx.cycle_graph(nagents)
swarm = Swarm(cycle_graph, agent_list)
# Output initial swarm configuration
swarm.write_output(basedir)

# Run PSO
nsteps = 50
plotter = spl.SwarmPlotter("PSO_ABdiblockTargetGYR", 2, nsteps, labels=[r'$f_A$', r'$\chi N$'], axislimits=boundaries, makemovie=True, OutputDir=basedir, fps=5)
for i in range(nsteps):
    plotter.updatePlot(agent_list, swarm.get_gbest().PBest, i, overlaycallbackfn=PlotABPhaseDiagram)
    # Update all the agents
    starttime = time.time()
    for agent in swarm.Agents:
        agent.update()
    swarm.write_output(basedir)
    print "Step {} Finished in {} hours\n\n".format(i, (time.time() - starttime) / 3600)
plotter.finishMovie()

