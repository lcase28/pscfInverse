#!/usr/bin/env python

# Test the PSO algorithm by optimizing both the Rosenbrock function and the Greiwank function
from context import psoinverse
from psoinverse.PSO.Swarm import *
import networkx as nx
import numpy as np
import psoinverse.PSO.SwarmPlot as spl
from psoinverse.PSO.Integrators import StandardIntegrator as stdInt
from psoinverse.PSO.SearchSpace import SearchBounds as Bounds
import time

# Function to be optimized
def optfn(x):
    # The Griewank function for arbitrary dimensions
    f = 1. + 1./4000 * np.sum(x*x) - np.product(np.cos(x))
    return -f

# A callback function to plot something onto the dynamic agent plot
def plotoptfn(mplaxis, **kwargs):
    N=150
    xlim = mplaxis.get_xlim()
    ylim = mplaxis.get_ylim()
    # Generate the mesh
    xlist = np.linspace(xlim[0], xlim[1], N)
    ylist = np.linspace(ylim[0], ylim[1], N)
    X, Y = np.meshgrid(xlist, ylist)
    # X, Y are 2D numpy arrays. optfn() expects to receive a 1D numpy array containing numplotcoords entries for each call (e.g., [x,y] of a single point)
    XY=np.reshape(np.array([X,Y]).T,(N,N,2))
    Z = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            Z[j][i] = kwargs['optfn'](XY[i][j])
    levels=np.linspace(np.min(Z),np.max(Z),100)
    #mplaxis.contour(X,Y,Z, cmap='cool',levels=levels,alpha=0.4)
    mplaxis.contourf(X,Y,Z, cmap='viridis',alpha=0.3,levels=levels)
    mplaxis.contour(X,Y,Z,15,alpha=0.3,colors='K')

# Set up the agents
# TODO: consider removing coordinate scaling entirely - it doesn't play a useful role and adds complexity
#  Ask Sean about why it is present and determine whether it helps in setting initial velocities.
nagents = 20
axislimit=100.0
ndim = 2
#boundaries = np.array([[-axislimit,axislimit],[-axislimit,axislimit]])
boundaries = Bounds(np.full(ndim,-axislimit),np.full(ndim,axislimit))
#initialvelocity = boundaries * 0.1;
agent_list=[]
for i in range(nagents):
    a = FunctionAgent(optfn, boundaries=boundaries, v0=None, spawnRange=None, useScale=np.ones(ndim))
    agent_list.append(a)
    print("Added agent {} to swarm".format(a.id))
# Initialize the swarm
cycle_graph = nx.cycle_graph(nagents)
integrator = stdInt(seekMax=True)
swarm = Swarm(cycle_graph, agent_list, integrator) #[0.729, 2.05, 2.05, True], True)

# Run PSO
nsteps = 200
#plotter = spl.SwarmPlotter("Griewank_PSO", 2, nsteps, labels=["x", "y"], axislimits=boundaries, OptimalCoord=[0,0], makemovie=True, fps=8)
for i in range(nsteps):
    #plotter.updatePlot(agent_list, swarm.get_gbest().PBest, i, overlaycallbackfn=plotoptfn, **{'optfn':optfn})
    # Update all the agents
    starttime = time.time()
    swarm.printState()
    swarm.step()
    #for agent in swarm.Agents:
    #    print("Agent {}: {}={}".format(agent.id, agent.get_coords(), agent.Location.Fitness))
    #    agent.update()
    print("Step {} Finished in {} hours".format(i, (time.time()-starttime)/3600))
#plotter.finishMovie()
