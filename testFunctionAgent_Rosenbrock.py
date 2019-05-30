#!/usr/bin/env python

# Test the PSO algorithm by optimizing both the Rosenbrock function and the Greiwank function
from PSO.Swarm import *
import networkx as nx
import numpy as np
from scipy.optimize import rosen # Rosenbrock function
import PSO.SwarmPlot as spl
import time

# Function to be optimized
def optfn(x):
    return -rosen(x)

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
nagents = 20
boundaries = np.array([[-2,2],[-1,3]])
#initialvelocity = boundaries * 0.1
agent_list=[]
for i in range(nagents):
    a = FunctionAgent(optfn, boundaries=boundaries)
    agent_list.append(a)
    print "Added agent {} to swarm".format(a.id)
# Initialize the swarm
cycle_graph = nx.cycle_graph(nagents)
swarm = Swarm(cycle_graph, agent_list)

# Run PSO
nsteps = 200
plotter = spl.SwarmPlotter("Rosenbrock_PSO", 2, nsteps, labels=["x", "y"], axislimits=boundaries, OptimalCoord=[1,1], makemovie=True, fps=8)
for i in range(nsteps):
    plotter.updatePlot(agent_list, swarm.get_gbest().PBest, i, overlaycallbackfn=plotoptfn, **{'optfn': optfn})
    # Update all the agents
    starttime = time.time()
    for agent in swarm.Agents:
        print "Agent {}: {}={}".format(agent.id, agent.get_coords(), agent.Location.Fitness)
        agent.update()
    print "Step {} Finished in {} hours".format(i, (time.time()-starttime)/3600)
plotter.finishMovie()
