# Dask Cluster Setup
from dask_mpi import initialize
initialize()

# Standard Library imports
import argparse
import distributed
import numpy as np
import networkx as nx
import os
import pathlib
import sys
import time

## Project Imports
from pscfinverse.pso.core import Velocity, OptimizationType, FITNESS_SELECTOR
from pscfinverse.pso.integrators import StandardIntegrator
from pscfinverse.pso.swarm import Swarm
from pscfinverse.util.parallelBatches import DaskCalculationManager
from pscfinverse.polymer.variables import (
    BlockRatioVariable,
    BlockLengthVariable,
    ChiVariable, 
    KuhnVariable,
    KuhnRatioVariable,
    PolymerVariableSet )
from pscfinverse.polymer.phases import MesophaseManager
from pscfinverse.polymer.pscf import PscfMesophase
from pscfinverse.polymer.scftAgent import ScftAgent
from pscfinverse.util.iotools import FileParser

client = DaskCalculationManager()

def wrongKey(expected,got):
    msg = "Expected key '{}'; got '{}'"
    raise(ValueError(msg.format(expected,got)))

def parseSearchSpace(words,key):
    """
    Read the input file to collate variables and constraints.
    """
    if not (key == "SearchSpace{"):
        raise(ValueError("Invalid key for SearchSpace block {}".format(key)))
    word = next(words)
    if not (key == "Variables{"):
        wrongKey("Variables{",word)
    
    
        

if __name__ == '__main__':
    ### Operation
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--file","-f", type=str,required=True)
    args = parser.parse_args()
    filepath = pathlib.Path(args.file)
    
    with FileParser(filepath) as words:
        
        # Check block opener
        word = next(words)
        if not word == "PscfInverse{":
            wrongKey("PscfInverse{",word)
        
        # Set core parameters
        FITNESS_SELECTOR.optimizationType = OptimizationType.maximize
        
        word = next(words)
        expected = "SearchSpace{"
        if not word == expected:
            wrongKey(expected,word)
        else:
            variables, constants = parseSearchSpace(words,word)
        
        # Create Variables and VariableSet
        variables = []
        variables.append( BlockRatioVariable((0,0),(0,1), value=0.0, lower=-2.5, upper=2.5 ) )
        print("Made {}".format(variables[0]))
        variables.append( ChiVariable( 0, 1, value=15, lower=10, upper=20 ) )
        print("Made {}".format(variables[1]))
        variables.append( KuhnVariable( 0, value=1.5, lower=1.0, upper=2.0 ) )
        print("Made {}".format(variables[2]))
        
        # Create Constants
        constants = []
        constants.append( BlockLengthVariable( [ (0,0),(0,1) ], value=1.0, upper=2.0 ) )
        
        varset = PolymerVariableSet(variables,constants)
        print("Made Variable Set")
        
        # Create Mesophases and Phase Managers
        pwd = pathlib.Path.cwd()
        inRoot = pwd / "inputs"
        
        p = "sigma"
        tgt = PscfMesophase.fromFieldGenFile(p,inRoot/p/"model")
        print("Made target {}".format(p))
        
        phases = []
        phaseNames = [ "a15", "fcc", "bcc", "gyr", "hex", "lam" ]
        for p in phaseNames:
            phases.append(PscfMesophase.fromFieldGenFile(p,inRoot/p/"model"))
            print("Made mesophase {}".format(p))
        
        pman = MesophaseManager(phases, tgt)
        print("Made MesophaseManager")
        
        # Create Integrator
        randGen = np.random.RandomState(300)
        integrator = StandardIntegrator(randGen)
        print("Created Integrator")
        
        # Set Run parameters
        nagent = 10
        nphases = 7
        nstep = 100
        
        # Create Agents
        velVal = np.ones(3)
        velCap = np.array([2.5, 6, 0.6])
        velsrc = Velocity(velVal,velCap)
        print("Created Velocity")
        agentList = ScftAgent.createSeveral(nagent,varset,pman,velsrc,client,pwd)
        print("Created Agents")
        agentList = ScftAgent.randomizeSeveral(agentList,randGen)
        print("Randomized Agents")
        
        # Create Swarm
        graph = nx.cycle_graph(nagent)
        sm = Swarm(graph, agentList, integrator, pwd)
        print("Created Swarm")
        
        ## Run iterations
        tryCooldown = 0.1
        print("StartingRun")
        runtime, proctime = sm.run(nstep,tryCooldown)
        print("Finished Running:\trun {}s;\tproc {}s".format(runtime,proctime))

client.close()
sys.stdout.flush()
