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

## Helper methods

def wrongKey(expected,got):
    msg = "Expected key '{}'; got '{}'"
    raise(ValueError(msg.format(expected,got)))

def checkKey(expected,got):
    if not (got == expected):
        wrongKey(expected,got)

def isblockend(word):
    return (word == "}")

def ensureblockend(word):
    if not isblockend(word):
        raise(ValueError("Expected end of block. Got {}.".format(word)))

def closeblock(words):
    word = next(words)
    ensureblockend(word)

## Search Space definition methods

def parseSearchSpace(words,key):
    """
    Read the input file to collate variables and constraints.
    """
    if not (key == "SearchSpace{"):
        raise(ValueError("Invalid key for SearchSpace block {}".format(key)))
    word = next(words) # Checking for start of variables block
    if not (word == "Variables{"):
        wrongKey("Variables{",word)
    variables = []
    velcap = []
    endblock = False
    while not endblock:
        word = next(words)
        if word == "}":
            endblock = True
        else:
            v1, v2 = parseVariable(words,word)
            variables.append(v1)
            velcap.append(v2)
    word = next(words) # checking for start of constraints block
    if not (word == "Constraints{"):
        wrongKey("Constraints{",word)
    constraints = []
    endblock = False
    while not endblock:
        word = next(words)
        if isblockend(word):
            endblock = True
        else:
            constraints.append(parseConstraint(words,word))
    closeblock(words)
    return variables, constraints, velcap

def parseVarBounds(words,word):
    checkKey("lower",word)
    lower = words.next_float()
    word = next(words)
    checkKey("upper",word)
    upper = words.next_float()
    word = next(words)
    ckeckKey("velocity_cap",word)
    velcap = words.next_float()
    initval = lower + ( (upper - lower) / 2 )
    closeblock(words)
    return lower, upper, velcap, initval

def parseVariable(words,key):
    """
    Read input data for search variable and return the variable object.
    """
    if key == "BlockLength{":
        blocks = []
        word = next(words)
        checkKey("block",word)
        while word == "block":
            pnum = words.next_int()
            bnum = words.next_int()
            blocks.append((pnum,bnum))
            word = next(words)
        lower, upper, velcap, initval = parseVarBounds(words, word)
        varobj = BlockLengthVariable(blocks=blocks,value=initval, lower=lower, upper=upper)
    elif key == "BlockRatio{":
        nblocks = []
        dblocks = []
        # numerator
        word = next(words)
        checkKey("Numerator{",word)
        word = next(words)
        checkKey("block",word)
        while word == "block":
            pnum = words.next_int()
            bnum = words.next_int()
            nblocks.append( (pnum,bnum) )
            word = next(words)
        ensureblockend(word)
        # denominator
        word = next(words)
        checkKey("Denominator{",word)
        word = next(words)
        checkKey("block",word)
        while word == "block":
            pnum = words.next_int()
            bnum = words.next_int()
            dblocks.append( (pnum,bnum) )
            word = next(words)
        ensureblockend(word)
        word = next(words)
        lower, upper, velcap, initval = parseVarBounds(words, word)
        varobj = BlockRatioVariable(nblocks, dblocks, initval, lower, upper)
    elif key == "Chi{":
        word = next(words)
        checkKey("monomers", word)
        mon1 = words.next_int()
        mon2 = words.next_int()
        word = next(words)
        lower, upper, velcap, initval = parseVarBounds(words, word)
        varobj = ChiVariable(mon1, mon2, initval, lower, upper)
    elif key == "KuhnLength{":
        word = next(words)
        checkKey("monomer", word)
        mon1 = words.next_int()
        word = next(words)
        lower, upper, velcap, initval = parseVarBounds(words, word)
        varobj = KuhnVariable(mon1, initval, lower, upper)
    elif key == "KuhnRatio{":
        word = next(words)
        checkKey("monomers", word)
        mon1 = words.next_int()
        mon2 = words.next_int()
        word = next(words)
        lower, upper, velcap, initval = parseVarBounds(words, word)
        varobj = KuhnRatioVariable(mon1, mon2, initval, lower, upper)
    else:
        raise(ValueError("Invalid key for search variable declaration, {}.".format(key)))
    return varobj, velcap

def parseConstraintValue(words,word):
    checkKey("value",word)
    value = words.next_float()
    lower = value - 0.1
    upper = value + 0.1
    closeblock(words)
    return value, lower, upper
    
def parseConstraint(words,key):
    """
    Read input data for search constraint and return the variable object.
    """
    if key == "BlockLength{":
        # block length
        blocks = []
        word = next(words)
        checkKey("block",word)
        while word == "block":
            pnum = words.next_int()
            bnum = words.next_int()
            blocks.append((pnum,bnum))
            word = next(words)
        initval, lower, upper = parseConstraintValue(words,word)
        varobj = BlockLengthVariable(blocks=blocks,value=initval, lower=lower, upper=upper)
    elif key == "BlockRatio{":
        # block ratio
        nblocks = []
        dblocks = []
        # numerator
        word = next(words)
        checkKey("Numerator{",word)
        word = next(words)
        checkKey("block",word)
        while word == "block":
            pnum = words.next_int()
            bnum = words.next_int()
            nblocks.append( (pnum,bnum) )
            word = next(words)
        ensureblockend(word)
        # denominator
        word = next(words)
        checkKey("Denominator{",word)
        word = next(words)
        checkKey("block",word)
        while word == "block":
            pnum = words.next_int()
            bnum = words.next_int()
            dblocks.append( (pnum,bnum) )
            word = next(words)
        ensureblockend(word)
        # value
        word = next(words)
        initval, lower, upper = parseConstraintValue(words,word)
        varobj = BlockRatioVariable(nblocks, dblocks, initval, lower, upper)
    elif key == "Chi{":
        # chi
        word = next(words)
        checkKey("monomers", word)
        mon1 = words.next_int()
        mon2 = words.next_int()
        word = next(words)
        initval, lower, upper = parseConstraintValue(words,word)
        varobj = ChiVariable(mon1, mon2, initval, lower, upper)
    elif key == "KuhnLength{":
        # Kuhn length
        word = next(words)
        checkKey("monomer", word)
        mon1 = words.next_int()
        word = next(words)
        initval, lower, upper = parseConstraintValue(words,word)
        varobj = KuhnVariable(mon1, initval, lower, upper)
    elif key == "KuhnRatio{":
        # Kuhn Ratio
        word = next(words)
        checkKey("monomers", word)
        mon1 = words.next_int()
        mon2 = words.next_int()
        word = next(words)
        initval, lower, upper = parseConstraintValue(words,word)
        varobj = KuhnRatioVariable(mon1, mon2, initval, lower, upper)
    else:
        raise(ValueError("Invalid key for search variable declaration, {}.".format(key)))
    return varobj

## Phase list definition methods

def parsePhases(words,key):
    checkKey("Phases{",key)
    key = next(words)
    checkKey("scft_solver",key)
    key = next(words)
    checkKey("pscf",key)
    key = next(words)
    if key == "input_root":
        root = pathlib.Path.cwd() / next(words)
    else:
        root = pathlib.Path.cwd()
    hasPhase = True
    tgts = []
    cmps = []
    while hasPhase:
        key = next(words)
        if key == "TargetPhase{":
            tgts.append(parsePhaseDetails(words,root))
        elif key == "CompetingPhase{":
            cmps.append(parsePhaseDetails(words,root))
        elif isblockend(key):
            hasPhase = False
        else:
            msg = "Unexpected key in Phases block, {}"
            raise(ValueError(msg.format(key)))
    return MesophaseManager(cmps, tgts)

def parsePhaseDetails(words,root):
    key = next(words)
    checkKey("name",key)
    name = next(words)
    key = next(words)
    checkKey("model_file",key)
    model = root / next(words)
    key = next(words)
    disfile = None
    distol = 0.001
    if key == "range_check":
        disfile = next(words)
        distol = words.next_float()
        key = next(words)
    coreOpts = None
    if key == "core_option":
        coreOpts = []
        while key == "core_option":
            coreOpts.append(words.next_int())
            key = next(words)
    ensureblockend(key) # At this point, key should hold block terminator "}"
    return PscfMesophase.fromFieldGenFile(name, model, coreOpts, disfile, distol)


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
        
        # Search Space definition
        word = next(words)
        expected = "SearchSpace{"
        checkKey(expected,word)
        variables, constants, velcap = parseSearchSpace(words,word)
        varset = PolymerVariableSet(variables,constants)
        
        # Create Mesophases and Phase Managers
        word = next(words)
        pman = parsePhases(words,word)
        
        # Check for start of Pso block
        key = next(words)
        checkKey("Pso{",key)
        
        # Check for specified random seed
        key = next(words)
        if key == "random_seed":
            seed = words.next_int()
            key = next(words)
        else:
            seed = time.time_ns() % (2**32)
        print("Random Seed Used: {}".format(seed))
        randGen = np.random.RandomState(seed)
        
        # Check for integrator specification
        if key == "StandardIntegrator{":
            # Read integrator specification
            key = next(words)
            if key == "constriction":
                chi = abs(words.next_float())
                key = next(words)
            else:
                chi = 0.729
            print("Using Constriction Factor: {}".format(chi))
            if key == "self_weight":
                eps1 = abs(words.next_float())
                key = next(words)
            else:
                eps1 = 2.05
            print("Using Self-Weight: {}".format(eps1))
            if key == "neighbor_weight":
                eps2 = abs(words.next_float())
                key = next(words)
            else:
                eps2 = 2.05
            print("Using Neighbor-Weight: {}".format(eps2))
            ensureblockend(key)
            integrator = StandardIntegrator(randGen, chi=chi, c1=eps1, c2=eps2
            key = next(words)
        else:
            integrator = StandardIntegrator(randGen)
        
        # Read Swarm Block
        checkKey("Swarm{",key)
        key = next(words)
        checkKey("n_agent",key)
        nagent = words.next_int()
        if nagent <= 0:
            raise(ValueError("Must have a positive number of agents; gave {}".format(nagent)))
        print("Using {} agents.".format(nagent))
        closeblock(words)
        key = next(words)
        checkKey("n_step",key)
        nstep = words.next_int()
        if nstep <= 0:
            raise(ValueError("Must have a positive number of steps; gave {}".format(nstep)))
        print("Running {} steps.".format(nstep))
        
        velcap = np.array(velcap)
        velval = np.zeros_like(velcap)
        velsrc = Velocity(velval,velcap)
        pwd = pathlib.Path.cwd()
        # Create Agents
        agentList = ScftAgent.createSeveral(nagent,varset,pman,velsrc,client,pwd)
        agentList = ScftAgent.randomizeSeveral(agentList,randGen)
        # Create Swarm
        graph = nx.cycle_graph(nagent)
        sm = Swarm(graph, agentList, integrator, pwd)
        
        ## Run iterations
        tryCooldown = 0.1
        runtime, proctime = sm.run(nstep,tryCooldown)
        print("Finished Running:\trun_time {}s;\tprocess_time {}s".format(runtime,proctime))

client.close()
sys.stdout.flush()
