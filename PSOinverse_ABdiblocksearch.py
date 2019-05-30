#!/usr/bin/python

# PSO imports
from PSO.Swarm import *
from PSO.morphology_agent import MorphologyAgent
from PSO.mesophase_factory import MesophaseFactory
from PSO.pso_chains import PSOChainSection
from ParamFactories import PolyFTSFactory

# utility imports
import threading
import cPickle as pickle
import time
import networkx as nx
from collections import OrderedDict


# ============= Helpers =======================

def run_threads(threads):
    for t in threads: t.start()
    for t in threads: t.join()


# ================= FTS_DB initialization =================

# REPLACE ME - need a new PolyFTSFactory that writes files locally
# and doesn't include database stuff
# Remove launcher throughout
pfac = PolyFTSFactory()
pfac.ParseParameterFile('params.in')
pfac.set(wall_time=[0, 45])

# ================= Set up polymer system =================
# Set the chains

chain_section = PSOChainSection(pfac)

AlphaKeys = []
AlphaKeys += chain_section.add_chain([1, 2],1)  # min/max defaults to [-4, 4]
AlphaKeys += chain_section.finish_chains()

# Restrictions on alpha and phi:
AlphaKeys[0][1] = [-4, 4]

pfac.set(chains=chain_section)

# Add the non-chain degrees of freedom you want to search with PSO here
d = OrderedDict()
#for k, v in AlphaKeys:
for k, v in [["InteractionsBK_chiN12", [10,30]]] +AlphaKeys:
    d[k] = v


def add_agent_to_swarm():
    global agent_list

    a = MorphologyAgent("o70", ["gyr", "hex", "dis", "lam"], launcher=launcher, parameters=d,
                        group=group)
    agent_list.append(a)
    output("Added agent {}".format(a.id))


def update_agent(agent):
    agent.update()


# Add 20 agents
n_steps = 500

t = []

try:
    with open("state.pkl", 'r') as f:
        swarm = pickle.load(f)
    s0 = int(swarm.Agents[0].Location.simulations.values()[0].get("step"))
except IOError:
    # ================= Swarm initialization 
    # Start fresh
    group.Empty(False)

    agent_list = []

    rows = 2
    cols = 4
    n_agents = rows * cols

    # Also can use nx.cycle_graph for ring topology
#    von_neumann_graph = nx.convert_node_labels_to_integers(nx.grid_2d_graph(rows, cols, periodic=True))
    cycle_graph = nx.cycle_graph(rows * cols)

    output("Starting PSO with n_agents={}, n_steps={}, agents/DOF = {:.2f}".format(n_agents, n_steps,
                                                                                   n_agents * 1. / len(d.items())))

    for key, Range in AlphaKeys:
        output("{}: {}".format(key, Range))

    # Generate the agents
    for i in range(n_agents):
        t.append(threading.Thread(target=add_agent_to_swarm))
    run_threads(t)

    swarm = Swarm(cycle_graph, agent_list, chi=0.80)
#    swarm.print_neighbor_graph()

    output("Finished adding all agents. Starting PSO update sweeps")

    s0 = 0

    # 
# Now Run
for step in range(s0, n_steps):
    print "\n=======================Step {}===================\n".format(step)
    # Print output
    swarm.write_output()
    s.commit()

    start = time.time()
    t = []

    # Update all the agents
    for agent in swarm.Agents:
        output("Agent {}: {}={}".format(agent.id, agent.get_coords(), agent.Location.Fitness))
        t.append(threading.Thread(target=update_agent, args=(agent,)))

    run_threads(t)

    output("\n Step {} Finished in {} hours".format(step, (time.time() - start) / 3600))
    output("\t Best = Agent: {}".format(swarm.Best[-1]))

#    sim = swarm.Best[-1].simulations[target]
#    XN = sim.get('InteractionsBK_chiN12')
#    alpha = sim.get('chain1Alpha0')
#    from Sean.Constraints import NSegmentDistribute as nsd
#    f_A = nsd.AlphaToParams(alpha)[0]
#
#    output("\t Best = Agent: f_A={},XN={},alpha={}; fitness ={}".format(f_A,XN,alpha,swarm.Best[-1].Fitness))
