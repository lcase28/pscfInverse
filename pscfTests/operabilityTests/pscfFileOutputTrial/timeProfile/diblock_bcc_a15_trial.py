from context import psoinverse

from psoinverse.PSO.Swarm import Swarm
from psoinverse.PSO.Integrators import StandardIntegrator

from psoinverse.mesophases.mesophaseVariables import BlockFractionVariable, \
                                                    ChiVariable, \
                                                    VariableSet
from psoinverse.mesophases.phaseManagement import MesophaseManager

from psoinverse.SCFT.PSCF.pscfMesophase import PSCFMesophase
from psoinverse.SCFT.scftAgent import ScftAgent

import networkx as nx
import numpy as np
import pathlib

bf = BlockFractionVariable(polyNum = 0, blockNum = 0, \
                        val = 0.25, lower = 0.0, upper = 0.5)
chi = ChiVariable(Monomer1 = 1, Monomer2 = 2, val = 2, \
                        lower = 0.0, upper = 4)
vs = VariableSet([bf,chi])

cwd = pathlib.Path.cwd()
bccrt = cwd.parent / "Inputs" / "BCC"
bccObj = PSCFMesophase.instanceFromFiles("BCC",bccrt/"rho_kgrid", \
                        bccrt/"param",bccrt/"model_in.txt")
a15rt = cwd.parent/"Inputs"/"A15"
a15Obj = PSCFMesophase.instanceFromFiles("A15",a15rt/"rho_kgrid", \
                        a15rt/"param",a15rt/"model_in.txt")
tgt = PSCFMesophase.instanceFromFiles("BCC2",bccrt/"rho_kgrid", \
                        bccrt/"param",bccrt/"model_in.txt")

candidates = {bccObj.phaseName : bccObj, a15Obj.phaseName : a15Obj}
mngr = MesophaseManager(candidates, tgt, vs)

rgen = np.random.RandomState(100)
agents = []
for i in range(4):
    rtpath = cwd/"agent{}".format(i)
    agents.append(ScftAgent(mngr, rgen, root=rtpath))

integr = StandardIntegrator(rgen, seekMax = True)
gr = nx.cycle_graph(4)
sm = Swarm(gr, agents, integr)
for i in range(2):
    sm.step()
