import numpy as np
from functools import total_ordering
import time
import traceback
import os
import shutil


# Helper function for debugging - just wraps the process of spitting out a string to a file
def debug(line):
    #with open("debug.out", 'a') as f:
    #    f.write("{}\n".format(line))
    print("DEBUG: {}".format(line))


def output(line):
    from datetime import datetime
    with open("runtime.out", 'a') as f:
        f.write("{} :: {}\n".format(datetime.now().strftime("%d/%m/%y %H:%M:%S"), line))

def log_exception():
    debug("Exception raised: {}".format(traceback.format_exc()))


# =========================================================
#
# Point Class - represents a point in phase space (block fraction, N_BCP, ChiN, ...)
#
# =========================================================
@total_ordering
class Point(object):
    def __init__(self):
        self.Fitness = None # This point's fitness when @ Coords
        self.Coords = None # Parameter values (current)
        self.Scale = None # Scale factor to apply to the coordinates
        self.simulations = {} # A record of all of the simulations used to determine fitness
        self.keys = None # Parameter names

    def fill_from(self, other_point):
        self.Coords = np.copy(other_point.Coords)
        self.Scale = other_point.Scale
        self.Fitness = other_point.Fitness
        self.simulations = {}
        self.simulations.update(other_point.simulations) # Update is a native dictionary operation
        self.keys = other_point.keys

#    def delete_simulations(self): # FTS_DB relete simulation records and resources
#        for k, simulation in self.simulations.items():
#            if simulation is not None:
#                FTS.DeleteSimulation(simulation)
#        self.clear_simulations()
#
#    def move_simulations_to(self, new_group): # FTS_DB - relocate the simulation records to a new group
#        for k, sim in self.simulations.items():
#            sim.Group = new_group
#        self.clear_simulations()

    def copy_simulations_to(self, newpath): # Copy local run files to a new path
        for k, sim in self.simulations.items():
            # Form destination and clean up existing
            newsim = newpath+k
            if os.path.isfile(newsim):
                os.remove(newsim)
            elif os.path.isdir(newsim):
                shutil.rmtree(newsim) # Remove any existing contents at newpath
            # This simulation may not have ran if not needed for computing fitness.
            # In that case skip the copy.
            if sim == None:
                continue
            if not os.path.isdir(sim):
                continue
            # Copy
            shutil.copytree(sim, newsim)

    def clear_simulations(self):
        for key in self.simulations.keys():
            self.simulations[key] = None

    def init_simulations(self, keys):
        for key in keys:
            self.simulations[key] = None

    def __gt__(self, other):
        if self.Fitness > other.Fitness:
            return True
        return False

    def __lt__(self, other):
        if self.Fitness < other.Fitness:
            return True
        return False

    def __eq__(self, other):
        if self.Fitness == other.Fitness and np.allclose(other.get_scaled_coords(), self.get_scaled_coords()):
            return True
        return False

    def get_dict(self):
        import OrderedDict
        return OrderedDict(zip(self.keys, self.get_scaled_coords()))

    def get_scaled_coords(self):
        return self.Coords * self.Scale

    # Allow friendly output when we print a Point() object
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        try:
            return "".join(["{}: {:.4f}\n".format(self.keys[i], self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)
        except ValueError:
            return "".join(["{}: {}\n".format(self.keys[i], self.get_scaled_coords()[i]) for i in
                            range(len(self.Coords))]) + "Fitness={:.4f}".format(self.Fitness)

#    # Support for pickling the Simulation record
#    def __getstate__(self):
#        try:
#            self.simulation_ids = {k: v.id for k, v in self.simulations.items()}
#        except:
#            pass
#
#        return self.__dict__
#
#    def __setstate__(self, state):
#        from table_def import Simulation
#
#        s = FTS.Globals['session']
#        for k, sim_id in state['simulation_ids'].items():
#            try:
#                state['simulations'][k] = s.query(Simulation).filter_by(id=state['simulation_ids'][k]).one()
#            except:
#                print "Failed to load sim_id {} unpickling Point".format(state['simulation_ids'][k])
#                state['simulations'][k] = None
#
#        self.__dict__.update(state)


# =========================================================
#
# Swarm Class
#
# =========================================================
class Swarm(object):
    def __init__(self, graph, agents, chi=None, c1=None, c2=None):
        assert len(graph) == len(agents), "Number of nodes in graph doesn't match number of agents"
        self.Agents = agents

        # For output
        self.History, self.Best = [], []

        # PSO Parameters - update from defaults if supplied
        if chi != None:
            Agent.chi = chi
            for a in self.Agents:
                a.chi = chi
        if c1 != None:
            Agent.c1 = c1
            for a in self.Agents:
                a.c1 = c1
        if c2 != None:
            Agent.c2 = c2
            for a in self.Agents:
                a.c2 = c2

        # Fill each agent's neighbor list according to the supplied graph
        import networkx as nx
        [self.Agents[node].connect(self.Agents[b]) for node in graph for b in graph.neighbors(node)]

    # Return global best of whole swarm
    def get_gbest(self):
        gbest_pt = self.Agents[0].PBest
        best_agent = self.Agents[0]

        for a in self.Agents:
            if a.PBest > gbest_pt:
                gbest_pt = a.PBest
                best_agent = a

        return best_agent

    def print_neighbor_graph(self, fname="SwarmNetwork.png"):
        import networkx as nx
        import pylab as py

        g = nx.Graph()
        [g.add_node(a.id) for a in self.Agents]
        [[g.add_edge(a.id, b.id) for b in a.neighbors] for a in self.Agents]

        # draw_graphviz doesn't work in networkx 1.11. Replace with 3 lines below.
        #nx.draw_graphviz(g)
        from networkx.drawing.nx_agraph import graphviz_layout
        pos = graphviz_layout(g)
        nx.draw(g, pos)

        py.savefig(fname)

    def write_output(self, outputbasedir):
        from copy import deepcopy

        # Maintain a history of the entire set of coordinates and fitnesses for all agents
        #self.History.append([list(agent.get_coords()) + [agent.Location.Fitness] for agent in self.Agents])

        # Maintain a list of the best agent at each step.
        # First, find the current best in the whole swarm.
        gbest_agent = self.get_gbest()
        # Only add to the list of Best coordinates if the gbest has been updated
        if gbest_agent.Location == gbest_agent.PBest or len(self.Best)<1:
            self.Best.append(deepcopy(gbest_agent.PBest))
            # Copy the simulation files to a gbest location
            gbest_agent.PBest.copy_simulations_to(outputbasedir+"GBest_step{}/".format(gbest_agent.steps))

        # Output all agent data
        for a in self.Agents:
            filename=outputbasedir+"History_Agent{}.dat".format(a.id)
            if a.steps <= 1:
                f=open(filename,'w')
                f.write("# ncoords = {}\n".format(len(a.Location.keys)))
                f.write("# step, fitness, {}, {}\n".format(", ".join(s for s in a.Location.keys), ", ".join("V_{}".format(s) for s in a.Location.keys)))
            else:
                f=open(filename,'a')
            f.write("{} {} {} {}\n".format(
                                    a.steps, a.Location.Fitness,
                                    " ".join(repr(s) for s in a.get_coords()),
                                    " ".join(repr(s) for s in a.Velocity)
                                    ))
            f.close()

        filename=outputbasedir+"History_GBest.dat"
        if gbest_agent.steps <= 1:
            f=open(filename,'w')
            f.write("# ncoords = {}\n".format(len(a.Location.keys)))
            f.write("# step, fitness, {}, AgentID\n".format(", ".join(s for s in a.Location.keys)))
        else:
            f=open(filename,'a')
        f.write("{} {} {} {}\n".format(
                                gbest_agent.steps, gbest_agent.PBest.Fitness,
                                " ".join(repr(s) for s in gbest_agent.PBest.get_scaled_coords()),
                                gbest_agent.id
                                ))
        f.close()

        #with open("Output.pkl", 'w') as f:
        #    pickle.dump(dict(history=self.History, best=self.Best), f)


# =========================================================
#
# Agent Class (Abstract)
#
# =========================================================
class Agent(object):

    # STATIC MEMBERS
    NextID = 0  # Static counter for generating unique agent IDs

    # These can be modified from their defaults in the Swarm object constructor, but
    # they should be shared by all agents and are hence static
    #
    # Influencer weights
    c1 = 2.05
    c2 = 2.05
    # Constriction factor
    chi = 0.729

    def __init__(self, boundaries, v0=None):
        """
            v0 is the initial velocity scale; all velocities will be randomized [-v0, v0]
              If v0 is None, the initial velocity will be set to 1/10 of the boundary range in each direction
            boundaries are the bounds of the search domain and should be a numpy array with shape (dims, 2).
        """
        # Initialize to zero steps taken
        self.steps = 0

        # Grab a unique ID for this agent
        self.id = Agent.NextID
        Agent.NextID += 1 # Update the static counter to deliver unique IDs

        #debug("Agent class init, ID = {}".format(self.id))

        # Store the seach boundaries and generate scale
        self.boundaries = boundaries
        scale = np.array([b[1] - b[0] for b in boundaries])

        #print "AGENT {}, boundaries = {}".format(self.id, scale)

        # The coords used in the update are scaled by the range of the search domain
        # Generate a Point object to store and manipulate the coordinate of the agent in phase space
        dims = len(boundaries)
        self.Location = Point()
        self.Location.Coords = (np.random.rand(dims) * (boundaries[:, 1] - boundaries[:, 0]) + boundaries[:, 0]) / scale
        self.Location.Scale = scale
        #print "AGENT {}, initial coordinates = {}, scaled coords = {} ".format(self.id, self.Location.Coords, self.Location.get_scaled_coords())
        self.Location.Fitness = -1e8

        # Randomize velocity uniformly on [-v0, v0]. Velocity is added to scaled coordinates.
        if v0 is None:
            self.Velocity = 2. * (np.random.rand(dims) - 0.5)
            self.Velocity *= 0.1 * np.array([b[1] - b[0] for b in boundaries]) / scale
        else:
            self.Velocity = 2 * v0 * (np.random.rand(dims) - 0.5) / scale

        #print "AGENT {}, initial velocity = {} ".format(self.id, self.Velocity)

        # Store PBest - best Point() location discovered in all history
        self.PBest = Point()
        self.PBest.fill_from(self.Location)
        #assert self.PBest.Fitness <  -1e7, "pbest fitness is too positive."

        # Initialize empty neighbor list and set fitness to zero. Will be maximized positive.
        self.neighbors = []
        #self.Location.Fitness = 0  # Do not do this - in some cases the fitness metric can be any sign and we don't want to pin the initial GBest/PBest to zero

        #debug("Agent {}: coords   = {}".format(self.id, self.get_coords()))
        #debug("Agent {}: velocity = {}".format(self.id, self.Velocity))
        #debug("Agent {}: fitness  = {}".format(self.id, self.Location.Fitness))


    def get_coords(self):
        return self.Location.get_scaled_coords()


    def get_nbest(self):
        # Find the best among the neighbors (including self)
        best = self.PBest
        for n in self.neighbors:
            if n.PBest > best:
                best = n.PBest
        return best


    def update(self, acceleration=None):
        """ Update the agent's Position (and Velocity) based on the current Position and neighbor information.

            ::Returns:: the agent's personal best fitness value (not Position)
        """

        # Find the best among the neighbors
        nbest = self.get_nbest()

        # Random variables for forcing terms
        e1 = np.random.rand(len(self.Location.Coords))
        e2 = np.random.rand(len(self.Location.Coords))

        # Inertia
        if acceleration is not None:
            acc = acceleration
        else:
            acc = np.zeros_like(e1)

        # Update the velocity
        new_velocity = (self.Velocity
                     + Agent.c1 * e1 * (self.PBest.Coords - self.Location.Coords)
                     + Agent.c2 * e2 * (nbest.Coords - self.Location.Coords))
        new_velocity = Agent.chi * new_velocity + acc

        # Update the position according to PSO dynamics. Retain in tmp array for boundary checks / constraints
        attempt = Point()
        attempt.fill_from(self.Location)
        attempt.Coords += new_velocity

        # Reflective boundaries
        if self.boundaries is not None:
            scaled_attempt = attempt.get_scaled_coords()
            for i, out_of_bounds in enumerate(
                    np.logical_or(scaled_attempt < self.boundaries[:, 0], scaled_attempt > self.boundaries[:, 1])):
                if out_of_bounds:
                    #print "REFLECT AGENT {}, V COMPONENT {}".format(self.id,i)
                    #print "OLD POSITION = {}".format(self.Location.get_scaled_coords())
                    #print "TRIAL POSITION = {}".format(scaled_attempt)
                    #print
                    #print
                    new_velocity[i] *= -1.0
                    attempt.Coords[i] = self.Location.Coords[i]

        # Replace internal Location and Velocity with new values
        move = True
        if move:
            from copy import deepcopy
            vtemp = deepcopy(self.Velocity)
            ltemp = deepcopy(self.Location.Coords)

            self.Velocity = new_velocity
            #output("Agent {}: V = {}".format(self.id, self.Velocity))
            self.Location.Coords = attempt.Coords

            # Fitness is incremented during Evaluate() since we use super() calls to generate the MRO
            self.Location.Fitness = 0
            # Use self.Evaluate as a validator
            if not self.evaluate():
                debug("Agent {} Resetting!".format(self.id))
                self.Velocity = vtemp
                self.Location.Coords = ltemp
                return False

        return True

    # Disconnect self and a neighbor 'a' from each other
    def disconnect(self, a):
        self.neighbors.remove(a)
        a.neighbors.remove(self)

    # Connect self and a neighbor 'a' to each other
    def connect(self, a):
        if a not in self.neighbors:
            self.neighbors.append(a)
        if self not in a.neighbors:
            a.neighbors.append(self)

    # Derived classes must override, update fitness, and call base using super calls to generate the MRO.
    def evaluate(self):
        # Now that fitness has been updated, compare to PBest
        if self.Location > self.PBest:
            # Re-use the old PBest simulation record
            #sim_tmp = {k: v for k, v in self.PBest.simulations.items()}

            self.PBest.fill_from(self.Location)

            # When we store the PBest point, it gets a reference to self.Location.sim
            # If sim_tmp = None (PBest never set), then a new record will be generated automatically in Agent.Evaluate()
            #self.Location.simulations = sim_tmp
        return True


# =========================================================
#
# FunctionAgent Class (Derives From: Agent)
#
# =========================================================
class FunctionAgent(Agent):
    def __init__(self, fn, **kwargs):
        self.Function = fn

        super(FunctionAgent, self).__init__(**kwargs)

        self.PBest.Fitness = -1

        self.evaluate()

    def evaluate(self):
        self.Location.Fitness = self.Function(self.get_coords())

        return super(FunctionAgent, self).evaluate()


# =========================================================
#
# SCFTAgent Class (Derives From: Agent)
#
# Send an SCFT calculation with the given Coords as parameter values
#
# =========================================================
class SCFTAgent(Agent):
    def __init__(self, simulation_keys, launcher, parameters, group, seed=None, **kwargs):

        self.Params = []
        boundaries = []
        dx = []

        self.seed = seed

        for k, v in parameters.items():
            self.Params.append(k)

            boundaries.append(v)
            dx.append(v[1] - v[0])

        dx = np.array(dx)
        boundaries = np.array(boundaries)

        # Initialize the Agent so we get its starting Coords
        super(SCFTAgent, self).__init__(v0=dx * .1, boundaries=boundaries, **kwargs)

        self.Location.init_simulations(simulation_keys)
        self.Location.keys = [k for k in self.Params]

        self.PBest.init_simulations(simulation_keys)
        self.PBest.keys = [k for k in self.Params]

        # pt = {k: v for k, v in zip(self.Params, self.Location.Coords)}

        self.launcher = launcher

        with FTS.Globals['lock']:
            self.group = group.Group("Agent_{}".format(self.id))

        self.evaluate()

    def RunSCFT_AutoDT(self, point, simulation_key, dt_lbound=1, dt_ubound=3, **kwargs):
        numsteps_block = self.launcher.pfac.get("SimulationBK_NumTimeStepsPerBlock")
        numblocks = self.launcher.pfac.get("SimulationBK_NumBlocks")

        # Run for 1 time step, scrape the eigenvalues, then run for the full period
        self.RunSCFT(point, simulation_key, SimulationBK_NumTimeStepsPerBlock=1, SimulationBK_NumBlocks=1, **kwargs)
        sim = point.simulations[simulation_key]

        # Get run.out
        with FTS.Globals['lock']:
            FTS.sshOpen(sim.cluster)
            nspecies = int(sim.get("monomersBK_nSpecies"))
            sed_cmd = ";".join(["n;p"] * nspecies)
            res, ign, ign = FTS.sshCMD("sed -n '/Eigen/{{{}}}' {}/run.out".format(sed_cmd, sim.getFullPath()))

        # Parse it for eigenvalues
        eigs = np.array([np.abs(float(l.split("->")[0])) for l in res])
        # Set the mobilities such that lambda / eigenvalue = 1
        DT = eigs
        # Band pass the time steps
        DT[eigs < dt_lbound] = dt_lbound
        DT[eigs > dt_ubound] = dt_ubound

        self.RunSCFT(point, simulation_key, SimulationBK_lambdaForceScale=DT,
                     SimulationBK_NumTimeStepsPerBlock=numsteps_block,
                     SimulationBK_NumBlocks=numblocks, **kwargs)

    def watch_scft(self, current_sim, **kwargs):
        start = time.time()
        while True:
            with FTS.Globals['lock']:
                # debug("{} Waiting for SCFT".format(self.id))
                pass

            time.sleep(2)
            with FTS.Globals['lock']:
                try:
                    FTS.UpdateSimulationStatus(current_sim)
                except:
                    log_exception()

                if current_sim.status in ["DONE", "TIMEDOUT"]:
                    # debug("{} Done".format(self.id))
                    return True

            # Check to see if this job has been queued for a while
            if time.time() - start > 60:

                with FTS.Globals['lock']:
                    pass
                    # debug("{} Checking sim".format(self.id))

                with FTS.Globals['lock']:
                    # This Freshen may not be necessary
                    FTS.Freshen(current_sim.experiment)

                    # If it failed, try reducing the time step
                    if current_sim.status == "FAILED":
                        dt = float(current_sim.get("SimulationBK_TimeStepDT"))
                        # Actually... don't keep re-trying failed sims. Probably in a pathological part of phase space
                        if dt < .5 or True:
                            return False

                        output("\t Simulation {} Failed, trying with dt={:.3f}".format(current_sim.id, dt / 2.0))
                        FTS.RerunSimulation(current_sim, SimulationBK_TimeStepDT=dt / 2.0, **kwargs)

                    # If the job isn't RUNNING
                    if current_sim.status not in ["RUNNING", "DONE", "TIMEDOUT"]:
                        if time.time() - start > 1200:
                            output(
                                "\t Taking a while... simulation {} [{}] on {}".format(current_sim.id, current_sim.status,
                                                                                       current_sim.cluster.name))
                        # Find the lowest load cluster
                        new_cluster = FTS.AutoSelectCluster(current_sim.code, cluster_list=self.launcher.cluster_list)
                        # If the current cluster isn't the lowest load, then move this sim to the lower load cluster
                        if new_cluster != current_sim.cluster or current_sim.status == None:
                            output(
                                "\t Rerunning simulation {} on {} to {}".format(current_sim, current_sim.cluster.name,
                                                                                new_cluster.name))
                            FTS.RerunSimulation(current_sim, cluster=new_cluster, **kwargs)
                    start = time.time()

    def RunSCFT(self, point, simulation_key, **kwargs):
        pt = {k: v for k, v in zip(point.keys, point.get_scaled_coords())}
        pt.update(kwargs)

        simulation_record = point.simulations[simulation_key]

        with FTS.Globals['lock']:
            # If the simulation record hasn't been created yet, do so
            if simulation_record is None:
                # debug("{} Launching new sim".format(self.id))
                try:
                    point.simulations[simulation_key] = self.launcher.launch(group=self.group, **pt)
                except:
                    log_exception()

            # If we have a record, just rerun using it - save on disk space
            else:
                #debug("{} Rerunning (for new) sim".format(self.id))
                try:
                    FTS.RerunSimulation(simulation_record, **pt)
                except:
                    debug("problem rerunning sim {}".format(simulation_record))
                    log_exception()
            simulation_record = point.simulations[simulation_key]
            simulation_record.set(step=self.steps)
            simulation_record.set(velocity=self.Velocity)
            simulation_record.set(coords=self.get_coords())

        if self.watch_scft(simulation_record, **pt):
            return True
        return False

    def evaluate(self):

        self.steps += 1

        from sqlalchemy.exc import SQLAlchemyError
        try:
            fitness = self.ComputeFitness()

            self.Location.Fitness += fitness

        except SQLAlchemyError, e:
            output("SQLAlchemy error: {}".format(e))
            log_exception()
            self.Location.Fitness = -2E4
            return False
        except Exception, e:
            output("\t Exception on agent {}: {}".format(self.id, e))
            log_exception()
            self.Location.Fitness = -2E4

            try:
                with FTS.Globals['lock']:
                    self.Location.delete_simulations()
            except:
                log_exception()

            self.Location.delete_simulations()

            return False

        return super(SCFTAgent, self).evaluate()

    def GetH(self, sim):
        FTS.WaitForFile(sim, 'operators.dat')

        with FTS.Globals['lock']:
            try:
                t, h = FTS.Load(sim, "operators.dat", cols=[0, 1])
            except Exception, e:
                print "simid = {}".format(sim.id)
                print "Failed to get operators.dat (for H)", sim
                print "Exception: ", e
                raise e

        return h[-1]

    # Default fitness metric is the hamiltonian. This is typically overloaded by subclassing SCFTAgent.
    def ComputeFitness(self):
        return -self.GetH(self.Location.simulations)

    # Enable friendly output upon print SCFTAgent() object
    def __repr__(self):
        return "ID: {}, pt: {}, v: {}".format(self.id, ", ".join(
            ["{}={:.2f}".format(k, v) for k, v in zip(self.Params, self.Location.Coords)]), ", ".join(
            ["{}={:.2f}".format(k, v) for k, v in zip(self.Params, self.Velocity)]))

# =========================================================
#
# LocalSCFTAgent Class (Derives From: Agent)
#   This is fully analogous to SCFTAgent, but all runs are
#   performed locally (without remote cluster job management)
#
#   This agent can accept multiple simulation keys, but it currently uses only the first's free energy
#   to determine the fitness.
#
# =========================================================
class LocalSCFTAgent(Agent):
    def __init__(self, simulation_keys, launcher, PSOparameters, SeedFields={}, **kwargs):
        # Initialize arrays
        self.Params = []
        boundaries = []
        v0 = []
        # Information for launching simulations
        # SeedFields is a dictionary of simulation_keys pointing to seed file paths.
        #  If key is not present, set to None
        self.SeedFields = SeedFields
        for k in simulation_keys:
            if k not in self.SeedFields:
                self.SeedFields[k] = None
        self.simulation_keys = simulation_keys

        # Store parameters involved in search
        # Passed as a dictionary; key = ParameterFactory keyword, value = two-entry list with bounds
        for k, v in PSOparameters.items():
            self.Params.append(k) # Parameter name
            boundaries.append(v)  # Parameter search bounds
            v0.append(0.1*(v[1] - v[0])) # Make initial velocity 1/10 of the search domain
        # Transform to numpy arrays
        v0 = np.array(v0)
        boundaries = np.array(boundaries)

        # Initialize the Agent so we get its starting Coords
        super(LocalSCFTAgent, self).__init__(v0=v0, boundaries=boundaries, **kwargs)

        self.agentrootdirectory="Agent_{}/".format(self.id)

        # The Point object will store a dictionary of simulation keys pointing to their paths
        self.Location.init_simulations(simulation_keys)
        # Set names of parameters in the Location Point() object
        self.Location.keys = [k for k in self.Params]

        self.PBest.init_simulations(simulation_keys)
        self.PBest.keys = [k for k in self.Params]

        # Store the launcher
        self.launcher = launcher

        # Evaluate my fitness
        self.evaluate()

    def RunSCFT_TuneMSEFieldMobilities(self, point, simulation_key, lbound=1, ubound=3, **kwargs):
        """ Run the SCFT job with tuning mobilities.
        """
        numsteps_block = self.launcher.pfac.get("SimulationBK_NumTimeStepsPerBlock")
        numblocks = self.launcher.pfac.get("SimulationBK_NumBlocks")

        # Run for 1 time step, scrape the eigenvalues, then run for the full period
        self.RunSCFT(point, simulation_key, SimulationBK_NumTimeStepsPerBlock=1, SimulationBK_NumBlocks=1, **kwargs)
        # Fetch simulation directory from the Location Point() object
        runpath = self.Location.simulations[simulation_key]

        # TODO: correctly obtain the MSE eigenvalues from the output
        #  e.g., if incompressible, one mode is eliminated. That mode gets mobility = 1 here, but there is potential for error in
        # capturing unintended output.
        nspecies = int(self.launcher.pfac.get("modelsBK_monomersBK_nSpecies"))
        evfound = False
        eigs = np.ones(nspecies)
        idx = 0
        with open(self.launcher.baserundir+runpath+"/STDOUT") as f:
            for line in f:
                if evfound and idx < nspecies-1:
                    if line.split()[1] != "->":
                        continue
                    eigs[idx] = line.split()[0]
                    idx += 1
                if "Eigenvalues" in line:
                    evfound = True
                if idx == nspecies:
                    break
        # Set the mobilities such that lambda / eigenvalue = 1
        mobility = eigs
        # Band pass the time steps
        mobility[eigs < lbound] = lbound
        mobility[eigs > ubound] = ubound

        # Now run the full SCFT with new mobilities
        self.RunSCFT(point, simulation_key, SimulationBK_lambdaForceScale=mobility,
                     SimulationBK_NumTimeStepsPerBlock=numsteps_block,
                     SimulationBK_NumBlocks=numblocks, **kwargs)


    def RunSCFT(self, point, simulation_key, **kwargs):
        # Create a dictionary with all settings from PSO search variables
        # in "point" and miscellaneous arguments in "kwargs"
        pt = {k: v for k, v in zip(point.keys, point.get_scaled_coords())}
        pt.update(kwargs)

        # Build location for running the SCFT calculation
        runpath = "step{}/".format(self.steps)+self.agentrootdirectory+"{}".format(simulation_key)
        point.simulations[simulation_key] = self.launcher.baserundir + runpath # Store for later retreival. Note that the launcher prepends a base directory to the runpath argument - we include it in the point record

        # Run the simulation
        try:
            if 'SeedField' in pt.keys(): # Override self.SeedFields[key] IF a SeedField has been passed in **pt.
                status = self.launcher.launch(runpath,**pt)
            else:
                status = self.launcher.launch(runpath,SeedField=self.SeedFields[simulation_key],**pt)

            if status == 'FAILED':
                return False
            # If fail or timeout here, could tweak dt etc.
            # Sean found that it's difficult to do this in an automated way - there is
            # probably some difficulty with the run due to entering an ill behaved region of the
            # parameter space.
        except:
            log_exception()

        return True


    def evaluate(self):
        self.steps += 1
        try:
            fitness = self.ComputeFitness()
            self.Location.Fitness += fitness
        except Exception, e:
            debug("\t Exception on agent {}: {}".format(self.id, e))
            log_exception()
            self.Location.Fitness = -2E4
            #self.Location.delete_simulations()
            return False

        #return super(LocalSCFTAgent, self).evaluate()

        # We update PBest here (instead of in Super's method)
        # and additionally keep a PBest copy of files and history
        if self.Location > self.PBest:
            oldPBfitness = self.PBest.Fitness
            newPBfitness = self.Location.Fitness
            self.PBest.fill_from(self.Location)
            newdir = self.launcher.baserundir+"PBests/Agent_{}/".format(self.id)
            self.PBest.copy_simulations_to(newdir)
            with open(newdir+"/history.out",'a') as f:
                f.write("PBest update on iteration {} : new fitness {}\n".format(self.steps, newPBfitness))


        return True


    # Default fitness metric is the hamiltonian. This is typically overloaded by subclassing SCFTAgent.
    def ComputeFitness(self):
        # Run all of the simulations that are used by this agent
        for simkey in self.simulation_keys:
            #self.RunSCFT_TuneMSEFieldMobilities(self.Location, simkey)
            self.RunSCFT(self.Location, simkey)
        # The fitness is just -H of the first simulation
        return -self.GetH(self.simulation_keys[0])


    def GetH(self, simulation_key=None):
        try:
            # Get just the simulation run path
            runpath = self.Location.simulations[simulation_key]
            # Load operators.dat
            opdata = np.loadtxt(runpath+"/operators.dat", usecols=[0,1])
        except Exception, e:
            debug("Agent ID = {}".format(self.id))
            debug("simulation_key = {}".format(simulation_key))
            debug("Exception: ", e)
            log_exception()
            raise e

        return opdata[-1][1] # Return the last H sample


    # Enable friendly output upon print LocalSCFTAgent() object
    def __repr__(self):
        return "ID: {}, pt: {}, v: {}".format(self.id, ", ".join(
            ["{}={:.2f}".format(k, v) for k, v in zip(self.Params, self.Location.Coords)]), ", ".join(
            ["{}={:.2f}".format(k, v) for k, v in zip(self.Params, self.Velocity)]))


# =========================================================
#
# LocalSCFTFieldAgent Class (Derives From: Agent)
#  This agent searches over field coefficients as the PSO variables.
#  WriteFields() transforms the PSO coords into an input field file.
#
# =========================================================
class LocalSCFTFieldAgent(Agent):
    def __init__(self, NFields, NPW, vmin, vmax):
        # Each mode has
        super(LocalSCFTFieldAgent, self).__init__(NFields * NPW, vmin)

        self.folder = "Agents/agent_{}".format(self.id)

        import os, shutil
        os.system("rm -rf {}".format(self.folder))
        os.makedirs(self.folder)

        # Copy in the params.in and initial field template
        shutil.copy("params.in", self.folder)
        shutil.copy("fields.dat", self.folder)
        shutil.copy("/home/sean/GIT/polyfts/bin/Release/PolyFTS_SP.x", self.folder)

    def evaluate(self):
        # I have a current field, so write it
        self.WriteFields()

        import os
        os.system("cd {};./PolyFTS_SP.x > run.out".format(self.folder))

        import PlotFTS as plt
        self.Location.Coords = plt.LoadDensity("{}/{}".format(self.folder, "fields.dat"), col=2)[1].flatten()

        self.Location.Fitness = -float(np.loadtxt("{}/operators.dat".format(self.folder), usecols=[1])[-1])
        super(LocalSCFTFieldAgent, self).evaluate()

    def WriteFields(self):
        data = np.loadtxt("{}/fields.dat".format(self.folder), usecols=range(3))

        import tempfile, shutil
        tmp = tempfile.NamedTemporaryFile()
        header = []
        header.append("# Format version 2")
        header.append("# nfields = 2")
        header.append("# NDim = 1")
        header.append("# PW grid = {}".format(len(self.Location.Coords)))
        header.append("# k-space data = 0 , complex fields = 0")
        header.append("# Columns: x field.real")

        with tmp:
            tmp.write("\n".join(header) + "\n")
            for field_val, line in zip(self.Location.Coords, data):
                line = list(line)
                line[2] = field_val
                tmp.write(" ".join(np.array(line, str)) + "\n")

            tmp.flush()

            shutil.copy(tmp.name, "{}/fields.in".format(self.folder))
