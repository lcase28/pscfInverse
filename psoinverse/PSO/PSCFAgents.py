
# third party imports
import numpy as np

# Package imports

# first two imports are relocations of classes
#   and should give the same behavior as in original writing
#from . import Swarm.Agent as Agent
#from . import SearchSpace.SimulationPoint as Point
from Swarm import Agent
from SearchSpace import SimulationPoint as Point

print("Successful imports")


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

        except SQLAlchemyError as e:
            output("SQLAlchemy error: {}".format(e))
            log_exception()
            self.Location.Fitness = -2E4
            return False
        except Exception as e:
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
            except Exception as e:
                print("simid = {}".format(sim.id))
                print("Failed to get operators.dat (for H)", sim)
                print("Exception: ", e)
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
        except Exception as e:
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
        except Exception as e:
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
            

