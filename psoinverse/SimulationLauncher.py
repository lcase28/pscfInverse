# Minimal launcher adapted from Sarean Paradiso's FTS_DB module by Kris Delaney

import os, sys
from subprocess import Popen, PIPE
from ParamFactories import *
from shutil import copyfile

STATUSDIC = {0: "RUNNING", 1: "FAILED", 2: "DONE", 3: "TIMEOUT"}
def getSimulationStatus(runpath):
    # Read status code and return error if not DONE
    statusfile = open(runpath+"/STATUS","r")
    status= STATUSDIC[int(statusfile.readline())]

    if status != "DONE":
        print "Run status : ",status

    return status


def LaunchSimulation(pfac, runpath, binary, SeedField=None, **passed_args):
    """
        pfac must be sent in ready to be printed (i.e. a 'primed' pfac)

        returns False if the run failed, and True if succeeded
    """

    print "Launching simulation in path ",runpath

    # Set this simulations's parameters in the pfac
    pfac.set(**passed_args)

    # Create the remote folders for this simulation (if required)
    if not os.path.exists(runpath):
        os.makedirs(runpath)
    cwd = os.getcwd()

    # Copy the seed field to the runpath
    if SeedField != None:
        copyfile(SeedField, runpath+"/"+os.path.basename(SeedField))

    # Change to run directory
    os.chdir(runpath)
    param_file = pfac.PrintCurrent('./params.in')

    # Run the simulation
    cmd = "{binary} {inputfile} > STDOUT".format(binary=binary, inputfile="./params.in")
    #print "Running command: ",cmd
    p = Popen(cmd, shell=True, stdout=None, stderr=PIPE)
    p.wait()

    errcode = p.returncode
    if p.returncode != 0:
        print "subprocess return code = ",p.returncode
        print "STDERR:"
        print p.stderr.readlines()

    status=getSimulationStatus("./")

    os.chdir(cwd)

    return status


# ============================================================================
# Class SimulationLauncher
# ============================================================================
class SimulationLauncher(object):
    """
    SimulationLauncher gets loaded with a parameterfactory, an executable,
    and a root directory for simulations. Each simulation launch specifies override parameters
    and a run subdirectory.
    """

    def __init__(self, pfac, code, baserundirectory):
        self.pfac = pfac
        self.code = code
        self.baserundir = baserundirectory

    def launch(self, rundirectory, code=None, SeedField=None, **params):
        """
            Launch a simulation.

            :param params: This kwargs list will be passed directly to the pfac record to set new parameter values.
            :type params: kwargs
        """
        # Set up the pfac for this run
        from copy import deepcopy
        localpfac = deepcopy(self.pfac)
        localpfac.set(**params)

        # If no code supplied, use default
        if code == None:
            code = self.code

        status = LaunchSimulation(localpfac, self.baserundir+rundirectory,code,SeedField)

        #return self.baserundir+rundirectory # Return the location of the run.
        return status
