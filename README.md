# Python PSO Tools
PSO tools for general global optimization problems with fixed dimensionality and continuous degrees of freedom (for problems with dynamic variable counts and/or discrete degrees of freedom, choose something else like a GA).  Specific derived utilities are targetted to the inverse problem of block polymer self assembly.

Most of the core software is located in the `PSO/` directory, with various application examples present in the root.

## Brief description of core source files:
- `PSO/Swarm.py` contains three main **base** classes:
    - `Point` stores information about the coordinate of a point in the search space.
    - `Agent` (abstract) defines the interface for handling evaluation and storage of the fitness of an agent.
    - `Swarm` - a collection of agents with an associated communication topology.
- `PSO/Swarm.py` **derived** classes:
    - `FunctionAgent` derives from Agent and permits the simplest evaluation method on an arbitrary function pointer.
    - `SCFTAgent` computes the fitness as the SCFT free energy. **Deprecated**. This class is intrinsically tied to the SimulationManager and is being replaced. This type of agent may also not be very useful - it computes the SCFT free energy of a list of simulations / phases, but then only uses the free energy of the first one as the objective function for the PSO search. Minimizing the absolute value of the free energy of a single phase is not a useful exercise.
    - `LocalSCFTAgent`: reimplementation of the above without using the SimulationManager.
    - `LocalSCFTFieldAgent`: This class expects field coefficients to be the PSO variables, and tries to minimize free energy to find the globally stable structure. This likely doesn't work very well due to the enormous dimensionality of the problem.
- `PSO/SwarmPlotter.py` uses matplotlib to help with 1D, 2D and 3D plots of the PSO agents' configurations within a potentially higher-dimensional search space. Includes movie generation for time dependence. 
- `PSO/morphology_agent.py` contains a subclass of `LocalSCFTAgent`. Each agent's fitness is evaluated by computing the SCFT free energy of a list of candidate phases and computing the free energy difference between the target phase and the next most stable. This class is central to the method used in [Mihir's paper](https://doi.org/10.1021/acs.macromol.7b01204). Uses `MesophaseFactory` and `Candidate` from the `PSO/mesophase_factory.py` listed below.
- `PSO/mesophase_factory.py`: handle the computation of the SCFT free energy of a single mesophase using `PolyFTS`, and specify the phase-specific input parameters. See also **`unittests/testMesophaseFactory.py`**, which tests this functionality in a non-PSO example of computing free-energy curves of all phases along a trajectory in the phase diagram.
- `PSO/mihirtools.py`: probably nothing useful left in here.
- `PSO/Constraints.py` **and** `PSO/pso_chains.py` handle the search over variables with constraints. This includes the block fractions on a chain summing to one, and volume fractions of chains summing to one. **These classes have not been re-integrated into the cleaned out code as of yet, and they need to be**. The test cases conducted so far are on diblock melts where only one block fraction is varied (the second is implicit as one-first). For anything more complicated (multi-block chains, or multi-chain assemblies) we must restore the constraint functionality. The method of enforcing constraints is described in the Supporting Info of [Sean's paper](https://doi.org/10.1021/acsmacrolett.6b00494)
- `ParamFactories.py` is a helper module for producing input files for a simulation code from a template by text substitution.
- `SimulationLauncher.py` combines a `ParameterFactory` with a simulation code specification and runs the simulation (potentially a sweep of simulations). This replaces the SimulationManager with a minimal local run agent that executes synchronously. See also **`unittests/testParamFactory.py`**.

## Description of example runs:
- `testFunctionAgent_Griewank.py` minimizes the Griewank function in arbitrary dimensions.
- `testFunctionAgent_Rosenbrock.py` minimizes the Rosenbrock function.
- `testPSO_LocalSCFTAgent_SinglePhase_LAMDstar.py` uses PSO to optimize the lattice parameter of the lamellar phase of a diblock polymer melt (obviously using a variable cell-shape SCFT algorithm is way more efficient, but this provides an easy 1D test case of spawning SCFT jobs).
- `testPSO_LocalSCFTAgent_SinglePhase_HEXortho.py` as above, but for the orthorhombic setting of the 2D hexagonally packed cylinder phase.
- `testPSO_LocalSCFTAgent_SinglePhase_FdddCellRelax.py` as above for 3D cell relaxation.
- `testPSO_MorhphologyAgent_diblockTargetGYR.py` a full run that aims to reproduce Fig. 2 of [Mihir's paper](https://doi.org/10.1021/acs.macromol.7b01204).


## Suggestions for what to do next:
- Add support for other SCFT packages (particularly PSCF). This involves writing a new implementation of `MorphologyAgent` (at a minimum).
- Check that the PSO fundamentals are working correctly (e.g., one Griewank and Rosenbrock functions), especially velcity values and position updates.
- Validate that the code works in its present state by comparing to [Mihir's paper](https://doi.org/10.1021/acs.macromol.7b01204).
- Simplify the code:
    - Point should be incorporated into Agent.
    - Influencer weights and constriction factor should reside only in Swarm, not Agents
    - Remove residual stuff from the Simulation Manager (e.g., FTS.Globals, locking, SQAlchemy stuff). This includes deleting SCFTAgent entirely when its functionality has been replicated and validated in LocalSCFTAgent.
    - Move LocalSCFTAgent out of Swarm.py file. Keep Swarm.py completely general.
    - Remove neighbor list (and associated connect/disconnect methods) from the Agent and relocate to the Swarm.
    - Determine whether domain scaling is required or useful, and potentially remove it from the coordinates and velocities.
- Reintegrate the constraints methods.

