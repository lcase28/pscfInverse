
******
Output
******

Several forms of output are produced by PscfInverse during a calculation.
All of this output is placed within the directory from which the program
was invoked. All of the output falls into one of three categories:
Runtime information, agent-specific data, and swarm-level data.

Runtime Information
===================

During program execution, some information is printed to standard output 
(*i.e.* printed to the console or piped to an output file).
The majority of this runtime information is produced while parsing
the parameter file to set up the problem. In particular, as the parameter
file is read, the program prints each "token" or "word" as it is read from
the file. If an error occurs during the setup of the problem, this output
can be used to identify when the error occurred and where, roughly, the
mistake may be in the parameter file (since the mistake is generally highlighted
by one of the last tokens read from the file). In addition to tracing the 
interpretation of the parameter file, values chosen for certain parameters
are also output. One key piece of data is the number used as the seed for 
the random number generator. If not specified, the program generates a seed
automatically based on the system clock, and this value is printed to 
standard output in a line of the form ``Random Seed Used: [seed]``, where
``[seed]`` would be replaced by the integer value of the random seed.

After the PSO calculation is complete, some time data is output.
This data is printed in a line similar to

::

    Finished Running:	run_time 3852.9761612415314s;	process_time 39.253919947s

In this line, the number following the label ``run_time`` indicates the wall (clock)
time that elapsed from the start to the end of the PSO search. 
The number following the label ``process_time`` is an estimate of the amount of time
spent actively performing operations in the PSO search, including checking for
completion of SCFT calculations, parsing SCFT results, and updating agent positions.
This ``process_time`` metric excludes time spent idling while waiting for SCFT calculations
to complete. In general, the core PSO operations are substantially lighter and faster than
the underlying SCFT calculations, which are the true bottle-neck in program operation.
Because of this, the process time for a PSO calculation is several orders of magnitude
smaller than the overall runtime.

Agent-Level Data
================

Each agent in the swarm outputs a significant amount of data throughout the calculation.
In order to organize this data, each agent generates their own dedicated directory within
the running directory. These directories are given the name ``agent[ID]`` where ``[ID]``
would be replaced with the agent's integer ID number. These ID numbers range in value
from 0 to *n_agent*.

Within each ``agent[ID]`` directory are three data files in Comma-Separated Value (CSV)
format, as well as a collection of further sub-directories with names following the
format ``step[##]`` where ``[##]`` would be replaced with an integer step number.
Here, we first detail the structure of the step sub-directories, and then describe the
CSV data files.

SCFT Calculation Output
-----------------------

Each ``step[##]`` sub-directory is used to keep each of the SCFT calculations performed
by an agent on a given step in one place, separate from other steps.
Numbering of the steps starts at 0, which are the SCFT calculations performed at the
Agent's initial position, before any PSO steps have been completed.
Subsequent step numbers indicate that the contained SCFT calculations were performed
at the Agent's position following the PSO step with the same number.
Thus, calculations within directory ``step1`` are done at the agent's position
following the first PSO position update (the first PSO step). Similarly is ``step2``
for the position following the second update during the second PSO step.
This continues out until ``step[n_step]`` which are calculations for the
agent's final position, after *n_step* PSO steps have been completed.
Altogether, each agent visits ( *n_step + 1* ) total positions.
Within each of these step directories, yet another sub-directory is generated for each
phase's SCFT calculation, which are named according to the ``name`` assigned to
each phase in the parameter file. These phase directories are the actual working
directories set for each of the SCFT calculations, and all output for each individual
SCFT calculation is stored in one of these phase sub-directories.

PSO Data File
-------------

Within each of the agent directories, the file ``psoData.csv`` contains
a tabulation of PSO-relevant data for the agent, tabulated in CSV format.
The file is summarized by the following generalized column headers:

::

    agent, step, bestStep, fitness, [ Variable_Positions ], [ Variable_Velocities ]

The first four headers are uniform across all calculations. 
Here, the ``agent`` column contains the ID number of the agent (the same
as that used in the agent's directory name).
The ``step`` column identifies the current step being detailed. These
numbers correspond directly to the numbers used in the ``step[##]`` directory
names, and indicate that the values detailed in that row are the current values
at the agent's position following the PSO step number given in this column.
The number in the ``bestStep`` column indicates the step number on which the
agent's current personal best position was seen. Rather than output duplicate data
for an agent's current and best position histories, the ``bestStep`` column is used
as an internally-referential tracker. For any step in the agent's history,
the personal best position seen by the agent up to that point in the search can be
found by taking the value in the ``bestStep`` column and tracing backward to the row
in which the same value is found in the ``step`` column. Thus, ``bestStep`` values
only change when the agent finds a new personal best position, and on the step when
that new best position is found, the values in ``step`` and ``bestStep`` will always
be equal. The value in the ``fitness`` column is the fitness of the agent at that step.

The remaining columns will depend on the search being performed. These columns specify
the position and velocity of the agent at each step according to the PSO search variables.
One column corresponding with each search variable is placed within the space designated
by ``[ Variable_Positions ]``. The headers placed here depend on the search variables
themselves, but all follow the general format of ``p_[label]`` where ``p_`` indicates
position data and ``[label]`` would be replaced by a unique label for the search variable.
Similar organization is used for velocity data, which is put in place of
``[ Variable_Velocities ]``, with column headers ``v_[label]`` where ``v_`` indicates 
velocity data and ``[label]`` is the same as in the corresponding position column.
These labels are algorithmically generated according to a few rules detailed below.
Depending on the detaild of the search, these labels can get very long, but will fully
specify the search variable.

 * The variable type is identified by the first part of the label.

    * Total Block Length : Length
    * Block Length Ratio : LengthRatio
    * Kuhn Length : Kunn
    * Kuhn Ratio : KuhnRatio
    * Chi parameter : ChiN
 
 * Specific polymer blocks are listed with ``_p[pid]b[bid]`` for block *[bid]* in polymer *[pid]*.
 * Monomers are identified with ``_m[mid]`` for monomer *[mid]*
 * For ratio variables, the components in the numerator are listed immediately following
   the variable type, followed by the key ``_r`` for ratio, followed by the components in the
   denominator.

An example of possible labels for each variable type are listed below:

 * Total Block Length:  ``Length_p0b0_p0b1`` for total length of a diblock copolymer.
 * Block Length Ratio:  ``LengthRatio_p0b0_r_p0b1`` for ratio between blocks of diblock polymer.
 * Kuhn Length : ``Kuhn_m0`` for statistical segment length of monomer 0.
 * Kuhn Ratio : ``KuhnRatio_m0_r_m1`` for ratio between segment lengths of monomers 0 and 1.
 * Chi :  ``ChiN_m0_m1`` for interaction parameter between monomers 0 and 1.

The Following is an example of the first four steps of psoData.csv for search in conformationally
symmetric diblock copolymers.

::

    agent,step,bestStep,fitness,p_LengthRatio_p0b0_r_p0b1,p_ChiN_m0_m1,v_LengthRatio_p0b0_r_p0b1,v_ChiN_m0_m1
    0,0,0,-inf,-1.5424027481055385,22.442175420796637,-0.06227226099288552,1.1414343348550768
    0,1,0,-inf,2.0836140548746216,18.86587413498185,3.0,-3.5763012858147856
    0,2,0,-inf,2.2586988939857373,28.28532777914046,0.17508483911111555,9.419453644158606
    0,3,0,-inf,0.1366481769710699,28.28532777914046,-2.1220507170146674,-6.689718874699529
    0,4,4,-0.14269907049999997,0.1366481769710699,20.989136121554353,2.871803070474223,-7.296191657586108

Phase Data File
---------------

Within the Agent directory, the file phaseData.csv contains a summary of polymer-specific
data includin phase energies and polymer parameters. The general format of the file, by
generalized column header outline, is given below.

::

    agent, step, [Target_Energies, ...], [Competitor_Energies, ...], [Polymer_Parameters, ...]

The first two columns match the columns of the same name from the PSO data files,
with agent ID numbers and step numbers. Each of the following three bracketed items
can vary in length. The first two ``[Target_Energies, ...]`` and ``[Competitor_Energies, ...]``
contain the Free energy, relative to disorder, of each target or competing phase.
The column labels are structured as ``F_tgt_[phaseName]`` for target phases
and ``F_cmp_[phaseName]`` for competitors, where ``[phaseName]`` would be replaced
with the name assigned to the phase in the parameter file.

The columns placed in the ``[Polymer_Parameter, ...]`` section at the end of the line
are the final values of the SCFT-relevant polymer parameters used in the calculation.
As with variable labels in the PSO data, these are algorithmically generated basaed on
the parameter type and identifying information. Presently, three types of Parameters
can appear here (as only three types of parameters are currently involed in variable
relationships). The first is block lengths for individual blocks which given a label
of the form ``Block_Len_p[pid]_b[bid]`` for block *[bid]* in polymer *[pid]* (both of
which are integer-indexed as they were specified in the SearchSpace section of the
parameter file). The second is the Statistical Segment Length of a monomer using the 
label ``Kuhn_m[mid]`` for monomer *[mid]*. Finally, interaction parameters labeled
with ``ChiN_m[mid_1]_m[mid_2]`` for the interaction between monomers [mid_1] and *[mid_2]*.
For both the segment length and the interaction parameters, the bracketed terms
would be the same integer-indexed monomer IDs used in declaration of the search space.

The following is an example of the first four steps for a search in conformationally
symmetric diblocks targeting BCC Spheres with cometitors Hexagonally-packed cylinders (hex),
double gyroid (gyr) and lamellae (lam).

::

    agent,step,F_tgt_bcc,F_cmp_gyr,F_cmp_hex,F_cmp_lam,Block_Len_p0_b0,Block_Len_p0_b1,ChiN_m0_m1
    0,0,inf,inf,inf,0.0,0.17618625733950488,0.8238137426604951,22.442175420796637
    0,1,inf,inf,inf,0.0,0.8893003207531411,0.11069967924685886,18.86587413498185
    0,2,inf,inf,inf,0.0,0.9053982471445895,0.09460175285541045,28.28532777914046
    0,3,inf,inf,-2.252871698,-2.4048796454,0.5341089851136952,0.46589101488630474,28.28532777914046
    0,4,-1.0064940203,-1.0877293566000001,-1.040020104,-1.1491930907999999,0.5341089851136952,0.46589101488630474,20.989136121554353

Swarm-Level Data
================

Swarm-level data is placed in a single CSV file in the main working directory for the program.
``swarmLog.csv`` contains data related to the asynchronicity of the agents, as well
as information on the best historical and current positions throughout the swarm.
The general format of the columns in this file are

::

    step, minStep, maxStep, histBestAgent, currBestAgent, [ agent[ID]step, ... ]

During operation, agent updates are allowed to proceed asynchronously to maximize 
use of computational resources and minimize idling time. Any agent is permitted to
begin the next PSO step as long as itself and all of its neighbors have completed
the current step. For example: Assume agent 1 has agents 0 and 2 as neighbors.
If each of agents 0, 1, and 2 have completed all required SCFT calculations for
their step 2 positions and have determined the fitness of their step 2 positions,
agent 1 is able to begin its update for step 3 and launch the required SCFT calculations
for step 3, even if agents 0 and 2 must wait for their own neighbors to complete step
2 calculations. Thus, for a period of time agent 1 is technically a step ahead of its
neighbors.

In ``swarmLog.csv``, the column ``step`` indicates the most recent step to have been
completed by all agents. The columns ``histBestAgent`` and ``currBestAgent`` both
specify information assuming that all agents had just completed the current step,
ignoring the asynchronicity. The ``histBestAgent`` value gives the ID number of 
the agent whose personal best seen position by the current step has the best fitness
among the personal best positions of all agents in the swarm; essentially, the agent
whose personal best-seen position is the best-seen position by the entire swarm
or the "global best" position at the end of that step.
The ``currBestAgent`` is the ID number of the agent whose current position has the
best fitness among all agent's current positions; this considers only the current
position of each agent at the end of the current step, independent of the personal
histories of the agents.

The remaining columns all relate to the asynchronicity of the agent updates.
``minStep`` indicate the lowest step number completed by any agent,
and should match the ``step`` value. ``maxStep`` indicates the highest step
number completed by any agent. The remaining columns (``agent[ID]step``) indicate
the most recent step completed by each agent, with one column per agent, and with *[ID]*
replaced by the integer ID of the agent in the header. Together, these columns give
a snapshot of how out-of-sync the swarm was at the moment when the last agent completed
the final calculation for the step identified in the ``step`` column.

An example of this file for the first four steps of a 5-agent search are shown below
(taken from the same calculation as the last two examples on this page).

::

    step,minStep,maxStep,histBestAgent,currBestAgent,agent0step,agent1step,agent2step,agent3step,agent4step
    0,0,1,3,3,0,0,1,1,0
    1,1,1,3,3,1,1,1,1,1
    2,2,3,3,3,2,2,2,3,3
    3,3,4,3,3,3,3,3,4,3
    4,4,4,3,4,4,4,4,4,4

Tools for Analysis of Output
============================

PscfInverse generates a lot of output files, and analysis of much of this data
requires cross-referencing several of these tabulated datasets. To assist in 
analysis of the PSO data, a Matlab script is included with the ``matlab`` directory
of the repository to assist in this cross-referencing.
This script, defining the Matlab function ``parsePsoData``, should be run from
the program's working directory (where ``swarmLog.csv`` is located) and takes 
the number of agents in the swarm as an argument. The function then collects data
from all of the output files and returns three Matlab table objects representing
collations of the data files. The first of these is ``agentdata``, which
combines the ``psoData.csv`` and ``phaseData.csv`` files from each agent into
a single table listing the data collected for each agent at each step.
The second table ``pbest`` uses the ``bestStep`` data to list the data from
each agent's personal best position at each step. Unlike in the original CSV
files or in the agentdata table, for each step that an agent does not find a new
best position, the data from the reigning best position is repeated in the row.
In this table all position and energy data listed is that at the best position seen
by the agent through the specified step, rather than the current position data at that
step.
The final table is ``gbest`` which lists data about the best position seen in the entire
swarm's history by the end of the step identified in the ``step`` column. This table
contains only one row per PSO step, identifying the best positions seen by the swarm
over time.

This function would be called with the line

::

    [agentdata, gbest, pbest] = parsePsoData(n_agent)

After these three tables are returned from the function, they may be saved to
.mat files in order to use them for later analyses without having to re-run the
parsing script.
