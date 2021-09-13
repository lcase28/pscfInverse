.. |name| replace:: ``Phases{}``
.. |tgt| replace:: ``TargetPhase{}``
.. |cmp| replace:: ``CompetingPhase{}``

.. _param-phases:

*****************
Phases Block
*****************

.. summary

The ``Phases{}`` block defines all phases considered in the search,
whether they are a target or competitor, and which SCFT solver
is to be used. The energies of each of these phases is used to
calculate the Fitness at each point in search space, :math:`{\Omega}`,
using

.. math:
    
    {\Omega} = \min_{m \in C} F_{m} - \min_{n \in T} F_{n}

for set :math:`T` of target phases and set :math:`C` of competitors,
and where :math:`F_{i}` is the free energy of phase :math:`i`.

.. summary

In order to determine the relative stability of a target phase,
and thus obtain the fitness at a point in the search space,
PscfInverse needs a pre-determined list of competing morphologies
to consider at each point in the search space. 
This set of morphologies is defined in the ``Phases{ ... }``
block of the parameter file. An example of such a block that could
be used in a simple search of conformationally symmetric diblock
copolymers seeking BCC spheres (bcc) against competing hexagonally packed
cylinders (hex), double gyroid (gyr), and lamellae (lam).

::

    Phases{
        scft_solver pscf
        input_root  inputs/
        TargetPhase{
            name    bcc
            model_file  bcc/model
        }
        CompetingPhase{
            name    hex
            model_file  hex/model
        }
        CompetingPhase{
            name    gyr
            model_file  gyr/model
        }
        CompetingPhase{
            name    lam
            model_file  lam/model
            range_check rho_rgrid.out   0.001
        }
    }

You may notice the absence of disorder among the competing morphologies
defined in the above example.
**Disorder is considered to be a competing phase implicitly 
and automatically accounted for by the program,
and does not need to be specified by the user as a competing phase.**
To do this, the program makes use of the fact that PSCF automatically
includes the homogeneous disordered free energy in the output file from
a calculation. When reading results from PSCF, PscfInverse first
determines the energy of each user-specified phase relative to disorder.
With the energy of each phase stored using homogeneous disorder as the
reference state, the search for the energy of the most stable competitor
can be initialized to :math:`0`, for the implicit disordered phase, and
proceed to check the user-specified competitors against that initial,
implicit disordered phase.

Within the ``Phases{ ... }`` block, the user first specifies the
scft solver to use for the calculation using the label ``scft_solver``.
Following this label, a string key is entered to identify the
desired SCFT solver. Presently, the only available solver is the
Fortran version of PSCF (using key *pscf*). As additional solvers
are added to the program, a list of available solvers and keys will
be added here.

Next, the user may (optionally) specify a shared root for the
phase input files. This entry uses the label ``input_root``
and is followed by a path to the phase input file root directory.
The path specified for this entry must be written relative to the
working directory of the program (i.e. the directory from which
the program was executed).
If this option is omitted, the input root is then taken to be
the program's working directory.

Finally, the user will specify the list of target and competing
phases. Each target phase is declared within its own 
``TargetPhase{ ... }`` block, while each competitor is 
defined in its own ``CompetingPhase{ ... }`` block.
These two declaration types differ only in the name given to
the block (*TargetPhase* or *CompetingPhase*), while their
internal syntax is identical. Each phase, whether a target
or competitor, must be given a unique name and a model file.
The name, speficied using the label ``name``, must be unique
among all target and competing phases, and must also represent
a valid unix-style directory name. The name will be used
by PscfInverse to label the phase in output data as well as to
create working directories for each of that phase's SCFT
calculations.
The model file, specified using the label ``model_file``,
must provide the path to a valid PscfFieldGen model file.
If an *input_root* was specified, the path specified here
must be written relative to that directory; otherwise,
it must be written relative to the program's working directory.
Although the model file specified here would contain its own
declaration of *scft_solver*, the declaration there should match
that given earlier in the *Phases* block (this is not checked by
the program, but a mismatch would likely cause the program to crash).
The model file identified here is used to set up the initial guess
field generator for the SCFT calculations, as well as to identify the
template parameter file being used for the calculations.
Specific requirements for the model file can be found in instructions
for PscfFieldGen.

Each phase block can also include the optional label ``range_check``.
This optional input should be placed immediately following the model
file input.
If present, this label sets the phase to check the total variation
in the final real-space density field to screen for a collapse
to disorder. In order to use this check, the template PSCF parameter file
must include a directive to convert the final density field output
by the ``ITERATE`` command (the *rho* file) into an *rgrid* density
field. In PSCF, this would be the ``FIELD_TO_RGRID`` command.
The format for this entry would be

::

    range_check [rgrid_field_file_name] [field_tolerance]

where ``[rgrid_field_file_name]`` is the name of the final
rgrid density field file, relative to the SCFT working directory
(i.e. it should match the *output_filename* value from the
*FIELD_TO_RGRID* command), and where ``[field_tolerance]``
is a floating point value identifying the minimum required 
variation in the density field to be considered "not disordered".

Finally, each phase block allows a set of core monomer options to be
to be specified (core monomer inputs in the model file are ignored).
These should be listed after the ``range_check``
input (if present) and should be the last inputs included for the phase.
For each monomer that should be considered an option for selecting the
structure cores, include the entry ``core_option    [monomerID]``.
If no core options are specified, all monomers in the system are considered
potential core monomers. **Note: If the phase being described is lamellar, 
these core options are ignored during field generation.**

Thus, a hypothetical phase declaration block with all optional details could
look like the following:

::

    TargetPhase{
        name    bcc
        model_file  bcc/model
        range_check out.rho_rgrid   0.001
        core_option 0
        core_option 2
        core_option 3
    }

In this example, a range check is used on the final rgrid field to
check for disorder. Additionally, monomer 1 is omitted from the core
selection process (as are monomers 4+, if they exist in the system.)
