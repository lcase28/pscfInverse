.. _param-main:

**************
Parameter File
**************

While users can write their own scripts to perform a PSO search -
perhaps to access custom features without the need to modify the
main program, for example - calculations are typically run
according to a parameter file, which both defines the pso
parameters and lists the actions to be executed.

The parameter file consists of a set of nested blocks enclosed by 
curly braces. At every level, block names use ``CamelCase`` labels - 
in which distinct words capitalize the first letter, while remaining
letters are lower-case - followed immediately by the opening bracket.
No space is placed between the BlockName and the opening bracket.
The end of a block is marked by a lone ``}`` with whitespace both
preceeding and following the bracket.
Here, "whitespace" includes such characters as spaces, tabs, line-breaks, etc.
Blocks may be written linearly::

    BlockName{ ... }

or they may be split over multiple lines::

    BlockName{
        ...
    }

according to user preference as long as required whitespace separates
all blocknames, data labels, data values, closing brackets, etc.
Blocks then contain combinations of
sub-blocks and parameter values.

.. _param-sections-sec:

Summary of Sections
===================

The entire parameter file is contained in one "parent" block with
label ``PscfInverse{`` which indicates the "active" portions 
of the parameter file. 
This parent block contains several sub-blocks, shown in skeletal 
form below, in order to group related information.
The order of these blocks is pre-defined and should be followed. 

::

    PscfInverse{
        SearchSpace{
            ...
        }
        Phases{
            ...
        }
        Pso{
            ...
        }
    }

Each section is briefly described below, with more detailed
information available on the block's own page (linked 
alongside the brief summary).

SearchSpace Block
-----------------

.. include:: searchspace.rst
    :start-after: summary
    :end-before: summary

More detail on the :ref:`SearchSpace Block Main Page <param-searchspace>`

Phases Block
------------

.. include:: phases.rst
    :start-after: summary
    :end-before: summary

More detail on the :ref:`Phases Block Main Page <param-phases>`

PSO Block
---------

.. include:: pso.rst
    :start-after: summary
    :end-before: summary

More detail on the :ref:`PSO Block Main Page <param-pso>`

.. _param-example-sec:

Example
========

Below is a complete example file for a 2-dimensional search in diblock
polymer phase space varying fA and ChiN and targeting the BCC sphere phase.

::

    PscfInverse{
        SearchSpace{
            Variables{
                BlockRatio{
                    Numerator{
                        block   0   0
                    }
                    Denominator{
                        block   0   1
                    }
                    lower   -2.5
                    upper   2.5
                    velocity_cap    3.0
                }
                Chi{
                    monomers    0   1
                    lower   10
                    upper   30
                    velocity_cap    10
                }
            }
            Constraints{
                BlockLength{
                    block   0   0
                    block   0   1
                    value   1.0
                }
            }
        }
        Phases{
            scft_solver pscf
            input_root  inputs
            TargetPhase{
                name    bcc
                model_file  bcc/model
            }
            CompetingPhase{
                name    gyr
                model_file  gyr/model
            }
            CompetingPhase{
                name    hex
                model_file  hex/model
            }
            CompetingPhase{
                name    lam
                model_file  lam/model
            }
        }
        Pso{
            random_seed 1234
            StandardIntegrator{
                constriction    0.7
                self_weight     2.05
                neighbor_weight 2.05
            }
            Swarm{
                n_agent 5
            }
            n_step  4
        }
    }

