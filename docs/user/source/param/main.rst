.. _param-main:

***********************
Parameter File Overview
***********************

While users can write their own scripts to perform a PSO search -
perhaps to access custom features without the need to modify the
main program, for example - calculations are typically run
according to a parameter file, which contains both defines the pso
parameters and lists the actions to be executed.

The parameter file consists of a set of nested blocks enclosed by 
curly braces. At every level, block names use ``CamelCase`` labels - 
in which distinct words capitalize the first letter, while remaining
letters are lower-case - followed immediately by the opening bracket.
No space is placed between the BlockName and the opening bracket.
The end of a block is marked by a lone ``}`` with whitespace both
preceeding and following the bracket. 
Blocks may be written linearly::

    BlockName{ ... }

or they may be split over multiple lines::

    BlockName{
        ...
    }

according to user preference. Blocks then contain combinations of
sub-blocks and parameter values.

.. _param-skeleton:

Primary Structure
=================

The entire parameter file is contained in one "parent" block with
label ``PscfInverse{`` which indicates the "active" portions 
of the parameter file. 
Content written before or after this block is ignored, and can
optionally be used for explanatory comments or other notes by the
user.

This parent block contains several sub-blocks, shown in skeletal 
form below, in order to group related information.

.. literalinclude:: skeleton_ex

The order of these sections is not pre-defined, although all
are required.
Each section is briefly described below, with more detailed
information available on the block's own page (linked 
alongside the brief summary.

:ref: `param-searchspace`
    .. include:: searchspace.rst
        :start-after: summary
        :end-before: summary
:ref: `param-phases`
    .. include:: phases.rst
        :start-after: summary
        :end-before: summary
:ref: `param-pso`
    .. include:: pso.rst
        :start-after: summary
        :end-before: summary
:ref: `param-commands`
    .. include:: commands.rst
        :start-after: summary
        :end-before: summary

