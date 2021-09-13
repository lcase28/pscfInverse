.. |name| replace:: ``SearchSpace{}``
.. |vars| replace:: ``Variables{}``
.. |cons| replace:: ``Constraints{}``

.. _param-searchspace:

*****************
SearchSpace Block
*****************

.. summary

The ``SearchSpace{}`` block defines all search variables
to use in the PSO calculation, as well as all constraints
required to uniquely define the values of all polymer
parameters impacted by those search variables.

.. summary

When running SCFT calculations directly using the PSCF software,
users define the polymer system, as well as other aspects
of their calculation, in the PSCF Parameter File.
The values defined in this file include, for example, 
lengths of individual blocks in a polymer molecule,
the statistical segment length of a particular monomer,
and the Flory-Huggins interaction parameter between sets of monomers.
These values directly define the polymer system at the coarse-grained
level considered in the mean-field calculation
and are referred to in the following text as "polymer parameters"
or simply "parameters".

The PSO search performed by ``PscfInverse`` operates on its own
set of search variables ("variables" in the following text).
These variables do not necessarily
map directly to a particular polymer parameter.
Such direct-mapping variables are available in ``PscfInverse``,
but generally, the variables seen by the PSO algorithm
represent mathematically linearizeable relationships between
two or more polymer parameter values.
Indeed, at a program level, even direct-map variables are treated
as a "relationship" involving only one parameter.
During program operation,
values of polymer parameters are determined from these relationships
via solution of the linear equations defined by these relationships.
While this construction does restricts the available PSO variables
to a pre-defined set, it also allows more flexible control of
the PSO search by allowing more complex relationships than a direct
map from variables to paramters. The set of available search variables
can also be expanded to include any linearizeable relationship of interest
via modest programming effort.

Depending on the search being performed, the search variables alone may not
be sufficient to fully define a solvable set of linear parameter equations.
To ensure that all systems are fully defined, any parameter relationship accessible
as a variable can also be used as a fixed constraint on the system.
In this case, the PSO algorithm does not modify the values in the relationship,
but the relationship is still used in determination of the parameter values.

Finally, the relationships defined in this section need only reference
the parameters involved in the search variable relationships.
For each phase considered in the calculation, the user must
provide a template PSCF Parameter file.
(For more information on this, see :ref:`Phases Block Main Page <param-phases>`).
Parameters defined in the template file but not referenced in the
variables or constraints retain the value defined in the template.
Thus, in the case of a 2D search of neat diblock polymers over 
:math:`f_{A}` and :math:`{\chi}N`, 
the user need not define any variables or constraints
related to statistical segment lengths, but can instead just set the
fixed values in the template file.

To differentiate variable relationships from constraint relationships,
each set of relationships is defined in a separate sub-block.
This gives the ``SearchSpace`` block the skeletal structure shown below.

::

        SearchSpace{
            Variables{
                ...
            }
            Constraints{
                ...
            }
        }

As shown above, |name| contains two sub-blocks: |vars|
and |cons|. Variables defined in the former represent
searchable dimensions, while those defined in the latter
represent fixed relationships.

.. _param-varconventions-sec:

Conventions
===========

Polymer and Block Indexing
--------------------------

While defining the search space, multiple available relationships
refer to the length of various blocks in the polymer system.
Within the PscfInverse parameter file, references to individual
blocks are made with pairs of integer ID or index numbers.
Each reference to a block in the polymer system is formatted as
shown below, with the keyword ``block`` followed by the pair of
ID numbers, represented by the labels ``[polymerID]`` and ``[blockID]``.

::

    block   [polymerID]     [blockID]

As with other aspects of the parameter file, specific formatting
of these references is flexible, as long as whitespace separates
all three elements, and whitespace preceeds and follows the whole
group.
Most available relationships will use this format. If an alternate
format is required for any relationships, they will be specified in
that relationship's description below. In any case, indexing will
be consistent.

The values of the ID numbers are determined by the PSCF parameter
file format. Both polymers and blocks are indexed starting at :math:`0`
with individual polymers and blocks indexed as they appear.
Thus, *polymer 0* will be the first polymer specified in the 
PSCF parameter file, *polymer 1* is the second polymer specified, and
so on. Block indexing is relative to each polymer ID. Thus,
for each polymer, the first block specified will be *block 0* and the
second block specified will be *block 1*.

This indexing is well-illustrated by considering the format of the
PSCF (Fortran) parameter file. In that format, block_monomer
and block_length values are specified in multi-row arrays,
with each row being potentially a different length, such that each
row contains the data for one polymer species and each entry within
a single row refer to different blocks in that polymer species.
Thus, the block indexing used here is such that [polymerID] identifies
the row of the block_length array, and [blockID] identifies the entry
within that row that is being considered in the relationship.

Using this convention, a blend of two diblock copolymers will have
the following four references as possibilities:

::

    block   0   0
    block   0   1
    block   1   0
    block   1   1

Monomer Indexing
----------------

As with polymers and blocks, monomers are indexed starting at :math:`0`.
The order of the monomer indexing likewise matches the order in which
the monomers are specified in the SCFT input file. In the case of
PSCF (Fortran), monomers are (internally) indexed starting at 1
in the order they are specified in the *kuhn* input array.
In this case, the index used here for a monomer is one less than the
index used in the PSCF parameter file, creating a "shift" in the indexing
scheme. 

.. _param-varvscon-sec:

Variables versus Constraints
============================

Internally, PscfInverse represents both variables and
constraints with the same classes. Thus, any parameter
relationship that can be used as a variable can also be
used as a constraint. Because of their differing roles,
the definition syntax for the same parameter relationship
varies based on its placement in |vars| or |cons|, however
the difference between the two uses is the same for all
relationship types. This section summarizes the entries
required for each relationship when used as a variable
or a constraint. The entries specified here are in addition
to the relationship-specific data outlined in the
:ref:`Relationship Types Section <param-vartypes-sec>`.

Relationships as Variables
--------------------------

When defined in |vars|, the parameter relationship requires
that a lower bound, upper bound, and maximum velocity magnitude
be given. Here, the lower bound and upper bound represent the
lowest and highest values the variable can assume during the
search. Altogether, the lower and upper bounds of all variables
in the search define full bounds of the search space.
The search bounds in PscfInverse are hard, reflective boundaries
such that, in any dimension, if a given PSO step would carry an
agent outside of these bounds, the position update in that dimension
is rejected and its velocity in that dimension is reversed
(as if the agent had "bounced" off a wall).
The maximum velocity, due to the unit time step applied in the PSO
algorithm, represents the maximum distance along that dimension
any agent is permitted to move in a single PSO step.
Any positive value may be chosen for this entry;
however, given the use of reflective bounds, the largest step
that can actually be taken is the full distance from lower to
upper bound. Allowing velocity values to exceed this logical limit
is permitted and, in theory, acceleration of the agent should
eventually lower the velocity to a smaller value as the swarm converges.
The effects of such large velocities on the performance of PscfInverse
has not been investigated, and choosing a velocity cap below
the full search space width is generally recommended.
The bounding and velocity entries are formatted uniformly across
relationship types and are summarized below.

    ==============  ================    =========================
    Variable        Label               Description
    ==============  ================    =========================
    Lower Bound     ``lower``           The lowest value allowed
                                        for the variable during 
                                        the search.
    Upper Bound     ``upper``           The highest value allowed
                                        for the variable during
                                        the search.
    Velocity Limit  ``velocity_cap``    *Optional*
                                        The fastest an agent is
                                        allowed to move in the 
                                        positive or negative
                                        direction along this 
                                        search dimension.
                                        (Symmetric about
                                        :math:`velocity = 0.0`).
                                        Default value is 
                                        :math:`(upper - lower)`.
    ==============  ================    =========================


Relationships as Constraints
----------------------------

When defined in |cons|, the parameter relationship requires
that its fixed value be defined. In this case, the bounds
and velocity limit defined for use as a variable are not
required. The constant value is specified with the label
``value`` followed by whitespace followed by the numerical value.

.. _param-vartypes-sec:

Relationship Types
==================

Each of the following subsections defines and describes one of
the available parameter relationships that can be used as variables
or constraints in PscfInverse.

    ===========================   ====================================
    Variable                      Description
    ===========================   ====================================
    Block Lengths                 Total length of one or more 
                                  polymer blocks.
    Block Length Ratios           Log of ratio between total 
                                  lengths of one or more polymer
                                  blocks.
    Statistical Segment Length    Statistical segment length of
                                  a monomer.
    Segment Length Ratio          Log of ratio between statistical
                                  segment lengths of two monomers.
    :math:`{\chi}` Parameter      Interaction parameter between two
                                  monomers.
    ===========================   ====================================

Block Lengths
-------------

A BlockLength variable relates one or more 
polymer blocks to a single, combined length
value. These blocks need not be adjacent to
each other, nor even in the same polymer chain.
Mathematically, the variable value is given by

.. math::
    BlockLength = \sum_{i=1}^{N} L_{p_i,b_i}

for a variable relating :math:`N` blocks,
when :math:`L_{p_i, b_i}` is the length of
block :math:`b_i` in polymer :math:`p_i`.

The relationship is defined within a block
named ``BlockLength{ ... }`` within either
the variables or constraints block.
Within the BlockLength block, the blocks to
be accounted in the total are specified as
a sequence of block references of the form
``block [polymerId] [blockID]``.
These block specifications should be placed
before either the variable bounds or constraint
values. Altogether, the general format would
resemble one of the following, depending on use
as a variable or constraint.

::

    BlockLength{
        block   0   0
        block   0   1
        lower   0.5
        upper   2.0
        velocity_cap    1.0
    }
    
    BlockLength{
        block   0   0
        block   0   1
        value   1.0
    }

If both examples above refer to a system in which
polymer 0 is a diblock copolymer, then the former
considers the total length of the diblock to be a
variable ranging from 0.5 to 2.0 (presumably relative
to another polymer of length unity), while the latter
fixes the total length of the diblock to unity
(which would perhaps allow the block fractions to be
varied at constant length).

Block Length Ratios
-------------------

The BlockRatio relationship considers the ratio
between the total lengths of two sets of blocks.
These blocks need not be adjacent to
each other, nor even in the same polymer chain.
Mathematically, the variable value is given by

.. math::
    BlockRatio = \ln \frac{\sum_{i=1}^{N} L_{p_i,b_i}}{\sum_{i=1}^{D} L_{p_i,b_i}}

for the ratio between blocks :math:`N`
and blocks :math:`D`,
when :math:`L_{p_i, b_i}` is the length of
block :math:`b_i` in polymer :math:`p_i`.
Here, logarithmic scaling is used to offer equal
numerical representation between a ratio and its
reciprocal.

Block Length Ratio relationships are defined in 
``BlockRatio`` blocks. The format of a ``BlockRatio``
block is similar to that used for ``BlockLength``
relationships, but the two sets of blocks
accounted in the relationship are placed in their
own sub-blocks. In this case,
the blocks accounted in the top half of the ratio
are specified within a ``Numerator{ ... }`` block,
while those in the bottom half of the ratio
are listed in the ``Denominator{ ... }`` block.
In both sets, block references are formatted
according to convention as 
``block [polymerID] [blockID]``.
These block specifications should be placed before
any variable bounds or constraint values.
Note that all variable bounds and constraint values
should refer directly to the logarithmic value.
Depending on use as a variable or constraint,
a BlockRatio relationship would thus be defined
similarly to one of the following:

::

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
    
    BlockRatio{
        Numerator{
            block   0   0
        }
        Denominator{
            block   0   1
        }
        value   0.0
    }

Consider a system in which polymer 0 is a diblock.
In this case, the first definition would allow the
log of the ratios of the
block fractions to vary between -2.5 and 2.5; this
is roughly equivallent to allowing :math:`f_{A}` to
vary from 0.075 to 0.925. The second definition
would force the diblock to remain compositionally
symmetric.


Statistical Segment Lengths
---------------------------

This simplified relationship offers a direct map
to the statistical segment length (or similar metric
used in the SCFT solver) from the PSO variable or
constraint. 
The monomer being referenced is specified with an
entry of the form ``monomer [monomerID]`` where
``[monomerID]`` follows the convention defined above.
The simple syntax for use as a variable
would resemble

::

    KuhnLength{
        monomer 0
        lower   1.0
        upper   2.0
        velocity_cap    0.6
    }

while use as a constraint (to fix the segment length)
would resemble

::

    KuhnLength{
        monomer 0
        value   1.0
    }



Statistical Segment Length Ratios
---------------------------------

This relationship references the ratio
between the statistical segment lengths
of two monomers. Similarly to the ratio of
block lengths, the log of the ratio is used
here as well. The relationship is defined in
a ``KuhnRatio{ ... }`` block. 
Since only two monomers are considered in one
relationship, both monomer IDs are specified
in a single entry of the form
``monomers  [numeratorID]   [denominatorID]``
where *[numeratorID]* is the ID of the monomer
in the numerator of the ratio and 
*[denominatorID]* is the ID of the monomer in
the denominator of the ratio.
As with the Block length ratio relationship,
the KuhnRatio variable bounds and constraint
values expect the logarithmic value.
When used as a variable, the relationship
would be specified with

::

    KuhnRatio{
        monomers    0   1
        lower   -1.0
        upper   1.0
        velocity_cap    1.0
    }

while use as a constraint would resemble

::

    KuhnRatio{
        monomers    0   1
        value   0.693
    }


Flory-Huggins :math:`{\chi}` Parameters
---------------------------------------

This relationship acts as a direct map to
the Flory-Huggins interaction parameter,
:math:`{\chi}_{i,j}` between monomers
:math:`i` and :math:`j`.
This relationship is specified within a
``Chi{ ... }`` block.
Monomers are identified in a similar format
to that used for segment length ratios,
``monomers  [monomerID_1]   [monomerID_2]``
where the two monomer IDs listed in brackets
would be the integer index for the interacting
monomers. The monomer IDs can be listed in any
order, but by convention the lower index 
is typically listed first. The format of this
relationship in the PscfInverse parameter file
would be similar to one of

::

    Chi{
        monomers    0   1
        lower   10.0
        upper   30.0
        velocity_cap    15.0
    }
    
    Chi{
        monomers    0   1
        value   25
    }

The latter format, that for use as a constraint, 
currently should see limited use. Given that
:math:`{\chi}` parameters are not presently used
in any relationships other than direct maps, and
that the program leaves unmodified any portion of
the PSCF parameter file not referenced from the
variables or constraints, :math:`{\chi}` parameters
not used as variables can simply be specified at their
constant values within the parameter file and excluded
from the PSO search space definition;
the value in the parameter file would then be left
as a constant throughout the search.
