.. _param_blocklen_sub

BlockLength
-----------

.. summary

A variable represents the sum of the lengths
of all constituent blocks.

.. summary

A BlockLength variable relates one or more 
polymer blocks to a single, combined length
value. These blocks need not be adjacent to
each other, nor even in the same polymer chain.
Mathematically, the variable value is given by

.. math::
    BlockLength = \sum_{i=1}^{N} Len_{p_i,b_i}

for a variable relating :math:`N` blocks,
when :math:`Len_{p_i, b_i}` is the length of
block :math:`b_i` in polymer :math:`p_i`.

General Syntax
..............

Regardless of its use as a search variable or
constraint, this relation requires the user
to define which blocks count toward the total
length. Within the variable's parameter file
block, each contributing polymer block is 
identified with a parameter labeled ``block``
followed by a polymer and block id, all 
separated by whitespace. If written on the same
line, this would appear as ::

    block   [polymer_id]    [block_id]

with ``[polymer_id]`` representing the ID
of the polymer containing the block and
``[block_id]`` representing the ID of the 
block within polymer *polymer_id*.

Indexing of polymers and blocks starts at 0,
with the order of indexing defined by the 
SCFT solver's interface. Generally speaking,
in instances when polymers and blocks are not
explicitly given ID numbers in the SCFT
solver's input files (such as the Fortan 
version of PSCF), the polymers and blocks
are indexed in the order they first appear in
the SCFT input file. When the SCFT input does
specify block or polymer ID numbers, the
given IDs would be kept in order, but shifted
to start indexing at 0. Thus if the SCFT input
file specified blocks *1*, *2*, and *3* in
polymer *1* (with no polymer *0*), these would
be accessed, respectively, with ::

    block   0   0
    block   0   1
    block   0   2

For example, assume you wish to specify
(as either a variable or constraint) the 
total length of blocks *0* and *1* of
polymer *0*. In this instance, the defining
block in the parameter file would be ::

    BlockLength{
        block   0   0
        block   0   1
        ...
    }

where ``...`` acts as a placeholder for 
required data to make this a variable or 
constraint.

Variable Syntax
...............

When using ``BlockLength`` as a search variable,
the user must specify search bounds and a maximum
velocity (magnitude) allowed for the agents in the
swarm. The following table summarizes these input
parameters.

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
                                        :math:`(upper - lower) / 10`.
    ==============  ================    =========================

.. note::
    Because |program| uses a unit time step in agent
    updates, the velocity limit is numerically
    equvallent to the farthest an agent is allowed 
    to travel in a single step.

Consider a researcher studying neat 
*ABA'C* tetrablock terpolymer systems.
If this researcher want to run a search
allowing the total fraction of *A* to
vary between 0.1 and 0.5. Assuming they
wish to allow velocities up to :math:`1/4`
of the total width. The following could 
be an excerpt from their parameter file::

    PscfInverse{
        SearchSpace{
            Variables{
                ...
                BlockLength{
                    block   0   0
                    block   0   2
                    lower   0.1
                    upper   0.5
                    velocity_cap    0.1
                }
                ...
            }
            ...
        }
        ...
    }


Constraint Syntax
.................

When using ``BlockLength`` as a search constraint
the user must specify the fixed total length.
This is specified with the label ``value`` followed
by a real number.

Suppose that polymer 1 is a diblock copolymer.
Suppose, also, that the user wishes to normalize it
to unit length (such that the total length of the
polymer is 1). The following could be an excerpt from
the required parameter file. ::

    PscfInverse{
        SearchSpace{
            ...
            Constraints{
                ...
                BlockLength{
                    block   1   0
                    block   1   1
                    value   1.0
                }
                ...
            }
        }
        ...
    }

