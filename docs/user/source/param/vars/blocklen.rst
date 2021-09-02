.. _param-blocklen-sub

BlockLength
-----------

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

General format (selected fields depend on usage):

::

    BlockLength{
        block   0   0
        block   0   1
        value   1.0
        lower   0.5
        upper   2.0
        velocity_cap    1.0
    }


.. summary

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

