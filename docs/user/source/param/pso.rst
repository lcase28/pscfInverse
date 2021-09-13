
.. _param-pso:

*****************
PSO Block
*****************

.. summary

The ``Pso{}`` block defines details of the PSO search.
This includes defining the PSO parameters, swarm size, and
communication graph, and length of the search.

.. summary

The ``Pso{ ... }`` block defines the details of the PSO algorithm itself.
The general format of this block is as shown below.

::

    Pso{
        random_seed     1234
        StandardIntegrator{
            constriction    0.729
            self_weight     2.05
            neighbor_weight 2.05
        }
        Swarm{
            n_agent 5
        }
        n_step  100
    }

The first entry in this block is the (optional) specification
of the random number generator seed. This entry uses the 
form ``random_seed  [seed]`` followed by an integer 
with 0 <= [seed] <= (2^32 - 1). If this option is omitted,
a seed will be automatically selected based on the system
clock. In either case (included or omitted), the random seed
used in the calculation is printed to standard output.

The next entry is an optional block specifying parameters for
the PSO integration. The ``StandardIntegrator{ ... }`` block
accepts three optional parameters corresponding to the
standard PSO update equation. If any parameter is omitted,
its standard value (shown in the above example) is selected.
If the block overall is omitted, all standard values are used
(giving the PSO "Standard" update equation).

The first required input in the ``Pso{ ... }`` block is the
``Swarm{ ... }`` block. At present, the swarm block requires
only the number of agents in the swarm (using the label *n_agent*
followed by an integer), and has no optional entries.
Several options will be available for the swarm block in future
updates (this is the reason we use a separate sub-block
for the swarm data), but these are still under development.

Finally, the number of steps to be completed in the PSO search
is given after the swarm block in the form ``n_step [step_count]``,
where *[step_count]* is a positive integer. When the search is complete
each agent will have visited *( [step_count] + 1)* positions,
specifically an initial position 0, and a position at the end of
each step (1 through *[step_count]*).
