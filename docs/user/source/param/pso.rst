
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
        random_seed 100
        StandardIntegrator{
            constriction    0.70
            self_weight     2.05
            neighbor_weight 2.05
        }
        Swarm{
            n_agent 5
        }
        n_step  100
    }



