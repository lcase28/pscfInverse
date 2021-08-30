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

The variables seen by the PSO algorithm do not necessarily
map identically to a particular polymer parameter or entry
in a PSCF input file. Instead, the variables seen by the
PSO algorithm may represent relationships between two or 
more polymer parameters. Similarly, a polymer parameter may
be related to two or more other parameters by two or more
variables. This means that the available PSO variables
are restricted to a pre-defined set; it also means that 
more complex relationships can be represented as variables
if separately implemented, as long as the relationship can
be expressed as a linear relationship of parameters.

The "parameter relationship" construction used to define 
the search variables creates the possibility that the
search variables alone may not adequately define the
polymer parameters involved. In this instance, constraints
must be defined in order to map variable values to polymer
parameters.

The |name| block of the parameter file defines
both the search variables, and any constraints necessary
to map those variables to polymer parameters. These two
roles ('variable' and 'constraint') are reflected in the
structure of the block.

.. .. literalinclude:: searchspace_ex

As shown above, |name| contains two sub-blocks: |vars|
and |cons|. Variables defined in the former represent
searchable dimensions, while those defined in the latter
represent fixed relationships.

.. _param-vartype-sec:

Variable Types
==============

Internally, |program| represents both variables and
constraints with the same classes. Thus, any parameter
relationship that can be used as a variable can also be
used as a constraint. 

Because of their differing roles,
the definition syntax for the same parameter relationship
varies based on its placement in |vars| or |cons|.
When defined in |vars|, the parameter relationship requires
that a lower bound, upper bound, and maximum velocity magnitude
be given.
When defined in |cons|, the parameter relationship requires
that its fixed value be defined.

Parameter Relationships:

..
    ===========================   ====================================
    Variable                      Description
    ===========================   ====================================
    :ref:`param-blocklen-sub`     Total length of one or more 
                                  polymer blocks.
    :ref:`param-blockratio-sub`   Log of ratio between total 
                                  lengths of one or more polymer
                                  blocks.
    :ref:`param-kuhnlen-sub`      Statistical segment length of
                                  a monomer.
    :ref:`param-kuhnratio-sub`    Log of ratio between statistical
                                  segment lengths of two monomers.
    :ref:`param-chi-sub`          Interaction parameter between two
                                  monomers.
    :ref:`param-blendfrac-sub`    Total volume fraction of one or 
                                  more species.
    :ref:`param-blendratio-sub`   Log of ratio between total volume
                                  fractions of one or more species.
    ===========================   ====================================
    
    .. include:: vars/blocklen.rst
    
    .. include:: vars/blockratio.rst
    
    .. include:: vars/kuhnlen.rst
    
    .. include:: vars/kuhnratio.rst
    
    .. include:: vars/chiinteraction.rst
    
    .. include:: vars/blendfrac.rst
    
    .. include:: vars/blendratio.rst
    
