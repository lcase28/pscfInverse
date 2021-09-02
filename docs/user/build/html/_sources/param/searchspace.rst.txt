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

.. _param-vartype-sec:

Relationship Types
==================

Internally, PscfInverse represents both variables and
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

.. include:: vars/blocklen.rst
    :start-after: summary
    :end-before: summary
    
More details on :ref:`Main Page <vars/param-blocklen-sub>`.

Block Length Ratios
-------------------

.. include:: vars/blockratio.rst
    :start-after: summary
    :end-before: summary

Statistical Segment Lengths
---------------------------

.. include:: vars/kuhnlen.rst
    :start-after: summary
    :end-before: summary

Statistical Segment Length Ratios
---------------------------------

.. include:: vars/kuhnratio.rst
    :start-after: summary
    :end-before: summary

Flory-Huggins :math:`{\chi}` Parameters
---------------------------------------

.. include:: vars/chi.rst
    :start-after: summary
    :end-before: summary
