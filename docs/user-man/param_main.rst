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
curly braces. At every level, block names use ``CamelCase{`` labels
followed immediately by the opening bracket.
