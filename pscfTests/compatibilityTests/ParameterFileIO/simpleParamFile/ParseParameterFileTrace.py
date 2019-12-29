#from ..ParamFactories import *
from context import psoinverse
from psoinverse import ParamFactories as pf

pfac = pf.PSCFPPFactory()

#fname = "/Users/logancase/Documents/Research/PSCF/CPP_Update/examples/fd1d/cylindrical/1/param"
#fname = "/Users/logancase/Documents/Research/PSO/Development/Primary/psoinverse/params_HEXortho_template.in"
fname = "param_example"

pfac.ParseParameterFile(fname)

print('\nAFTER PARAMETER READ:\n')
print(pfac)
