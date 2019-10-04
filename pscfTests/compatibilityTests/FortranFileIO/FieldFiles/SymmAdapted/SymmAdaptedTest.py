from context import psoinverse
from psoinverse.SCFT.PSCF.FileManagers.fieldfile import SymFieldFile as sff

fname = 'omega'

field = sff(fname)

print(field)

field.write('omega_out')

