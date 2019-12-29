from context import psoinverse
from psoinverse.SCFT.PSCF.FileManagers.fieldfile import CoordFieldFile as sff

fname = 'rho_rgrid'

field = sff(fname)

print(field)

field.write('rho_rgrid_out')

