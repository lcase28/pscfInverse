from context import psoinverse
from psoinverse.SCFT.PSCF.FileManagers.fieldfile import WaveVectFieldFile as sff

fname = 'rho_kgrid'

field = sff(fname)

print(field)

field.write('rho_kgrid_out')

