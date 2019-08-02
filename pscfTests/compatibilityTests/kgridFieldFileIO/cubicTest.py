# Imports
from context import psoinverse
from psoinverse.mesophases.FieldGenerators import FieldGenerator

fg = FieldGenerator.from_file("model_in.txt")

print(fg)

print(fg.to_kgrid([0.75, 0.25]))

fg.to_file([0.75, 0.25])
