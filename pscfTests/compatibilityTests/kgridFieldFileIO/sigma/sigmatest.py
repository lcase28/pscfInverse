# Imports
from context import psoinverse
from psoinverse.mesophases.FieldGenerators import FieldGenerator

fg = FieldGenerator.from_file("model")

print(fg)

print(fg.to_kgrid([0.25, 0.75]))

fg.to_file([0.25, 0.75])
