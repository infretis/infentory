from ase import Atoms
from ase.io import write

# Create one atom (use any element, e.g. He as placeholder)
left_well= Atoms(
    symbols="He",
    positions=[(-1.0, 0.0, 0.0)],
    pbc=False
)
right_well= Atoms(
    symbols="He",
    positions=[(0.4, 0.0, 0.0)],
    pbc=False
)

# Write configuration to file
left = "left_well.traj"
right = "right_well.traj"
write(left, left_well)
write(right, right_well)
print("created", left, right, "!")
