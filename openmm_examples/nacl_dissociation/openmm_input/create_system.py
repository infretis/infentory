import openmm
import openmm.app
import openmm.unit
from openff.toolkit import ForceField, Molecule, unit

from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box

from ase.io import read

Na = Molecule.from_smiles("[Na+]")
Cl = Molecule.from_smiles("[Cl-]")
water = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")

Nwat = 135 # 383
box_size = 1.6 #2.3
### Naming the residue is not needed to parameterize the system or run the simulation, but makes visualization easier
for atom in water.atoms:
    atom.metadata["residue_name"] = "HOH"

topology = pack_box(
    molecules=[Na, Cl, water],
    number_of_copies=[1, 1, Nwat],
    box_vectors=box_size * UNIT_CUBE * unit.nanometer,
)

positions = topology.get_positions()
topology = topology.to_openmm()
forcefield = openmm.app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
system = forcefield.createSystem(topology, nonbondedMethod = openmm.app.CutoffPeriodic, constraints=openmm.app.HBonds, rigidWater=True, nonbondedCutoff=0.8*openmm.unit.nanometer)

with open ("topology.pdb", "w") as pdbfile:
    openmm.app.PDBFile.writeFile(topology, positions.magnitude*openmm.unit.nanometer, pdbfile)

with open("system.xml", "w") as output:
    output.write(openmm.XmlSerializer.serialize(system))

# finally write ase atoms
atoms = read("topology.pdb")
atoms.write("initial.traj")
