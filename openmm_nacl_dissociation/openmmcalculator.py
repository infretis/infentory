import openmm
import openmm.app
import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase import units as ase_unit
import openmm.unit as omm_unit

class OpenMMCalculator(Calculator):
    """
    Minimalistic OpenMM calculator.

    To control the number of cpus, run
        export OPENMM_CPU_THREADS=1

    To choose the platform, run
        export OPENMM_DEFAULT_PLATFORM={CPU, CUDA,...,}
    """
    def __init__(self):
        super().__init__()
        system_xml = "openmm_input/system.xml"
        topology_pdb = "openmm_input/topology.pdb"

        self.subcycles = 5
        self.timestep = 0.5*omm_unit.femtoseconds
        integrator = openmm.VerletIntegrator(self.timestep)

        # setup system, topology and create the openmm simulation object
        with open(system_xml) as rfile:
            system = openmm.XmlSerializer.deserialize(rfile.read())

        pdb = openmm.app.PDBFile(topology_pdb)
        topology = pdb.topology
        self.simulation = openmm.app.Simulation(topology, system, integrator)
        self.context = self.simulation.context
        # for outputting the potential energy at each subcycles step
        self.implemented_properties = ["energy"]

        # print some Hardware info
        pf = openmm.Platform.getPlatformByName("CPU")
        print("Using N_threads =", pf.getPropertyValue(self.context, "Threads"))

    def calculate(self, atoms = None, properties = None, system_changes = all_changes):
        # only set context positions once when propagating from a phasepoint
        if atoms.info.get("not_setup", False):
            print("setting positions, velocities and cell")
            self.context.setPositions(atoms.positions/10)
            self.context.setVelocities(atoms.get_velocities()/10*1000*ase_unit.fs)
            self.context.setPeriodicBoxVectors(*atoms.cell.array/10)
            atoms.info["not_setup"] = False
        else:
            self.simulation.step(self.subcycles)

        state = self.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
        energy = state.getPotentialEnergy()._value
        atoms.positions = state.getPositions(asNumpy=True)*10
        atoms.set_velocities(state.getVelocities(asNumpy=True)*10/1000)
        self.results["energy"] = energy * 0.0103641 # eV
        self.results["forces"] = None
        self.results["stress"] = None
