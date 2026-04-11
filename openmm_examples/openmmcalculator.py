import openmm
import openmm.app
import openmm.unit as omm_units
import numpy as np
from sys import stdout
from openmm.app import StateDataReporter

from ase.calculators.calculator import Calculator, all_changes
from ase import units as ase_units


class OpenMMCalculator(Calculator):
    """
    """

    def __init__(self, subcycles, temperature, timestep_ps, friction_per_ps):
        super().__init__()
        system_xml = "openmm_input/system.xml"
        topology_pdb = "openmm_input/topology.pdb"

        self.subcycles = subcycles
        self.timestep = timestep_ps * omm_units.picoseconds
        self.temperature = temperature*omm_units.kelvin
        self.friction_coef = friction_per_ps/omm_units.picoseconds

        integrator = openmm.LangevinIntegrator(self.temperature, self.friction_coef, self.timestep)

        # setup system, topology and create the openmm simulation object
        with open(system_xml) as rfile:
            system = openmm.XmlSerializer.deserialize(rfile.read())
        pdb = openmm.app.PDBFile(topology_pdb)
        topology = pdb.topology
        self.simulation = openmm.app.Simulation(topology, system, integrator)
        #self.simulation.reporters.append(StateDataReporter(stdout, 1, step=True, potentialEnergy=True, kineticEnergy=True,temperature=True))

        # get degrees of freedom for calculating the temperature
        self.dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > 0*omm_units.dalton:
                self.dof += 3
        # Subtract constrained DOFs
        self.dof -= system.getNumConstraints()

        # output printing for testing
        self.context = self.simulation.context
        # for outputting the potential energy at each subcycles step
        self.implemented_properties = ["energy"]

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        # only set context positions once when propagating from a new phasepoint
        if atoms.info.get("not_setup", False):
            self.update_omm_context(atoms)
            if atoms.info.get("reverse_vel", False):
                self.simulation.step(1)
                state = self.context.getState(getVelocities=True)
                vel = state.getVelocities(asNumpy=True)/omm_units.angstrom*omm_units.femtosecond/ase_units.fs
                atoms.set_velocities(-vel)
                self.update_omm_context(atoms)
                atoms.info["reverse_vel"] = False
            atoms.info["not_setup"] = False
        else:
            self.simulation.step(self.subcycles)

        state = self.context.getState(
            getEnergy=True, getPositions=True, getVelocities=True
        )

        # set positions and velocities in ASE units
        atoms.positions = state.getPositions(asNumpy=True) / omm_units.angstrom
        atoms.set_velocities(
            state.getVelocities(asNumpy=True) / omm_units.angstrom * omm_units.femtosecond / ase_units.fs
        )
        atoms.set_cell(state.getPeriodicBoxVectors(asNumpy=True) / omm_units.angstrom)

        # set potential, temperature, and kinetic energy here, as ASE might not be
        # aware of constraints to compute the correct temperature
        vpot = state.getPotentialEnergy()
        ekin = state.getKineticEnergy()
        # convert omm units to eV
        self.results["energy"] = (
            vpot / omm_units.kilojoule * omm_units.mole * ase_units.kJ / ase_units.mol
        )
        self.results["ekin"] = (
            ekin / omm_units.kilojoule * omm_units.mole * ase_units.kJ / ase_units.mol
        )
        self.results["temp"] = (
            2
            * ekin
            / (
                self.dof
                * omm_units.BOLTZMANN_CONSTANT_kB
                * omm_units.AVOGADRO_CONSTANT_NA
            )
            / omm_units.kelvin
        )
        self.results["forces"] = np.array([[]])

    def update_omm_context(self, atoms):
        """Update positions and velocities of the openmm context with the ASE
        positions, velocities and box vectors.
        """
        self.context.setPositions(atoms.positions / ase_units.nm)
        self.context.setVelocities(
            atoms.get_velocities() / ase_units.nm * (1000 * ase_units.fs)
        )
        self.context.setPeriodicBoxVectors(*atoms.cell.array / ase_units.nm)
        # ensure velocity constraints are applied before propagation, such that the
        # velocites generated with ASE have their constraint components removed
        self.context.applyVelocityConstraints(1e-12)
