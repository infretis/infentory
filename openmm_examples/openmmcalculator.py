import openmm
import openmm.app
import openmm.unit as omm_units
import numpy as np
from sys import stdout
from openmm.app import StateDataReporter

from ase.calculators.calculator import Calculator, all_changes
from ase import units as ase_units
from ase.io import read


class OpenMMCalculator(Calculator):
    """
    """

    def __init__(self, subcycles, temperature, timestep_ps, friction_per_ps, halfstep_vel):
        super().__init__()
        print("="*10,"init OpenMMCalculator")
        system_xml = "openmm_input/system.xml"
        topology_pdb = "openmm_input/topology.pdb"

        self.subcycles = subcycles
        self.timestep = timestep_ps * omm_units.picoseconds
        self.temperature = temperature*omm_units.kelvin
        self.friction_coef = friction_per_ps/omm_units.picoseconds
        self.halfstep_vel = halfstep_vel

        integrator = openmm.LangevinIntegrator(self.temperature, self.friction_coef, self.timestep)
        #integrator = openmm.VerletIntegrator(self.timestep)
        #integrator = openmm.CustomIntegrator(self.timestep)
        #integrator.addPerDofVariable("x1", 0)
        #integrator.addUpdateContextState()
        #integrator.addComputePerDof("v", "v+0.5*dt*f/m")
        #integrator.addComputePerDof("x", "x+dt*v")
        #integrator.addComputePerDof("x1", "x")
        #integrator.addConstrainPositions()
        #integrator.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        #integrator.addConstrainVelocities()

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
        with open("calculator_output","a") as wfile:
            wfile.write(f"Active Platform: {self.simulation.context.getPlatform().getName()}\n")

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        # only set context positions once when propagating from a new phasepoint
        if atoms.info.get("not_setup", False):
            #print("Setting up atoms")
            atoms.info["not_setup"] = False
            self.update_omm_context(atoms)
            if self.halfstep_vel and atoms.info.get("reverse_vel", False):
                #print("Reversing velocities with single step MD propagation")
                self.simulation.step(1)
                state = self.context.getState(getVelocities=True)
                vel = state.getVelocities(asNumpy=True)/omm_units.angstrom*omm_units.femtosecond/ase_units.fs
                atoms.set_velocities(-vel)
                self.update_omm_context(atoms)
                atoms.info["reverse_vel"] = False
            elif atoms.info.get("reverse_vel", False):
                #print("Reversing velocities")
                vel = atoms.get_velocities()
                atoms.set_velocities(-vel)
                self.update_omm_context(atoms)
                atoms.info["reverse_vel"] = False
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

    def generate_velocities(self, conf_in, conf_out):
            #print("Generating velocities")
            atoms = read(conf_in)
            self.update_omm_context(atoms)
            if atoms.info.get("reverse_vel", False):
                print("Generating velocities and reversing at the same time should not happen.")
            self.context.setVelocitiesToTemperature(self.temperature)
            state = self.context.getState(getVelocities=True)
            vel = state.getVelocities(asNumpy=True)/omm_units.angstrom*omm_units.femtosecond/ase_units.fs
            atoms.set_velocities(vel)
            atoms.write(conf_out)
