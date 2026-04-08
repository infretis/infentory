import openmm
import openmm.app
import openmm.unit as omm_units
import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase import units as ase_units


class OpenMMCalculator(Calculator):
    """
    Minimalistic *external* OpenMM calculator.

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
        self.timestep = 2 * omm_units.femtoseconds
        # the temperature has to be the same as in the .toml
        self.temperature = 300*omm_units*kelvin
        slef.friction_coef = 1.0/omm_units.picoseconds
        print("Temperature in openmmcalculator.py: T = ", self.temperature, "K")

        #   set up the velocity-verlet integrator from:
        #     https://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.openmm.CustomIntegrator.html
        #   which is reversible wrt. the momenta. The default VelocityVerlet would not work out of the box
        #   when reversing velocities to propagate backwards in time
        #   integrator = openmm.CustomIntegrator(self.timestep)
        #   integrator.addPerDofVariable("x1", 0)
        #   integrator.addUpdateContextState()
        #   integrator.addComputePerDof("v", "v+0.5*dt*f/m")
        #   integrator.addComputePerDof("x", "x+dt*v")
        #   integrator.addComputePerDof("x1", "x")
        #   integrator.addConstrainPositions()
        #   integrator.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        #   integrator.addConstrainVelocities()

        # the LangevinIntegrator stores velocites on half-steps
        self.half_step_velocities = True
        integrator = openmm.LangevinIntegrator(self.temperature, self.friction_coef, self.emperature, self.friction_coef, self.timestep)

        # setup system, topology and create the openmm simulation object
        with open(system_xml) as rfile:
            system = openmm.XmlSerializer.deserialize(rfile.read())
        pdb = openmm.app.PDBFile(topology_pdb)
        topology = pdb.topology
        self.simulation = openmm.app.Simulation(topology, system, integrator)

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
        else:
            self.simulation.step(self.subcycles)

        state = self.context.getState(
            getEnergy=True, getPositions=True, getVelocities=True
        )

        # set positions and velocities in ASE units
        atoms.positions = state.getPositions(asNumpy=True) / omm_units.angstrom
        atoms.set_velocities(
            state.getVelocities(asNumpy=True)
            / omm_units.angstrom
            * omm_units.femtosecond
            * ase_units.fs
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

    def _reverse_velocities(self, filename: str, outfile: str) -> None:
        """Reverse velcoities using the ASE method if we store on-step velocities,
        else reverse the half-step velocities of the next step by performing a single
        integration step.
        """
        if not self.half_step_velocities:
            self._reverse_velocities(filename, outfile)
        else:
            atoms = read(filename)
            if isinstance(atoms, list):
                atoms = atoms[0]
            self.update_omm_context(atoms)
            # perform 1 MD step
            self.simulation.step(1)
            state = self.context.getState(getVelocities=True)
            # update atoms with half-step velocities and reverse
            vel = state.getVelocities(asNumpy=True)
                / omm_units.angstrom
                * omm_units.femtosecond
                * ase_units.fs
            atoms.set_velocities(-vel)
            write(outfile, atoms)

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
        atoms.info["not_setup"] = False
