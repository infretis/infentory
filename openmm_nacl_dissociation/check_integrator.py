import openmm
import openmm.app
from sys import stdout
import openmm.unit as omm_unit

system_xml = "openmm_input/system.xml"
topology_pdb = "openmm_input/topology.pdb"


timestep = 2 * omm_unit.femtosecond
temperature = 300

# non-reversible (out of the box) integrator
integrator = openmm.LangevinIntegrator(temperature, 1 / omm_unit.picosecond, timestep)

# reversible integrator
# integrator = openmm.CustomIntegrator(timestep)
# integrator.addPerDofVariable("x1", 0)
# integrator.addUpdateContextState()
# integrator.addComputePerDof("v", "v+0.5*dt*f/m")
# integrator.addComputePerDof("x", "x+dt*v")
# integrator.addComputePerDof("x1", "x")
# integrator.addConstrainPositions()
# integrator.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
# integrator.addConstrainVelocities()

# set up the simulation
with open(system_xml) as rfile:
    system = openmm.XmlSerializer.deserialize(rfile.read())

pdb = openmm.app.PDBFile(topology_pdb)
topology = pdb.topology
simulation = openmm.app.Simulation(topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
print("Equilibrating...")
simulation.step(1000)
simulation.reporters.append(
    openmm.app.StateDataReporter(stdout, 1, step=True, temperature=True)
)

print(f"Setting velocities to {temperature} K")
simulation.context.setVelocitiesToTemperature(temperature)


# equilibration
print("Forward propagation:")
simulation.step(1)
simulation.step(1)
simulation.step(1)
simulation.step(1)
state = simulation.context.getState(getVelocities=True)
vel = state.getVelocities(asNumpy=True)
simulation.context.setVelocities(-vel)
print("Backward propagation:")
simulation.step(1)
simulation.step(1)
simulation.step(1)
