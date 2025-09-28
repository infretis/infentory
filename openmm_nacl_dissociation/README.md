# Simple NaCl dissociation example (WIP)
Here we illustrate the flexibility of the ASE engine, where we essentially only use ASE to handle the storage and modification of phasepoints (writing trajectories, extracting trajectory frames, randomizing velocities and reversing velocities). Everything else is handled by the ASE calculator.

Currently, this is only possible on the `external_ase` branch of infretis.

Note that this does not equilibrate the system, this tutorial merely illustrates how to use OpenMM in combination with infretis.

Instead of using an ASE calculator and send positions and forces back and forth between ASE and OpenMM runs `subcycles` steps, and then simply passes positions, velocities and box information to ASE, which is used to calculate the orderparameter and write trajectories in .traj format.

To use openmm, one needs a `system.xml` and `topology.pdb` file, which is used to set up the openmm simulation. See `openmmcalculator.py` for how this is done. Here, one can change `subcycles`, `timestep` and all other system settings. The files `openmm_input/` were generated with the `openmm_input/create_system.py` script.

Constraints are treated by first generating velocites with ASE wihout constraints. The velocity components of the constraints are then enforced by openmm before propagating the dynamics.

Note that OpenMM does not wrap coordinates back into the box. For a permeation simulation where the membrane stays in place, this can be advantegous as one can use the absolute z position as the orderparameter. One can wrap coordinates by adding `atoms.wrap()` somewhere in the `infretis/classes/engines/propagator.py` script, or somewhere in the openmm calculator.

# A note on integrators
The default integrators in OpenMM seemingly all output velocities at half-steps. As such, one can not use these out of the box. To see if your integrator is viable, set up a simple simulation with an openmm reporter that reports e.g. the temperature every timestep.

Then run some single MD steps and inspect the temperture output. Then reverse the velocites and continue the propagation, which should now run in reverse. The temperature should now be the same as the previous steps. Running the `check_integrator.py` script, the output using the CustomIntegrator is:

```
Equilibrating...
Setting velocities to 300 K
Forward propagation:
#"Step","Temperature (K)"
1001,301.9886790207829
1002,289.81835670673865
1003,275.13610867320745
1004,260.9088850729737
Backward propagation:
1005,275.13610922108796
1006,289.81835977189746
1007,301.98868302614306
```
where we see that the tempearture at step 1005 is the same as step 1003, step 1006 the same as 1002, etc.

For the LangevinIntegrator we get:
```
Equilibrating...
Setting velocities to 300 K
Forward propagation:
#"Step","Temperature (K)"
1001,293.25347385263774
1002,288.497064447396
1003,284.3419835632597
1004,279.0789574478106
Backward propagation:
1005,305.33070545807885
1006,318.9661593495396
1007,326.3560979474691
```
where we see large temperature changes in the system.


# infretis.toml
To use the external ase\_calculator you have to set `engine = "ase_external"`.

The external ase engine can be used in two modes. We will use a custom mode where OpenMM does all the hard work - like updating the positions and velocities every time step. As such, we have to set `subcycles = 0`, `timestep = 0` and `integrator = "external"`.

The timestep and the subcycles is then set to the desired value in openmmcalculator.py.

Note that the temperature has to still be set in the .toml, and also in the openmm calculator thermostat.

# intial paths

`cstep=-1` does not work with infinit as of now due to some folder issues. This can be sidestepped with:
```bash
export OPENMM_CPU_THREADS=1
inft generate_zero_paths -conf openmm_input/initial.traj -toml infretis0.toml &
python -m infretis.classes.engines.propagator infretis0.toml temporary_load
```

Then kill the propagator process, and run
```
python -m infretis.classes.engines.propagator infretis0.toml worker0 &
inft infinit -toml infretis0.toml &
```

To use srun and multiple workers, something like this might be used:
```
for i in {0..7}; do
	srun --nodes 1 --ntasks 1 --cores 7 --gpus=1 python -m infretis.classes.engines.propagator infretis0.toml worker$i &
done
inf infinit -toml infretis0.toml
# in some cases the master process above has to also be launched with srun as well:
#	srun --ntasks=1 --cpus-per-task=2 --gpus=0 inft infinit -toml infretis0.toml &
```
