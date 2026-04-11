This example must be run on the `external_ase` infretis GitHub branch.

This study uses the double-well settings (and reproduces the results) of
*Exchanging Replicas with Unequal Cost, Infinitely and Permanently J. Phys. Chem. A 2022, 126, 47, 8878–8886*


Generate initial paths:
```
inft generate_zero_paths -conf openmm_input/initial.traj -toml infretis0.toml &
python -m infretis.classes.engines.propagator infretis0.toml temporary_load

```

Run short inf-init:

```
python -m infretis.classes.engines.propagator infretis0.toml worker0 &
inft infinit -toml infretis0.toml
```

Remember to kill the background processes when finished (pkill python).
Remember to also remove the worker folders before restarting/continuing a simulation
