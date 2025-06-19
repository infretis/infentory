# Introduction
This tutorial illustrates the use of the `infinit` functionality, which generates initial paths and places the interfaces in an optimal manner, starting out from just a single configuration! The process we will study is the pyramidal inversion of the NH3 molecule

This tutorial can easily be adapted to a large number of different systems with minimal modifications.

# Required packages
Be sure to be on the latest versions of the main branches of `infretis` and `inftools`.

We use the xtb hamiltonian to describe NH3, so we need the `xtb-python` package, which can be installed with conda

```bash
conda install xtb-python
```

# Step 0: Generating the intial configuration
Ideally, we would start `infinit` from a multitude of independent equlibrated initial configurations, but as of now, this option is not implemented yet to do this in an automated fashion. We here start from a single configuration

```python
from ase.build import molecule
atoms = molecule("NH3")
atoms.write("conf.traj")
```

# Step 1: Setting up the .toml
Usually when setting up a new system, we copy the infretis.toml from one of the example systems of the engines since few of the settings there change.

The engine section:

```toml
[engine]
class = "ase"
engine = "ase"
temperature = 300
input_path = "."
timestep = 0.5
subcycles = 2
integrator = "langevin"
langevin_fixcm = true
langevin_friction = 0.005

[engine.calculator_settings]
module = "xtbcalc.py"
class = "XTBCalculator"
```
We use a timestep of 0.5 fs with the langevin integrator. With the ASE engine, we need to define a calculator in an external python file, which you can inspect in `xtbcalc.py`.


We here the dihedral angle between the 4 atoms of NH3 to describe the progress of the reaction, this is already defined in the file `infretis0.toml`.

```toml
[orderparameter]
class = "Dihedral"
index = [0,1,2,3]
```

We give here the name `infretis0.toml` so that we have a backup of the toml, as infinit will create a multitude of `infretis.toml` and `infretis_X.toml` files, where X is a number.

# Step 2: The [infinit] section
In `infretis0.toml` you should see the following in the [infinit] section.
```toml
[infinit]
cstep = -1
initial_conf = "conf.traj"
steps_per_iter = [ 40, 80, 150, 150],
pL = 0.3
skip = 0.05
lamres = 0.005
```

* `cstep` is the current infinit itreation we are on. `cstep = -1` lets infinit know that we do not have initial paths, and that a `load/` folder is absent. Infinit will therefore first generate a [0-] and a [0+] path from the initial configuration (under the hood it uses `inft generate_zero_paths` and then copyies the [0+] path N worker times). If we had a `load/` folder containing some initial paths from e.g. a long MD simulation, we could pass that to infinit as well, but then having `cstep = 0`.
* `initial_conf` is the initial configuration we will generate the [0-] and [0+] paths from by propagating forwards and backwards until we hit the interface and have 1 valid path. Then, the last point is extended and another path created, giving a valid [0-] and [0+] path.
* `steps_per_iter` tells how many infretis steps we should run before updating the interfaces. In our case, `[40, 80, ...]` means __after__ generating the [0-] and [0+] path, we will run 40 infretis steps (cstep = 0), then update the interfaces, fill these with new initial paths from the previous simulation, and the do another infretis simulation with 80 steps (cstep = 1).
* `pL` is the local crossing probability between the interfaces infinit will place. So during the interface updates, new interfaces are placed based on the crossing probability estimate using all data from the previous infretis simulations. Often we would like pL = 0.3, but it could also be higher if we have available a large number of workers.
* `skip = 0.05` means that 5% of the first `infretis_data.txt` entries are not used in the estimation of the crossing probability, so the data of the first 5% paths are assumed to be discarded for equilibration purposes.
* `lamres = 0.005` means that after the interface estimation, the interfaces are rounded to a precision of 0.005. This is handy for later WHAM analaysis of the data.

# Running infinit
We should now have everything set up to run the simulation

```bash
inft infinit -toml infretis0.toml
```

# Restarting infinit or continuing the simulation
