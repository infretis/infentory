# Introduction
This tutorial illustrates the use of the `infinit` functionality, which generates initial paths and optimizes the interfaces. All we need to start is an initial configuration and an orderparameter!

The process we will study is the pyramidal inversion of the NH3 molecule, but this tutorial can easily be adapted to a large number of different systems with minimal modifications.

# Required packages
Be sure to be on the latest versions of the main branches of `infretis` and `inftools`.

We use the XTB Hamiltonian to describe NH3, so we need the `xtb-python` package, which can be installed with conda

```bash
conda install xtb-python
```

# Step 0: The initial configuration and order parameter
Ideally, we would start `infinit` from a multitude of independent equlibrated initial configurations, but as of now, this option is not implemented yet to do this in an automated fashion. We start here from a single configuration

```python
from ase.build import molecule
atoms = molecule("NH3")
atoms.write("conf.traj")
```

The orderparameter we are using is just the dihedral angle between the 4 atoms.

```toml
[orderparameter]
class = "Dihedral"
index = [ 0, 3, 2, 1]
periodic = false
```

We give here the name `infretis0.toml` so that we have a backup of the toml, as infinit will create a multitude of `infretis.toml` and `infretis_X.toml` files, where X is a number.

# Step 2: The [infinit] section
In `infretis0.toml`, you should see the following in the [infinit] section.
```toml
[infinit]
cstep = -1
initial_conf = "conf.traj"
steps_per_iter = [ 40, 80, 150, 150]
pL = 0.3
skip = 0.05
lamres = 0.005
```

* `cstep` is the current infinit itreation we are on. `cstep = -1` lets infinit know that we do not have initial paths, and that a `load/` folder is absent. Infinit will therefore first generate a [0-] and a [0+] path from the initial configuration (under the hood it uses `inft generate_zero_paths` and then copies the [0+] path N worker times). If we had a `load/` folder containing some initial paths from e.g. a long MD simulation, we could pass that to infinit as well, but then having `cstep = 0`.
* `initial_conf` is the initial configuration we will generate the [0-] and [0+] paths from by propagating forwards and backwards until we hit the interface and have 1 valid path. Then, the last point is extended and another path created, giving a valid [0-] and [0+] path.
* `steps_per_iter` tells how many infretis steps we should run before updating the interfaces. In our case, `[40, 80, ...]` means __after__ generating the [0-] and [0+] path, we will run 40 infretis steps (cstep = 0), then update the interfaces, fill these with new initial paths from the previous simulation, and the do another infretis simulation with 80 steps (cstep = 1).
* `pL` is the local crossing probability between the interfaces infinit will place. So during the interface updates, new interfaces are placed based on the crossing probability estimate using all data from the previous infretis simulations. Often we would like pL = 0.3, but it could also be higher if we have available a large number of workers.
* `skip = 0.05` means that 5% of the first `infretis_data.txt` entries are not used in the estimation of the crossing probability, so the data of the first 5% paths are assumed to be discarded for equilibration purposes.
* `lamres = 0.005` means that after the interface estimation, the interfaces are rounded to a precision of 0.005. This is handy for later WHAM analysis of the data.

# Running infinit
We should now have everything set up to run the simulation, and you can run infinit with the following command.

```bash
export OMP_NUM_THREADS=1 # use only 1 OpenMP thread for this small system for XTB
inft infinit -toml infretis0.toml
```
The simulation should complete in approximately one minute.

# Restarting infinit or continuing the simulation
If the simulation crashes at any point, you can restart the simulation by runnining
```bash
inft infinit -toml infretis.toml
```
Alternatively, you can change or add steps to the `steps_per_iter` list in `infretis.toml` to add more steps.

Infinit should be able to figure out on its own where to pick up simulations. Infinit should also be able to figure out if the `restart.toml` is usable to restart the simulation.


# Output files
<details>
<summary>
:eyes: The output may give you some hints of what infinit is doing under the hood :eyes: </summary>
</summary>

  
conf.traj  
xtbcalc.py  
infretis0.toml  - _orignal .toml file, not changed or overwritten if not called infretis.toml_  
zero_paths.toml  - _.toml file that was used to generate the [0-] and [0+] paths_  
infretis_data.txt  - _empty data file after generating zero paths_  
**temporary_load** - _the [0-] and [0+] trajectories were generated in here_  
**run0** - _this was the first load/ folder, now renamed to run0_  
infretis_data_1.txt  - _first data file from the paths resent in run0/_  
combo_0.txt  - _a combined infretis_data.txt file with all data generated up til now, with 5% skipped (skip=0.05 in [infinit])_  
combo_0.toml  - _a combined .toml file, having all combined interfaces from all simulations til now_  
infretis_1.toml  - _.toml file that was used for the first infretis simulation (for paths in run0/)_  
**run1**  - _the directory containing paths of the second infretis simulation_  
infretis_data_2.txt  - _first data file from the paths resent in run1/_  
combo_1.txt  - _combined data from infretis_data_1.txt and infretis_data_2.txt, with 5% skipped from each file_  
combo_1.toml  - _combined interfaces from infretis_1.toml and infretis_2.toml_  
infretis_2.toml  - _.toml used to run the second infretis simulation_  
**run2**  - _paths from third infretis simulation_  
infretis_data_3.txt  - _data from third simulation_  
combo_2.toml  - _combined interfaces from sim 1, 2 and 3_  
combo_2.txt  - _combined data from sim 1, 2 and 3_  
infretis_3.toml  - _toml used for sim 3_  
worker1.log  
**worker0**  
worker0.log  
**run3**  
**worker1**  
infretis_data_4.txt  
sim.log  
restart.toml  
combo_3.toml  
combo_3.txt  
last_infretis_pcross.txt  - _estimate of crossing probability using all data that has been generated up til now, calculated after each infinit iteration_  
last_infretis_path_weigths.txt  - _path weights, not used atm_  
infretis_init.log  - _a basic logger containing some un-informative prints_  
infretis_4.toml  
infretis.toml  - _new infretis.toml with updated interfaces, ready to be used for production with infreisrun by changing `steps`, or continuing with infinit by adding to `steps_per_iter`_  
**load** - _current load/ folder, ready to be run with infretis.toml_  


</details>
