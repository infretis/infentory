# Introduction
This tutorial illustrates the use of the `infinit` functionality, which generates initial paths and optimizes the interfaces.

All we need to start is an initial configuration and an order parameter!

The process we will study is the pyramidal inversion of the NH3 molecule, but this tutorial can easily be adapted to a large number of different systems with minimal modifications.

# Required packages
<details>


Be sure to be on the latest versions of the main branches of `infretis` and `inftools`.

We use the XTB Hamiltonian to describe NH3, so we need the `xtb-python` package, which can be installed with conda

```bash
conda install xtb-python
```

</details>

# The initial configuration and order parameter
<details>

Ideally, we would start `infinit` from a multitude of independent equilibrated initial configurations, but as of now, this option is not implemented in an automated fashion yet. We start here from a single configuration

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

</details>

# The [infinit] section
<details>

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

* `cstep` is the current infinit itreation we are on. `cstep = -1` lets infinit know that we do not have initial paths, and that a `load/` folder is absent. Infinit will therefore first generate a [0-] and a [0+] path from the initial configuration (under the hood, it uses `inft generate_zero_paths` and then copies the [0+] path N worker times). If we had a `load/` folder containing some initial paths from e.g. a long MD simulation, we could pass that to infinit as well, but then having `cstep = 0`.
* `initial_conf` is the initial configuration we will generate the [0-] and [0+] paths from by propagating forwards and backwards until we hit the interface and have 1 valid path. Then, the last point is extended and another path created, giving a valid [0-] and [0+] path.
* `steps_per_iter` tells how many infretis steps we should run before updating the interfaces. In our case, `[40, 80, ...]` means __after__ generating the [0-] and [0+] path, we will run 40 infretis steps (cstep = 0), then update the interfaces, fill these with new initial paths from the previous simulation, and the do another infretis simulation with 80 steps (cstep = 1).
* `pL` is the local crossing probability between the interfaces infinit will place. So during the interface updates, new interfaces are placed based on the crossing probability estimate using all data from the previous infretis simulations. Often we would like pL = 0.3, but it could also be higher if we have available a large number of workers.
* `skip = 0.05` means that 5% of the first `infretis_data.txt` entries are not used in the estimation of the crossing probability, so the data of the first 5% paths are assumed to be discarded for equilibration purposes.
* `lamres = 0.005` means that after the interface estimation, the interfaces are rounded to a precision of 0.005. This is handy for later WHAM analysis of the data.

</details>


# Running infinit
<details>


We should now have everything set up to run the simulation, and you can run infinit with the following command.

```bash
export OMP_NUM_THREADS=1 # use only 1 OpenMP thread for this small system for XTB
inft infinit -toml infretis0.toml
```
The simulation should complete in approximately one minute.

</details>

# Restarting an interrupted infinit simulation
<details>

If the simulation crashes at any point, you can restart the simulation by running
```bash
inft infinit -toml infretis.toml
```
Alternatively, you can change or add steps to the `steps_per_iter` list in `infretis.toml` to add more steps. This is illustrated further down.

Infinit should be able to figure out on its own where to pick up simulations. Infinit should also be able to figure out if the `restart.toml` is usable to restart the simulation.

</details>

# Output files created by infinit
<details>

The output may give you some hints of what infinit is doing under the hood
  
* infretis0.toml  - _original .toml file, not changed or overwritten if not called infretis.toml_  
* zero_paths.toml  - _.toml file that was used to generate the [0-] and [0+] paths_  
* infretis_data.txt  - _empty data file after generating zero paths_  
* **temporary_load** - _the [0-] and [0+] trajectories were generated in here_  
* **run0** - _this was the first load/ folder, now renamed to run0_  
* infretis_data_1.txt  - _first data file from the paths resent in run0/_  
* combo_0.txt  - _a combined infretis_data.txt file with all data generated up til now, with 5% skipped (skip=0.05 in [infinit])_  
* combo_0.toml  - _a combined .toml file, having all combined interfaces from all simulations til now_  
* infretis_1.toml  - _.toml file that was used for the first infretis simulation (for paths in run0/)_  
* **run1**  - _the directory containing paths of the second infretis simulation_  
* infretis_data_2.txt  - _first data file from the paths resent in run1/_  
* combo_1.txt  - _combined data from infretis_data_1.txt and infretis_data_2.txt, with 5% skipped from each file_  
* combo_1.toml  - _combined interfaces from infretis_1.toml and infretis_2.toml_  
* infretis_2.toml  - _.toml used to run the second infretis simulation_  
* ...
* last_infretis_pcross.txt  - _estimate of crossing probability using all data that has been generated up til now, calculated after each infinit iteration_  
* last_infretis_path_weigths.txt  - _path weights, not used atm_  
* infretis_init.log  - _a basic logger containing some uninformative prints_  
* infretis_4.toml  
* infretis.toml  - _new infretis.toml with updated interfaces, ready to be used for production with infreisrun by changing `steps`, or continuing with infinit by adding to `steps_per_iter`_  
* **load** - _current load/ folder, ready to be run with infretis.toml_  

</details>

# Analysis: Are the interfaces reasonable?
<details>

 You can plot the order parameter of the previous simulations with the previous interfaces:
 ```bash
inft plot_order -traj run3 -toml infretis_4.toml
```
or the current paths and interfaces estimated til now
```bash
inft plot_order -traj load -toml infretis.toml
```
![tmp](https://github.com/user-attachments/assets/e3a5b5bc-ad16-4530-ba90-ff65c67fd5c3)

We see that there are reactive paths in the current **load/** folder!

We also see that the interfaces seem smoothly spaced and placed, but there are some irregularities (the distance between interface 4 and 5). Most converged crossing probability curves change slowly (they are quite smooth), so we also expect a very regular interval between interface locations. We will now investigate further whether the interfaces are placed well enough or if we need more simulations.


To do this, we WHAM all of the combined data up til now, meaning using the latest `combo.txt` files (which contain the combined infretis data) and the `combo.toml`, containing the combined interfaces.


```bash
inft wham -data combo_3.txt -toml combo_3.toml -nskip 0 -lamres 0.005 -folder wham_combo
```
* `nskip = 0` because the lines are already trimmed in the combo.txt files wrt skip in the [infinit] section
* `lamres` should be the same or less than specified in the [infinit] section.

The crossing probability `wham_combo/Pcross.txt` from the WHAM analysis with the most recent interfaces (infretis.tom) is shown below:

![tmp2](https://github.com/user-attachments/assets/a424b75e-2d53-4369-ae87-0eaf7a8398d0)

We see that the crossing probability looks smooth-ish, but there are some blocky segments. So what do we do now - should we run a long infretis simulation with those interfaces, or should we run some more steps with infinit to get the probability? Of course, this also depends on how expensive the simulations are, but adding more steps with infinit might be the wiser choice, as the data either way can be used in the final rate estimate.

It is possible that the estimated interfaces are not placed optimally, and if we run a single long infertis simulation, the efficiency might be suboptimal. There might be more to gain if we add one or more steps of infinit.


</details>

# Continuing the simulation with more steps
<details>

Here, we decide on adding one additional shorter infinit step, and then one longer step, which we take as our production run.

Open the file `infretis.toml`. You should see something like

```toml
[infinit]
cstep = 4
initial_conf = "conf.traj"
steps_per_iter = [
    40,
    80,
    150,
    150,
]
```
We see that `cstep = 4`, but the 4th element `steps_per_iter` does not exist (counting from 0). Add two more infinit steps by changing the keyword to `steps_per_iter = [40, 80, 150, 150, 250, 750]`.

Now, continue the infinit simulation by running
```bash
export OMP_NUM_THREADS=1 # XTB related
inft infinit -toml infretis.toml
```
which continues the infinit loop from `cstep = 4`.
</details>

# Analysis part 2
<details>
Depending on how much data you can generate in the final production run, you might also want to exclude the infinit simulations as equilibration. We can choose only to analyze the data from the final infretis run using

```bash
mv load new_load # rename the newly created load folder temporarily
mv run5 load # required so that 'inft wham' finds the order.txt files
inft wham -data infretis_data_6.txt -toml infretis_6.toml -nskip 75 -lamres 0.005 -fener -xcol 1 -nbx 30 -minx -1 -maxx 1
```

If you want to increase the accuracy of the rate estimate, you might want to continue the simulation. Below we plot `Pcross.txt`  from the combined data, in addition to the previous and current interfaces (interfaces from `infretis_6.toml` and `infretis.toml`):

![tmp3](https://github.com/user-attachments/assets/94f6d71d-2272-41eb-92d0-2dc9e1da4a65)

We see only a slight change in the last 3 interfaces. The spacing between the latest interfaces (red curve, but also in the blue curve) seems regular, and the black crossing proability curve looks smooth

🏁🏁🏁

Since the previous interfaces (blue lines) didn't change much as compared to the updated interfaces (red lines), one could also just have continued from the last 750 steps with the `restart.toml` and regular infretis. We will showcase this here as well. 

In the `restart.toml`, change the `steps` keyword from 750 to e.g. 2500.

We already renamed the run5/ folder back to load/, so we can now simply run

```bash
infretisrun -i restart.toml
```
and infretis continues the simulation from step 750 until 2500.

</details>
