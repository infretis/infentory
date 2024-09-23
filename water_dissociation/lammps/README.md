<h1 align="center">
Water autoionization
</h1>

### TODO
* [x] reduce load path size
* [x] upload to infentory
* [ ] adapt inft plot_order
* [x] make scripts for visualization
* [x] fix lammpstrj processor
* [x] fix xyz processor?
* [ ] also adapt puckering exercise for cosy; exactly the same as here but with gromacs. Only production md run and path sampling, then visualize
* [x] move puckering and water dissociation README.md to infentory
* [ ] fix installation after moving to infentory
## Aloha 👋
In this session, we will study the autoionization of water using **path sampling**. The main outcomes of this simulation allow us to

* calculate 🖥️ exactly how often water dissociates into H3O+ and OH-
* and visualize 👀 how this chemical reaction actually happens

An essential ingredient of path sampling is using molecular dynamics (MD) to make the atoms and molecules wiggle around and react. As an introduction to the &infin;RETIS method you will therefore perform the following steps:

* 1️⃣ An MD simulation using [LAMMPS](https://www.lammps.org/#nogo)
* 2️⃣ A path sampling simulation on this system with [&infin;RETIS](https://github.com/infretis/infretis) + LAMMPS
* 3️⃣ See and learn how water dissociates at the molecular scale 🔎

## Step 0: Installation
Open a terminal 💻

If you don't already have conda or mamba:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b
```

Close and then re-open the terminal. You should now see **(base)** in the lower left of your screen.

Download python and lammps:

```bash
mamba create --name cosy_24 python==3.11 lammps
```

Install infretis, and the exercise files.
```bash
mamba activate cosy_24
python -m pip install git+https://github.com/infretis/infretis.git@cosy_24
python -m pip install git+https://github.com/infretis/inftools.git@main
git clone https://github.com/infretis/infentory.git
cd infentory/water_dissociation/lammps/
echo ========== We will perform the exercise from this folder ===============
```

## Step 1: MD with LAMMPS
Familiarize yourself with the files in the directory `lammps_input/`. Can you explain what these files contain?

Now, change to the `step1_md_run` directory and modify `lammp.input` to run an MD simulation at 300K for around 1 picosecond with a 0.5 fs timestep. We want to analyze some of the output, so write output with reasonable frequency, e.g. every 5 steps.

LAMMPS can be run with the command `lmp -i lammps.input`

📈 Does the system reach the desired temperature?

Animate the trajectory by opening the .dump file in Avogadro.

```bash
# Download Avogadro
wget https://github.com/OpenChemistry/avogadrolibs/releases/download/1.99.0/Avogadro2-x86_64.AppImage ~/
chmod +x ~/Avogadro2-x86_64.AppImage
# Run Avogadro
~/Avogadro2-x86_64.AppImage
```
Use the `animation tool` and click the `dynamic bonding` checkbox. You can rotate and zoom with the `navigation tool`. Do you see anything interesting? 🔎 

We will now investigate the value of the order parameter.

```bash
inft recalculate_order -toml infretis.toml -traj md_run.dump -format 'lammpsdump'
```
Plot the order parameter with gnuplot
```bash
gnuplot
plot 'order_rec.txt' with linespoint
```

Values up to around 1.5 mean we have only water present (the largest O-H bond length of all water molecules). Values greater than 1.5 tell us we have one OH- and one H3O+ present (the shortest distance between the OH- oxygen and H3O+ hydrogen).

<p align="center">
<img src="https://github.com/infretis/infentory/blob/main/water_dissociation/orderp.png" width="50%" height="50%">
</p>

Does the value of the order parameter during the simulation make sense with your conclusions from visualizing the trajectory? 


## Step 2: Path sampling with &infin;RETIS + LAMMPS
You now know how to run an MD simulation and calculate the order parameter. This is what &infin;RETIS does under the hood; a single Monte Carlo (MC) step with &infin;RETIS will run a LAMMPS simulation given some initial configuration and calculate the order parameter. If this trajectory meets the ensemble criterion we may accept and add it to our sampled states. If not we resample the old trajectory. So in path sampling, we combine both MC and MD in a hybrid approach. We need to do these steps repeatedly, which can take some time. Therefore we start the infretis simulation now. 

Navigate to the `step2_infretis` folder and fire off `infretisrun -i infretis.toml`.

At this point, reviewing the main outcomes of a path sampling simulation may be useful to remind yourself why we are doing this ⏮️

Then move on to the next step.

## Step 3: Visualization and analysis
Open a new terminal, run `mamba activate cosy_24`, and navigate to the same directory.

We may now start visualizing the results as infretis produces them. Of interest are paths (trajectories) with large order parameter values. These may be reactive and contain information on how the water deprotonation occurs 👀 The next task is therefore to identify a path with a large order parameter.

Open the `sim.log` and look for accepted MC moves by searching `'ACC'`. These lines give you the length of the path `len` and the min/max order parameter value `op: [min, max]`. Find a path with a large OP value (above 4.0). Identify the `new_path_nr` by looking at the line above for `old_path_nr -> new_path_nr`.

The path is stored in `load/new_path_nr`. Gnuplot the order parameter value `order.txt`. Do you see large jumps in the values🐇? What do you think they mean?

The jumps mean that a proton jumps from one water molecule to another. Therefore, to visualize the path nicely in Avogadro, we want to center the view on the oxygen the proton jumps away from to become OH-. Open `order.txt` and look at the 3rd column. The value of this column is the index of the oxygen in OH-. Take note of this number. 

Now, in the `load/path` folder, center the trajectory on the atom index you found:

```bash
inft trjcat -centersel "index 72"  -out traj.pdb -traj traj.txt -topology ../../../../cp2k/cp2k_data/initial.xyz -format lammpsdump 
```
but replace 72 with the number you found. 

Now, visualize `traj.pdb` in Avogadro. Do you see anything interesting?

### Optional: Dissociation rate and waiting time

We can also calculate how often this event happens after constructing the crossing probability curves and calculating the flux. This can be achieved by

```bash
inft wham -data infretis_data.txt -nskip 0
```
If you get an error you may not have enough data yet, or infretis wrote to another file `infretis_data_X.txt` where X is some number.

Gnuplot the crossing probability in `wham/Pcross.txt`. Use `set logscale y; set xrange [1:8]` before plotting to get a nicer view. This is the probability of reaching an order parameter value of $\lambda$ given that we start with pure water. Write down the lowest y-value, which we call $P_{tot}$.

Also, plot the `wham/runav_flux.txt` and write down the last value, which we call $f$. To get the rate $v$ in units of nanoseconds $^{-1}$, use the formula

$$v = 2'000'000 \cdot f \cdot P_{tot}$$

The interpretation is that we see a dissociation event every $1/v$ nanoseconds.

We can also calculate how many days we would have to wait to observe a single event in a regular MD simulation ⏳

Open the `log.lammps` file in `step1_md_run` and search for `'ns/day`. So, given that you can run X nanoseconds per day, and have to wait $1/v$ nanoseconds for an event. How many days would you have to wait to observe it in a simulation?

## 🏁 Further information

If you are interested in using &infin;RETIS in your work, feel free to contact the infretis team to help you get started 🤝

* titus.van.erp@ntnu.no
* anders.lervik@ntnu.no

If you are interested to learn more, this is a reproduction of the work studied in [this](http://www.pnas.org/cgi/doi/10.1073/pnas.1714070115) paper. It was also studied with infretis [here](https://doi.org/10.1073/pnas.2318731121). However, we used _reaxff_ while those papers used density function theory.
