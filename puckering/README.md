<h1 align="center">
The Ring Flip Enigma:

Unveiling Molecular Secrets with Path Sampling
</h1>

# Motivation
The motivation for this assignment is to introduce students to software designed for conducting path sampling. The software employed here is an in-house-developed Python code for running &infin;RETIS that interfaces among others with Gromacs to perform the molecular dynamics (MD) steps. Through this exercise, we aim to demonstrate the capabilities of path sampling, where we can study transition processes that can be hard or impossible to investigate using conventional brute-force MD methods due to their rare event nature. The algorithms and software utilized in this assignment are the result of very recent active developments within the research group of Theoretical Chemistry.

If path sampling and software development sound interesting to you, to the extent that you would like to study them in more detail, please don't hesitate to get in touch with Titus and Anders to explore potential master projects. You can contact them at titus.van.erp@ntnu.no and anders.lervik@ntnu.no.

# Goals
In this exercise, you'll journey into the heart of molecular mysteries. Your primary goal is to gain hands-on experience by simulating the [ring flip](https://en.wikipedia.org/wiki/Ring_flip), an intriguing phenomenon often referred to as puckering. This transition, a rare occurrence at the molecular timescale, has puzzled scientists for ages. With your newfound knowledge in path sampling, you now hold the key to understanding its mechanisms. Your quest? To reveal the secrets hidden within the molecular world.

# The system
<p align="center">
<img src="https://github.com/infretis/infretis/blob/molmod_exercise5/examples/gromacs/puckering/graphics/puckering.gif" width="30%" height="30%">
</p>

This transition occurs very rarely at the molecular time scale, making it extremely challenging to study with standard molecular dynamics simulations. On the macroscale, these systems are awfully small and the transition happens exceedingly fast, making it almost impossible to study experimentally. Truly, this process remains hidden within the world of molecules! However, <ins>we would like to know exactly how often this transition occurs and the mechanism behind it </ins>. We can obtain this information by performing a path-sampling simulation.


<details>
<summary> ... even more about the system 🤓 </summary>

## Even more about the system

6-rings play a vital role in the world of chemistry and biology, impacting systems as diverse as carbohydrates being broken down by enzymes within your very body. The physical and chemical properties of 6-rings are intimately linked to their shapes, and their conformational landscape is a puzzle to be unraveled, with **C**hair, **H**alf-chair, **B**oat, **S**kew-boat, and **E**nvelope conformations. The conformations of 6-rings can be projected onto the surface of a sphere, where each conformer is uniquely specified by the angles $\theta$ and $\phi$.

<img src="http://enzyme13.bt.a.u-tokyo.ac.jp/CP/sugarconf.png" width="90%" height="90%">

These angles should not be viewed as regular angles between atoms, but rather as a coordinate transformation of the atoms that can be mapped onto the surface of a sphere [[1](https://doi.org/10.1021/ja00839a011)]. But the "hows" aren't important right now. The essential thing you need to know for now is that there is a high energy barrier between the north pole and the equator, and again between the equator and the south pole. We will study the transition over the first barrier; _starting at the north pole and ending at any of the structures on the equator_. By the end of this exercise, you will be able to say exactly how often this transition happens, and the mechanism behind it.

### Can you answer these?
* Given that the 6-ring in the animation above starts as $^4\text{C}_1$, can you see that the ending structure is $^{3,O}B$? Hint: The super- and subscripts refer to which atoms are above and below the mean plane of the ring, respectively.

* What is the initial value of the angle $\theta$, and what are the final values of the angles $\phi$ and $\theta$?

* Can you suggest an order parameter for this transition?

</details>



# Step 0: Installing the required packages
We first need to install the required programs to run this exercise. This includes a program that generates the parameters of a modern force field ([OpenFF 2.1](https://openforcefield.org/force-fields/force-fields/)) for your molecule and the ∞RETIS software developed at the theoretical chemistry group at NTNU.

Download and install mamba with the following commands (if you don't already have conda installed). Click the copy button on the box below and paste it into a terminal, and then do what is asked in the output on your screen (on Ubuntu, pressing down the mouse-wheel-button often works better for pasting than ctrl+V).
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

```
Now close the terminal.

If everything went successfully, you should see `(base)` in the left of your terminal window after reopening.

Then download and install the required python packages to run this exercise. Again copy-paste the code and do what is asked of you in the output.
```bash
mamba create --name molmod python==3.11
```
```bash
mamba activate molmod
mkdir software
cd software
git clone https://github.com/infretis/infretis.git
cd infretis
python -m pip install -e .
cd -
git clone https://github.com/infretis/inftools.git
cd inftools
python -m pip install -e .
cd ~
git clone https://github.com/infretis/infentory.git
cd infentory/puckering
echo "All done! We will perform the exercise from this folder."
```

You should now see `(molmod)` in the left of your terminal. Whenever you open a new terminal, write `mamba activate molmod` to activate the required Python packages.

We will perform the exercise from the directory `~/infentory/puckering/`. Get an overview of the folder structures from the terminal.

# Step 1: Equilibration
Run the following command:
```bash
cd gromacs_input
gmx solvate -cs spc216.gro -cp mol.gro -p topol.top -o conf.g96
```
#### 🤔 Question 1:
* What does the above command do?


Navigate to the `step1_equilibration` folder and get an overview of the directory structure. Perform an energy minimization (EM) and an NVT and NPT equilibration in the provided directories. Here are some commands to speed up the process.


```bash
gmx grompp -f em.mdp -p ../../gromacs_input/topol.top -c ../../gromacs_input/conf.g96 -o em.tpr
gmx mdrun -deffnm em -ntomp 2 -ntmpi 1 -pin on -v
```
```bash
gmx grompp -f nvt.mdp -p ../../gromacs_input/topol.top -c ../em/em.gro -o nvt.tpr
gmx mdrun -deffnm nvt -ntomp 2 -ntmpi 1 -pin on -v

```
```bash
gmx grompp -f npt.mdp -p ../../gromacs_input/topol.top -c ../nvt/nvt.gro -t ../nvt/nvt.cpt -o npt.tpr
gmx mdrun -deffnm npt -ntomp 2 -ntmpi 1 -pin on -v

```
#### 🤔 Question 2:
* Has the temperature and density reached the expected values during the NPT equilibration? The properties are accessible using `gmx energy -f npt.edr`. (Hint: retaw yltsom si metsys ruoY. Hint2: The letters of the previous hint are reversed to avoid spoilers.)

# Step 2: MD run
We have now equilibrated our system, and are now going to perform a slightly longer MD run. Navigate to the `step2_md_run` folder and run an MD run with the NPT equilibrated structure.

Use this command for the mdrun:
```bash
gmx grompp -f md.mdp -p ../gromacs_input/topol.top -c ../step1_equilibration/npt/npt.gro -t ../step1_equilibration/npt/npt.cpt -o md-run.tpr
gmx mdrun -deffnm md-run -ntomp 2 -ntmpi 1 -pin on -v -c confout.g96
```

This run should take a couple of minutes. You can use the time to answer the following question.

#### 🤔 Question 3:
* What is an order parameter in path sampling, and why do we need it?

Calculate the order parameter for each frame in the trajectory by using:
```bash
inft recalculate_order -traj md.trr -toml infretis.toml -out md-order.txt

```
Plot the order parameter values (column 1) vs time (column 0) from the MD run using e.g. gnuplot.

#### 🤔 Question 4:
* Given that the product state of your molecule is defined by $\lambda=90$, are you optimistic that you could observe a spontaneous transition during a plain MD simulation?

It is always a good idea to visualize trajectories to ensure everything is running as expected, and that our molecules haven't blown up 💥

We will use the popular visual molecular dynamics ([VMD](https://www.ks.uiuc.edu/Research/vmd/)) software for this:

```bash
vmd md-run.trr md-run.gro -e ../graphics/vmd-script.tcl
```

#### 🤔 Question 4 - 5:
* Do you see any interesting conformational changes when visualizing the trajectory?
* How can path sampling help us here?

# Step 3: ∞RETIS
In this section, we will finally perform the path simulation. However, before we can do that, we need to provide the ∞RETIS program with a set of interfaces and an initial path in each of the path ensembles defined by the interfaces. We can use the ∞RETIS initial path generator `infinit` for this. The way it works is illustrated below.

<img src="https://github.com/infretis/infretis/blob/molmod_exercise5/examples/gromacs/puckering/graphics/initial-paths.gif" width="45%" height="45%">


Navigate to the `step3_infretis` directory. The file `infretis.toml` defines all the path sampling setup. 

In the [simulation] section, define the initial state and final state by specifying two interfaces at $\lambda=10$ and $\lambda=90$ in `infretis.toml`.



# Step 4: Analysis
The following analysis is performed within the `step3_infretis` folder.
## The transition mechanism
We can say something about the mechanism of the complete $^4\text{C}_1 \rightarrow ^1\text{C}_4$ transition of your molecule if we assume that the second barrier from the equator to the south pole is negligible. The final configuration of your reactive paths would then be the transition state of the whole $^4\text{C}_1 \rightarrow ^1\text{C}_4$ transition. This may be a crude approximation, and we could test it by running another path simulation.

Plot the $\phi$ vs. $\theta$ values of the trajectories using the `-xy 2 1` option in `plot_order`. Looking at the reactive trajectories, what is/are the preferred route(s) from $^4\text{C}_1$ to $^1\text{C}_4$?

If you want, you can confirm this by visualizing some of the reactive trajectories. The following command removes the solvent, centers your molecule, and reorders the trajectories output from ∞RETIS:

```bash
# replace 'nr' with the path number of some trajectory you want to visualize
nr=46
inft concatenate -path load/${nr} -tpr ../gromacs_input/topol.tpr -out path${nr}.xyz
```
Now you get a file `path${nr}.xyz`that you can visualize in Avogadro.

## The transition rate

When you approach a reasonable number of paths in your simulation you can start analyzing the output. The following script calculates the rate, along with some other properties such as the crossing probability and error estimates.

```bash
inft wham -toml infretis.toml -data infretis_data.txt
```
The running average of the rate is written to the `runav_rate.txt` file, with the value in the fourth column giving the best estimate for the rate.
You can plot it in `gnuplot`

```bash
# in gnuplot
set logscale y
plot 'runav_rate.txt' using 1:4 with linespoints title 'rate'
```

The last line/point in this file is the estimated transition rate using all paths. To get this into units of $\text{ps}^{-1}$, divide the rate by $c$ where

$$c=\text{subcycles}\cdot \text{timestep}$$

which is found in the `infretis.toml` file.

Other files you may want to plot are the `Pcross.txt` for the crossing probability as a function of $\theta$, the `runav_flux` and `runav_Pcross.txt` for the running average of the flux and the crossing probability, and the `errRATE.txt`, `errFLUX.txt`, and `errPtot.txt` files for estimates of the relative error in the corresponding properties.

## Questions
* **10:** What is/are the preferred transition structures of your molecule on the equator?
* **11:** What is the rate in units of $\text{ns}^{-1}$?
* **12:** What is the interpretation of the inverse of the rate (1/rate)? (Hint: noitisnart rep emit ni era stinu ehT).
* **13:** Inspect the last part of the `md.log` file from `step2_md_run` and write down the Performance in ns/day. This number says how many nanoseconds of simulation you generate in one day on your machine. From the value of the inverse rate, how many days would you have to wait to observe a single transition in a standard MD simulation?


# How to pass this exercise
Answer all of the 13 questions and show/discuss them with the teaching assistants.
