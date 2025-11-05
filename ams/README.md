  <h1 align="center">
Unveiling Molecular Secrets with Path Sampling in AMS
</h1>

Through this tutorial, we aim to demonstrate the capabilities of path sampling. It allows us to study rare events that can be hard or impossible to investigate using conventional molecular dynamics methods. The algorithms and software utilized in this assignment are the result of very recent active developments. 
The following points shall guide your way from an initial structure towards running a RETIS simulation and what tools are there to make the most out of your simulation data. Contributions, questions and suggestions are always welcome!

The prerequisites to follow the tutorial are a recent version of AMS2025 or newer, including PLAMS. 

We will guide you step by step using the example of NaCl dissociation in water. This simple test case was chosen to showcase the latest improvements in our approach, which now enables us to extract both rate constants and mechanistic insights from a RETIS simulation (see [J. Chem. Theory Comput. 2023, 19, 7, 2222–2236](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00054) for details). The subfolders include all necessary scripts and files to run the simulations. All other systems can be set up and analyzed in the same way. 

# Step 0: Molecular Dynamics

The purpose of running MD as a preparation step is not only, to properly equilibrate your system to get an initial configuration to start your simulations from! We also need to learn about our systems state A, to place the first interface λ<sub>A</sub> in the most efficient way. Therefore, we run a long MD run and calculate the order parameter λ along the path. The optimal placement is now roughly, where 20% of the order parameter values are above.
In the [MD](./MD/) folder you can find the AMS input scripts to run those simulations. In this simple example, you can do this by visualizing the trajectory and use the 'Distance' tool of AMSMovie, but you can also recalculate it with inftools, a package that includes many helper functions for infRETIS. The required toml file is located in [2_explore_A](./0_MD/2_explore_A/).

# Step 1: Setup - Finding initial trajectories and interface placement

This step is as crucial, as it can be cumbersome. InfRETIS needs at least one path, physically meaningful or not, that connects state A and state B. Those paths can be obtained through various methods. Two convenient ones are presented in the following and only require one Geometry and the order parameter description:    
  1. InfInit, a tool contained in the [inftools](https://github.com/infretis/inftools) package and InfRETIS native. 
  2. PLUMED, as contained in AMS

Their comparison in an overview:

|           | InfInit                                                                                                                                                  | PLUMED                                                                                                                                       |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| Principle | <ul><li> Unbiased MD </li> <li> Multiple short RETIS runs</li> <li> Interfaces are updated at each step based on local crossing probabilities </li></ul> | <ul><li> Biased MD </li> <li> Exploration of free energy surface</li> <li> Interface placement estimated based on energy gradient </li></ul> |
| Pros      | + Same order parameter as production run <br> + Precise interface placement <br> + Minimal supervision                                                   | + Included in AMS <br> + Faster, if you know what you are doing <br> + Free energy surface obtained                                          |
| Cons      | - Separate python environment for InfTools <br> - early stages of development                                                                            | - Interface placement depends on convergence <br> - No python for order parameter <br> - Large number of RETIS steps to 'forget' the bias    |

Find examples [here](./setup/)


## InfInit
InfInit is part of the inftools package, a toolbox of different features to visualize and analyse a RETIS run. It is actively developed and not all of its functions are accessible with AMS and rkf-trajectories. We are always happy about contributions of new features! 

<details>

  <summary>Click to expand</summary>

### (Very detailed) Installation
<details>

  <summary>Click to expand</summary>
Download and install mamba with the following commands (if you don't already have conda installed).
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
Press enter and type yes every time when prompted. 
You should see `(base)` in the left of your terminal window after reopening if everything went successfully.

Then download and install the required python packages.

```bash
mamba create --name molmod python==3.11 MDAnalysis tqdm
```
```bash
mamba activate molmod
cd ~/opt
git clone https://github.com/infretis/infretis.git
cd infretis
python -m pip install -e .
cd ../
git clone https://github.com/infretis/inftools.git
```

#### Install PLAMS
AMS and ∞RETIS run in different python versions. We need to apply some tricks to work with both at the same time. First we need an external version of [PLAMS](https://www.scm.com/doc/Scripting/PLAMS/PLAMS.html), which is the AMS-internal python package. We install that by copying following lines to our terminal. Make sure, that the molmod environment is active.

```bash
cd ~/opt
git clone https://github.com/SCM-NV/PLAMS.git
cd PLAMS
git checkout oldparser
```
Now you probably get a `vim`-prompt. Just type `:wq` and hit `Enter`. Now you can continue with the next modifications: 
```bash
pip install py-ubjson
python -m pip uninstall plams
python -m pip install .
```
Now you have plams in your python3.11 stack. To make it work, you need to find the location of your `site-packages/scm/` folder. If you followed the instructions until here, you should be able to run the following lines:
```bash
cd ~/miniforge3/envs/molmod/lib/python3.11/site-packages/scm/ #if miniforge path is left standard
cp -r $AMSHOME/scripting/scm/input_parser/ .
```
The setup is done. You should test your install by running the following example simulation.

#### Run the test

For test purpose, a small and very fast system has been set up [here](./examples/water_dimer/). 
Run the script with the command `bash runner`.
Your folder should now contain `amsworker`- and `worker`-directories. If you want to run on slurm systems, check out [srunner](./examples/water_dimer/srunner).

</details>

### Initializing the run
#### [Prerequisites](./setup/NaCl_infinit/) 
- ams_inp-folder with the [AMS-engine settings](./setup/NaCl_infinit/ams_inp/ams.inp) and an initial [geometry](./setup/NaCl_infinit/ams_inp/initial.rkf) in rkf format
- [infinit.toml](./setup/NaCl_infinit/infinit.toml) - the InfRETIS input format. The key sections are closer explained [here](../chignolin/infretis_settings.md)
  
The key settings, we have to take care of, are defined in the [infinit]-section:

```
[infinit]
pL = 0.3 # the local crossing probability between interfaces
initial_conf = "ams_inp/initial.rkf" # absolute or relative path to the initial geometry
lamres = 0.001 # Resolution for the interface placement
nskip = 10  # Skipped paths for Crossing probability calculation
steps_per_iter = [
    100,
    200,
    500,
    1000,
    2000
] # Steps per short RETIS run
cstep = -1" # Current RETIS run -> -1, as first step is only a [0-] and [0+] path
```

With the WHAM (Weighted Histogram Analysis Method) feature, you can obtain detailed report of the convergence.
```
inft wham -data {infretis_data_x.txt} -lamres 0.001
```
 The running average of the rate should be reasonably smooth, as well as the overall crossing probability curve. For more detailed information, check out [this](../chignolin/infretis_output.md).

You can start the production run from here!
</details>

## PLUMED
<details>
  <summary>Click to expand</summary>
The exemplar input script is given [here](./setup/NaCl_metaD/).

For detailed explanation on how to use PLUMED, check out the [AMS Documentation](https://www.scm.com/doc/plams/examples/AMSPlumedMD/AMSPlumedMD.html) and the [PLUMED website](https://www.plumed.org/).

When we have a converged simulation, we can use [metadynminer](https://github.com/spiwokv/metadynminer) to plot the free energy surface and place the interfaces. A [python script](./setup/NaCl_metaD/get_energy_profile.py) for that can be found in the respective folder.

</details>

# Step 2: Running InfRETIS

Running infRETIS after the setup requires:
1. infretis.toml file (cf. )
2. load/ folder with a valid trajectory per ensemble (can also just be the same reactive path)
3. ams_inp folder with the [settings](./examples/NaCl/ams_inp/ams.inp) for ams and an initial [geometry](./examples/NaCl/ams_inp/initial.rkf)
4. A couple of environment variables: `NSCM` is the number of cores used per worker, unless `OMP_NUM_THREADS` is set, then `NSCM=1`.

```
export NSCM=1 # cores used
export OMP_NUM_THREADS=2 # N CPUs, only when using ReaxFF
export AMS_COPYSTATE_MOVES_RESULTS=1
infretisrun -i infretis.toml
```

To run it on SLURM clusters, we have to specify:

```
#SBATCH --nodes=X                        # Number of nodes
#SBATCH --ntasks-per-node=N_RETISWorkers/X            # Number of processes per node
#SBATCH --cpus-per-task=N_CPUs per worker

export SCM_TMPDIR={Location of temp directory}
```


# Step 3: Analyze the data
Besides the setup scripts, [inftools](https://github.com/infretis/inftools) also includes many more useful tools, to get the most out of the RETIS simulation.
You can:
- Calculate the rate (based on WHAM - publication about the aplication in InfRETIS is coming!)
- Post calculate parameters to check their predictive power for the reaction outcome 
  - This can be used to screen a large amount of parameters, to learn about reaction paths and mechanisms, as shown for [NaCl](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00054)
  - The function is soon to be added in inftools
- Obtain a free energy profile along the order parameter
