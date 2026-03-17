
# 1D single-particle systems (ASE)

This folder contains the input files for two 1D single particle systems simulated with *infretis*. Both systems adhere to a potential $V(x)$ symmetric around $x=0$. 

## Systems

### 1. Flat potential (`flat/`)
A simple 1D system with flat potential:
$V(x)=0$
<img src="./flat/flat_pot.png" width="350">

### 2. Cosine bump (`cos_bump/`)
A 1D system featuring a cosine bump potential barrier centered at $x=0$.
<img src="./cos_bump/cos_bump_pot.png" width="350">

## Simulation Details

The infretis simulation defines the order parameter $\lambda$ as $x$ (logically), with 3 interfaces: $\lambda=[-0.1,0,0.1]$, along with a $\lambda_{-1}$ at $-0.2$.
The Langevin dynamics parameters in reduced units (cf. PyRETIS definition) are:
- Friction parameter $\gamma=20$
- Mass $m=1.0$
- Timestep $\Delta t=0.002$
- Temperature $T=1.0$

### NOTE: unit conversion
For these systems, the parameters in reduced units (PyRETIS convention) were converted to ASE internal units, while preserving the distance and energy scale:
- $\sigma = 1.0 \text{ Å}$ (so distance units match)
- $\epsilon = 1.0 \text{ eV}$ (energy scale)

Resulting in the following conversion factors:

| Parameter | Red. unit value | Formula | ASE value |
| :--- | :--- | :--- | :--- |
| **Friction** ($\gamma$) | 20.0 | $\gamma_{red} \sqrt{\text{eV} / \text{amu}} / \text{Å}$ | 1.965 |
| **Timestep** ($\Delta t$) | 0.002 | $\Delta t_{red} \cdot \text{Å} \sqrt{\text{amu} / \text{eV}}$ | 0.0203 |
| **Temperature** ($T$) | 1.0 | $T_{red} \cdot \text{eV} / k_B$ | 11604.5 K |


## Contents

The directory structure is as follows:
```text
.
├── flat/
│   ├── load_copy/
│   └── wham_example/
├── cos_bump/
│   ├── load_copy/
│   └── wham_example/
└── results_plot/
```
In every system's folder there is an `infretis.toml` file with the settings, a `runner.sh` script from which you can launch an infretis simulation, and a `load_copy/` folder containing initial paths.
`wham_example/` contains the output of the `inft wham` command for a simulation with $100\;000$ MC moves.

In `results_plot`, the analysis plot for varying $\gamma$ and $m$ is included for both systems.

## Usage
Simply run the `runner.sh` script. You'll need the infretis package installed.
```
chmod -x runner.sh
./runner.sh
```
After runnning, the kinetic properties and conditional free energy can be computed with the following command:
```
inft wham -lamres 0.005 -fener
```