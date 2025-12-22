# Introduction

This tutoral utilizes `infinit` to create a metastable pore from an intact bilayer consisting of 128 DMPC lipids.
The input files used provided here are representative of those used in Ref [1].

The order parameter function `orderp.py` included here utilizes a linear combination of the pore expansion CV [2] and local leaflet proximity [1].

Generating 1 reactive trajectory may require around 50-200 GB of space and 12-48h of compute using a node containing a GPU.


## Required packages

Be sure to be on the latest versions of the main branches of `infretis` and `inftools`. For example in an environment, download and install,

```bash
git clone https://github.com/infretis/infretis.git
cd infretis
pip install -e .
```

and

```bash
git clone https://github.com/infretis/inftools.git
cd infretis
pip install -e .
```

To utilize the pore expansion CV imported from `dztools` in `orderp.py`, install also

```bash
pip install git+https://github.com/dz24/dztools.git@v1.0.1
pip install MDAnalysis
```

## Tuning GROMACS in infretis0.toml

In `infretis0.toml`, the default setting is to run 4 workers in parallel on 1 GPU node, modify accordingly to your hardware specs:

```toml
[runner]
workers = 4
wmdrun = [
    "gmx_mpi mdrun -ntomp 5 -pinstride 1 -pinoffset 0  -pin on -notunepme -nb gpu -bonded gpu -pme gpu",
    "gmx_mpi mdrun -ntomp 5 -pinstride 1 -pinoffset 5  -pin on -notunepme -nb gpu -bonded gpu -pme gpu",
    "gmx_mpi mdrun -ntomp 5 -pinstride 1 -pinoffset 10 -pin on -notunepme -nb gpu -bonded gpu -pme gpu",
    "gmx_mpi mdrun -ntomp 5 -pinstride 1 -pinoffset 15 -pin on -notunepme -nb gpu -bonded gpu -pme gpu",
]
```

Not only there, but also in the following section, specify if `gmx` or `gmx_mpi` is to be utilized for running the GROMACS engine.

```toml
[engine]
gmx = "gmx_mpi"
```

## Running infinit

In the root folder (`infentory/pore_formation`), we suggest to run infinit in a separate folder,

```bash
mkdir init_1
cp -r scripts gromacs_input infretis0.toml orderp.py init_1
cd init_1
```

Now in `infentory/pore_formation/init_1`, run

```bash
infinit -toml infretis0.toml
```

to run infinit.

## Analysis

The progress of `infinit` can be monitored by running `inft plot_msg` or `xip_llp.py` inside the scripts folder,

```bash
inft plot_msg
cd scripts
python3 xip_llp.py
```

We attach here the result of an `infinit` simulation that has generated 1 reactive trajectory,

![image](https://github.com/infretis/infentory/blob/pf/pore_formation/scripts/progress.png)


[1] To be added

[2] J. Chem. Theory Comput. 2021, 17, 1229−1239
