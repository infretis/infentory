This tutorial illustrates the use of the `infinit` functionality, which generates initial paths and places the interfaces in an optimal manner, starting out from just a single configuration!

This tutorial can easily be adapted to a large number of different systems with minimal modifications!

The process we will study is the pyramidal inversion of the NH3 molecule using the `xtb-python` package, which can be install with conda

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

# Step 1: Setting up the order parameter
We here the dihedral angle between the 4 atoms of NH3 to describe the progress of the reaction.

```toml
[orderparameter]
class = "Dihedral"
index = [0,1,2,3]
```
