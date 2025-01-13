
# mini wip pp tutorial

### requirements

1. Require infpp from fork https://github.com/dz24/infretis/

Maybe something this works (did not test myself)


```
git remote add dz24 https://github.com/dz24/infretis.git
git pull
git checkout infpp
```

then install.

2. pyretis for analysis

see [https://pyretis.org/current/user/install.html](https://pyretis.org/current/user/install.html)

### Example run

To see if it works you can run the example (peek ``runner.sh`` to see the actual commands being ran)

```bash runner.sh```

and from inftools, run in terminal

```inft plot_ens -pp```

to see that sampled trajectories stop between neighbouring interfaces

### To analyze results using pyretisanalyse -i retis.rst

`python3 conv_inf_py.py`

for converting `infretis_data.txt` to `0**/pathensemble.txt` files. then 

```pyretisanalyse -i fakeretis.rst```

### For creating initial paths for any other system

use `init_paths_pp.py`
which require a initial trajectory and its assoicated order parameter values in a text file.


### NB

Only works with `sh` moves or you'll get `KeyError: 'must_cross_M'`.
