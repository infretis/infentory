An example system that reproduces [RETIS in a 2D potential](https://pyretis.org/current/examples/examples-2d-hysteresis.html), the only difference being that we use a mass of 1.008 and not 1.0

Also, note the bug in `ase_engine.py`. If your not on the main branch, line 80 should say

```python
"friction": langevin_friction / units.fs,
```

To convert the infretis flux to pyretis units we can use $f_\text{pyretis} = f_\text{infretis}/0.002$.



### infretisrun

A 1'000 step simulation with 4 workers can be run with:

```bash
infretisrun -i infretis.toml
```

### 2D Conditional Free Energy

The free energy for phase points lying on paths that visited state A more recently than state B can be calculated using WHAM.

The wham script uses additional variables that are stored alongside the order parameter. This is achieved by returning multiple elements from the orderparameter calculation. In our case, the returned values from [orderx.py](orderx.py) are `[order, x-coord, y-coord, potential]`. We will now calculate the conditional free energy in the x-y plane. This is done by specifying `-xcol 2` and  `-ycol 3` , which selects columns 2 and 3 from the order.txt files.

```bash
# run 'inft wham -h' for a detailed description of the keywords
inft wham -data infretis_data.txt -toml infretis.toml -lamres 0.001 -nskip 100 -fener \
-nbx 100 -nby 100 -xcol 2 -ycol 3 -minx -1.2 -maxx 1.2 -miny -1.2 -maxy 1.2
```

The WHAM script needs for each accepted path a corresponding order.txt file, as output by infretis in the load/ folder. The order.txt files may also be recalculated from each accepted trajectory post-simulation, given that all trajectories are kept.

Note that the `minx, miny, maxx, maxy` values should span the whole range of the order parameter values, i.e. no value in order.txt should be outside the range spanned by these 4 values.



A minimal plotting script is available in `inftools/tistools/fener_histo.py` and can be adapted to your needs.

```bash
cd wham/
inft plot_hist
```
Below are the results after running 10'000 infretis steps with a temperature of 7'000:
![wham_fener](https://github.com/user-attachments/assets/1e94bccc-85ce-4580-98fa-72deeabf3a72)


