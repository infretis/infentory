An example system that reproduces [RETIS in a 2D potential](https://pyretis.org/current/examples/examples-2d-hysteresis.html), the only difference being that we use a mass of 1.008 and not 1.0

Also note the bug in `ase_engine.py`. Line 80 should say

```python
"friction": langevin_friction / units.fs,
```

To convert the infretis flux to pyretis units we can use $f_\text{pyretis} = f_\text{infretis}/0.002$.


