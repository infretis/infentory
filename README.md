# infentory
An inventory of example systems for the infretis software.

These examples are not directly runnable out of the box, as initial paths are not provided.


<p>
<img src="/chignolin/movie.gif" width="45%" height="45%">
</p>

## Initial paths
One strategy to generate initial paths is via the [inftools](https://github.com/infretis/inftools/tree/main) functionality: 

```bash
# running in the the chignolin folder
inft generate_zero_paths -toml infretis.toml -conf initial.xyz -maxlen 100
```

This generates paths in the [0-] and [0+] ensembles. You can then follow the initial-path procedure in the [puckering tutorial](https://github.com/infretis/infretis/tree/main/examples/gromacs/puckering#step-3-retis) to push the system over the barrier using path-sampling.

