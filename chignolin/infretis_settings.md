Create a file `infretis0.toml`, and modify the number of workers and each of the workers gromacs `mdrun` settings to the settings found during the performance test.

```toml
# asyncrunner.py specific settings
[runner]
# The number of infretis workers used for the simulation.
workers = 8

# A list with commands on how to run gromacs. Each worker uses one of these commands.
wmdrun = [ 
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        'gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -nb gpu -bonded cpu -pme gpu -notunepme',
        ]
# Specific settings for the inftools software infinit.
# Infinit can create initial paths and optimize the interfaces given a single configuration file.
[infinit]
# Place the interfaces such that the local crossing probabilites are **at least** pL. The number of
# interfaces is determined by the number of workers and the current estimate of the total crossing
# probability. We place n_worker*2 interface with interfaces spaced such that the local crossing
# probabilites are equal. If they are lower than pL, we add additional interfaces such that it is
# **at least** pL. We can add a maximum number of interfaces by also specifying the variable 'num_ens'.
# For example,
#     num_ens = 14
# overwrites the behaiviour by pL and places **at most** 14 interfaces.
pL = 0.3
# The number of infretis steps to run between interfaces updates. In this case, we run 80 infretis steps,
# then we estimate the crossing proabilites using WHAM and place the interfaces as mentioned in the 'pL'
# settings section. The we pick out valid paths from the simulation we jut ran and give them as initial
# paths to the next infretis run, which will run for 160 steps.
steps_per_iter = [80, 160, 320, 640]
# Between infretis simulations, skip this many initial paths when calculating the crossing probabilities
# using WHAM
nskip = 10
# The resolution along the orderparameter when using the WHAM procedure. Should be lower than the spacing
# between interfaces by a factor of around 10
lamres = 0.001values
# The current infinit step we are on. If `cstep = -1`, we tell infinit that it should start the 
# simulation from a single configuration file. If `cstep >= 0', we can supply load/ paths and a
# set of interfaces.
cstep = -1
# Since 'cstep = -1' we have to give infinit the configuration to start the simulation from. If we
# give it a configuration with a rather large orderparameter value, we should also pump up the numbers
# in 'steps_per_iter', as we need to cover the whole orderparamter space when sampling such that we 
# can construct adequate estimates of the crossing probability
initial_conf = "gromacs_input/conf.g96"


[simulation]
# When using infinit, we have to supply the state A and state B system definitions by setting
# the first and last interfaces. infinit will then put interfaces inbetween these two states
# such that they are optimally placed
interfaces = [0.6, 6.0]
# the number of infretis steps. This will be overwritte by infinit if 'cstep = -1'
steps = 100000
# the seed infretis uses
seed = 0
# paths are created in the load/ folder
load_dir = 'load'
# In [0-] and [0+] we use the shooting method to generate new paths. infinit will then add wire-fencing
# moves to this list for all except the [0-] and [0+] ensembles as more interfaces are placed
shooting_moves = ['sh', 'sh']

[simulation.tis_set]
# The maximum path-length, which is a number that should never be reached in practice to achieve correct
# sampling. However, we keep it here because in case of a crash, a zombie process might keep running and
# never terminate
maxlength = 160000
# When using 'sh' moves, paths can be rejected based on the early rejection scheme. This means that when
# a path that we currently are generating is longer than this value, the path generation is stopped and
# we reject the path. If 'allowmaxlength = true', we don't terminate and accept the path even if it is
# too long
allowmaxlength = false
# After generating velocities, subtract the total center-of-mass motion such that it is zero. This is
# identical to setting the total momentum of the system to zero, thereby the name. This option should
# be false for small dimensional systems like e.g. double-wells or vacuum simulations
zero_momentum = true # momentum true
# number of subtrajectories to generate
n_jumps = 3 
# interface cap, meaning we never shoot from phasepoints having orderparameter values larger than
# 'interface_cap'. The interface cap should be placed to avoid B to B paths, and must be between the
# second-to-last and last interface.
interface_cap = 3.0 

[engine]
class = 'gromacs'
engine = 'gmx'
timestep = 0.002
gmx_format = 'g96'
input_path = 'gromacs_input'
gmx = 'gmx_2024.4_gpu'
subcycles = 500 
temperature = 340 

[orderparameter]
class = 'NeuralNetDiffmap'
module = '/cluster/work/lukasba/rerun-chig/orderp.py'

[output]
data_dir = "./"
screen = 1 
pattern = false
# Specifies if we delete trajectories that are no longer "active", that is, thet have been pushed
# out of the sampling because a new trajetory has been accepted that was generated by this "old" 
# trajectory. This option needs to be set to false with infinit, because we may need these new paths
# again later when filling up new ensembles with initial paths. When running infretis, this option can
# be set to true, in which case the trajectory files are deleted. This saves alot of space. 
delete_old = false
# If we want to start a main simulation with infretis, it could be wise to set 'delete_old = true'.
# Then, via the gromacs .mdp options, we can tell gromacs to also write .xtc files every N-steps.This 
# can save some space. So .trr files are deleted because they correspond to the trajectory format of
# gromacs (.trr) if 'delte_old = true', but .xtc files are not deleted with the below option. We can also
# add e.g. '.edr' to the list.
keep_traj_fnames = [".xtc"]
```

```sh-session
RuntimeError: Execution of external program (GROMACS engine zamn) failed with command:
 gmx_2024.4_gpu mdrun -ntomp 2 -ntmpi 1 -notunepme -nb gpu -bonded cpu -pme gpu -notunepme -s genvel.tpr -deffnm genvel -c genvel.g96.
Return code: 1
```

We can locate the `stderr.txt` and `stdout.txt` files to see what went wrong.

### Starting from a previous MD simulation (GROMACS only)

In some cases we might want to use an MD simulation as a starting point for our simulation instead of a single configuration. This allows us to imediately use N workers, instead of generating zero-paths with 1 worker, which can be tedious, and doesn't immediately use all of the allocated hardware.

```bash
inft initial_path_from_md -trr md_run.trr -order order_rec.txt -toml infretis0.toml
```

The `order_rec` should correspond to an orderparameter file, which can be calculated with

```bash
inft recalculate_order -traj md_run.trr -toml infretis0.toml -out order_rec.txt
```

### Running an infinit simulation

```bash
# first time we call infinit
inft infinit -toml infretis0.toml
# to continue the simulation after changing some options (which have to be changed in both restart.toml
# and infretis.toml)
inft infinit -toml infretis.toml
```
