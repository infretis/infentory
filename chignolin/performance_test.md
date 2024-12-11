### Performance test

In the below, we assume gromacs has been installed with the script `install-gromacs-slurm.sh`, as we need to load some specific modules.

Before starting an infretis simulation, you have to do a performance test so you know that all the allocated hardware is used efficiently.

These settings can be rather difficult to set because we now have a situation where we have to optimize the performance of not a single simulation, but rather the total simulation time.

In the usual situation, with a single MD simulation, we often just pump up the numbers of the CPUs used, even if we get a slight performance increase from e.g. 8 to 12 CPUs. With infretis however, it might be more beneficial to keep the number of CPUs low, and instead increase the number of workers. For example, using 2 workers with 6 CPUs each, or maybe even 6 workers with 2 CPUs each could be more beneficial than using just 1 worker with 12 CPUs, depending on the scaling with the number of CPUs.

A way to determine the settings could be to start with 1 worker and 1 CPU and increase the number of workers to what is reasonable for that node. Then, again start with 1 worker and now 2 CPUs and increase it to what is possible on that node, and so on with 4 CPUs, etc. Then, one would compute the total simulation throughput and pick the settings that maximize this number.

Depending on your situation and how fast you want your results, you could also then increase the number of nodes used, and as such reach a large number of total simulation throughput and minimize the time you have to wait to get your results. 

In our case, we have decided beforehand that each worker will use 2 CPUs and a single GPU that is shared between all workers, and we don't want to run on more than 1 node. Here, we then only optimize the number of workers, since the hardware was specified beforehand.

We run the script with 1, 2, 4, and 8 workers, and then check the performance.

```bash
#!/bin/bash
#SBATCH --partition=GPUQ
#SBATCH --gres=gpu:1
#SBATCH --account=nv-ikj
#SBATCH --job-name=gmx-test
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1

module purge
module load CUDA/12.6.0
module load OpenMPI/5.0.3-GCC-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0

source ~/software/gromacs-2024.4_gpu/bin/GMXRC

# only need a runnable .tpr file run.tpr
mkdir w4
cd w4
cp ../run.tpr .

export CUDA_MPS_LOG_DIRECTORY=mps_logs
mkdir -p $CUDA_MPS_LOG_DIRECTORY
nvidia-cuda-mps-control -d # when using more than 1 worker with gromacs

gmx_2024.4_gpu mdrun -deffnm md0 -s run.tpr -ntomp 2 -ntmpi 1 -resethway -nsteps 5000 -pme gpu -nb gpu -bonded cpu -notunepme &
gmx_2024.4_gpu mdrun -deffnm md1 -s run.tpr -ntomp 2 -ntmpi 1 -resethway -nsteps 5000 -pme gpu -nb gpu -bonded cpu -notunepme &
gmx_2024.4_gpu mdrun -deffnm md2 -s run.tpr -ntomp 2 -ntmpi 1 -resethway -nsteps 5000 -pme gpu -nb gpu -bonded cpu -notunepme &
gmx_2024.4_gpu mdrun -deffnm md3 -s run.tpr -ntomp 2 -ntmpi 1 -resethway -nsteps 5000 -pme gpu -nb gpu -bonded cpu -notunepme &

wait

echo quit | nvidia-cuda-mps-control # when using more than 1 worker with gromacs
```

Note that it is wise to test other settings, e.g. varying the number of CPUs with the -ntomp flag, etc.

```bash
grep 'Performance' w*/*
# output:
w1/md0.log:Performance:      343.885        0.070

w2/md0.log:Performance:      351.684        0.068
w2/md1.log:Performance:      344.305        0.070

w4/md0.log:Performance:      286.784        0.084
w4/md1.log:Performance:      288.329        0.083
w4/md2.log:Performance:      296.531        0.081
w4/md3.log:Performance:      290.186        0.083

w8/md0.log:Performance:      230.124        0.104
w8/md1.log:Performance:      225.502        0.106
w8/md2.log:Performance:      226.892        0.106
w8/md3.log:Performance:      220.539        0.109
w8/md4.log:Performance:      226.064        0.106
w8/md5.log:Performance:      226.791        0.106
w8/md6.log:Performance:      224.986        0.107
w8/md7.log:Performance:      227.983        0.105
```

We stick to 8 workers in our case, since this gives around $8 \times 230 \text{ ns/day}$.

### Installing Infretis (and inftools)

```bash
# first load required packages to run gromacs
# in my case i need to load the following
module purge
module load CUDA/12.6.0
module load OpenMPI/5.0.3-GCC-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0

source ~/software/gromacs-2024.4_gpu/bin/GMXRC

python -m venv ~/iretis # create virtual environment
source ~/iretis/bin/activate # activate environment
git clone https://github.com/infretis/infretis.git
cd infretis/
python -m pip install -e .
cd -

git clone https://github.com/infretis/inftools.git
cd inftools
python -m pip install -e .
cd -
```

