#!/bin/bash
#SBATCH --partition=GPUQ
#SBATCH --gres=gpu:1
#SBATCH --account=nv-ikj
#SBATCH --job-name=gmx-install
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=22G

module purge
module load CUDA/12.6.0
module load OpenMPI/5.0.3-GCC-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0

version=2024.4
suffix_g="_${version}_gpu"
install_path=/cluster/home/$USER/software/gromacs-${version}_gpu/

mkdir /home/$USER/software_slurm
cd /home/$USER/software_slurm
wget ftp://ftp.gromacs.org/gromacs/gromacs-$version.tar.gz
tar zxvf gromacs-$version.tar.gz
cd gromacs-$version

# Build GPU version
mkdir build_slurm
cd build_slurm
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$install_path -DGMX_DEFAULT_SUFFIX=OFF -DGMX_DOUBLE=OFF -DGMX_BINARY_SUFFIX=$suffix_g -DGMX_LIBS_SUFFIX=$suffix_g -DGMX_GPU=CUDA
make -j 12
make install -j 12
