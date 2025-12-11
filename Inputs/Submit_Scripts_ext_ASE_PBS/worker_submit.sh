#!/bin/bash
#PBS -j oe

### replace with your conda module and env
module load Miniconda/3.1
conda activate /path/to/your/env

cd $PBS_O_WORKDIR

python -m infretis.classes.engines.propagator infretis.toml worker${WORKER_ID}

