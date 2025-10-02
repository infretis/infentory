#!/bin/bash
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=71:59:00
#PBS -r n
#PBS -N "Main_Infinit"
#PBS -A "reaction-network"

### replace with your conda module and env
module load Miniconda/3.1
conda activate /path/to/your/env

export NUM_WORKERS=4

### optinal HPC-specific job variables (remove from L21 if not used)
export ACCEL_TYPE="gtx1080ti"

cd $PBS_O_WORKDIR

job_ids=()
for ((i=0; i<NUM_WORKERS; i++)); do
    jid=$(qsub -l select=1:ncpus=2:mem=8gb:ngpus=1:accelerator_model=$ACCEL_TYPE \    ### change as needed
              -l walltime=71:59:00 \
              -N "WORKER_${i}" \
              -v WORKER_ID=$i \
              worker_submit.sh)

    echo "Submitted worker $i with job ID $jid" 
    job_ids+=($jid)
done

inft infinit -toml infretis.toml

echo "Killing worker jobs..." 
for jid in "${job_ids[@]}"
do
   echo "Killing $jid" 
   qdel $jid
done


