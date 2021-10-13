#!/bin/bash
#SBATCH -c 1 # Number of cores
#SBATCH -p arguelles_delgado # Partition
#SBATCH --mem 1000 # Memory request (4Gb)
#SBATCH -t 0-2:00 # Maximum execution time (D-HH:MM)
python3 dummy.pkl ${SLURM_ARRAY_TASK_ID} foo.py
