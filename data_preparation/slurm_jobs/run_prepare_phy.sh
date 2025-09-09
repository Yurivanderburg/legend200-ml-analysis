#!/bin/bash
#SBATCH -q shared
#SBATCH -C cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH -J Prepare
#SBATCH -A XXYY

export HDF5_USE_FILE_LOCKING=FALSE 
cd ..

module load python

python -m data_preparation.scripts.prepare_pht_dataset --config data_preparation/configs/example_config_pht.py