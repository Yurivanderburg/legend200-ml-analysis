#!/bin/bash
#SBATCH -q shared
#SBATCH -C cpu
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --mem=64GB
#SBATCH -J pkl2hdf5
#SBATCH -A XXYY


cd ..

module load python

python -m data_preparation.scripts.pickle_to_hdf5.py input.pkl output.hdf5
