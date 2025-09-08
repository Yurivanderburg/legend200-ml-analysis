#!/bin/sh
#SBATCH -q shared
#SBATCH -C cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=12:00:00
#SBATCH -J ReMaGe
#SBATCH -A XXYY


# Parameters
N_THREADS=128
GEOMETRY="path/to/geometry_file.gdml"
OUTPUT_FILE="path/to/output.lh5" 
MACRO_NAME="path/to/settings.mac" 


# Run remage -> Create $N_THREADS .lh5 files containing GEANT4 output
cenv remage remage --threads $N_THREADS -g $GEOMETRY -o $OUTPUT_FILE $MACRO_NAME
