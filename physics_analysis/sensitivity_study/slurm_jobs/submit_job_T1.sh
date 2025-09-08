#!/bin/bash                                                                                                                                                 
#SBATCH -q regular                                                                                                                                       
#SBATCH --constraint=cpu                                                                                                                                    
#SBATCH -t 35:00:00
#SBATCH -J sens_T                                                                                                                                         
#SBATCH --output parallel_T_fix.log                                                     
#SBATCH --error parallel_T_fix.err  

# Go to parent directory
cd ..

module load parallel 
module load julia
srun="srun -N 1"
parallel="parallel --delay 1 -j 128"


$srun $parallel "julia sensitivity.jl -c config/fake_config_T1.json -i {1} -f true -b 2" ::: {1..10000} &

wait