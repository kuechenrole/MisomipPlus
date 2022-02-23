#!/bin/bash

#SBATCH --job-name=io0038UA
#SBATCH -p matlab
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=30GB
##SBATCH -w fat-0004
#SBATCH -w prod-0304
#SBATCH -t 00:30:00
#SBATCH --qos=normal
#SBATCH -o /work/ollie/orichter/MisomipPlus/io0038/uarun/slurm-%j.out


module unload matlab
module load matlab

##srun matlab.x -nodisplay -r "Ua; exit;" > ice.log 2>ice.err
srun matlab.x -nodisplay -r "Ua; exit;"
#srun matlab.x -nodisplay -r "run('/home/ollie/cwekerle/matlab/write_mesh_netcdf.m')"

