#!/bin/bash
#SBATCH --job-name=u2m41 
#SBATCH -p matlab
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=30GB
###SBATCH -w prod-0306
#SBATCH -t 12:00:00
#SBATCH --qos=normal
#SBATCH -o /work/ollie/orichter/MisomipPlus/io0041/postprocessing/slurm-%j.out

module unload matlab
module load matlab/R2020b

#matlab.sh -s -S"-wprod-0304" -M"-nojvm -r ua2miso('/work/ollie/orichter/MisomipPlus/io0041/uadata','/work/ollie/orichter/MisomipPlus/io0041/IceOcean1_TYP_ice_uaFesom_SSATsai.nc',1,2400);exit"

srun matlab.x "-nojvm -r ua2miso('/work/ollie/orichter/MisomipPlus/io0041/uadata','/work/ollie/orichter/MisomipPlus/io0041/IceOcean1_TYP_ice_uaFesom_SSATsai.nc',1,2400); exit"
#srun matlab.x "-nojvm -r ver;exit"

#srun matlab.x -nodisplay -r "run('/home/ollie/cwekerle/matlab/write_mesh_netcdf.m')"
