#!/bin/bash
##SBATCH --time=17:30:00 ##fuer 100 jahre
#SBATCH --time=12:00:00
#SBATCH --mem=5000M
#SBATCH -o /work/ollie/orichter/MisomipPlus/io0039/postprocessing/slurm-%j.out
#SBATCH --job-name=f2mio0039

module unload python
module load python/2.7.14-intel2018

date

cd /home/ollie/orichter/tools/
srun python fesom2miso.py io0039 0 2401 IceOcean1_COM_ocean_uaFesom_SSATsai.nc

date
