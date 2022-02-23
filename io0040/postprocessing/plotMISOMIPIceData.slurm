#!/bin/bash
##SBATCH --time=17:30:00 ##fuer 100 jahre
#SBATCH --time=12:00:00
#SBATCH --mem=5GB
#SBATCH -o /work/ollie/orichter/MisomipPlus/io0040/postprocessing/slurm-%j.out
#SBATCH --job-name=f2mio0040

module unload python
module load python/2.7.14-intel2018

date

cd /home/ollie/orichter/XylarSupplement/misomip1
srun python plotMISOMIPIceData.py /work/ollie/orichter/MisomipPlus/io0040/postprocessing/IceOcean1_COM_ice_uaFesom_SSAWeertman.nc /work/ollie/orichter/MisomipPlus/io0040/postprocessing/figures/iceData/.

date
