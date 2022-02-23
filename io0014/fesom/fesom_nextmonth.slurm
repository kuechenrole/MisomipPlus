#!/bin/sh
#SBATCH --nodes=16 
#SBATCH --ntasks-per-node=36
#SBATCH --time=01:00:00
#SBATCH --job-name="RG47906"
#SBATCH -p mpp

module load intel.compiler
module load intel.mpi
module load netcdf
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/AWIsoft/metis/5.1.0_intel/lib

#set | grep SLURM_

ulimit -s unlimited 

rm goodfile
rm goonfile

rm namelist.oce
sleep 10
cp namelist.oce.norm namelist.oce


EXE=${PWD}/fesom.x
echo now run $EXE
srun --mpi=pmi2 -n 576 --tasks-per-node=36  ${EXE}

#if [ -f goonfile ]; then      # use this for ice sheetr/ocean coupled runs
if [ -f goodfile ]; then      # use this for standalone FESOM runs
 echo "goodfile exists, launch next job"
 sbatch fesom_nextmonth.slurm
else 
 echo "goodfile does not exist, we do break the chain here"
fi 





