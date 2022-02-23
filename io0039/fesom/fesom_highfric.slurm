#!/bin/sh
#SBATCH --nodes=16 
#SBATCH --ntasks-per-node=36
#SBATCH --time=04:00:00
#SBATCH --job-name="RG45981"
#SBATCH -p mpp

module load intel.compiler
module load intel.mpi
module load netcdf
module list

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/AWIsoft/metis/5.1.0_intel/lib

#set | grep SLURM_

ulimit -s unlimited 

rm goonfile

rm namelist.oce
sleep 10
cp namelist.oce.highfric namelist.oce

EXE=${PWD}/fesom.x
echo now run $EXE
srun --mpi=pmi2 -n 576 --tasks-per-node=36  ${EXE}

if [ -f goonfile ]; then
 echo "goonfile exists, launch next job"
 sbatch fesom_nextmonth.slurm
else 
 echo "goonfile does not exist, we do break the chain here"
fi 





