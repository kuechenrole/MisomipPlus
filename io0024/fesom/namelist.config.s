! This is the namelist file for model general configuration

&modelname
runid='RG47911'
case_initial='warm'
case_forcing='warm'
variable_GammaTS=.false.
GammaT=0.01
/

&timestep
step_per_day= 120!360        ! 960
run_length=5
run_length_unit='y' 		! y, m, d, s
/

&clockinit			! the model starts at
timenew=0.0
daynew=1
yearnew=2021
/

&paths
MeshPath='/work/ollie/orichter/mesh/oce0_s/'
OpbndPath=' '
ClimateDataPath= ''
ForcingDataPath=''
TideForcingPath= ''
ResultPath     = '/work/ollie/orichter/data/oce0_s/'


/

&initialization
use_prepared_init_ice=.false.   !how to init. ice; runid.initial_ice.nc
OceClimaDataName=''
/

&inout
restartflag='last'		!restart from which saved record,'last,'#'
output_length=1			!only required for m,d,h,s cases,  y takes 1
output_length_unit='m'   	!output period: y, m, d, h, s 
logfile_outfreq=1  	        !in logfile info. output frequency, # steps
/

&mesh_def
grid_type=2			!1 z-level, 2 sigma, 3 z+sigma
/

&geometry
domain_length=360.    	        ![degree]
cartesian=.false.
fplane=.true.
betaplane=.false.
f_fplane=-1.405e-4        	![1/s]
beta_betaplane=2.0e-11  	![1/s/m]
rotated_grid=.false. 	  	!option only valid for coupled model case now
alphaEuler=50.			![degree] Euler angles, convention:
betaEuler=15.	 		![degree] first around z, then around new x,
gammaEuler=-90.			![degree] then around new z.
/

&calendar
include_fleapyear=.false.
/

