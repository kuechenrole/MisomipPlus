! module of model configuration parameters

module g_config
  implicit none
  save

  ! *** Modelname ***
  !character(9)             	:: runid='iceOcean1'                ! a model/setup name
  !character(7)             	:: runid='RG47911'                ! a model/setup name
  character(6)             	:: runid='io0001'                ! a model/setup name
    ! ISOMIP+ settings
  character(4)                 :: case_initial='warm'           ! 'warm', 'cold'
  character(4)                 :: case_forcing='warm'           ! 'warm', 'cold'
  logical                      :: variable_GammaTS=.true.
  real(kind=8)                 :: GammaT=0.01             !

  namelist /modelname/ runid, case_initial, case_forcing, variable_GammaTS, GammaT

  ! *** time step ***
  integer                  	:: step_per_day=12           	!number of steps per day
  integer                  	:: run_length=1	                !run length
  character                     :: run_length_unit='y'          !unit: y, d, s

  namelist /timestep/ step_per_day, run_length, run_length_unit

  ! *** Paths for all in and out ***
  character(100)                :: MeshPath='./mesh/'
  character(100)                :: OpbndPath='./opbnd/'
  character(100)                :: ClimateDataPath='./hydrography/'
  character(100)                :: ForcingDataPath='./forcing/'
  character(100)                :: TideForcingPath='./tide_forcing/'
  character(100)                :: ResultPath='./result/'

  namelist /paths/  MeshPath, OpbndPath, ClimateDataPath, ForcingDataPath, &
       TideForcingPath, ResultPath

  ! *** ocean climatology data name ***
  character(100)                :: OceClimaDataName='annual_woa01_ts.out'
  logical                       :: use_prepared_init_ice=.false.     !how to initial. ice at the beginning 

  namelist /initialization/ OceClimaDataName, use_prepared_init_ice

  ! *** in out ***
  character*4               	:: restartflag='last'  	             !restart from which saved record,'#','last'
  integer                       :: output_length=1                   !valid for d,h,s
  character                	:: output_length_unit='m'      	     !output period: y, m, d, h, s 
  integer                       :: logfile_outfreq=1                 !in logfile info. output frequency, # steps

  namelist /inout/ restartflag, output_length, output_length_unit, logfile_outfreq

  ! *** mesh ***
  integer                       :: grid_type=1              	! z-level, 2 sigma, 3 sigma + z-level

  namelist /mesh_def/ grid_type

  ! *** model geometry
  logical                  	:: cartesian=.false.
  logical                  	:: fplane=.false.
  logical                  	:: betaplane=.false.
  real(kind=8)             	:: f_fplane=-1.4e-4        	![1/s]
  real(kind=8)             	:: beta_betaplane=2.0e-11  	![1/s/m]
  real(kind=8)             	:: domain_length=360.    	![degree]
  !
  logical                  	:: rotated_grid=.true.    	!option only valid for coupled model case now
  real(kind=8)             	:: alphaEuler=50. 		![degree] Euler angles, convention:
  real(kind=8)             	:: betaEuler=15.  		![degree] first around z, then around new x,
  real(kind=8)			:: gammaEuler=-90.		![degree] then around new z.

  namelist /geometry/  cartesian, fplane, betaplane, f_fplane, beta_betaplane, &
       domain_length, rotated_grid, alphaEuler, betaEuler, gammaEuler

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  
  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                       :: system=1                     ! XD1 2(byte), HLRN 1(word)
  
  namelist /machine/ system


  ! *** others ***
  real(kind=8)             	:: dt, dt_inv
  integer                  	:: istep, nsteps
  integer                       :: save_count
  logical                       :: r_restart

end module g_config


module g_forcing_param
  implicit none
  save

  ! *** exchange coefficients ***
  real*8    :: Ce_atm_oce=1.75e-3 ! exchange coeff. of latent heat over open water
  real*8    :: Ch_atm_oce=1.75e-3 ! exchange coeff. of sensible heat over open water
  real*8    :: Cd_atm_oce=1.0e-3  ! drag coefficient between atmosphere and water

  real*8    :: Ce_atm_ice=1.75e-3 ! exchange coeff. of latent heat over ice
  real*8    :: Ch_atm_ice=1.75e-3 ! exchange coeff. of sensible heat over ice
  real*8    :: Cd_atm_ice=1.32e-3 ! drag coefficient between atmosphere and ice

  namelist /forcing_exchange_coeff/ Ce_atm_oce, Ch_atm_oce, Cd_atm_oce, &
       Ce_atm_ice, Ch_atm_ice, Cd_atm_ice


  ! *** forcing source and type ***
  character(14)                 :: wind_data_source='CFSR'   !cfsr 6hourly
  character(14)                 :: rad_data_source='CFSR'
  character(14)                 :: precip_data_source='CFSR'
  character(14)                 :: runoff_data_source='none'
  character(14)                 :: sss_data_source='none'
  integer                       :: wind_ttp_ind=1
  integer                       :: rad_ttp_ind=1
  integer                       :: precip_ttp_ind=1
  integer                       :: runoff_ttp_ind=0
  integer                       :: sss_ttp_ind=0


  namelist /forcing_source/ wind_data_source, rad_data_source, precip_data_source, &
       runoff_data_source, sss_data_source, wind_ttp_ind, rad_ttp_ind, precip_ttp_ind, &
       runoff_ttp_ind, sss_ttp_ind

  ! *** coefficients in bulk formulae ***
  logical                       :: AOMIP_drag_coeff=.false.
  logical                       :: ncar_bulk_formulae=.false.

  namelist /forcing_bulk/ AOMIP_drag_coeff, ncar_bulk_formulae

  ! *** add land ice melt water ***
  logical                       :: use_landice_water=.false.
  integer                       :: landice_start_mon=1
  integer                       :: landice_end_mon=12

  namelist /land_ice/ use_landice_water, landice_start_mon, landice_end_mon

end module g_forcing_param

