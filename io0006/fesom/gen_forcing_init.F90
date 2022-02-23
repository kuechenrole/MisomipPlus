! allocate the surface forcing arrays

subroutine forcing_array_setup
  !inializing forcing fields 
  ! 
  ! Coded by Qiang Wang
  !
  ! modified for more universal use (on the expense of increased memory usage)
  ! by Ralph Timmermann, 09.06.12
  !------------------------------------------------------------------
  
  use o_param
  use o_mesh
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe
  use g_config
  implicit none
  
  integer    :: n2

  n2=ToDim_nod2D  

  ! Allocate memory for atmospheric forcing
  allocate(shortwave(n2), longwave(n2))
  allocate(prec_rain(n2), prec_snow(n2))
  allocate(u_wind(n2), v_wind(n2), Pair(n2))
  allocate(Tair(n2), shum(n2), Tdew(n2))
  allocate(runoff(n2), evaporation(n2))
  allocate(shortwave_t(2,n2), longwave_t(2,n2))
  allocate(prec_rain_t(2,n2), prec_snow_t(2,n2), prec_net_t(2,n2))
  allocate(u_wind_t(2,n2), v_wind_t(2,n2))
  allocate(Tair_t(2,n2), shum_t(2,n2), Tdew_t(2,n2))
  allocate(runoff_t(2,n2), cloud_t(2,n2))
  allocate(Pair_t(2,n2),evaporation_t(2,n2), e_vapor_t(2,n2))
  shortwave=0.
  longwave=0.
  prec_rain=0.
  prec_snow=0.
  u_wind=0.
  v_wind=0.
  Tair=0.
  Tdew=0.
  Pair=0.
  shum=0.
  runoff=0.
  shortwave_t=0.
  longwave_t=0.
  prec_rain_t=0.
  prec_snow_t=0.
  u_wind_t=0.
  v_wind_t=0.
  Tair_t=0.
  Tdew_t=0.
  Pair_t=0.
  shum_t=0.
  runoff_t=0.


  if(use_landice_water) then
     allocate(runoff_landice(n2))
     runoff_landice=0.0
  end if
 
  ! shortwave penetration
#ifdef use_sw_pene
  allocate(chl(n2))
  allocate(sw_3d(myDim_nod3d+eDim_nod3D))
  chl=0.0
#endif

  !for ice diagnose
#ifdef use_ice
  allocate(thdgr(n2), thdgrsn(n2), flice(n2))
  allocate(olat_heat(n2), osen_heat(n2), olwout(n2))
  thdgr=0.
  thdgrsn=0.
  flice=0.
  olat_heat=0.
  osen_heat=0.
  olwout=0.
#endif 

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  allocate(Ce_atm_oce_arr(n2))
  allocate(Ch_atm_oce_arr(n2))
  Cd_atm_oce_arr=Cd_atm_oce
  Ce_atm_oce_arr=Ce_atm_oce 
  Ch_atm_oce_arr=Ch_atm_oce
#ifdef use_ice
  allocate(Cd_atm_ice_arr(n2)) 
  Cd_atm_ice_arr=Cd_atm_ice   
#endif

  if(mype==0) write(*,*) 'forcing arrays have been set up'   

end subroutine forcing_array_setup
!
!----------------------------------------------------------------------
!
subroutine forcing_array_setup_OnlyOcean
  !inializing forcing fields for an ocean-alone case
  !currently only wind is applied.
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_param
  use o_mesh
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer    :: n2

  n2=myDim_nod2D+eDim_nod2D  

  ! Allocate memory for atmospheric forcing 
  allocate(u_wind(n2), v_wind(n2))
  u_wind=0.
  v_wind=0.

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  Cd_atm_oce_arr=Cd_atm_oce

  if(mype==0) write(*,*) 'forcing arrays (for an ocean-alone case) have been set up'   

end subroutine forcing_array_setup_OnlyOcean
