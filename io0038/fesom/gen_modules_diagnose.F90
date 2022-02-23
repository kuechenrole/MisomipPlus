! modules of diagnose control paramers and diagnose arrays

module g_diag
  ! diagnose flag
  implicit none
  save

  ! run state
  logical                                  :: check_run_state=.true.       ! use salinity to check blowup

  ! ocean
  logical                                  :: diag_oce=.true.
  logical                                  :: diag_oce_KE=.true.
  logical                                  :: diag_oce_energy_conv=.true.
  logical                                  :: diag_oce_mix_layer=.true.
  logical                                  :: diag_oce_transp=.true.
  logical                                  :: diag_oce_GM_vel=.true.
  logical                                  :: diag_oce_SGS_transp=.true.
  logical                                  :: diag_oce_Kv=.true.
  
  ! ice
  logical                                  :: diag_ice=.true.
  ! forcing
  logical                                  :: diag_forcing=.true.
  ! mesh
  logical                                  :: diag_mesh=.true.


  namelist /diag_flag/ check_run_state, &
       diag_oce, diag_oce_KE, diag_oce_energy_conv, &
       diag_oce_mix_layer, diag_oce_transp, &
       diag_oce_GM_vel, diag_oce_SGS_transp, diag_oce_Kv, &
       diag_ice, diag_forcing, diag_mesh

end module g_diag
!
!------------------------------------------------------------------------------
!
module g_meanarrays
  ! mean and (mean)diagnose arrays
  implicit none
  save

  ! counter
  integer                                  :: meancounter

  ! mean of prediction variables

  ! ocean
  real(kind=8), allocatable, dimension(:,:):: tracermean
  real(kind=8), allocatable, dimension(:)  :: ufmean, sshmean
#ifndef use_non_hydrostatic
  real(kind=8), allocatable, dimension(:)  :: wrhsmean
#endif

  ! ice
  real(kind=8), allocatable, dimension(:)  :: a_ice_mean, m_ice_mean, m_snow_mean
  real(kind=8), allocatable, dimension(:)  :: u_ice_mean, v_ice_mean


  ! (mean) diagnose variables

  ! ocean
  real(kind=8), allocatable, dimension(:)  :: uTFmean, vTFmean
  real(kind=8), allocatable, dimension(:)  :: uSFmean, vSFmean
  real(kind=8), allocatable, dimension(:)  :: sgs_u, sgs_v
  real(kind=8), allocatable, dimension(:)  :: sgs_ut, sgs_vt
  real(kind=8), allocatable, dimension(:)  :: sgs_us, sgs_vs
  real(kind=8), allocatable, dimension(:)  :: mixlay_dep_mean
  real(kind=8), allocatable, dimension(:)  :: uumean, vvmean
  real(kind=8), allocatable, dimension(:)  :: rhomean, urhomean
  real(kind=8), allocatable, dimension(:)  :: vrhomean, uvmean

  ! ice
  real(kind=8), allocatable, dimension(:)  :: thdgr_mean, thdgrsn_mean
  real(kind=8), allocatable, dimension(:)  :: uhice_mean, vhice_mean
  real(kind=8), allocatable, dimension(:)  :: uhsnow_mean, vhsnow_mean
  real(kind=8), allocatable, dimension(:)  :: flice_mean

  ! forcing
  real(kind=8), allocatable, dimension(:)  :: tair_mean, tdew_mean, shum_mean
  real(kind=8), allocatable, dimension(:)  :: uwind_mean, vwind_mean
  real(kind=8), allocatable, dimension(:)  :: rain_mean, snow_mean
  real(kind=8), allocatable, dimension(:)  :: runoff_mean, evap_mean
  real(kind=8), allocatable, dimension(:)  :: lwrd_mean, swrd_mean
  real(kind=8), allocatable, dimension(:)  :: qnet_mean, wnet_mean
  real(kind=8), allocatable, dimension(:)  :: olat_mean, osen_mean
  real(kind=8), allocatable, dimension(:)  :: olwout_mean
  real(kind=8), allocatable, dimension(:)  :: virtual_salt_mean, relax_salt_mean
  real(kind=8), allocatable, dimension(:)  :: stress_x_mean, stress_y_mean
  real(kind=8), allocatable, dimension(:)  :: Tsurfmean, Ssurfmean
end module g_meanarrays
!
!----------------------------------------------------------------------------

