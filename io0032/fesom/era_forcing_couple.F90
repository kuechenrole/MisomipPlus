subroutine init_atm_forcing
  ! initialize the atmospheric forcing data
  use o_param
  use o_mesh
  use o_array
  use i_therm_parms
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_eraint_NetCDF
  use g_clock
  use g_parfe
  
  implicit none
  !
!  integer, parameter        		:: nci=96, ncj=72     ! hadcm3 grid
!  integer, parameter        		:: nci=96, ncj=48     ! mpeh5c grid
  integer, parameter        		:: nci=480, ncj=241     ! eraint grid
  integer                   		:: itime, i, k, n2
  integer                               :: readtype
  character(180)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy

  n2=myDim_nod2D+eDim_nod2D       

  ! predefinition/correction
  ! for the mpeh5c case:
  
  wind_ttp_ind   = 1
  rad_ttp_ind    = 1
  precip_ttp_ind = 1
  runoff_ttp_ind = 0
  sss_ttp_ind    = 0


  ! compute forcing index
  call forcing_index


  !==========================================================================
  ! wind u and v, Tair, and Tdew

  if (wind_data_source(1:10).eq.'ERAINT_20C') then
   if (yearnew.ge.1979.and.yearnew.lt.2016) then   ! RTnew 07.08.2018
    itime=(yearnew-1979)*365 + int(float(yearnew-1977)/4.) + daynew  
    ForcingDataPath=trim(ForcingDataPath)//'eraint_new/'
!    write(*,*) yearnew,daynew,'pre-2016 init itime=',itime, ForcingDataPath
   else if (yearnew.ge.2016 .and. yearnew .lt. 2018) then 
    itime= daynew
    ForcingDataPath=trim(ForcingDataPath)//'eraint_years/'//cyearnew//'/'
    write(*,*) yearnew,daynew,'post-2016 init itime=',itime, ForcingDataPath
   else
    write(*,*)'STOP: yearnew error'
    stop
   endif
  else
   write(*,*)'STOP: please define what to do with the wind'
   stop
  endif

  do i=1,2  
!   write(*,*)'init wind i,itime=',i,itime
   if (i.eq.1) then
#include "read_eraint_data1_1.incF90" ! 1 to 1_1 ha read itime ga chigau
!  write(*,*)'back from read_eraint_data1_1, now i=',i
   end if
!  write(*,*)'init wind between ifs. i=',i

   if (i.eq.2) then
#include "read_eraint_data2.incF90"
   end if
!   write(*,*)'init wind after ifs. i=',i

   if(update_forcing_flag(wind_ttp_ind)==1) then
    ! updating will take place in update_forcing in the first iteration
    u_wind_t(2,:)=u_wind_t(1,:)
    v_wind_t(2,:)=v_wind_t(1,:)
    Tair_t(2,:)=Tair_t(1,:)
    Tdew_t(2,:)=Tdew_t(1,:)
    shortwave_t(2,:)=shortwave_t(1,:)
    longwave_t(2,:)=longwave_t(1,:)
    prec_rain_t(2,:)=prec_rain_t(1,:)
    evaporation_t(2,:)=evaporation_t(1,:)
!    Pair_t(2,:)=Pair_t(1,:)
    Pair_t(2,:)=1013.e2  ! [Pa]
    exit
   end if
!   write(*,*)'init wind end of do loop , now i=',i
  end do


  !==========================================================================
  ! runoff
  ! is not provided
  !==========================================================================
  ! sss restoring
  ! is not provided

end subroutine init_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine update_atm_forcing
  ! update atmospheric forcing data
  
  
  use o_PARAM
  use o_MESH
  use o_array
  use i_array
  use i_dyn_parms
  use i_therm_parms
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_parfe
  use g_clock
  implicit none

  integer		:: i
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy
  real              	:: t1, t2

  t1=MPI_Wtime()  

  ! first, read forcing data
  call read_new_atm_forcing

  ! second, do time interpolation 

  ! wind, Tair, Tdew
  i_coef=interp_coef(wind_ttp_ind)
  do i=1,myDim_nod2d+eDim_nod2d                                        
     u_wind(i)=u_wind_t(1,i)+i_coef*(u_wind_t(2,i)-u_wind_t(1,i))
     v_wind(i)=v_wind_t(1,i)+i_coef*(v_wind_t(2,i)-v_wind_t(1,i))
     Tair(i)=Tair_t(1,i)+i_coef*(Tair_t(2,i)-Tair_t(1,i))
     Tdew(i)=Tdew_t(1,i)+i_coef*(Tdew_t(2,i)-Tdew_t(1,i))
!     Pair(i)=Pair_t(1,i)+i_coef*(Pair_t(2,i)-Pair_t(1,i))
     Pair(i)=1013.e2  ! [Pa]

  end do
  if (mype==0) write(*,*) 'update1',u_wind(100),v_wind(100),Tair(100),Tdew(100)
  
  ! radiation
  i_coef=interp_coef(rad_ttp_ind)
  do i=1,myDim_nod2d+eDim_nod2d    
     shortwave(i)=shortwave_t(1,i)!+i_coef*(shortwave_t(2,i)-shortwave_t(1,i))
     longwave(i)=longwave_t(1,i)!+i_coef*(longwave_t(2,i)-longwave_t(1,i))
  end do

  ! precipitation 
  if(precip_ttp_ind>0) then 
     i_coef=interp_coef(precip_ttp_ind)
     do i=1,myDim_nod2d+eDim_nod2d   
        prec_rain(i)=prec_rain_t(1,i)!+i_coef*(prec_rain_t(2,i)-prec_rain_t(1,i))
        evaporation(i)=evaporation_t(1,i)!+i_coef*(evaporation_t(2,i)-evaporation_t(1,i))
     end do
  end if
  if (mype==0) write(*,*)'update2', shortwave(100),longwave(100),prec_rain(100),evaporation(100)


 ! no runoff, no SSS, no chlorophyll in this setup
 

  ! third, compute exchange coefficients
  ! 1) drag coefficient 
  if(AOMIP_drag_coeff) then
     call cal_wind_drag_coeff
  end if
  ! 2) drag coeff. and heat exchange coeff. over ocean in case using ncar formulae
  if(ncar_bulk_formulae) then
     call ncar_ocean_fluxes_mode
  elseif(AOMIP_drag_coeff) then
     cd_atm_oce_arr=cd_atm_ice_arr
  end if

  ! forth, compute wind stress
  do i=1,myDim_nod2d+eDim_nod2d     
     dux=u_wind(i)-u_w(i) 
     dvy=v_wind(i)-v_w(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmoce_x(i) = Cd_atm_oce_arr(i)*aux*dux
     stress_atmoce_y(i) = Cd_atm_oce_arr(i)*aux*dvy
     dux=u_wind(i)-u_ice(i) 
     dvy=v_wind(i)-v_ice(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmice_x(i) = Cd_atm_ice_arr(i)*aux*dux
     stress_atmice_y(i) = Cd_atm_ice_arr(i)*aux*dvy
  end do

  ! heat and fresh water fluxes are treated in i_therm and ice2ocean


  t2=MPI_Wtime()

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if

end subroutine update_atm_forcing
!
!------------------------------------------------------------------------------------
!
subroutine read_new_atm_forcing
! read the second record of HadCM3 atmospheric forcing data 

use o_PARAM
use o_MESH
use o_array
use i_therm_parms
use g_forcing_param
use g_forcing_arrays
use g_forcing_index
use g_forcing_interp
use g_read_eraint_NetCDF
!use g_read_mpeh5c_NetCDF
use g_read_other_NetCDF
use g_clock
use g_parfe
implicit none
!
!integer, parameter        		:: nci=96, ncj=72 ! hadCM3 grid
!integer, parameter        		:: nci=96, ncj=48     ! mpeh5c grid
integer, parameter        		:: nci=480, ncj=241     ! mpeh5c grid
integer                   		:: itime, m, i, k, n2
integer                               :: readtype
character(180)             		:: file
character(15)             		:: vari, filevari
character(4)				:: fileyear
real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
real(kind=8), dimension(nod2D)    	:: array_fe
logical                               :: check_dummy
real(kind=8), allocatable             :: aux(:)   
real                                  :: th_era    

!==========================================================================
! wind u and v, Tair, and Tdew
n2=myDim_nod2D+eDim_nod2D                 

th_era=mod(nint(timeold),86400)/86400. !modified SB

! write(*,*) th_era

if(update_forcing_flag(wind_ttp_ind)==1) then

 !save the second record to the first record
 do i=1,myDim_nod2d+eDim_nod2d       
  u_wind_t(1,i)=u_wind_t(2,i)
  v_wind_t(1,i)=v_wind_t(2,i)
  Tair_t(1,i)=Tair_t(2,i)
  Tdew_t(1,i)=Tdew_t(2,i)
  shortwave_t(1,i)=shortwave_t(2,i)
  longwave_t(1,i)=longwave_t(2,i)
  prec_rain_t(1,i)=prec_rain_t(2,i)
  evaporation_t(1,i)=evaporation_t(2,i)
!  Pair_t(1,:)=Pair_t(2,:)
  Pair_t(1,:)=1013.e2  ! [Pa]

 end do
 
if (wind_data_source(1:10).eq.'ERAINT_20C') then
 if (yearnew.ge.1979.and.yearnew.lt.2016) then   ! RTnew 07.08.2018
  itime=(yearnew-1979)*365 + int(float(yearnew-1977)/4.) + daynew  
!  write(*,*) yearnew,daynew,'pre-2016 readnew itime=',itime, ForcingDataPath
 else if (yearnew.ge.2016) then 
  itime= daynew
  write(*,*) yearnew,daynew,'post-2016 readnew itime=',itime, ForcingDataPath
  if (yearnew.eq.2017 .and. itime.eq.365) then
   write(*,*) 'last record of 2017, set itime=:364'
   itime=364
  endif
 else
  write(*,*)'STOP: yearnew error'
  stop
 endif
else
 write(*,*)'STOP: please define your forcing data set'
 stop
endif
 

!! itime=itime+1 !this has to be removed for ERA!  !LS 20.01.2017
 
i=2  ! this is important ! new data are to be read into time level 2

if (th_era.eq.0.) then
! write(*,*) th_era, 'data2'
#include "read_eraint_data2.incF90"
elseif (th_era.eq.0.25) then
! write(*,*) th_era, 'data3' 
#include "read_eraint_data3.incF90"
elseif (th_era.eq.0.5) then
! write(*,*) th_era, 'data4'
#include "read_eraint_data4.incF90"
elseif (th_era.eq.0.75) then
! write(*,*) th_era, 'data1'
#include "read_eraint_data1.incF90"
 endif

endif

end subroutine read_new_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_eraint_wind(xarray, yarray)
  ! rotate wind on mpeh5c grid from geographical coord. to rotated coordinates.
  use o_param
  use g_rotate_grid
  implicit none

!  integer, parameter 	:: ni=96, nj=72  ! HadCM3 grid
!  integer, parameter 	:: ni=96, nj=48  ! MPEH5C grid
  integer, parameter 	:: ni=480, nj=241  ! MPEH5C grid
  integer               :: i, j
  real(kind=8)      	:: cx(ni), cy(nj), xarray(ni,nj), yarray(ni,nj) 

  
!  ! NCEP/CORE latitude
!  cy=(/-88.542, -86.6531, -84.7532, -82.8508, -80.9473, -79.0435, &  
!       -77.1394, -75.2351, -73.3307, -71.4262, -69.5217, -67.6171, &  
!       -65.7125, -63.8079, -61.9033, -59.9986, -58.0939, -56.1893, &  
!       -54.2846, -52.3799, -50.4752, -48.5705, -46.6658, -44.7611,&  
!       -42.8564, -40.9517, -39.0470, -37.1422, -35.2375, -33.3328, &  
!       -31.4281, -29.5234, -27.6186, -25.7139, -23.8092, -21.9044, &  
!       -19.9997, -18.0950, -16.1902, -14.2855, -12.3808, -10.47604, &  
!      -8.57131, -6.66657, -4.76184, -2.8571, -0.952368, 0.952368, &  
!       2.8571, 4.76184, 6.66657, 8.57131, 10.47604, 12.3808, &  
!       14.2855, 16.1902, 18.095, 19.9997, 21.9044, 23.8092, &  
!       25.7139, 27.6186, 29.5234, 31.4281, 33.3328, 35.2375,&  
!       37.1422, 39.047,  40.9517, 42.8564, 44.7611, 46.6658,&  
!       48.5705, 50.4752, 52.3799, 54.2846, 56.1893, 58.0939,&  
!       59.9986, 61.9033, 63.8079, 65.7125, 67.6171, 69.5217, &  
!       71.4262, 73.3307, 75.2351, 77.1394, 79.0435, 80.9473, &  
!       82.8508, 84.7532, 86.6531, 88.542 /)*rad

! HadCM3 latitude
!   cy=(/ 88.75, 86.25, 83.75, 81.25, 78.75, 76.25, 73.75, 71.25, 68.75, &
!    66.25, 63.75, 61.25, 58.75, 56.25, 53.75, 51.25, 48.75, 46.25, 43.75, & 
!    41.25, 38.75, 36.25, 33.75, 31.25, 28.75, 26.25, 23.75, 21.25, 18.75,  &
!    16.25, 13.75, 11.25, 8.75, 6.25, 3.75, 1.25, -1.25, -3.75, -6.25, -8.75, & 
!    -11.25, -13.75, -16.25, -18.75, -21.25, -23.75, -26.25, -28.75, -31.25,  &
!    -33.75, -36.25, -38.75, -41.25, -43.75, -46.25, -48.75, -51.25, -53.75,  &
!    -56.25, -58.75, -61.25, -63.75, -66.25, -68.75, -71.25, -73.75, -76.25, &
!    -78.75, -81.25, -83.75, -86.25, -88.75 /)*rad*(-1.)

! MPEH5C latitude
! cy = (/ 87.1590945558629, 83.4789366693172, 79.7770456548256, &
!    76.0702444625451, 72.3615810293448, 68.6520167895175, 64.9419494887575, &
!    61.2315731880771, 57.52099379797, 53.8102740319414, 50.0994534129868, &
!    46.3885581116054, 42.6776061726049, 38.966610469454, 35.2555804613682, &
!    31.5445232840217, 27.8334444519932, 24.122348326088, 20.4112384335678, &
!    16.7001176938427, 12.9889885820881, 9.27785325150786, 5.56671362791359, &
!    1.85557148599326, -1.85557148599326, -5.56671362791359, &
!    -9.27785325150786, -12.9889885820881, -16.7001176938427, &
!    -20.4112384335678, -24.122348326088, -27.8334444519932, &
!    -31.5445232840217, -35.2555804613682, -38.966610469454, &
!    -42.6776061726049, -46.3885581116054, -50.0994534129868, &
!    -53.8102740319414, -57.52099379797, -61.2315731880771, -64.9419494887575, &
!    -68.6520167895175, -72.3615810293448, -76.0702444625451, &
!    -79.7770456548256, -83.4789366693172, -87.1590945558629 /) * rad * (-1.)

! ERAINT
  cy(1)=-90.0*rad
  do i=2,nj
    cy(i)=cy(i-1)+0.75*rad
  end do

  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+0.75*rad
  enddo
    
!! Remember we turned latitudes upside down here.

!  ! NCEP/CORE longitude
!  cx(1)=0.0
!  do i=2,ni
!     cx(i)=cx(i-1)+1.875*rad
!  enddo

! HadCM3 longitude
! cx =(/ 1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, &
!    31.875, 35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, &
!    65.625, 69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, &
!    99.375, 103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, &
!    129.375, 133.125, 136.875, 140.625, 144.375, 148.125, 151.875, 155.625, &
!    159.375, 163.125, 166.875, 170.625, 174.375, 178.125, 181.875, 185.625, &
!    189.375, 193.125, 196.875, 200.625, 204.375, 208.125, 211.875, 215.625, &
!    219.375, 223.125, 226.875, 230.625, 234.375, 238.125, 241.875, 245.625, &
!    249.375, 253.125, 256.875, 260.625, 264.375, 268.125, 271.875, 275.625, &
!    279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875, 305.625, &
!    309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625, &
!    339.375, 343.125, 346.875, 350.625, 354.375, 358.125 /)*rad

! MPEH5C longitude
! cx = (/ 0., 3.75, 7.5, 11.25, 15., 18.75, 22.5, 26.25, 30., 33.75, 37.5, 41.25, &
!    45., 48.75, 52.5, 56.25, 60., 63.75, 67.5, 71.25, 75., 78.75, 82.5, 86.25, &
!    90., 93.75, 97.5, 101.25, 105., 108.75, 112.5, 116.25, 120., 123.75, 127.5, &
!    131.25, 135., 138.75, 142.5, 146.25, 150., 153.75, 157.5, 161.25, 165., &
!    168.75, 172.5, 176.25, 180., 183.75, 187.5, 191.25, 195., 198.75, 202.5, &
!    206.25, 210., 213.75, 217.5, 221.25, 225., 228.75, 232.5, 236.25, 240., &
!    243.75, 247.5, 251.25, 255., 258.75, 262.5, 266.25, 270., 273.75, 277.5, &
!    281.25, 285., 288.75, 292.5, 296.25, 300., 303.75, 307.5, 311.25, 315., &
!    318.75, 322.5, 326.25, 330., 333.75, 337.5, 341.25, 345., 348.75, 352.5, &
!    356.25 /) * rad



  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do
end subroutine rotate_eraint_wind
!
!---------------------------------------------------------------------------------------------------
