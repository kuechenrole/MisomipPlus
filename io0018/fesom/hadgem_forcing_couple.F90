subroutine init_atm_forcing
  ! initialize the atmospheric forcing data 
  ! assume forcing data from hadgem2
  use o_param
  use o_mesh
  use o_array
  use i_therm_parms
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_hadgem2_NetCDF
!  use g_read_other_NetCDF
  use g_clock
  use g_parfe
  
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=144     ! hadgem2 grid
  integer                   		:: itime, i, k, n2, itime_evspsbl
  integer                               :: readtype
  character(110)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
!rt  real(kind=8), allocatable             :: aux(:) 

  n2=myDim_nod2D+eDim_nod2D       

  ! predefinition/correction
  ! for the hadgem2 case:
  
  wind_ttp_ind   = 2
  rad_ttp_ind    = 2
  precip_ttp_ind = 2
  runoff_ttp_ind = 0
  sss_ttp_ind    = 0


  ! compute forcing index
  call forcing_index


  !==========================================================================
  ! wind u and v, Tair, and shum

  if (wind_data_source(1:7).ne.'HadGem2') then
   write(*,*)'STOP: This subroutine is designed only for use of HadGem2 data.'
   stop
  endif

  if (wind_data_source.eq.'HadGem2_20C') then
   itime=(yearnew-1950)*360+30+daynew ! HadGem2 20C data start 01.12.1949
  elseif (wind_data_source.eq.'HadGem2_RCP85') then
   itime=(yearnew-2005)*360+daynew ! for HadGem2 RCP data start 01.01.2005
  else
   write(*,*)'STOP: please define what to do in this case'
   stop
  endif
  

  do i=1,2
  
     write(*,*)'i,itime=',i,itime

#include "read_hadgem2_data.incF90"

   if(update_forcing_flag(wind_ttp_ind)==1) then
    ! updating will take place in update_forcing in the first iteration
    u_wind_t(2,:)=u_wind_t(1,:)
    v_wind_t(2,:)=v_wind_t(1,:)
    Tair_t(2,:)=Tair_t(1,:)
    shum_t(2,:)=shum_t(1,:)
    shortwave_t(2,:)=shortwave_t(1,:)
    longwave_t(2,:)=longwave_t(1,:)
    prec_rain_t(2,:)=prec_rain_t(1,:)
    evaporation_t(2,:)=evaporation_t(1,:)
    Pair_t(2,:)=Pair_t(1,:)
    exit
   end if

   itime=itime+1
   if (itime>34200) itime=34200  ! HadCM3 files contain 100 years: very last record is used twice (dirty)


  enddo
  write(*,*)'init_atm_forcing completed'


  !==========================================================================
  ! runoff
  ! is not provided for the HadCM3 and HadGem2 experiments


  !==========================================================================
  ! sss restoring
  ! is not provided for the HadCM3 and HadGem2 experiments

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

  ! wind, Tair, shum
  i_coef=interp_coef(wind_ttp_ind)
  do i=1,myDim_nod2d+eDim_nod2d                                        
     u_wind(i)=u_wind_t(1,i)+i_coef*(u_wind_t(2,i)-u_wind_t(1,i))
     v_wind(i)=v_wind_t(1,i)+i_coef*(v_wind_t(2,i)-v_wind_t(1,i))
     Tair(i)=Tair_t(1,i)+i_coef*(Tair_t(2,i)-Tair_t(1,i))
     shum(i)=shum_t(1,i)+i_coef*(shum_t(2,i)-shum_t(1,i))
     Pair(i)=Pair_t(1,i)+i_coef*(Pair_t(2,i)-Pair_t(1,i))
  end do
!  write(*,*) 'update1',u_wind(100),v_wind(100),Tair(100),shum(100),Pair(100)
  
  ! radiation
  i_coef=interp_coef(rad_ttp_ind)
  do i=1,myDim_nod2d+eDim_nod2d    
     shortwave(i)=shortwave_t(1,i)+i_coef*(shortwave_t(2,i)-shortwave_t(1,i))
     longwave(i)=longwave_t(1,i)+i_coef*(longwave_t(2,i)-longwave_t(1,i))
  end do

  ! precipitation 
  if(precip_ttp_ind>0) then 
     i_coef=interp_coef(precip_ttp_ind)
     do i=1,myDim_nod2d+eDim_nod2d   
        prec_rain(i)=prec_rain_t(1,i)+i_coef*(prec_rain_t(2,i)-prec_rain_t(1,i))
        evaporation(i)=evaporation_t(1,i)+i_coef*(evaporation_t(2,i)-evaporation_t(1,i))
     end do
  end if
!  write(*,*)'update2', shortwave(100),longwave(100),prec_rain(100),evaporation(100)


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
use g_read_hadgem2_NetCDF
use g_read_other_NetCDF
use g_clock
use g_parfe
implicit none
!
integer, parameter        		:: nci=192, ncj=144 ! hadgem2 grid
integer                   		:: itime, m, i, k, n2, itime_evspsbl
integer                               :: readtype
character(110)             		:: file
character(15)             		:: vari, filevari
character(4)				:: fileyear
real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
real(kind=8), dimension(nod2D)    	:: array_fe
logical                               :: check_dummy
real(kind=8), allocatable             :: aux(:)       

!==========================================================================
! wind u and v, Tair, and shum
n2=myDim_nod2D+eDim_nod2D                 

if(update_forcing_flag(wind_ttp_ind)==1) then

 !save the second record to the first record
 do i=1,myDim_nod2d+eDim_nod2d       
  u_wind_t(1,i)=u_wind_t(2,i)
  v_wind_t(1,i)=v_wind_t(2,i)
  Tair_t(1,i)=Tair_t(2,i)
  shum_t(1,i)=shum_t(2,i)
  shortwave_t(1,i)=shortwave_t(2,i)
  longwave_t(1,i)=longwave_t(2,i)
  prec_rain_t(1,i)=prec_rain_t(2,i)
  evaporation_t(1,i)=evaporation_t(2,i)
  Pair_t(1,:)=Pair_t(2,:)
 end do
 
 if (wind_data_source.eq.'HadGem2_20C') then
  itime=(yearnew-1950)*360+30+daynew ! HadGem2 20C data start 01.12.1949
 elseif (wind_data_source.eq.'HadGem2_RCP85') then
  itime=(yearnew-2005)*360+daynew ! for HadGem2 RCP data start 01.01.2005
 else
  write(*,*)'STOP: please define what to do in this case'
  stop
 endif
 itime=itime+1

 if (itime>34200) then
  write(*,*)'limit itime to 34200'
  itime=34200  ! HadCM3 files contain 95 years: very last record is used twice (dirty)
 endif
 
 i=2  ! this is important ! new data are supposed to be read into time level 2

#include "read_hadgem2_data.incF90"

endif

end subroutine read_new_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_hadgem2_wind(xarray, yarray)
  ! rotate wind on hadgem2 grid from geographical coord. to rotated coordinates.
  use o_param
  use g_rotate_grid
  implicit none

  integer, parameter 	:: ni=192, nj=144  ! HadCM3 grid
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
!! Remember the *(-1) turns latitudes upside down here.

! HadGem2 latitudes

cy=(/ -89.375, -88.125, -86.875, -85.625, -84.375, -83.125, -81.875,  &
    -80.625, -79.375, -78.125, -76.875, -75.625, -74.375, -73.125, -71.875, &
    -70.625, -69.375, -68.125, -66.875, -65.625, -64.375, -63.125, -61.875, &
    -60.625, -59.375, -58.125, -56.875, -55.625, -54.375, -53.125, -51.875, &
    -50.625, -49.375, -48.125, -46.875, -45.625, -44.375, -43.125, -41.875, &
    -40.625, -39.375, -38.125, -36.875, -35.625, -34.375, -33.125, -31.875, &
    -30.625, -29.375, -28.125, -26.875, -25.625, -24.375, -23.125, -21.875, &
    -20.625, -19.375, -18.125, -16.875, -15.625, -14.375, -13.125, -11.875, &
    -10.625, -9.375, -8.125, -6.875, -5.625, -4.375, -3.125, -1.875, -0.625, &
    0.625, 1.875, 3.125, 4.375, 5.625, 6.875, 8.125, 9.375, 10.625, 11.875, &
    13.125, 14.375, 15.625, 16.875, 18.125, 19.375, 20.625, 21.875, 23.125, &
    24.375, 25.625, 26.875, 28.125, 29.375, 30.625, 31.875, 33.125, 34.375, &
    35.625, 36.875, 38.125, 39.375, 40.625, 41.875, 43.125, 44.375, 45.625, &
    46.875, 48.125, 49.375, 50.625, 51.875, 53.125, 54.375, 55.625, 56.875, &
    58.125, 59.375, 60.625, 61.875, 63.125, 64.375, 65.625, 66.875, 68.125, &
    69.375, 70.625, 71.875, 73.125, 74.375, 75.625, 76.875, 78.125, 79.375, &
    80.625, 81.875, 83.125, 84.375, 85.625, 86.875, 88.125, 89.375 /)*rad



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


! HadGem2 longitudes
 cx=(/0.9375, 2.8125, 4.6875, 6.5625, 8.4375, 10.3125, 12.1875, 14.0625,  &
    15.9375, 17.8125, 19.6875, 21.5625, 23.4375, 25.3125, 27.1875, 29.0625,  &
    30.9375, 32.8125, 34.6875, 36.5625, 38.4375, 40.3125, 42.1875, 44.0625,  &
    45.9375, 47.8125, 49.6875, 51.5625, 53.4375, 55.3125, 57.1875, 59.0625,  &
    60.9375, 62.8125, 64.6875, 66.5625, 68.4375, 70.3125, 72.1875, 74.0625,  &
    75.9375, 77.8125, 79.6875, 81.5625, 83.4375, 85.3125, 87.1875, 89.0625,  &
    90.9375, 92.8125, 94.6875, 96.5625, 98.4375, 100.3125, 102.1875,  &
    104.0625, 105.9375, 107.8125, 109.6875, 111.5625, 113.4375, 115.3125,  &
    117.1875, 119.0625, 120.9375, 122.8125, 124.6875, 126.5625, 128.4375,  &
    130.3125, 132.1875, 134.0625, 135.9375, 137.8125, 139.6875, 141.5625,  &
    143.4375, 145.3125, 147.1875, 149.0625, 150.9375, 152.8125, 154.6875,  &
    156.5625, 158.4375, 160.3125, 162.1875, 164.0625, 165.9375, 167.8125,  &
    169.6875, 171.5625, 173.4375, 175.3125, 177.1875, 179.0625, 180.9375,  &
    182.8125, 184.6875, 186.5625, 188.4375, 190.3125, 192.1875, 194.0625,  &
    195.9375, 197.8125, 199.6875, 201.5625, 203.4375, 205.3125, 207.1875,  &
    209.0625, 210.9375, 212.8125, 214.6875, 216.5625, 218.4375, 220.3125,  &
    222.1875, 224.0625, 225.9375, 227.8125, 229.6875, 231.5625, 233.4375,  &
    235.3125, 237.1875, 239.0625, 240.9375, 242.8125, 244.6875, 246.5625,  &
    248.4375, 250.3125, 252.1875, 254.0625, 255.9375, 257.8125, 259.6875,  &
    261.5625, 263.4375, 265.3125, 267.1875, 269.0625, 270.9375, 272.8125,  &
    274.6875, 276.5625, 278.4375, 280.3125, 282.1875, 284.0625, 285.9375,  &
    287.8125, 289.6875, 291.5625, 293.4375, 295.3125, 297.1875, 299.0625,  &
    300.9375, 302.8125, 304.6875, 306.5625, 308.4375, 310.3125, 312.1875,  &
    314.0625, 315.9375, 317.8125, 319.6875, 321.5625, 323.4375, 325.3125,  &
    327.1875, 329.0625, 330.9375, 332.8125, 334.6875, 336.5625, 338.4375,  &
    340.3125, 342.1875, 344.0625, 345.9375, 347.8125, 349.6875, 351.5625,  &
    353.4375, 355.3125, 357.1875, 359.0625 /)*rad

  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do
end subroutine rotate_hadgem2_wind
!
!---------------------------------------------------------------------------------------------------
