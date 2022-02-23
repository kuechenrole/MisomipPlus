module g_forcing_param_notused
  implicit none
  save   

  ! *** exchange coefficients ***
  real*8    :: Ce_atm_oce=1.5e-3 ! exchange coeff. of latent heat over open water
  real*8    :: Ch_atm_oce=1.2e-3 ! exchange coeff. of sensible heat over open water
  real*8    :: Cd_atm_oce=1.3e-3 ! drag coefficient between atmosphere and water

  real*8    :: Ce_atm_ice=1.5e-3 ! exchange coeff. of latent heat over ice
  real*8    :: Ch_atm_ice=1.2e-3 ! exchange coeff. of sensible heat over ice
  real*8    :: Cd_atm_ice=1.3e-3 ! drag coefficient between atmosphere and ice !1.63

  namelist /forcing_exchange_coeff/ Ce_atm_oce, Ch_atm_oce, Cd_atm_oce, &
       Ce_atm_ice, Ch_atm_ice, Cd_atm_ice


  ! *** forcing source and type ***
  character(10)                 :: wind_data_source
  character(10)                 :: rad_data_source
  character(10)                 :: precip_data_source
  character(10)                 :: runoff_data_source
  character(10)                 :: sss_data_source
  integer                       :: wind_ttp_ind
  integer                       :: rad_ttp_ind
  integer                       :: precip_ttp_ind
  integer                       :: runoff_ttp_ind
  integer                       :: sss_ttp_ind


  namelist /forcing_source/ wind_data_source, rad_data_source, precip_data_source, &
       runoff_data_source, sss_data_source, wind_ttp_ind, rad_ttp_ind, precip_ttp_ind, &
       runoff_ttp_ind, sss_ttp_ind, ncar_bulk_formulae

  ! *** coefficients in bulk formulae ***
  logical                       :: AOMIP_drag_coeff=.false.
  logical                       :: ncar_bulk_formulae=.false.

  namelist /forcing_bulk/ AOMIP_drag_coeff, ncar_bulk_formulae

 ! *** add land ice melt water ***
  logical                       :: use_landice_water=.false.
  integer                       :: landice_start_mon=1
  integer                       :: landice_end_mon=12

  namelist /land_ice/ use_landice_water, landice_start_mon, landice_end_mon


end module g_forcing_param_notused
!
!----------------------------------------------------------------------------
!
module g_forcing_arrays
  implicit none
  save    

  ! forcing arrays
  real(kind=8), allocatable, dimension(:)         :: u_wind, v_wind, Pair
  real(kind=8), allocatable, dimension(:)         :: Tair, shum, tdew
  real(kind=8), allocatable, dimension(:)         :: shortwave, longwave
  real(kind=8), allocatable, dimension(:)         :: prec_rain, prec_snow
  real(kind=8), allocatable, dimension(:)         :: runoff, evaporation
  real(kind=8), allocatable, dimension(:,:)       :: u_wind_t, v_wind_t, tdew_t
  real(kind=8), allocatable, dimension(:,:)       :: Tair_t, shum_t, evaporation_t
  real(kind=8), allocatable, dimension(:,:)       :: shortwave_t, longwave_t, e_vapor_t
  real(kind=8), allocatable, dimension(:,:)       :: shortwave_e, longwave_e
  real(kind=8), allocatable, dimension(:,:)       :: prec_rain_t, prec_snow_t, prec_net_t
  real(kind=8), allocatable, dimension(:,:)	  :: runoff_t, Pair_t, cloud_t

  real(kind=8), allocatable, dimension(:)         :: runoff_landice
  real(kind=8)                                    :: landice_season(12)


  ! shortwave penetration
  real(kind=8), allocatable, dimension(:)         :: chl, sw_3d
  real(kind=8), allocatable, dimension(:,:)       :: chl_t

  real(kind=8), allocatable, dimension(:)         :: thdgr, thdgrsn, flice
  real(kind=8), allocatable, dimension(:)         :: olat_heat, osen_heat, olwout
  real(kind=8), allocatable, dimension(:)         :: fwat_down_only

  ! drag coefficient Cd_atm_oce and transfer coefficients for evaporation
  ! Ce_atm_oce and sensible heat Ch_atm_oce between atmosphere and ocean
  real(kind=8), allocatable, dimension(:)	  :: Cd_atm_oce_arr
  real(kind=8), allocatable, dimension(:)	  :: Ch_atm_oce_arr
  real(kind=8), allocatable, dimension(:)	  :: Ce_atm_oce_arr

  ! drag coefficient Cd_atm_oce between atmosphere and ice
  real(kind=8), allocatable, dimension(:)	  :: Cd_atm_ice_arr

end module g_forcing_arrays
!
!----------------------------------------------------------------------------
!
module g_forcing_index
  
  use g_parfe
  
  implicit none
  save

  ! arrays for temporal interpolation
  integer                                         :: update_forcing_flag(6)
  integer                                         :: forcing_rec(6)
  real(kind=8)                                    :: interp_coef(6)

contains

  subroutine forcing_index
    use g_clock
    implicit none

    real(kind=8)          :: sixhour_sec, onehour_sec
    real(kind=8)          :: oneday_sec
    real(kind=8)          :: modtimeold

    data sixhour_sec /21600.0/, oneday_sec /86400.0/, onehour_sec /3600.0/

    modtimeold=mod(timeold,oneday_sec)

    ! if update forcing or not
    update_forcing_flag=0
    if(mod(timeold, sixhour_sec)==0.0)        update_forcing_flag(1)=1
    if(modtimeold==0.0)                       update_forcing_flag(2)=1
    if(day_in_month==1 .and. modtimeold==0.0) update_forcing_flag(3:4)=1
    if(modtimeold==0.0)                       update_forcing_flag(5)=1
    if(mod(timeold, onehour_sec)==0.0)        update_forcing_flag(6)=1
    
    ! which record should be used as the first one in interpolation
    forcing_rec(1) = 1+int(modtimeold/sixhour_sec)+4*(daynew-1)
    forcing_rec(2) = daynew
    forcing_rec(3) = month
    forcing_rec(4) = month
    forcing_rec(5) = daynew
    forcing_rec(6) = 1+int(modtimeold/onehour_sec)+24*(daynew-1)
    
    !if(mype==0) write(*,*) 'forcing_rec ',forcing_rec
    !if(mype==0) write(*,*) modtimeold,daynew,istep

    ! interpolation coefficients
    interp_coef(1)=mod(timeold, sixhour_sec)/sixhour_sec
    interp_coef(2)=modtimeold/oneday_sec
    interp_coef(3)=(day_in_month-1.0+modtimeold/oneday_sec) &
                   /real(num_day_in_month(fleapyear,month))
    interp_coef(4)=interp_coef(3)
    interp_coef(5)=interp_coef(2)
    interp_coef(6)=mod(timeold, onehour_sec)/onehour_sec


    if(any(interp_coef>1.) .or. any(interp_coef<0.)) then
       write(*,*) 'error in interp_coef'
       stop
    end if

  end subroutine forcing_index

end module g_forcing_index
!
!----------------------------------------------------------------------------
!
module g_forcing_interp
  !This module prepare the weights for interpolating 
  !forcing data Tair, humidity, wind velocities,
  !precipitation, radiations, etc.
  !Based on assumption forcing data are on the T62 NCEP/NCAR grid

    use o_MESH
    use o_PARAM
    use g_rotate_grid
    use g_parfe

    implicit none

  integer, allocatable      :: lint_ind(:,:,:)
  real(kind=8), allocatable :: lint_weight(:,:)
  real(kind=8),allocatable  	:: xmod(:), ymod(:)


contains

  !------------------------------------------------------------------
  subroutine  init_forcing_interp

    ! routine to calculate neighbours and weights for linear interpolation
    !
    ! required information
    !    xmod(nmp)  longitudes of model point on geographical grid in Bogenmass
    !    ymod(nmp)  latitudes of model points on geographical grid in Bogenmass
    !         nmp   number of model points where data are needed
    !    cx(ni)     longitudes of data points on regular geographical grid
    !               by now must be in range[0,..360] in ascending order
    !    cy(nj)     latitudes of data points on regular geographical grid 
    !               by now must be in range [-90,..90] in ascending order
    !
    ! OUTPUT
    !    lint_ind(4,2,nmp)   /i,j/ adress for four neighbors of each model node
    !    lint_weight(4,nmp)  interpolation weights
  
  use g_forcing_param

    integer             :: ni, nj, nic, njc  ! NCEP and CORE are on the same grid. BUT GME IS NOT. und erst recht nicht COSMO.  VH 
                                             !                                     Und HadCM3 auch nicht.  RT
    integer     	:: i, ii, j, n, row, n2
    real(kind=8)      	:: rlon, rlat, aux
    real(kind=8), allocatable, dimension(:)           :: cx, cy, cxc, cyc
    real(kind=8)        :: wt(4)

    n2=myDim_nod2D+eDim_nod2D     
    
    if (wind_data_source=='COSMO') then
      nic=600
      njc=180 ! COSMO grid                  !   baustelle
      allocate(cxc(nic), cyc(njc))
    elseif (wind_data_source=='GME' .or. wind_data_source=='COSMO') then
       ni=720
       nj=361
    elseif (wind_data_source=='NCEP' .or. wind_data_source(1:4)=='CORE') then
       ni=192
       nj=94
    elseif (wind_data_source(1:6)=='HadCM3') then
       ni=96
       nj=72
    elseif (wind_data_source(1:6)=='MPEH5C') then
       ni=96
       nj=48
    elseif (wind_data_source(1:6)=='ERAINT') then
       ni=480
       nj=241
    elseif (wind_data_source(1:4)=='CFSR') then
       ni=1152
       nj=576
    else
       write(*,*)'STOP: undefined wind_data_source'
       stop
    endif
    
    
    ! define grids
    
    allocate(cx(ni), cy(nj))
    allocate(lint_ind(4,2,n2))   
    allocate(lint_weight(4,n2))  
    allocate(xmod(myDim_nod2D+eDim_nod2D), ymod(myDim_nod2D+eDim_nod2D))

    if (wind_data_source=='COSMO') then
      ! COSMO longitude
      cxc(1)=285.0
      do i=2,nic
         cxc(i)=cxc(i-1)+0.125
      enddo
                                                ! baustelle
      ! COSMO latitude
      cyc(1)=-80.0
      do i=2,njc
         cyc(i)=cyc(i-1)+0.125
      enddo  
       
        ! CHECK
      if (cyc(njc) .ne. -57.625 .or. cxc(nic) .ne. 359.875) then
        write(*,*) 'ERROR: COSMO lon and/or lat calculated wrongly in init_forcing_interp!'
	call abort
      endif
      cxc=cxc*rad
      cyc=cyc*rad

    elseif (wind_data_source=='GME' .or. wind_data_source=='COSMO') then

      ! GME longitude
      cx(1)=0.0
      do i=2,ni
         cx(i)=cx(i-1)+0.5
      enddo
    
      ! GME latitude
      cy(1)=-90.0
      do i=2,nj
         cy(i)=cy(i-1)+0.5
      enddo   
     
    elseif (wind_data_source=='NCEP' .or. wind_data_source(1:4)=='CORE') then
    
      ! NCEP/CORE latitude
      cy =(/-88.542, -86.6531, -84.7532, -82.8508, -80.9473, -79.0435, &  
         -77.1394, -75.2351, -73.3307, -71.4262, -69.5217, -67.6171, &  
         -65.7125, -63.8079, -61.9033, -59.9986, -58.0939, -56.1893, &  
         -54.2846, -52.3799, -50.4752, -48.5705, -46.6658, -44.7611,&  
         -42.8564, -40.9517, -39.0470, -37.1422, -35.2375, -33.3328, &  
         -31.4281, -29.5234, -27.6186, -25.7139, -23.8092, -21.9044, &  
         -19.9997, -18.0950, -16.1902, -14.2855, -12.3808, -10.47604, &  
         -8.57131, -6.66657, -4.76184, -2.8571, -0.952368, 0.952368, &  
         2.8571, 4.76184, 6.66657, 8.57131, 10.47604, 12.3808, &  
         14.2855, 16.1902, 18.095, 19.9997, 21.9044, 23.8092, &  
         25.7139, 27.6186, 29.5234, 31.4281, 33.3328, 35.2375,&  
         37.1422, 39.047,  40.9517, 42.8564, 44.7611, 46.6658,&  
         48.5705, 50.4752, 52.3799, 54.2846, 56.1893, 58.0939,&  
         59.9986, 61.9033, 63.8079, 65.7125, 67.6171, 69.5217, &  
         71.4262, 73.3307, 75.2351, 77.1394, 79.0435, 80.9473, &  
         82.8508, 84.7532, 86.6531, 88.542 /)

      ! NCEP/CORE longitude
      cx(1)=0.0
      do i=2,ni
         cx(i)=cx(i-1)+1.875
      enddo
     
     elseif (wind_data_source(1:4)=='CFSR' ) then

   cy=(/-89.7609978067649, -89.4513839813147, -89.1399444898031, -88.8280918544348, -88.5160824601939, &
        -88.2039972821210, -87.8918697257972, -87.5797161067180, -87.2675453237151, -86.9553626407769, &
        -86.6431713700617, -86.3309736994237, -86.0187711319782, -85.7065647344594, -85.3943552848229, &
        -85.0821433637477, -84.7699294134275, -84.4577137765236, -84.1454967226726, -83.8332784669523, &
        -83.5210591830149, -83.2088390126008, -82.8966180725465, -82.5843964600232, -82.2721742565071, &
        -81.9599515308245, -81.6477283415154, -81.3355047386851, -81.0232807654711, -80.7110564592133, &
        -80.3988318523959, -80.0866069734125, -79.7743818471896, -79.4621564957006, -79.1499309383903, &
        -78.8377051925286, -78.5254792735062, -78.2132531950833, -77.9010269696002, -77.5888006081557, &
        -77.2765741207594, -76.9643475164622, -76.6521208034683, -76.3398939892319, -76.0276670805412, &
        -75.7154400835910, -75.4032130040461, -75.0909858470969, -74.7787586175082, -74.4665313196615, &
        -74.1543039575933, -73.8420765350281, -73.5298490554079, -73.2176215219189, -72.9053939375144, &
        -72.5931663049358, -72.2809386267315, -71.9687109052734, -71.6564831427719, -71.3442553412894, &
        -71.0320275027528, -70.7197996289641, -70.4075717216104, -70.0953437822730, -69.7831158124356, &
        -69.4708878134917, -69.1586597867514, -68.8464317334473, -68.5342036547404, -68.2219755517253, &
        -67.9097474254345, -67.5975192768432, -67.2852911068728, -66.9730629163950, -66.6608347062347, &
        -66.3486064771734, -66.0363782299518, -65.7241499652726, -65.4119216838094, -65.0996933861819, &
	-64.7874650729996, -64.4752367448354, -64.1630084022344, -63.8507800457158, -63.5385516757741, &
	-63.2263232928810, -62.9140948974863, -62.6018664900196, -62.2896380708909, -61.9774096404920, &
	-61.6651811991977, -61.3529527473660, -61.0407242853399, -60.7284958134478, -60.4162673320040, & 
	-60.1040388413099, -59.7918103416544, -59.4795818333146, -59.1673533165564, -58.8551247916351, &
	-58.5428962587959, -58.2306677182743, -57.9184391702968, -57.6062106150810, -57.2939820528364, &
	-56.9817534837646, -56.6695249080595, -56.3572963259080, -56.0450677374901, -55.7328391429793, &
	-55.4206105425426, -55.1083819363415, -54.7961533245315, -54.4839247072627, -54.1716960846800, &
	-53.8594674569233, -53.5472388241279, -53.2350101864244, -52.9227815439389, -52.6105528967935, &
	-52.2983242451062, -51.9860955889913, -51.6738669285591, -51.3616382639167, -51.0494095951675, &
	-50.7371809224118, -50.4249522457468, -50.1127235652666, -49.8004948810624, -49.4882661932226, &
	-49.1760375018330, -48.8638088069768, -48.5515801087346, -48.2393514071848, -47.9271227024034, &
	-47.6148939944642, -47.3026652834388, -46.9904365693970, -46.6782078524062, -46.3659791325323, &	
	-46.0537504098392, -45.7415216843889, -45.4292929562419, -45.1170642254568, -44.8048354920910, &
	-44.4926067561998, -44.1803780178376, -43.8681492770570, -43.5559205339092, -43.2436917884443, &
	-42.9314630407108, -42.6192342907563, -42.3070055386268, -41.9947767843674, -41.6825480280221, &
	-41.3703192696335, -41.0580905092434, -40.7458617468925, -40.4336329826204, -40.1214042164660, &
	-39.8091754484670, -39.4969466786603, -39.1847179070820, -38.8724891337671, -38.5602603587501, &
	-38.2480315820645, -37.9358028037430, -37.6235740238177, -37.3113452423198, -36.9991164592800, &
	-36.6868876747282, -36.3746588886935, -36.0624301012046, -35.7502013122894, -35.4379725219752, &
	-35.1257437302889, -34.8135149372567, -34.5012861429041, -34.1890573472562, -33.8768285503376, &
	-33.5645997521724, -33.2523709527841, -32.9401421521959, -32.6279133504301, -32.3156845475092, &
	-32.0034557434547, -31.6912269382879, -31.3789981320297, -31.0667693247005, -30.7545405163204, &
	-30.4423117069090, -30.1300828964857, -29.8178540850693, -29.5056252726785, -29.1933964593314, &
	-28.8811676450460, -28.5689388298397, -28.2567100137299, -27.9444811967335, -27.6322523788671, &	
	-27.3200235601471, -27.0077947405894, -26.6955659202099, -26.3833370990241, -26.0711082770472, &
	-25.7588794542941, -25.4466506307797, -25.1344218065184, -24.8221929815244, -24.5099641558118, &
	-24.1977353293943, -23.8855065022855, -23.5732776744988, -23.2610488460473, -22.9488200169440, &
	-22.6365911872016, -22.3243623568326, -22.0121335258495, -21.6999046942644, -21.3876758620893, &
	-21.0754470293361, -20.7632181960163, -20.4509893621416, -20.1387605277232, -19.8265316927724, &
	-19.5143028573001, -19.2020740213172, -18.8898451848345, -18.5776163478626, -18.2653875104119, &
	-17.9531586724928, -17.6409298341154, -17.3287009952899, -17.0164721560262, -16.7042433163340, &
	-16.3920144762232, -16.0797856357033, -15.7675567947838, -15.4553279534741, -15.1430991117835, &
	-14.8308702697211, -14.5186414272961, -14.2064125845175, -13.8941837413940, -13.5819548979346, &
	-13.2697260541479, -12.9574972100426, -12.6452683656272, -12.3330395209102, -12.0208106759000, &
	-11.7085818306049, -11.3963529850331, -11.0841241391930, -10.7718952930924, -10.4596664467395, &
	-10.1474376001423, -9.83520875330866, -9.52297990624647, -9.21075105896353, -8.89852221146760, &
	-8.58629336376639, -8.27406451586750, -7.96183566777854, -7.64960681950711, -7.33737797106064, &
	-7.02514912244664, -6.71292027367251, -6.40069142474572, -6.08846257567353, -5.77623372646320, &
	-5.46400487712210, -5.15177602765742, -4.83954717807647, -4.52731832838632, -4.21508947859429, &
	-3.90286062870728, -3.59063177873268, -3.27840292867733, -2.96617407854850, -2.65394522835317, &
	-2.34171637809846, -2.02948752779123, -1.71725867743865, -1.40502982704764, -1.09280097662518, &
	-0.780572126178202, -0.468343275713578, -0.156114425239382, 0.156114425239382, 0.468343275713578, &
	0.780572126178202, 1.09280097662518, 1.40502982704764, 1.71725867743865, 2.02948752779123, &
	2.34171637809846, 2.65394522835317, 2.96617407854850, 3.27840292867733, 3.59063177873268, &
	3.90286062870728, 4.21508947859429, 4.52731832838632, 4.83954717807647, 5.15177602765742, &
	5.46400487712210, 5.77623372646320, 6.08846257567353, 6.40069142474572, 6.71292027367251, &
	7.02514912244664, 7.33737797106064, 7.64960681950711, 7.96183566777854, 8.27406451586750, &
	8.58629336376639, 8.89852221146760, 9.21075105896353, 9.52297990624647, 9.83520875330866, &
	10.1474376001423, 10.4596664467395, 10.7718952930924, 11.0841241391930, 11.3963529850331, &
	11.7085818306049, 12.0208106759000, 12.3330395209102, 12.6452683656272, 12.9574972100426, &
	13.2697260541479, 13.5819548979346, 13.8941837413940, 14.2064125845175, 14.5186414272961, &
	14.8308702697211, 15.1430991117835, 15.4553279534741, 15.7675567947838, 16.0797856357033, &
	16.3920144762232, 16.7042433163340, 17.0164721560262, 17.3287009952899, 17.6409298341154, &
	17.9531586724928, 18.2653875104119, 18.5776163478626, 18.8898451848345, 19.2020740213172, &
	19.5143028573001, 19.8265316927724, 20.1387605277232, 20.4509893621416, 20.7632181960163, &
	21.0754470293361, 21.3876758620893, 21.6999046942644, 22.0121335258495, 22.3243623568326, &
	22.6365911872016, 22.9488200169440, 23.2610488460473, 23.5732776744988, 23.8855065022855, &
	24.1977353293943, 24.5099641558118, 24.8221929815244, 25.1344218065184, 25.4466506307797, &
	25.7588794542941, 26.0711082770472, 26.3833370990241, 26.6955659202099, 27.0077947405894, &
	27.3200235601471, 27.6322523788671, 27.9444811967335, 28.2567100137299, 28.5689388298397, &
	28.8811676450460, 29.1933964593314, 29.5056252726785, 29.8178540850693, 30.1300828964857, &
	30.4423117069090, 30.7545405163204, 31.0667693247005, 31.3789981320297, 31.6912269382879, &
	32.0034557434547, 32.3156845475092, 32.6279133504301, 32.9401421521959, 33.2523709527841, &
	33.5645997521724, 33.8768285503376, 34.1890573472562, 34.5012861429041, 34.8135149372567, &
	35.1257437302889, 35.4379725219752, 35.7502013122894, 36.0624301012046, 36.3746588886935, & 
	36.6868876747282, 36.9991164592800, 37.3113452423198, 37.6235740238177, 37.9358028037430, &
	38.2480315820645, 38.5602603587501, 38.8724891337671, 39.1847179070820, 39.4969466786603, &
	39.8091754484670, 40.1214042164660, 40.4336329826204, 40.7458617468925, 41.0580905092434, &
	41.3703192696335, 41.6825480280221, 41.9947767843674, 42.3070055386268, 42.6192342907563, &
	42.9314630407108, 43.2436917884443, 43.5559205339092, 43.8681492770570, 44.1803780178376, &
	44.4926067561998, 44.8048354920910, 45.1170642254568, 45.4292929562419, 45.7415216843889, &
	46.0537504098392, 46.3659791325323, 46.6782078524062, 46.9904365693970, 47.3026652834388, &
	47.6148939944642, 47.9271227024034, 48.2393514071848, 48.5515801087346, 48.8638088069768, &
	49.1760375018330, 49.4882661932226, 49.8004948810624, 50.1127235652666, 50.4249522457468, &
	50.7371809224118, 51.0494095951675, 51.3616382639167, 51.6738669285591, 51.9860955889913, &
	52.2983242451062, 52.6105528967935, 52.9227815439389, 53.2350101864244, 53.5472388241279, &
	53.8594674569233, 54.1716960846800, 54.4839247072627, 54.7961533245315, 55.1083819363415, &
	55.4206105425426, 55.7328391429793, 56.0450677374901, 56.3572963259080, 56.6695249080595, &
	56.9817534837646, 57.2939820528364, 57.6062106150810, 57.9184391702968, 58.2306677182743, &
	58.5428962587959, 58.8551247916351, 59.1673533165564, 59.4795818333146, 59.7918103416544, &
	60.1040388413099, 60.4162673320040, 60.7284958134478, 61.0407242853399, 61.3529527473660, &
	61.6651811991977, 61.9774096404920, 62.2896380708909, 62.6018664900196, 62.9140948974863, &
	63.2263232928810, 63.5385516757741, 63.8507800457158, 64.1630084022344, 64.4752367448354, &
	64.7874650729996, 65.0996933861819, 65.4119216838094, 65.7241499652726, 66.0363782299518, & 
	66.3486064771734, 66.6608347062347, 66.9730629163950, 67.2852911068728, 67.5975192768432, & 
	67.9097474254345, 68.2219755517253, 68.5342036547404, 68.8464317334473, 69.1586597867514, &
	69.4708878134917, 69.7831158124356, 70.0953437822730, 70.4075717216104, 70.7197996289641, &
	71.0320275027528, 71.3442553412894, 71.6564831427719, 71.9687109052734, 72.2809386267315, &
	72.5931663049358, 72.9053939375144, 73.2176215219189, 73.5298490554079, 73.8420765350281, &
	74.1543039575933, 74.4665313196615, 74.7787586175082, 75.0909858470969, 75.4032130040461, &
	75.7154400835910, 76.0276670805412, 76.3398939892319, 76.6521208034683, 76.9643475164622, &
	77.2765741207594, 77.5888006081557, 77.9010269696002, 78.2132531950833, 78.5254792735062, &
	78.8377051925286, 79.1499309383903, 79.4621564957006, 79.7743818471896, 80.0866069734125, &
	80.3988318523959, 80.7110564592133, 81.0232807654711, 81.3355047386851, 81.6477283415154, &
	81.9599515308245, 82.2721742565071, 82.5843964600232, 82.8966180725465, 83.2088390126008, &
	83.5210591830149, 83.8332784669523, 84.1454967226726, 84.4577137765236, 84.7699294134275, &
	85.0821433637477, 85.3943552848229, 85.7065647344594, 86.0187711319782, 86.3309736994237, &
	86.6431713700617, 86.9553626407769, 87.2675453237151, 87.5797161067180, 87.8918697257972, &
	88.2039972821210, 88.5160824601939, 88.8280918544348, 89.1399444898031, 89.4513839813147, &
	89.7609978067649 /)
    ! NCEP/CORE longitude
  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+0.3125
  enddo


    elseif (wind_data_source(1:6)=='HadCM3') then

! HadCM3 longitude
   cx =(/ 1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, &
    31.875, 35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, &
    65.625, 69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, &
    99.375, 103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, &
    129.375, 133.125, 136.875, 140.625, 144.375, 148.125, 151.875, 155.625, &
    159.375, 163.125, 166.875, 170.625, 174.375, 178.125, 181.875, 185.625, &
    189.375, 193.125, 196.875, 200.625, 204.375, 208.125, 211.875, 215.625, &
    219.375, 223.125, 226.875, 230.625, 234.375, 238.125, 241.875, 245.625, &
    249.375, 253.125, 256.875, 260.625, 264.375, 268.125, 271.875, 275.625, &
    279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875, 305.625, &
    309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625, &
    339.375, 343.125, 346.875, 350.625, 354.375, 358.125 /)

 ! HadCM3 latitude
   cy=(/ 88.75, 86.25, 83.75, 81.25, 78.75, 76.25, 73.75, 71.25, 68.75, &
    66.25, 63.75, 61.25, 58.75, 56.25, 53.75, 51.25, 48.75, 46.25, 43.75, & 
    41.25, 38.75, 36.25, 33.75, 31.25, 28.75, 26.25, 23.75, 21.25, 18.75,  &
    16.25, 13.75, 11.25, 8.75, 6.25, 3.75, 1.25, -1.25, -3.75, -6.25, -8.75, & 
    -11.25, -13.75, -16.25, -18.75, -21.25, -23.75, -26.25, -28.75, -31.25,  &
    -33.75, -36.25, -38.75, -41.25, -43.75, -46.25, -48.75, -51.25, -53.75,  &
    -56.25, -58.75, -61.25, -63.75, -66.25, -68.75, -71.25, -73.75, -76.25, &
    -78.75, -81.25, -83.75, -86.25, -88.75 /)*(-1.)

    elseif (wind_data_source(1:6)=='MPEH5C') then

 ! MPEH5C longitude
 cx = (/ 0., 3.75, 7.5, 11.25, 15., 18.75, 22.5, 26.25, 30., 33.75, 37.5, 41.25, &
    45., 48.75, 52.5, 56.25, 60., 63.75, 67.5, 71.25, 75., 78.75, 82.5, 86.25, &
    90., 93.75, 97.5, 101.25, 105., 108.75, 112.5, 116.25, 120., 123.75, 127.5, &
    131.25, 135., 138.75, 142.5, 146.25, 150., 153.75, 157.5, 161.25, 165., &
    168.75, 172.5, 176.25, 180., 183.75, 187.5, 191.25, 195., 198.75, 202.5, &
    206.25, 210., 213.75, 217.5, 221.25, 225., 228.75, 232.5, 236.25, 240., &
    243.75, 247.5, 251.25, 255., 258.75, 262.5, 266.25, 270., 273.75, 277.5, &
    281.25, 285., 288.75, 292.5, 296.25, 300., 303.75, 307.5, 311.25, 315., &
    318.75, 322.5, 326.25, 330., 333.75, 337.5, 341.25, 345., 348.75, 352.5, &
    356.25 /)


! MPEH5C latitude
 cy = (/ 87.1590945558629, 83.4789366693172, 79.7770456548256, &
    76.0702444625451, 72.3615810293448, 68.6520167895175, 64.9419494887575, &
    61.2315731880771, 57.52099379797, 53.8102740319414, 50.0994534129868, &
    46.3885581116054, 42.6776061726049, 38.966610469454, 35.2555804613682, &
    31.5445232840217, 27.8334444519932, 24.122348326088, 20.4112384335678, &
    16.7001176938427, 12.9889885820881, 9.27785325150786, 5.56671362791359, &
    1.85557148599326, -1.85557148599326, -5.56671362791359, &
    -9.27785325150786, -12.9889885820881, -16.7001176938427, &
    -20.4112384335678, -24.122348326088, -27.8334444519932, &
    -31.5445232840217, -35.2555804613682, -38.966610469454, &
    -42.6776061726049, -46.3885581116054, -50.0994534129868, &
    -53.8102740319414, -57.52099379797, -61.2315731880771, -64.9419494887575, &
    -68.6520167895175, -72.3615810293448, -76.0702444625451, &
    -79.7770456548256, -83.4789366693172, -87.1590945558629 /) * (-1.)

    elseif (wind_data_source(1:6)=='ERAINT') then
       ! ERAINT

  cy(1)=-90.0
  do i=2,nj
    cy(i)=cy(i-1)+0.75
  end do

  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+0.75
  enddo

    else
     write(*,*)'STOP: Please define the lat/lon grid of your forcing data'
     write(*,*)'Remember they should be degrees (not rad) here!'
     stop
    endif

    ! some checks, range of cx and cy
    if(cx(ni)-cx(1).gt.360.) then
       write(*,*) 'gauss_init: x-interval gt 360'
       call abort
    endif
    if(cy(nj)-cy(1).gt.180.) then
       write(*,*) 'gauss_init: y-interval gt 180'
       call abort
    endif
    if(cx(1).ge.360.0) then
       aux=int(cx(1)/360.)*360.
       do i=1,ni
          cx(i)=cx(i)-aux
       enddo
    endif
    if(cx(ni).gt.360.) then
     write(*,*)'STOP: cx(ni).gt.360'
     stop
    endif

    ! in the following we need cx and cy in unit radian
    cx=cx*rad
    cy=cy*rad

    ! model grid coordinate (in radian, between 0 and 2*pi)
    do row=1,n2                    
       if(rotated_grid) then
          rlon=coord_nod2D(1,row)
          rlat=coord_nod2D(2,row)
          call r2g(xmod(row), ymod(row), rlon, rlat)
          if (xmod(row)<0.0) xmod(row)=2*pi+xmod(row)
       else
          xmod(row)=coord_nod2D(1,row)
          ymod(row)=coord_nod2D(2,row)
          if (xmod(row)<0.0) xmod(row)=2*pi+xmod(row)	
       endif
    enddo

 !    baustelle ###########################################   Gebietsabfrage??

    ! linear interpolation: nodes and weight
    
  if (wind_data_source=='COSMO') then
      
    do row=1,n2
     !  if innerhalb gebiet                            !    baustelle 
      if (xmod(row) > 285.*rad .and. xmod(row) < 360.*rad .and. ymod(row) < -57.5*rad) then 
       if(xmod(row)<cxc(nic)) then
          do i=1,nic
             if(xmod(row)<cxc(i)) then
                lint_ind(1,1,row)=i-1
                lint_ind(2,1,row)=i-1
                lint_ind(3,1,row)=i
                lint_ind(4,1,row)=i
                aux=(cxc(i)-xmod(row))/(cxc(i)-cxc(i-1))
                exit
             end if
          end do
       else
          lint_ind(1,1,row)=nic
          lint_ind(2,1,row)=nic
          lint_ind(3,1,row)=1
          lint_ind(4,1,row)=1
          aux=(360.0_8-xmod(row))/(360.0_8-cxc(nic))
       end if
       wt(1)=aux
       wt(2)=aux
       wt(3)=1.0_8-aux
       wt(4)=1.0_8-aux

       if(ymod(row)<cyc(njc)) then
          do j=1,njc
             if(ymod(row)<cyc(j)) then
                lint_ind(1,2,row)=j-1
                lint_ind(2,2,row)=j
                lint_ind(3,2,row)=j
                lint_ind(4,2,row)=j-1
                aux=(cyc(j)-ymod(row))/(cyc(j)-cyc(j-1))
                exit
             end if
          end do
       else
          lint_ind(1,2,row)=njc
          lint_ind(2,2,row)=njc
          lint_ind(3,2,row)=njc
          lint_ind(4,2,row)=njc
          aux=1.0_8
       end if
       lint_weight(1,row)=wt(1)*aux
       lint_weight(2,row)=wt(2)*(1.0_8-aux)
       lint_weight(3,row)=wt(3)*(1.0_8-aux)
       lint_weight(4,row)=wt(4)*aux
       
      else 
      
       if(xmod(row)<cx(ni)) then
          do i=1,ni
             if(xmod(row)<cx(i)) then
                lint_ind(1,1,row)=i-1
                lint_ind(2,1,row)=i-1
                lint_ind(3,1,row)=i
                lint_ind(4,1,row)=i
                aux=(cx(i)-xmod(row))/(cx(i)-cx(i-1))
                exit
             end if
          end do
       else
          lint_ind(1,1,row)=ni
          lint_ind(2,1,row)=ni
          lint_ind(3,1,row)=1
          lint_ind(4,1,row)=1
          aux=(360.0_8-xmod(row))/(360.0_8-cx(ni))
       end if
       wt(1)=aux
       wt(2)=aux
       wt(3)=1.0_8-aux
       wt(4)=1.0_8-aux

       if(ymod(row)<cy(nj)) then
          do j=1,nj
             if(ymod(row)<cy(j)) then
                lint_ind(1,2,row)=j-1
                lint_ind(2,2,row)=j
                lint_ind(3,2,row)=j
                lint_ind(4,2,row)=j-1
                aux=(cy(j)-ymod(row))/(cy(j)-cy(j-1))
                exit
             end if
          end do
       else
          lint_ind(1,2,row)=nj
          lint_ind(2,2,row)=nj
          lint_ind(3,2,row)=nj
          lint_ind(4,2,row)=nj
          aux=1.0_8
       end if
       lint_weight(1,row)=wt(1)*aux
       lint_weight(2,row)=wt(2)*(1.0_8-aux)
       lint_weight(3,row)=wt(3)*(1.0_8-aux)
       lint_weight(4,row)=wt(4)*aux
       
      end if
    end do
    
  else   ! if COSMO 
  
    do row=1,n2
    
      if(xmod(row)<cx(ni)) then
          do i=1,ni
             if(xmod(row)<cx(i)) then
                lint_ind(1,1,row)=i-1
                lint_ind(2,1,row)=i-1
                lint_ind(3,1,row)=i
                lint_ind(4,1,row)=i
                aux=(cx(i)-xmod(row))/(cx(i)-cx(i-1))
                exit
             end if
          end do
       else
          lint_ind(1,1,row)=ni
          lint_ind(2,1,row)=ni
          lint_ind(3,1,row)=1
          lint_ind(4,1,row)=1
          aux=(360.0_8-xmod(row))/(360.0_8-cx(ni))
       end if
       wt(1)=aux
       wt(2)=aux
       wt(3)=1.0_8-aux
       wt(4)=1.0_8-aux

       if(ymod(row)<cy(nj)) then
          do j=1,nj
             if(ymod(row)<cy(j)) then
                lint_ind(1,2,row)=j-1
                lint_ind(2,2,row)=j
                lint_ind(3,2,row)=j
                lint_ind(4,2,row)=j-1
                aux=(cy(j)-ymod(row))/(cy(j)-cy(j-1))
                exit
             end if
          end do
       else
          lint_ind(1,2,row)=nj
          lint_ind(2,2,row)=nj
          lint_ind(3,2,row)=nj
          lint_ind(4,2,row)=nj
          aux=1.0_8
       end if
       lint_weight(1,row)=wt(1)*aux
       lint_weight(2,row)=wt(2)*(1.0_8-aux)
       lint_weight(3,row)=wt(3)*(1.0_8-aux)
       lint_weight(4,row)=wt(4)*aux
      
    end do
    
  endif  ! if COSMO  else
    
    

    if(mype==0)  write(*,*) 'weights for interpolating forcing / 2D fields computed'
    return     
  end subroutine init_forcing_interp


  !------------------------------------------------------------------------
  subroutine forcing_linear_ip(zd,idim,jdim,ind,weights,zi,nmpt)
    !  this subroutine interpolates data using prepared weights- 
    !  see subroutine init_forcing_interp
    !
    !  INPUT
    !        zd(idim,jdim)         available data set
    !        nmpt                  number of model positions, where data are wanted
    !        indx(4,2,nmpt)        i,j index of four neighbors of each node
    !        weights(4,nmpt)       interpolation weights
    !        
    !  OUTPUT
    !        zi(nmpt)              array of interpolated values for model points

    use g_parfe
    implicit none                                             

    integer      :: idim, jdim, nmpt                          
    integer      :: ind(4,2,nmpt)
    integer      :: i, n                      
    real(kind=8) :: zd(idim,jdim), zi(nmpt)
    real(kind=8) :: weights(4,nmpt)
    real(kind=8) :: fx

    do n=1,nmpt
       fx=0.0
       do i=1,4
          fx=fx+weights(i,n)*zd(ind(i,1,n),ind(i,2,n))
       enddo
       zi(n)=fx
    enddo

    return
  end subroutine forcing_linear_ip
  !------------------------------------------------------------------------
  subroutine forcing_linear_ip_COSMO(zd,zd2,idim,jdim,idim2,jdim2,ind,weights,zi,nmpt)
    !  this subroutine interpolates data using prepared weights- 
    !  see subroutine init_forcing_interp
    !
    !  INPUT
    !        zd(idim,jdim)         available data set
    !        zd2(idim2,jdim2)      another available data set
    !        nmpt                  number of model positions, where data are wanted
    !        indx(4,2,nmpt)        i,j index of four neighbors of each node
    !        weights(4,nmpt)       interpolation weights
    !        
    !  OUTPUT
    !        zi(nmpt)              array of interpolated values for model points

    use g_parfe
    implicit none                                             

    integer      :: idim, jdim, idim2, jdim2, nmpt                          
    integer      :: ind(4,2,nmpt)
    integer      :: i, n                      
    real(kind=8) :: zd(idim,jdim),zd2(idim2,jdim2), zi(nmpt)
    real(kind=8) :: weights(4,nmpt)
    real(kind=8) :: fx

    do n=1,nmpt
     ! gebietsabfrage
     if (xmod(n) > 285.*rad .and. xmod(n) < 360.*rad .and. ymod(n) < -57.5*rad) then   ! baustelle
       fx=0.0
       do i=1,4
          fx=fx+weights(i,n)*zd2(ind(i,1,n),ind(i,2,n))
       enddo
       zi(n)=fx
     else
       fx=0.0
       do i=1,4
          fx=fx+weights(i,n)*zd(ind(i,1,n),ind(i,2,n))
       enddo
       zi(n)=fx
     endif
    enddo

    return
  end subroutine forcing_linear_ip_COSMO
  !------------------------------------------------------------------------
  subroutine forcing_linear_ip_update(zd,idim,jdim,ind,weights,zi,nmpt)
    !  this subroutine interpolates data using prepared weights- 
    !  see subroutine init_forcing_interp
    !
    !  INPUT
    !        zd(idim,jdim)         available data set
    !        nmpt                  number of model positions, where data are wanted
    !        indx(4,2,nmpt)        i,j index of four neighbors of each node
    !        weights(4,nmpt)       interpolation weights
    !        
    !  OUTPUT
    !        zi(nmpt)              array of interpolated values for model points

    use g_parfe
    implicit none                                             

    integer      :: idim, jdim, nmpt                          
    integer      :: ind(4,2,nmpt)
    integer      :: i, n                      
    real(kind=8) :: zd(idim,jdim), zi(nmpt)
    real(kind=8) :: weights(4,nmpt)
    real(kind=8) :: fx

    do n=1,nmpt
     ! gebietsabfrage
     if (xmod(n) > 285.*rad .and. xmod(n) < 360.*rad .and. ymod(n) < -57.5*rad) then   ! baustelle
       fx=0.0
       do i=1,4
          fx=fx+weights(i,n)*zd(ind(i,1,n),ind(i,2,n))
       enddo
       zi(n)=fx
     endif
    enddo

    return
  end subroutine forcing_linear_ip_update

end module g_forcing_interp
