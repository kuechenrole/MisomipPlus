subroutine init_atm_forcing
! cfsr, 6hourly

  use o_PARAM
  use o_MESH
  use o_array
  use i_therm_parms
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_config
  use g_parfe
  implicit none
  
  integer, parameter			:: nci=1152, ncj=576 !cfsr
  integer                   		:: itime, i, k, n2
  integer                               :: readtype
  character(180)             		:: file
  character(25)             		:: vari
  character(25)             		:: filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy

  n2=myDim_nod2D+eDim_nod2D  

  ! predefinition 6hourly
  
 wind_ttp_ind   = 1
 rad_ttp_ind    = 1
 precip_ttp_ind = 1
 runoff_ttp_ind = 0
 sss_ttp_ind    = 0


  ! compute forcing index
  call forcing_index


  !==========================================================================
  ! wind u and v, Tair, and shum

  if(wind_data_source=='CFSR') then

     if(wind_ttp_ind==1) then !6hourly

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind)
        do i=1,2

           ! 10-m wind m/s ----------------------------------------

	   write(*,*) 'read_U_start'
           vari='U_GRD_L103'
           file=trim(ForcingDataPath)//'CFSR/wnd10m.gdas.'//fileyear//'.nc'
           write(*,*) 'read_U_start2'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           write(*,*) 'read_U_start3'
           call upside_down(array_nc,nci,ncj)
           write(*,*) 'read_U_start4'

           vari='V_GRD_L103'
           file=trim(ForcingDataPath)//'CFSR/wnd10m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)
           write(*,*) 'read_U_start5'

           ! rotate wind
           if(rotated_grid) call rotate_T382_wind(array_nc, array_nc2)
           write(*,*) 'read_U_start6'

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind_t(i,:),n2)   
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind_t(i,:),n2) 
	   write(*,*) 'check_interpolation', array_nc(10,10), u_wind_t(1,10),v_wind_t(1,10)
	    
           ! 2-m temperature --------------------------------------
           write(*,*) 'read_U_start7'
           vari='TMP_L103'
           file=trim(ForcingDataPath)//'CFSR/tmp2mq2m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair_t(i,:),n2)   
           Tair_t(i,:)=Tair_t(i,:)-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           write(*,*) 'read_U_start8'
           vari='SPF_H_L103'
           file=trim(ForcingDataPath)//'CFSR/tmp2mq2m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum_t(i,:),n2)   

           if(rad_ttp_ind==1) then
              vari='PRES_L1'
              file=trim(ForcingDataPath)//'CFSR/pressfc.gdas.'//fileyear//'.nc'
              call read_CFSR_NetCDF(file, vari, itime, array_nc)
              call upside_down(array_nc,nci,ncj)
              call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Pair_t(i,:),n2)
           end if

           if(update_forcing_flag(wind_ttp_ind)==1) then

              ! updating will take place in update_forcing in the first iteration

              u_wind_t(2,:)=u_wind_t(1,:)
              v_wind_t(2,:)=v_wind_t(1,:)
              Tair_t(2,:)=Tair_t(1,:)
              shum_t(2,:)=shum_t(1,:)
	      if(rad_ttp_ind==5) Pair_t(2,:)=Pair_t(1,:)

              exit
           end if

           itime=itime+1
           if (itime>ndpyr) then	! to read the first record of the next year
              itime=1
              write(fileyear,'(i4)') yearnew+1
           endif

        end do

     end if

  else ! other data source

     write(*,*) 'The wind data source ', wind_data_source, ' is not supported.'
     write(*,*) 'Please check and update the code!'
     stop

  endif

  !==========================================================================
  ! precipitation
   
   if(precip_data_source=='CFSR') then
     write(*,*)'in CFSR precip, precip_ttp_ind=',precip_ttp_ind

     if(precip_ttp_ind==1) then  !6hourly

        fileyear=cyearnew
        itime=forcing_rec(precip_ttp_ind)

        do i=1,2

           ! total precip mm/s ------------------------------------

           write(*,*) 'read_U_start9'
           vari='PRATE_L1_Avg_1'
           file=trim(ForcingDataPath)//'CFSR/prate_lhtfl.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain_t(i,:),n2) 
           prec_rain_t(i,:)=prec_rain_t(i,:)/1000.  ! mm/s --> m/s

           vari='LHTFL_L1_Avg_1'
           file=trim(ForcingDataPath)//'CFSR/prate_lhtfl.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,evaporation_t(i,:),n2)
           evaporation_t(i,:)=evaporation_t(i,:)/2.5e6/1000. ! F_w [m/s]= Q_l / l_evap / rho_w

           ! net precipitation (p-e)
           prec_net_t(i,:)=prec_rain_t(i,:)-evaporation_t(i,:)
   
           if(update_forcing_flag(precip_ttp_ind)==1) then

              ! updating will take place in update_forcing in the first iteration

              prec_rain_t(2,:)=prec_rain_t(1,:)
	      cloud_t(2,:)=cloud_t(1,:)
	      evaporation_t(2,:)=evaporation_t(1,:)

              exit
           end if

           itime=itime+1
           if (itime>ndpyr) then	! to read the first record of the next year
              itime=1
              write(fileyear,'(i4)') yearnew+1
           endif

        end do
     end if

  else ! other precip data source

     write(*,*) 'The precipitation data source ', precip_data_source, ' is not supported.'
     write(*,*) 'Please check and update the code!'
     stop

  end if

  !==========================================================================
  ! radiation 

  if(rad_data_source=='CFSR') then
     write(*,*)'in CFSR radiation, rad_ttp_ind=',rad_ttp_ind

     if(rad_ttp_ind==1) then !6hourly

        fileyear=cyearnew
        itime=forcing_rec(rad_ttp_ind)
        
        do i=1,2

           ! short wave W/m2 --------------------------------------
           write(*,*) 'read_U_start10'

           vari='DSWRF_L1'
	   file=trim(ForcingDataPath)//'CFSR/dlwrf_dswrf.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave_t(i,:),n2)   

           ! long wave W/m2 ---------------------------------------

           vari='DLWRF_L1'
	   file=trim(ForcingDataPath)//'CFSR/dlwrf_dswrf.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave_t(i,:),n2)   

           if(update_forcing_flag(rad_ttp_ind)==1) then

              ! updating will take place in update_forcing in the first iteration

              shortwave_t(2,:)=shortwave_t(1,:)
              longwave_t(2,:)=longwave_t(1,:)

              exit
           end if

           itime=itime+1
           if (itime>4*ndpyr) then	! to read the first record of the next year
              itime=1
              write(fileyear,'(i4)') yearnew+1
           endif

        end do

     else !other temporal type
        write(*,*) 'The radiation-flux temporal type ', rad_ttp_ind, ' is not supported.'
        write(*,*) 'Please check and update the code!'
        stop
     end if

  else ! other rad data source

     write(*,*) 'The radiation data source ', rad_data_source, ' is not supported.'
     write(*,*) 'Please check and update the code!'
     stop

  end if


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
  ! cfsr, 6hourly
  use o_PARAM
  use o_MESH
  use o_array
  use i_array
  use i_dyn_parms
  use i_therm_parms
  use g_forcing_param
  use g_forcing_interp
  use g_forcing_arrays
  use g_forcing_index
  use g_parfe
  use g_clock
  use g_config
  implicit none

  integer		:: i
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy
  real              	:: t1, t2

  t1=MPI_Wtime()  

  ! first, read forcing data
  call read_new_atm_forcing

  ! second, do time interpolation 

 if(wind_data_source=='COSMO') then
   
 else  ! if COSMO

  ! wind, Tair, shum
  i_coef=interp_coef(wind_ttp_ind)
  do i=1,myDim_nod2d+eDim_nod2d                                        
     u_wind(i)=u_wind_t(1,i)+i_coef*(u_wind_t(2,i)-u_wind_t(1,i))
     v_wind(i)=v_wind_t(1,i)+i_coef*(v_wind_t(2,i)-v_wind_t(1,i))
     Tair(i)=Tair_t(1,i)+i_coef*(Tair_t(2,i)-Tair_t(1,i))
     shum(i)=shum_t(1,i)+i_coef*(shum_t(2,i)-shum_t(1,i))

  end do

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
        
        if(precip_data_source=='CFSR') then
           evaporation(i)=evaporation_t(1,i)+i_coef*(evaporation_t(2,i)-evaporation_t(1,i))
        end if
     end do
  end if
  
 endif  ! if COSMO
 
 
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

  ! fourth, compute wind stress
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
  ! read the second record of atmospheric forcing data
  ! cfsr, 6hourly
  use o_PARAM
  use o_MESH
  use o_array
  use i_therm_parms
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parfe
  use g_config
  implicit none
  !
  integer, parameter			:: nci=1152, ncj=576 !cfsr
  integer                   		:: itime,m, i, k, n2
  integer                               :: readtype
  character(180)             		:: file
  character(25)             		:: vari
  character(25)             		:: filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy

  n2=myDim_nod2D+eDim_nod2D                 

  !==========================================================================
  ! wind u and v, Tair, and shum
 if(wind_data_source=='CFSR') then
     if(update_forcing_flag(wind_ttp_ind)==1) then
        if(mype .eq. 0) write(*,*) 'read new wind record'
      
        !save the second record to the first record
        do i=1,myDim_nod2d+eDim_nod2d       
           u_wind_t(1,i)=u_wind_t(2,i)
           v_wind_t(1,i)=v_wind_t(2,i)
           Tair_t(1,i)=Tair_t(2,i)
           shum_t(1,i)=shum_t(2,i)
	   if(rad_ttp_ind==1) Pair_t(1,:)=Pair_t(2,:)
        end do

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind)
        itime=itime+1

        ! three temporal types (6 hourly, daily and monthly) are possible 

        if(wind_ttp_ind==1) then ! 6 hourly data

           if (itime>4*ndpyr) then	! to read the first record of the next year
              itime=1
              write(fileyear,'(i4)') yearnew+1
           endif

           ! 10-m wind m/s ----------------------------------------

           vari='U_GRD_L103'
           file=trim(ForcingDataPath)//'CFSR/wnd10m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='V_GRD_L103'
           file=trim(ForcingDataPath)//'CFSR/wnd10m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T382_wind(array_nc, array_nc2) !T62

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind_t(2,:),n2)   
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind_t(2,:),n2) 

           ! 2-m temperature --------------------------------------
           vari='TMP_L103'
           file=trim(ForcingDataPath)//'CFSR/tmp2mq2m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair_t(2,:),n2)   
           Tair_t(2,:)=Tair_t(2,:)-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------
           vari='SPF_H_L103'
           file=trim(ForcingDataPath)//'CFSR/tmp2mq2m.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum_t(2,:),n2)    
           
	   if(rad_ttp_ind==1) then
              vari='PRES_L1'
              file=trim(ForcingDataPath)//'CFSR/pressfc.gdas.'//fileyear//'.nc'
              call read_CFSR_NetCDF(file, vari, itime, array_nc)
              call upside_down(array_nc,nci,ncj)
              call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Pair_t(2,:),n2)
           end if
        end if
     end if

  endif


  !==========================================================================
  ! precipitation

  if(precip_data_source=='CFSR') then
     if(update_forcing_flag(precip_ttp_ind)==1) then
        if(mype .eq. 0) write(*,*) 'read new precip record'

        !save the second record to the first record
        do i=1,myDim_nod2d+eDim_nod2d              
           prec_rain_t(1,i)=prec_rain_t(2,i)
	   cloud_t(1,i)=cloud_t(2,i)
	   evaporation_t(1,i)=evaporation_t(2,i)
        end do

        fileyear=cyearnew
        itime=forcing_rec(precip_ttp_ind)
        itime=itime+1

        if(precip_ttp_ind==1) then ! 6 hourly data

           if (itime>4*ndpyr) then	! to read the first record of the next year
              itime=1
              write(fileyear,'(i4)') yearnew+1
           endif

           ! total precip mm/s ------------------------------------
           vari='PRATE_L1_Avg_1'
           file=trim(ForcingDataPath)//'CFSR/prate_lhtfl.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain_t(2,:),n2)  
           prec_rain_t(2,:)=prec_rain_t(2,:)/1000.  ! mm/s --> m/s
           
           vari='LHTFL_L1_Avg_1'
           file=trim(ForcingDataPath)//'CFSR/prate_lhtfl.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)	     
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,evaporation_t(2,:),n2)
           evaporation_t(2,:)=evaporation_t(2,:)/2.5e6/1000. ! F_w [m/s]= Q_l / l_evap / rho_w

        end if
     end if
  end if


  !==========================================================================
  ! radiation 

  if(rad_data_source=='CFSR') then
     if(update_forcing_flag(rad_ttp_ind)==1) then
        if(mype .eq. 0) write(*,*) 'read new rad record'

        !save the second record to the first record
        do i=1,myDim_nod2d+eDim_nod2d          
           shortwave_t(1,i)=shortwave_t(2,i)
           longwave_t(1,i)=longwave_t(2,i)
        end do

        fileyear=cyearnew
        itime=forcing_rec(rad_ttp_ind)
        itime=itime+1

         if(rad_ttp_ind==1) then ! 6 hourly data

           if (itime>4*ndpyr) then	! to read the first record of the next year
              itime=1
              write(fileyear,'(i4)') yearnew+1
           endif

           ! short wave W/m2 --------------------------------------
           vari='DSWRF_L1'
	   file=trim(ForcingDataPath)//'CFSR/dlwrf_dswrf.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave_t(2,:),n2) 

           ! long wave W/m2 ---------------------------------------

           vari='DLWRF_L1'
	   file=trim(ForcingDataPath)//'CFSR/dlwrf_dswrf.gdas.'//fileyear//'.nc'
           call read_CFSR_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave_t(2,:),n2)  

        end if

     end if
  end if


  !==========================================================================
  ! runoff not provided
  ! sss restoring not provided 

end subroutine read_new_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_T62_wind(xarray, yarray)
  ! rotate wind on T62 grid from geographical coord. to rotated coordinates.
  use o_param
  use g_rotate_grid
  use g_config
  implicit none

  integer, parameter 	:: ni=192, nj=94  ! NCEP and CORE are on the same grid.
  integer               :: i, j
  real(kind=8)      	:: cx(ni), cy(nj), xarray(ni,nj), yarray(ni,nj) 

  ! NCEP/CORE latitude
  cy=(/-88.542, -86.6531, -84.7532, -82.8508, -80.9473, -79.0435, &  
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
       82.8508, 84.7532, 86.6531, 88.542 /)*rad

  ! NCEP/CORE longitude
  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+1.875*rad
  enddo

  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do
end subroutine rotate_T62_wind
!
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_T382_wind(xarray, yarray)
  ! rotate wind on T382 grid from geographical coord. to rotated coordinates.
  use o_param
  use g_rotate_grid
  use g_config
  implicit none

  integer, parameter 	:: ni=1152, nj=576  ! NCEP and CORE are on the same grid.
  integer               :: i, j
  real(kind=8)      	:: cx(ni), cy(nj), xarray(ni,nj), yarray(ni,nj) 

  ! NCEP/CORE latitude
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
	89.7609978067649 /)*rad

  ! NCEP/CORE longitude
  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+0.3125*rad
  enddo

  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do
end subroutine rotate_T382_wind
!


!---------------------------------------------------------------------------------------------------
!
subroutine rotate_GME_wind(xarray, yarray)

  ! rotate wind on GME grid from geographical coord. to rotated coordinates.
  use o_param
  use g_parfe
  use g_rotate_grid
  use g_config
  implicit none

  integer, parameter    :: ni=720, nj=361 ! GME grid.
  integer               :: i, j
  real(kind=8)          :: cx(ni), cy(nj), xarray(ni,nj), yarray(ni,nj)

  ! GME latitude
  cy(1)=-90.0*rad
  do i=2,nj
     cy(i)=cy(i-1)+0.5*rad
  enddo

  ! GME longitude
  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+0.5*rad
  enddo

  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do

end subroutine rotate_GME_wind
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_COSMO_wind(xarray, yarray)

  ! rotate wind on COSMO interpolated grid from geographical coord. to rotated coordinates.
  use o_param
  use g_parfe
  use g_rotate_grid
  use g_config
  implicit none

  integer, parameter    :: nic=600, njc=180 ! COSMO interpolated grid.
  integer               :: i, j
  real(kind=8)          :: cxc(nic), cyc(njc), xarray(nic,njc), yarray(nic,njc)

  ! COSMO latitude
  cyc(1)=-80.0*rad
  do i=2,njc
     cyc(i)=cyc(i-1)+0.125*rad
  enddo

  ! COSMO longitude
  cxc(1)=285.0*rad
  do i=2,nic
     cxc(i)=cxc(i-1)+0.125*rad
  enddo
  
  ! CHECK
  if (abs(cyc(njc)+57.625*rad) .gt. 1e-5 .or. abs(cxc(nic)-359.875*rad) .gt. 1e-5) then
    write(*,*) 'x: ',cxc(1)/rad,' to ',cxc(nic)/rad
    write(*,*) 'y: ',cyc(1)/rad,' to ',cyc(njc)/rad
    write(*,*) 'ERROR: COSMO lon and/or lat calculated wrongly in rotate_COSMO_wind!'
  endif
  
  !rotate wind
  !cxc cyc are in radian
  do i=1,nic
     do j=1,njc
        call vector_g2r(xarray(i,j), yarray(i,j), cxc(i), cyc(j), 1)
     end do
  end do

end subroutine rotate_COSMO_wind
!

!=========================================================================
subroutine shortwave_radiation(daynew,nhi,nhf,i,n2)
!
!calculates the incoming shortwave radiation
!
! input data:
! ndoyr               day of year  (1..365 or 366)
! nhi                 initial hour (0...23)
! nhf                 final hour   (0...23)
!                     (last hour to be completely considered)
!
! For instance, for a 6h-timestep:     nhi    nhf
!                                       0      5
!                                      6     11
!                                      12     17
!                                      18     23
!
!            for a daily timestep:      0     23
!
!
! derived from subroutine swr of BRIOS sea ice model by
! Ralph Timmermann, 23.9.2004
!=========================================================================
USE o_MESH
USE o_PARAM
USE g_forcing_arrays
USE g_forcing_interp
USE g_config

IMPLICIT NONE

INTEGER               ::      ndoyrd,daynew, nhi, nhf
INTEGER               ::      it, iday, itd, i, n2
real*8, dimension(0:23) ::      cosha
real*8, dimension(n2)   ::      sinlat, coslat, cosz
real*8, dimension(n2)   ::    cf, fsh
real*8                  ::      dec, sindec, cosdec

        ! write (*,*) 'shortwave radiation',daynew,':',nhi,'-',nhf
! ----------------------------------------------------------------------
! Cosine of hour angle for every hour
! [Parkinson and Washington, 1979, eq. (3)]
! ----------------------------------------------------------------------
do it=0,23
 cosha(it)=cos(2*pi*(12-it)/24.)                       ! correct for GM
enddo

! ----------------------------------------------------------------------
! Declination of the sun for every day
! -----------------------------------------------------------------------
ndoyrd=min(daynew,365)    !modified SEB
dec=rad*23.44*cos((174.-(ndoyrd-1))*2*pi/365.) ! 365 days/year
sindec=sin(dec)
cosdec=cos(dec)

!-----------------------------------------------------------------------
! Sine and cosine of geographical latitude
!-----------------------------------------------------------------------
sinlat=sin(ymod)
coslat=cos(ymod)

!-----------------------------------------------------------------------
!       fsh             shortwave radiation
!       cf              cloud factor
!-----------------------------------------------------------------------
fsh  = 0.
!e_vapor = shum_t*Pair_t/0.622*1.e-3    !rt ausgelagert
!tdew = 1.e20    ! dummy
!ecmwf e_vapor(i,j)=611.e-5*exp(19.*tdew(i,j)/(tdew(i,j)+250.))
cf(:)   = 1.-0.6*cloud_t(i,:)**3


!-----------------------------------------------------------------------
!LOOP for every hour
!-----------------------------------------------------------------------
do it=nhi,nhf
 if (it.lt.0) then
  itd=it+24
 else
  itd=it
 endif
!-----------------------------------------------------------------------
! Cosine of zenith distance
! Parkinson and Washington, 1979, eq. (2)
!-----------------------------------------------------------------------
 cosz=sinlat*sindec+coslat*cosdec*cosha(itd)
!-----------------------------------------------------------------------
! At night, when no sun is present and cosz < 0.,
! the incoming solar radiation is zero (and not negative).
!-----------------------------------------------------------------------
 cosz=max(cosz,0.)
!-----------------------------------------------------------------------
! Add up incoming shortwave radiation for each and every hour
! Parkinson and Washington, 1979, eq. (1)
!
! Solar constant = 1353 [W/m**2]
!         56.375 = 1353 / 24 for contributions from 24 hours
!         1353/float(nhf-nhi+1) for contributions from nhf-nhi hours
!-----------------------------------------------------------------------
 fsh(:)=fsh(:)+1353./float(nhf-nhi+1)*cosz*cosz*cf(:) &
              /((cosz+2.7)*e_vapor_t(i,:)+1.085*cosz+0.1)
!-----------------------------------------------------------------------
enddo

shortwave_t(i,:)=fsh

end subroutine shortwave_radiation
!=========================================================================

