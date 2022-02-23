program main
  !=============================================================================!
  !
  !                 Finite Element Sea-ice Ocean Model
  !
  !=============================================================================!
  !                      The main driving routine
  !=============================================================================!         

  use g_config
  use o_param
  use o_array          
  use o_solver
  use o_mixing_kpp_mod
  use g_parfe
  use g_clock
  use g_forcing_index
  use g_forcing_param
  use g_forcing_arrays
  use g_diag
  use g_forcing_interp
#ifdef use_ice
  use i_array
#endif
  implicit none
  integer yearfromicemod


  ! MPI initialization
  call par_init                 ! initializes MPI
  if(mype==0) write(*, *) 'Running on ', npes, ' PEs'


  if(mype==0) write(*,*) '*************************************************************'
  ! read namelist, initialize clock, prepare basic configuration etc.

  call setup_model              ! setup basic config, do it before clock_init
  call clock_init               ! read the clock file
  call get_run_steps            ! compute total steps to run 
  call config_remind_warning_info

!!RT RG45901: check that clock entry is consistent with yearfromicemod.dat
!  open(99,file=trim(meshpath)//'yearfromicemod.dat',status='old')
!  read(99,*) yearfromicemod
!  close(99)
!  if (yearnew.ne.yearfromicemod) then
!   write(*,*)'ERROR: yearnew ne yearfromicemod',yearnew,yearfromicemod
!   stop
!  endif
!RT-



  if(mype==0) write(*,*) '*************************************************************'
  ! mesh and communication buffers

  call ocean_mesh_setup         ! setup the 2D/3D meshes
  call set_par_support
  call mesh_cluster_setup       ! build cluster area and volume and save to file

  if(mype==0) write(*,*) '*************************************************************'
  ! ocean: matrices, arrays, initialization, buffer zone, tide etc.

  call ocean_matrices_setup     ! Builds matrices and call partitioning
  call ocean_array_setup        ! allocate ocean arrays 
  call ocean_init               ! initialize the oce or read restart files

  if(use_ref_density) then
     call compute_ref_density   ! Fills in ocean reference density 
  endif

#ifdef use_opbnd_restoring
  call init_restoring_vel
#endif

  if(buffer_zone) then
     call init_restoring_bufferzone_isomip
  end if

#ifdef use_opbnd_tide
  call init_tidal_opbnd         ! initialize tidal ocean open boundary
#endif

#ifdef use_ice
  if(mype==0) write(*,*) '*************************************************************'
  ! ice: matrices, arrays, initialization

  call ice_matrices_setup       ! Build ice matrices
  call ice_array_setup          ! allocate ice arrays, setup ice adv matrix
  call ice_init                 ! initialize the ice or read restart files
#endif


  if(mype==0) write(*,*) '*************************************************************'
  ! forcing: arrays, initialization, interpolation preparation  

#if defined(use_ice)
  call forcing_array_setup
  !call init_forcing_interp !calculates the forcing interpolation weights
  !call init_atm_forcing         ! initialize forcing fields
#else
#ifndef toy_ocean
  call forcing_array_setup_OnlyOcean
  call init_forcing_interp 
  call init_atm_forcing_OnlyOcean 
#endif
#endif 

  if(use_landice_water) call landice_water_init


  if(mype==0) write(*,*) '*************************************************************'
  ! init mean arrays, adjust clock and create output files if required  

#if defined(allow_calcmeans) || defined(allow_diag)
  call init_meanarrays          ! allocate arrays for mean fields
#endif
  call clock_newyear  		! check if it is a new year
  call init_output              ! create new output files


#ifdef use_fullfreesurf
  if(mype==0) write(*,*) '*************************************************************'
  ! updating mesh and matrices in case of full free surface setup 

  if(any(abs(ssh)>small)) then
     call update_mesh
     call update_matrices
     call update_mesh
  endif
#endif

  ! set some flags for solvers
  iter_first=.true.             ! iter_first & iteruv_first should be 'true' at start
  iteruv_first=.true.


  if(mype==0) then
     write(*,*) '*************************************************************'
     write(*,*) 'iteration starts ...'
     write(*,*) '*************************************************************'
     write(*,*)
  end if

!RT:
!  if (daynew.eq.1) then
!   write(*,*)'new year, write initial data into file'
!   call output(1)
!  endif
!RT-

  ! preparation done
  !----------------------------------------------------------------------------
  ! start iteration

  do istep=1, nsteps   
     call clock
     call init_output

     call forcing_index 

#ifdef use_ice    
     call ocean2ice
#ifndef isomip
     call update_atm_forcing
     call ice_step
     call ice2ocean
#endif
     if(use_landice_water) call add_landice_water
#else
#ifndef toy_ocean
     call update_atm_forcing_OnlyOcean
#endif
#endif
#ifdef use_cavity
  call cavity_momentum_fluxes
  call cavity_heat_water_fluxes_3eq
#endif

#ifdef use_fullfreesurf 
     if(balance_salt_water) call check_imb_freshwater
#endif

#ifdef use_sw_pene
     call cal_shortwave_rad
#endif

     call ocean_step

#if defined(allow_calcmeans) || defined(allow_diag)
     call add2meanarrays      
#endif

     ! save (NetCDF)
     call output(0) 

!!$     ! save (ascii, for an easy debugging during new setup test phase)
!!$     if(mod(istep,step_per_day*5)==0) then
!!$        call oce_out     
!!$#ifdef use_ice
!!$        call ice_out
!!$#endif
!!$     end if

     if(mod(istep,logfile_outfreq)==0) then	
        !log file output (for debugging)
!        write(*,*) 'uf max', maxval(abs(uf)), 's min', minval(tracer(:,2)), istep

        if (mype==0) then
	 !  write(*,*) 'uf max', maxval(abs(uf)), istep
           write(*,*) 'completed step', istep, '  day', daynew, '  year', yearnew
  	   write(*,*)
        endif
     end if

     if(check_run_state) call check_blowup	! check if the program blows up

  end do

  ! iteration done
  !--------------------------------------------------------------------
  ! some finishing-up routines follow

  if(mix_scheme=='MY2p5') call save_MY_vara  ! save MY2.5 variables for next restart

  ! save (ascii, for an easy debugging during new setup test phase)
  !call oce_out
#ifdef use_ice
  !call ice_outfesom_main.F90
#endif

  call clock_finish		! save clock

  if(mype==0) then
     open(unit=50, file='goodfile')
     write(50,*)'go on'
     close(50)
  end if

 
  if(mype==0) then 
    if (yearold.eq.yearnew) then 
     open(unit=50, file='goonfile')
     write(50,*)'go on'
     close(50)
    endif
  end if


  call par_ex 			! finalizes MPI
  write(*,*) 'Experiment '//runid//' successfully completed'

end program main


