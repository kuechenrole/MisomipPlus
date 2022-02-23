!=============================================================================
!  Performs a time step of the ocean model
!=============================================================================

subroutine ocean_step
  ! driving routine for ocean stepping
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !----------------------------------------------------------

  use o_param
  use o_array
  use o_mixing_kpp_mod
  use o_mixing_pp_mod
  use o_mixing_my2p5_mod
  use o_mixing_tidal_mod
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use o_mesh
  use o_elements
  use o_solver
  use g_config
  use g_parfe

  implicit none
  integer      :: i, m, row, row2, row3, n3
  real(kind=8) :: t0,t1,t2,t3,t4,t5,t6,t7, t8,t9,t10, t11
  real(kind=8) :: t111, t112, t113, t114, t115, t116, t117, t118, t119, t120

  t0=MPI_Wtime() 

  n3=ToDim_nod3d

  do row=1,n3                           
     row2=row+n3                   
     uf0(row)=uf(row)                 ! uf0 & uf: u^n
     uf0(row2)=uf(row2) 
#ifdef use_non_hydrostatic
     row3=row2+n3                     
     uf0(row3)=uf(row3)
#endif
  end do

  !-----------------------------------------------------------

  call compute_density
  call compute_bfsq
  call compute_dbsfc

  if(grid_type/=2) then
     call compute_pressure            ! compute hpressure
  end if
  if(grid_type/=1) then
     call compute_pressure_force      ! compute pressure gradient
  end if
  t2=MPI_Wtime()    

  !-----------------------------------------------------------

  if(trim(mix_scheme)=='KPP') then
     call oce_mixing_kpp(Av, Kv)
     call convect_adjust
  elseif(trim(mix_scheme)=='PP') then
     call oce_mixing_pp
  elseif(trim(mix_scheme)=='MY2p5') then
     call oce_mixing_MY2p5
  elseif(trim(mix_scheme)=='no') then
     call convect_adjust
  end if
  
  if(tidal_mixing) call oce_mixing_tidal(Av, Kv)
  t3=MPI_Wtime()

  !------------------------------------------------------------

  call velocity_rhs                 ! rhs for u*
  if(use_vertvisc_impl) call uv_sfc_bott_bc
  t4=MPI_Wtime()    

  if(lump_uv_matrix) then           ! solver for du*
     call uv_solve
  else
     call solve(solve_u)   
     iteruv_first=.false.
     call solve(solve_v)
     call com_3d(duf(1:n3))              
     call com_3d(duf(1+n3:2*n3))           
  endif

#ifdef use_non_hydrostatic 
  call solve(solve_w)
  call com_3d(duf(1+2*n3:3*n3))          
#endif 
  do row=1,n3                            
     row2=row+n3                   
     uf(row)=uf(row)+duf(row)      ! uf: u*
     uf(row2)=uf(row2)+duf(row2)
#ifdef use_non_hydrostatic
     row3=row2+n3                         
     uf(row3)=uf(row3)+duf(row3)
#endif
  end do

  if(use_vertvisc_impl) then      ! apply implicit vertical viscosity
     call impl_vertvisc
  end if
  t5=MPI_Wtime() 

  !--------------------------------------------------------------

  ssh0=ssh                        ! ssh & ssh0: ssh^n  

#ifdef use_opbnd_tide
  call update_tidal_opbnd         ! update tidal open boundary
#endif

  call compute_ssh_rhs            ! ssh rhs

  call solve(solve_ssh)           ! solve dssh
  ssh=ssh0+dssh                   ! ssh: ssh^n+1
  call com_2D(ssh)                       

#ifdef use_non_hydrostatic
  nhp0=nhp
  call compute_nhp_rhs            ! nhp rhs
  call solve(solve_nhp)           ! solve nhp
  call com_3D(nhp)                       
#endif   
  t6=MPI_Wtime()      

  do row=1,n3                              
     row2=row+n3                        
     uf0(row)=uf(row)             ! uf0: u*     
     uf0(row2)=uf(row2) 
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uf0(row3)=uf(row3)
#endif	
  end do

  !---------------------------------------------------------------

  call velocity_rhs_update        ! Update rhs: contribution from ssh/nhp
  t7=MPI_Wtime()

  if(lump_uv_matrix) then         ! solve for full du 
     call uv_solve
  else
     call solve(solve_u)                          
     call solve(solve_v)
     call com_3d(duf(1:n3))             
     call com_3d(duf(1+n3:2*n3))        
  endif
#ifdef use_non_hydrostatic
  call solve(solve_w)
  call com_3d(duf(1+2*n3:3*n3))         
#endif
  do row=1,n3                            
     row2=row+n3                         
     uf(row)=uf(row)+duf(row)     ! uf: u^n+1  
     uf(row2)=uf(row2)+duf(row2)
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uf(row3)=uf(row3)+duf(row3)
#endif	
  end do
  t8=MPI_Wtime()

  !---------------------------------------------------------

#ifndef use_non_hydrostatic
  call compute_vvel_rhs           ! vertical rhs
  call solve_wpot                 ! solve for w potential 
#endif
  t9=MPI_Wtime()  

  !----------------------------------------------------------

#ifdef use_fullfreesurf
  call update_mesh
  call update_matrices
#endif 
  t10=MPI_Wtime() 

  !-----------------------------------------------------------

  if(Redi_GM) call compute_neutral_slope   ! calc. neutral slope

  !-----------------------------------------------------------

  if(brine_rejection_param) call cal_brine_rejection

#ifdef use_tracer_gls
  call tsstiff_fill_gls          ! tracer matrix/rhs
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  do i=1,num_tracer
     call solve(solve_tra+i-1)   ! solve for tracer
     call com_3D(dtracer(:,i))          
  end do
  do i=1,num_tracer                     
     tracer(:,i)=tracer(:,i)+dtracer(:,i)
  end do

#else
 
!if use_tracer_gls is undefined
#ifdef use_tracer_fct
  !if (mype.eq.0) write(*,*) 'before call tracer_rhs_tg'
  t111=MPI_Wtime()
  call tracer_rhs_tg
  t112=MPI_Wtime()
  !if (mype.eq.0) write(*,*) 'before call ts_sfc_bc'
  call ts_sfc_bc
  t113=MPI_Wtime()
  !if (mype.eq.0) write(*,*) 'before call ptr_sfc_bc'
  if(use_passive_tracer) call ptr_sfc_bc
  t114=MPI_Wtime()
  if(use_age_tracer) call age_tracer_tendency
  t115=MPI_Wtime()
  t11=MPI_Wtime()
  !if (mype.eq.0) write(*,*) 'before call fct_tracer_solve'
  call fct_tracer_solve 
  t116=MPI_Wtime()
#else
  call tracer_rhs_tg
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  if(lump_ts_matrix) then
     call tracer_solve
  else
     do i=1,num_tracer
        call solve(solve_tra+i-1)
	call com_3D(dtracer(:,i))      
     end do
  endif
  do i=1,num_tracer
     do row=1,n3                         
        tracer(row,i)=tracer(row,i)+dtracer(row,i)
     end do
  end do
#endif

  t117=MPI_Wtime()
  if(brine_rejection_param) call apply_brine_rejection
  t118=MPI_Wtime()
  if(use_vertdiff_impl) then
     call impl_vertdiff        ! apply implicit vertical diff.
  end if
  t119=MPI_Wtime()
#endif

  if(use_passive_tracer) call ptr_cutoff_restore
  t120=MPI_Wtime()
  if(use_age_tracer) call age_tracer_cutoff_restore

  !--------------------------------------------------------------

  t1=MPI_Wtime()
  iter_first = .false.

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then      
     write(*,*)
     write(*,*) 'ocean  took   ', t1-t0
     write(*,*) 'Dens/pressure  took', t2-t0
     write(*,*) 'mixing scheme  took', t3-t2
     write(*,*) 'v rhs          took', t4-t3
     write(*,*) 'solve v_star   took', t5-t4
     write(*,*) 'pressure_solve took', t6-t5
     write(*,*) 'rhs_update     took', t7-t6
     write(*,*) 'solve full v   took', t8-t7
#ifndef use_non_hydrostatic
     write(*,*) 'vert_vel_solve took', t9-t8
#endif
#ifdef use_fullfreesurf
     write(*,*) 'update mesh    took', t10-t9
#endif
     write(*,*) 'tra_assemble   took', t11-t10
     write(*,*) 'tra_solve etc. took', t1-t11

     write(*,*) 'tracer_rhs_tg  took', t112-t111
     write(*,*) 'ts_sfc_bc      took', t113-t112
     write(*,*) 'ptr_sfc_bc     took', t114-t113
     write(*,*) 'age_tracer_tendency took', t115-t114
     write(*,*) 'fct_tracer_solve    took', t116-t115
     write(*,*) 'apply_brine_rejection took', t118-t118
     write(*,*) 'impl_vertdiff  took', t119-t118
     write(*,*) 'ptr_cutoff_restore  took', t120-t119
  endif

end subroutine ocean_step
!=================================================================
