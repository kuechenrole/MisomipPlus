! initialize mean and diagnose arrays, compute and clean these variables.
! Coded by Ralph Timmermann
! Modified by Qiang Wang to include more diagnose variables
! Reviewed by ??
!---------------------------------------------------------------

subroutine init_meanarrays
  ! allocates and initializes fields used for mean arrays
  use o_mesh
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe

  implicit none
  integer       :: n2, n3

  n2=myDim_nod2D+eDim_nod3D        
  n3=myDim_nod3D+eDim_nod3D       

  ! for mean fields

#ifdef allow_calcmeans

  ! ocean part                    
  allocate(sshmean(n2))
#ifndef use_non_hydrostatic
  allocate(ufmean(2*n3)) 
  allocate(wrhsmean(n3))
  
#else
  allocate(ufmean(3*n3)) 
#endif
  allocate(tracermean(n3,num_tracer))

  ! ice part
#ifdef use_ice
  allocate(m_ice_mean(n2), a_ice_mean(n2), m_snow_mean(n2))
  allocate(u_ice_mean(n2), v_ice_mean(n2))
#endif

#endif

  ! for diagnostics

#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(diag_oce_KE) then
        allocate(uumean(n3), vvmean(n3))
     end if
     if(diag_oce_energy_conv) then
        allocate(rhomean(n3),urhomean(n3))
        allocate(vrhomean(n3),uvmean(n3))    
     end if
     if(diag_oce_transp) then
        allocate(uTFmean(n3),vTFmean(n3))
        allocate(uSFmean(n3),vSFmean(n3))
     end if
     if(Redi_GM .and. diag_oce_GM_vel) then
        allocate(sgs_u(myDim_elem3d), sgs_v(myDim_elem3d))
     end if
     if(diag_oce_SGS_transp) then
        allocate(sgs_ut(myDim_elem3d), sgs_vt(myDim_elem3d))
        allocate(sgs_us(myDim_elem3d), sgs_vs(myDim_elem3d))
     end if
     if(diag_oce_mix_layer) then
        allocate(mixlay_dep_mean(n2))
     end if
  end if

#ifdef use_ice
  ! ice
  if(diag_ice) then
     allocate(thdgr_mean(n2), thdgrsn_mean(n2))
     allocate(uhice_mean(n2), vhice_mean(n2))
     allocate(uhsnow_mean(n2), vhsnow_mean(n2))
     allocate(flice_mean(n2))
  end if

  ! forcing
  if(diag_forcing) then
     allocate(tair_mean(n2), tdew_mean(n2), shum_mean(n2))
     allocate(uwind_mean(n2), vwind_mean(n2))
     allocate(rain_mean(n2), snow_mean(n2))
     allocate(runoff_mean(n2), evap_mean(n2))
     allocate(lwrd_mean(n2), swrd_mean(n2))
     allocate(qnet_mean(n2), wnet_mean(n2))
     allocate(olat_mean(n2), osen_mean(n2))
     allocate(olwout_mean(n2))
     allocate(virtual_salt_mean(n2), relax_salt_mean(n2))
     allocate(stress_x_mean(n2), stress_y_mean(n2))
     allocate(Tsurfmean(n2), Ssurfmean(n2))
  end if
#endif

#endif

  call clean_meanarrays

  if(mype==0) write(*,*) 'Mean arrays have been set up'

  return
end subroutine init_meanarrays
!=============================================================================!


!=============================================================================!
subroutine add2meanarrays
  ! adds values to the mean-arrays
  use o_mesh
  use o_param
  use o_array
  use i_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_forcing_arrays
  use g_parfe
  implicit none
  !
  integer       :: m, row, row2, row3, j

  meancounter=meancounter+1

#ifdef allow_calcmeans
#ifndef use_non_hydrostatic
  call vvel_nodes 
#endif
#endif

  do row=1,myDim_nod3d                      

     row2=row+myDim_nod3d+eDim_nod3D       

#ifdef allow_calcmeans

     ! ocean
     ufmean(row)        = ufmean(row) + uf(row)
     ufmean(row2)       = ufmean(row2) + uf(row2)
#ifdef use_non_hydrostatic
     row3=row2+myDim_nod3D+eDim_nod3D       
     ufmean(row3)       = ufmean(row3) + uf(row3)
#else
     wrhsmean(row)      = wrhsmean(row) + wrhs(row)
#endif
     do j=1,num_tracer
        tracermean(row,j)  = tracermean(row,j) + tracer(row,j)
     end do
#endif

#ifdef allow_diag
     
     if(diag_oce) then
        if(diag_oce_KE) then
           uumean(row)   = uumean(row)   + uf(row)*uf(row)
           vvmean(row)   = vvmean(row)   + uf(row2)*uf(row2)
        end if

        if(diag_oce_energy_conv) then
           rhomean(row)  = rhomean(row)  + density_insitu(row)
           urhomean(row) = urhomean(row) + uf(row)*density_insitu(row)
           vrhomean(row) = vrhomean(row) + uf(row2)*density_insitu(row)   
           uvmean(row)   = uvmean(row)   + uf(row)*uf(row2)  
        end if
        
        if(diag_oce_transp) then
           uTFmean(row)  = uTFmean(row) + uf(row)*tracer(row,1)
           vTFmean(row)  = vTFmean(row) + uf(row2)*tracer(row,1)
           uSFmean(row)  = uSFmean(row) + uf(row)*tracer(row,2)
           vSFmean(row)  = vSFmean(row) + uf(row2)*tracer(row,2)
        endif
     end if
#endif
  end do

  !-----------------------------------------------------------

  do row=1,myDim_nod2d  

     ! ocean
#ifdef allow_calcmeans 
     sshmean(row)       = sshmean(row) + ssh(row)
#endif

     ! ice
#ifdef use_ice
#ifdef allow_calcmeans
     a_ice_mean(row)    = a_ice_mean(row) + a_ice(row)
     m_ice_mean(row)    = m_ice_mean(row) + m_ice(row)     
     m_snow_mean(row)   = m_snow_mean(row) + m_snow(row)
     u_ice_mean(row)    = u_ice_mean(row)  + u_ice(row)
     v_ice_mean(row)    = v_ice_mean(row)  + v_ice(row)
#endif
#ifdef allow_diag
     if(diag_ice) then
        thdgr_mean(row)    = thdgr_mean(row)  + thdgr(row)
        thdgrsn_mean(row)  = thdgrsn_mean(row) + thdgrsn(row)
        uhice_mean(row)    = uhice_mean(row) + u_ice(row)*m_ice(row)
        vhice_mean(row)    = vhice_mean(row) + v_ice(row)*m_ice(row)
        uhsnow_mean(row)   = uhsnow_mean(row) + u_ice(row)*m_snow(row)
        vhsnow_mean(row)   = vhsnow_mean(row) + v_ice(row)*m_snow(row)
        flice_mean(row)    = flice_mean(row) + flice(row)
     endif
#endif
#endif

     ! forcing
#ifdef allow_diag
#ifdef use_ice
     if(diag_forcing) then
        tair_mean(row)     = tair_mean(row) + Tair(row)
        tdew_mean(row)     = tdew_mean(row) + Tdew(row)
        shum_mean(row)     = shum_mean(row) + shum(row)
        uwind_mean(row)    = uwind_mean(row) + u_wind(row)
        vwind_mean(row)    = vwind_mean(row) + v_wind(row)   
        rain_mean(row)     = rain_mean(row) + prec_rain(row)
        snow_mean(row)     = snow_mean(row) + prec_snow(row)
        runoff_mean(row)   = runoff_mean(row) + runoff(row)
        evap_mean(row)     = evap_mean(row) + evaporation(row)
        lwrd_mean(row)     = lwrd_mean(row) + longwave(row)
        swrd_mean(row)     = swrd_mean(row) + shortwave(row)
        qnet_mean(row)     = qnet_mean(row) + net_heat_flux(row)
        wnet_mean(row)     = wnet_mean(row) + fresh_wa_flux(row)
        olat_mean(row)     = olat_mean(row) + olat_heat(row)
        osen_mean(row)     = osen_mean(row) + osen_heat(row)
        olwout_mean(row)   = olwout_mean(row) + olwout(row)
        virtual_salt_mean(row) = virtual_salt_mean(row) + virtual_salt(row)
        relax_salt_mean(row)   = relax_salt_mean(row) + relax_salt(row)
        stress_x_mean(row) = stress_x_mean(row) + stress_x(row)
        stress_y_mean(row) = stress_y_mean(row) + stress_y(row)
        Tsurfmean(row)     = Tsurfmean(row) + Tsurf(row)
        Ssurfmean(row)     = Ssurfmean(row) + Ssurf(row)
     endif
#endif
#endif

  enddo

  ! other diagnostics

  ! mixed layer thickness
#ifdef allow_diag
  if(diag_oce .and. diag_oce_mix_layer) then
     call compute_mixlay
  end if
#endif

  ! diagnose for SGS parameterization is done in the ts_rhs routine

  return
end subroutine add2meanarrays
!=============================================================================!


!=============================================================================!
subroutine compute_means
  !computes the mean values for output
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  implicit none

  real(kind=8)   :: cnt

  cnt=float(max(meancounter,1))

#ifdef allow_calcmeans

  ! ocean
  sshmean               = sshmean            /cnt
  ufmean                = ufmean             /cnt
#ifndef use_non_hydrostatic
  wrhsmean              = wrhsmean           /cnt
#endif
  tracermean            = tracermean         /cnt

  ! ice
#ifdef use_ice
  a_ice_mean            = a_ice_mean         /cnt
  m_ice_mean            = m_ice_mean         /cnt
  m_snow_mean           = m_snow_mean        /cnt
  u_ice_mean            = u_ice_mean         /cnt
  v_ice_mean            = v_ice_mean         /cnt
#endif

#endif

  ! diagnose
#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(diag_oce_KE) then
        uumean                = uumean             /cnt
        vvmean                = vvmean             /cnt
     end if
     if(diag_oce_energy_conv) then
        rhomean               = rhomean            /cnt
        urhomean              = urhomean           /cnt
        vrhomean              = vrhomean           /cnt
        uvmean                = uvmean             /cnt
     end if
     if(diag_oce_transp) then
        uTFmean               = uTFmean            /cnt
        vTFmean               = vTFmean            /cnt
        uSFmean               = uSFmean            /cnt
        vSFmean               = vSFmean            /cnt
     endif
     if(Redi_GM .and. diag_oce_GM_vel) then
        sgs_u                 = sgs_u              /cnt
        sgs_v                 = sgs_v              /cnt
     endif
     if(diag_oce_SGS_transp) then
        sgs_ut                = sgs_ut             /cnt
        sgs_vt                = sgs_vt             /cnt
        sgs_us                = sgs_us             /cnt
        sgs_vs                = sgs_vs             /cnt
     endif
     if(diag_oce_mix_layer) then
        mixlay_dep_mean       = mixlay_dep_mean    /cnt
     endif    
  endif

  ! ice
#ifdef use_ice
  if(diag_ice) then
     thdgr_mean            = thdgr_mean         /cnt
     thdgrsn_mean          = thdgrsn_mean       /cnt
     uhice_mean            = uhice_mean         /cnt
     vhice_mean            = vhice_mean         /cnt
     uhsnow_mean           = uhsnow_mean        /cnt
     vhsnow_mean           = vhsnow_mean        /cnt
     flice_mean            = flice_mean         /cnt
  endif
#endif

  ! forcing
#ifdef use_ice
  if(diag_forcing) then
     tair_mean             = tair_mean          /cnt
     tdew_mean             = tdew_mean          /cnt
     shum_mean             = shum_mean          /cnt
     uwind_mean            = uwind_mean         /cnt
     vwind_mean            = vwind_mean         /cnt
     rain_mean             = rain_mean          /cnt
     snow_mean             = snow_mean          /cnt
     runoff_mean           = runoff_mean        /cnt
     evap_mean             = evap_mean          /cnt
     lwrd_mean             = lwrd_mean          /cnt
     swrd_mean             = swrd_mean          /cnt
     qnet_mean             = qnet_mean          /cnt
     wnet_mean             = wnet_mean          /cnt
     olat_mean             = olat_mean          /cnt
     osen_mean             = osen_mean          /cnt
     olwout_mean           = olwout_mean        /cnt
     virtual_salt_mean     = virtual_salt_mean  /cnt
     relax_salt_mean       = relax_salt_mean    /cnt
     stress_x_mean         = stress_x_mean      /cnt
     stress_y_mean         = stress_y_mean      /cnt
     Tsurfmean             = Tsurfmean          /cnt
     Ssurfmean             = Ssurfmean          /cnt    
  endif
#endif

#endif

end subroutine compute_means
!=============================================================================!


!=============================================================================!
subroutine clean_meanarrays
  ! puts zeros into the mean-arrays
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  implicit none

  meancounter=0

#ifdef allow_calcmeans

  ! coean
  sshmean=0.
  ufmean=0.
#ifndef use_non_hydrostatic
  wrhsmean=0.
#endif
  tracermean=0.

  ! ice
#ifdef use_ice
  a_ice_mean=0.
  m_ice_mean=0.
  m_snow_mean=0.
  u_ice_mean=0.
  v_ice_mean=0.
  thdgr_mean=0.
  thdgrsn_mean=0.
  uhice_mean=0.
  vhice_mean=0.
  uhsnow_mean=0.
  vhsnow_mean=0.
  flice_mean=0.
#endif

#endif

  ! diagnose
#ifdef allow_diag

  if(diag_oce) then
     ! ocean
     if(diag_oce_KE) then
        uumean=0.
        vvmean=0.
     end if
     if(diag_oce_energy_conv) then
        rhomean=0.
        urhomean=0.
        vrhomean=0.
        uvmean=0.
     end if
     if(diag_oce_transp) then
        uTFmean=0.
        vTFmean=0.
        uSFmean=0.
        vSFmean=0.
     endif
     if(Redi_GM .and. diag_oce_GM_vel) then
        sgs_u=0.
        sgs_v=0.
     endif
     if(diag_oce_SGS_transp) then
        sgs_ut=0.
        sgs_vt=0.
        sgs_us=0.
        sgs_vs=0.
     endif
     if(diag_oce_mix_layer) then
        mixlay_dep_mean=0.
     endif
  endif

  ! ice
#ifdef use_ice
  if(diag_ice) then
     thdgr_mean=0.
     thdgrsn_mean=0.
     uhice_mean=0.
     vhice_mean=0.
     uhsnow_mean=0.
     vhsnow_mean=0.
     flice_mean=0.
  endif
#endif

  ! forcing
#ifdef use_ice
  if(diag_forcing) then
     tair_mean=0.
     tdew_mean=0.
     shum_mean=0.
     uwind_mean=0.
     vwind_mean=0.
     rain_mean=0.
     snow_mean=0.
     runoff_mean=0.
     evap_mean=0.
     lwrd_mean=0.
     swrd_mean=0.
     qnet_mean=0.
     wnet_mean=0.
     olat_mean=0.
     osen_mean=0.
     olwout_mean=0.
     virtual_salt_mean=0.
     relax_salt_mean=0.
     stress_x_mean=0.
     stress_y_mean=0.
     Tsurfmean=0.
     Ssurfmean=0.
  endif
#endif

#endif

  return
end subroutine clean_meanarrays
!=============================================================================!
