! Sea-ice thermodynamics routines
!
! Coded by N. Yakovlev and S. Danilov.
! Adjusted for upgraded model physics (NCEP forcing data; parameterization
! of ice-ocean heat flux considering friction velocity, etc) 
! by Ralph Timmermann.
! Adjusted for general forcing data and NlFs option, cleaned up, bug fixing,
! by Qiang Wang, 13.01.2009
! Upgraded for fail-save operation in a world with more than just two forcing datasets
! by Ralph Timmermann, 08.06.2012.
!----------------------------------------------------------------------------

subroutine thermodynamics
  !
  ! For every surface node, this subroutine extracts the information
  ! needed for computation of thermodydnamics, calls the relevant
  ! subroutine, and returns the information to the vectors of prognostic
  ! variables.
  !------------------------------------------------------------------------
  
  use o_param
  use o_mesh
  use o_elements 
  use o_array,  only: real_salt_flux
  use i_therm_parms
  use i_dyn_parms
  use i_array
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe

  implicit none
  real*8  :: h,hsn,A,fsh,flo,Ta,Td,p_sl,qa,rain,snow,runo,rsss,rsf,evap_in
  real*8  :: ug,ustar,T_oc,S_oc,h_ml,t,ch,ce,ch_i,ce_i,fw,ehf,evap
  real*8  :: ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout
  real*8  :: lat
  integer :: i, j

  rsss=ref_sss

  do i=1, ToDim_nod2d   
#ifdef use_cavity
     if(cavity_flag_nod2d(i)==1) cycle
#endif
     h       = m_ice(i)
     hsn     = m_snow(i)
     A       = a_ice(i)
     fsh     = shortwave(i)
     flo     = longwave(i)
     Ta      = Tair(i)

    if(wind_data_source=='GME' .or. wind_data_source=='COSMO' .or. wind_data_source(1:6)=='MPEH5C'.or. &
       wind_data_source(1:6)=='ERAINT') then
       Td   = Tdew(i)
       qa   = 1.e20
     elseif (wind_data_source=='NCEP'.or.wind_data_source(1:4)=='CORE'.or.wind_data_source(1:6)=='HadCM3'.or.wind_data_source(1:4)=='CFSR' &
          .or.wind_data_source(1:7)=='HadGem2'.or.wind_data_source(1:7)=='mpi-esm') then
       Td   = 1.e20
       qa   = shum(i) 
     else
      write(*,*)'STOP: please adjust treatment of humidity-vs-dew-point issue for your forcing fields'
      stop
     endif

    if(wind_data_source=='GME' .or. wind_data_source=='COSMO' .or. wind_data_source(1:6)=='HadCM3') then
       p_sl = Pair(i)
     elseif (wind_data_source=='NCEP'.or.wind_data_source(1:4)=='CORE'.or.wind_data_source(1:6)=='MPEH5C'.or.wind_data_source(1:4)=='CFSR'&
         .or.wind_data_source(1:6)=='ERAINT') then
       p_sl = 980.e2   ! [Pa] 
     else
      write(*,*)'STOP: please adjust treatment of surface pressure for your forcing fields'
      stop
     endif


    if(precip_data_source=='NCEP' .or. precip_data_source=='GME' .or. precip_data_source=='COSMO' &
                                   .or. precip_data_source(1:6)=='HadCM3' &
                                   .or. precip_data_source(1:7)=='HadGem2' &
                                   .or. precip_data_source(1:6)=='MPEH5C' &
                                   .or. precip_data_source(1:6)=='ERAINT' &
                                   .or.precip_data_source(1:4)=='CFSR' &
                                   .or. precip_data_source(1:7)=='mpi-esm') then
        if(Ta>=0.0) then
           rain=prec_rain(i)
           snow=0.0
        else
           rain=0.0
           snow=prec_rain(i)
        endif

        if (precip_data_source=='NCEP'.or.precip_data_source(1:6)=='HadCM3'&
        .or.precip_data_source(1:7)=='HadGem2'.or.precip_data_source(1:6)=='MPEH5C' &
        .or.precip_data_source(1:7)=='mpi-esm'.or.precip_data_source(1:6)=='ERAINT'.or.precip_data_source(1:4)=='CFSR') then 
           evap    = -evaporation(i)    ! evap < 0 from here on
           evap_in = -evaporation(i)    ! (positive down)
        elseif(precip_data_source=='COSMO' .or. precip_data_source=='GME') then
           evap    = 0.0
           evap_in = 0.0
        else
           write(*,*)'STOP: Please think about treatment of evaporation here.'
           stop
        endif
     elseif (precip_data_source(1:4)=='CORE') then
        rain    = prec_rain(i)
        snow    = prec_snow(i)
        evap    = 0.0
        evap_in = 0.0
     else
        write(*,*)'STOP: please adjust treatment of rain vs snow to your forcing fields'
        stop
     end if


     runo    = runoff(i)
     ug      = sqrt(u_wind(i)**2+v_wind(i)**2)
     ustar   = sqrt(Cd_oce_ice)*sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
     T_oc    = T_oc_array(i)      
     S_oc    = S_oc_array(i)
     if(ref_sss_local) rsss = S_oc
     t       = t_skin(i)   
     ch	     = Ch_atm_oce_arr(i)
     ce	     = Ce_atm_oce_arr(i)
     ch_i    = Ch_atm_ice
     ce_i    = Ce_atm_ice
!rt     h_ml    = 10.       	         ! 10.0 or 30. used previously
     h_ml    = coord_nod3d(3,nod3d_below_nod2d(1,i))-coord_nod3d(3,nod3d_below_nod2d(2,i)) !RT: RG44006
     fw      = 0.
     ehf     = 0.

     call therm_ice(h,hsn,A,fsh,flo,Ta,qa,Td,p_sl,rain,snow,runo,rsss, &
          ug,ustar,T_oc,S_oc,h_ml,t,dt,ch,ce,ch_i,ce_i,evap_in,fw,ehf,evap, &
          rsf, ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout)

     m_ice(i)         = h
     m_snow(i)        = hsn
     a_ice(i)         = A
     t_skin(i)        = t
     fresh_wa_flux(i) = fw      !positive down
     net_heat_flux(i) = ehf     !positive down
     evaporation(i)   = -evap   !positive up
#ifdef use_fullfreesurf
     real_salt_flux(i)= rsf
#endif         

     thdgr(i)         = ithdgr
     thdgrsn(i)       = ithdgrsn
     flice(i)         = iflice
     olat_heat(i)     = hflatow
     osen_heat(i)     = hfsenow
     olwout(i)        = hflwrdout

  end do

end subroutine thermodynamics
!
!===================================================================
!
subroutine therm_ice(h,hsn,A,fsh,flo,Ta,qa,Td,p_sl,rain,snow,runo,rsss, &
     ug,ustar,T_oc,S_oc,H_ML,t,dt,ch,ce,ch_i,ce_i,evap_in,fw,ehf,evap, &
     rsf, dhgrowth, dhsngrowth, iflice, hflatow, hfsenow, hflwrdout)
  ! Ice Thermodynamic growth model     
  !
  ! Input parameters:
  !------------------
  ! h - ice mass [m]
  ! hsn - snow mass [m]
  ! A - ice compactness
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! Ta - air temperature
  ! qa - specific humidity
  ! rain - precipitation rain
  ! snow - precipitation snow
  ! runo - runoff
  ! ug - wind speed
  ! ustar - friction velocity
  ! T_oc, S_oc - ocean temperature and salinity beneath the ice (mixed layer)
  ! H_ML - mixed layer depth - should be specified.
  ! t - temperature of snow/ice top surface
  ! dt - time step [s]
  ! ch - transfer coefficient for sensible heat (for open ocean)
  ! ce - transfer coefficient for evaporation   (for open ocean)
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice)  

  ! Output parameters:
  !-------------------
  ! h - ice mass
  ! hsn - snow mass
  ! A - ice compactness
  ! t - temperature of snow/ice top surface
  ! fw - freshwater flux due to ice melting [m water/dt]
  ! ehf - net heat flux at the ocean surface [W/m2]        !RTnew

  use i_therm_parms
  use g_forcing_param,  only: precip_data_source
  implicit none

  integer k
  real*8  h,hsn,A,fsh,flo,Ta,Td,p_sl,qa,rain,snow,runo,rsss,evap_in
  real*8  ug,ustar,T_oc,S_oc,H_ML,t,dt,ch,ce,ch_i,ce_i,fw,ehf
  real*8  dhgrowth,dhsngrowth,ahf,prec,subli,subli_i,rsf
  real*8  rhow,show,rhice,shice,sh,thick,thact,lat
  real*8  rh,rA,qhst,sn,hsntmp,o2ihf,evap
  real*8  iflice,hflatow,hfsenow,hflwrdout
  real*8, external  :: TFrez  ! Sea water freeze temperature.

  ! Store ice thickness at start of growth routine
  dhgrowth=h  	  

  ! determine h(i,j)/a(i,j) = actual ice thickness.
  ! if snow layer is present, add hsn weighted with quotient
  ! of conductivities of ice and snow, according to 0-layer approach
  ! of Semtner (1976).   	    
  ! thickness at the ice covered part
  thick=hsn*(con/consn)/max(A,Armin)    ! Effective snow thickness
  thick=thick+h/max(A,Armin)            ! Effective total snow-ice thickness

  ! Growth rate for ice in open ocean
  rhow=0.0
  evap=0.0
!  call obudget(qa,fsh,flo,T_oc,ug,ta,ch,ce,rhow,evap,hflatow,hfsenow,hflwrdout) 
  call obudget_rt(T_oc,ug,Ta,Td,qa,p_sl,fsh,flo,A,ch,ce,rhow,evap,hflatow,hfsenow,hflwrdout)

  hflatow=hflatow*(1.0-A)
  hfsenow=hfsenow*(1.0-A)
  hflwrdout=hflwrdout*(1.0-A)

  ! add heat loss at open ocean due to melting snow fall
  !rhow=rhow+snow*1000.0/rhoice !qiang
  ! dt and (1-A) will be multiplied afterwards

  ! growth rate of ice in ice covered part
  ! following Hibler 1984
  ! assuming ice thickness has an euqal, 7-level distribution from zero to two times h 
  rhice=0.0                      
  subli=0.0
  if (thick.gt.hmin) then
     do k=1,iclasses			  
        thact = (2*k-1)*thick/float(iclasses)  	! Thicknesses of actual ice class
        call budget(thact,hsn,t,Ta,qa,Td,p_sl,fsh,flo,ug,S_oc,ch_i,ce_i,shice,subli_i) 
        !Thick ice K-class growth rate
        rhice=rhice+shice/float(iclasses)      	! Add to average heat flux
        subli=subli+subli_i/float(iclasses)
     end do
  end if

  ! Convert growth rates [m ice/sec] into growth per time step DT.
  rhow=rhow*dt
  rhice=rhice*dt

  ! Multiply ice growth of open water and ice
  ! with the corresponding areal fractions of grid cell
  show =rhow*(1.0-A)
  shice=rhice*A
  sh   =show+shice

  ! Store atmospheric heat flux, average over grid cell [W/m**2]
  ahf=-cl*sh/dt   

  ! precipitation (into the ocean)
  prec=rain+runo+snow*(1.0-A)  	        ! m water/s

  ! snow fall above ice
  hsn=hsn+snow*dt*A*rhowat/rhosno	! Add snow fall to temporary snow thickness    !!!
  dhsngrowth=hsn   		        ! Store snow thickness after snow fall 

  ! evaporation/sublimation water flux  (m water /s)
  if (precip_data_source=='NCEP'.or.precip_data_source(1:6)=='HadCM3'&
        .or.precip_data_source(1:7)=='HadGem2'.or.precip_data_source(1:7)=='mpi-esm'.or.precip_data_source(1:6)=='ERAINT'.or.precip_data_source(1:4)=='CFSR') then
     evap  = evap_in*(1.0-A)    ! use prescribed evaporation for evap water flux
     subli = evap_in*A
!RG44001     evap  = evap_in                   ! use prescribed evaporation for evap water flux
!RG44001     subli = 0.
  else
     write(*,*)'STOP. Are you sure you want to compute evaporation from the latent heat flux?'
     stop
     evap  = evap  * (1.0-A)    		! use latent heat flux computed in model
     subli = subli * A
  end if

  ! If there is atmospheric melting, first melt any snow that is present.
  ! Atmospheric heat flux available for melting
  ! sh = MINUS atm. heat flux / specific latent heat of sea ice
  ! Note: (sh<0) for melting, (sh>0) for freezing
  hsntmp= -min(sh,0.0)*rhoice/rhosno

  ! hsntmp is the decrease in snow thickness due to atmospheric melting
  ! [m/DT]. Do not melt more snow than available
  hsntmp=min(hsntmp,hsn)
  hsn=hsn-hsntmp  ! Update snow thickness after atmospheric snow melt

  ! Negative atmospheric heat flux left after melting of snow
  ! Note: (sh<0) and (hsntmp>0) for melting conditions
  ! hsntmp=0 for non-snow-melting conditions
  rh=sh+hsntmp*rhosno/rhoice
  h=max(h,0.)

  ! Compute heat content qhst of mixed layer - sea ice system
  !
  ! Total heat content is the sum of
  !	h	ice thickness after calculation of dynamic effects
  !	rh	change in ice thickness due to atmospheric forcing
  ! and heat available in mixed layer, with
  !	T_oc	temperature of ocean surface layer
  !	Tfrez	freezing point of sea water
  !	H_ML	thickness of uppermost layer
  !
  !RT:
  ! There are three possibilities to do this.
  ! 1.: Assume an instantaneous adjustment of mixed layer heat content.
  !     Any heat available is then instantaneously used to melt ice.
  !     (so-called ice-bath approach)
  !     This is what used to be used in the Lemke sea ice-mixed layer model.
  ! rh=rh-(T_oc-TFrez(S_oc))*H_ML*cc/cl
  ! qhst=h+rh 
  !
  ! 2.: Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference. For a first step 
  !     we can assume a constant exchange coefficient gamma_t:
  ! o2ihf= (T_oc-TFrez(S_oc))*gamma_t*cc*A     &
  !        +(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A) ! [W/m2]
  ! rh=rh-o2ihf*dt/cl
  ! qhst=h+rh		                      	! [m]
  !
  ! 3.  Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference and the
  !     friction velocity:
  o2ihf= (T_oc-TFrez(S_oc))*0.012*ustar*cc*A  &
       +(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A)  	! [W/m2]
  rh=rh-o2ihf*dt/cl
  qhst=h+rh		              		! [m]

  ! Melt snow if there is any ML heat content left (qhst<0).
  ! This may be the case if advection moves ice (with snow) to regions
  ! with a warm mixed layer.
  sn=hsn+min(qhst,0.)*rhoice/rhosno

  ! New temporary snow thickness must not be negative:
  sn=max(sn,0.)

  ! Update snow and ice depth
  hsn=sn
  h=max(qhst,0.)     
  if (h.lt.1E-6) h=0.        ! Avoid very small ice thicknesses

  ! heat and fresh water fluxes
  dhgrowth=h-dhgrowth        ! Change in ice thickness due to thermodynamic effects
  dhsngrowth=hsn-dhsngrowth  ! Change in snow thickness due to thermodynamic melting
  ! (without snow fall). This is a negative value (MINUS snow melt)

  dhgrowth=dhgrowth/dt       ! Conversion: 'per time step' -> 'per second'
  dhsngrowth=dhsngrowth/dt   ! Conversion: 'per time step' -> 'per second'

  ! (radiation+turbulent) + freezing(-melting) sea-ice&snow 
  ehf = ahf + cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth) 
  ! (prec+runoff)+evap - freezing(+melting) ice&snow
#ifdef use_fullfreesurf
  fw= prec+evap - dhgrowth*rhoice/rhowat - dhsngrowth*rhosno/rhowat 
  rsf= -dhgrowth*rhoice/rhowat*Sice
#else
  fw= prec+evap - dhgrowth*rhoice/rhowat*(rsss-Sice)/rsss - dhsngrowth*rhosno/rhowat 
#endif

  ! Changes in compactnesses (equation 16 of Hibler 1979)
  rh=-min(h,-rh)   ! Make sure we do not try to melt more ice than is available
  rA= rhow - o2ihf*dt/cl !Qiang: it was -(T_oc-TFrez(S_oc))*H_ML*cc/cl, changed in June 2010
  !rA= rhow - (T_oc-TFrez(S_oc))*H_ML*cc/cl*(1.0-A)
!rt  A=A + 0.5*min(rh,0.)*A/max(h,hmin) + max(rA,0.)*(1.-A)/h0   
  A=A + 0.1*min(rh,0.)*A/max(h,hmin) + max(rA,0.)*(1.-A)/h0     ! RT
  !meaning:           melting                         freezing
 
  A=min(A,h*1.e6)     ! A -> 0 for h -> 0
  A=min(max(A,0.),1.) ! A >= 0, A <= 1

  ! Flooding (snow to ice conversion)
  iflice=h
  call flooding(h,hsn)     
  iflice=(h-iflice)/dt

  ! to maintain salt conservation for the current model version
  !(a way to avoid producing net salt from snow-type-ice) 
#ifdef use_fullfreesurf
  rsf=rsf-iflice*rhoice/rhowat*Sice
#else
  fw=fw+iflice*rhoice/rhowat*Sice/rsss
#endif   

  return
end subroutine therm_ice
!
!=====================================================================================
!
subroutine budget (hice,hsn,t,Ta,qa,Td,p_sl,fsh,flo,ug,S_oc,ch_i,ce_i,fh,subli)
! from call budget(thact,hsn,t,Ta,qa,Td,p_sl,fsh,flo,ug,S_oc,ch_i,ce_i,shice,subli_i)

  ! Thick ice growth rate [m ice/sec]
  !
  ! INPUT:
  ! hice - actual ice thickness [m]
  ! hsn - snow thickness, used for albedo parameterization [m]
  ! t - temperature of snow/ice surface [C]
  ! ta - air temperature [C]
  ! qa - specific humidity [Kg/Kg]
  ! fsh - shortwave radiation [W/m**2]
  ! flo - longwave radiation  [W/m**2]
  ! ug - wind speed [m/sec]
  ! S_oc - ocean salinity for the temperature of the ice base calculation [ppt]
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice) 
  !
  ! OUTPUT: fh - growth rate
  !
  ! qiang: The formular for saturated humidity was modified according to Large/Yeager2004
  ! to allow direct comparison with the CORE results (Griffies et al. 2009). The new
  ! formular does not require sea level pressure.
  ! A similar change was also made for the obudget routine.
  ! It was found through experiments that the results are quite similar to that from the
  ! original code, and the simulated ice volume is only slightly larger after modification. 
  
  use i_therm_parms
  use g_forcing_param
  
  implicit none

  integer iter, imax      ! Number of iterations
  real*8  hice,hsn,t,ta,td,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh,p_sl
  real*8  hfsen,hfrad,hflat,hftot,subli         
  real*8  alb             ! Albedo of sea ice
  real*8  q1, q2	  ! coefficients for saturated specific humidity
  real*8  A1,A2,A3,B,C, d1, d2, d3   
  real*8, external :: TFrez

  data q1 /11637800.0/, q2 /-5897.8/ 
  data imax /5/

  if(wind_data_source=='GME' .or. wind_data_source=='COSMO' .or. wind_data_source(1:6)=='MPEH5C' .or. wind_data_source(1:6)=='ERAINT') then
  ! prescribe dew point temperature and
  ! calculate specific humidities in the air and over water [kg/kg]
   qa=(0.622/(p_sl*0.01))*6.11*exp(21.8746*td/(td+265.5))
  elseif (wind_data_source=='NCEP'.or.wind_data_source(1:4)=='CORE'.or.wind_data_source(1:6)=='HadCM3'.or.wind_data_source(1:4)=='CFSR' & 
     .or.wind_data_source(1:7)=='HadGem2'.or.wind_data_source(1:7)=='mpi-esm') then
   continue   ! specific humidity qa already exists
  else
   write(*,*)'STOP: please adjust treatment of humidity-vs-dew-point issue before calling obudget'
   stop
  endif



  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished.
  if (t<0.0) then	        ! freezing condition    
     if (hsn.gt.0.0) then	!   snow cover present  
        alb=albsn         	
     else              		!   no snow cover       
        alb=albi       	
     endif
  else			        ! melting condition     
     if (hsn.gt.0.0) then	!   snow cover present  
        alb=albsnm	    	
     else			!   no snow cover       
        alb=albim		
     endif
  endif

  d1=rhoair*cpair*Ch_i
  d2=rhoair*Ce_i
  d3=d2*clhi

  ! total incoming atmospheric heat flux
  A1=(1.0-alb)*fsh + flo + d1*ug*ta + d3*ug*qa   ! in LY2004 emiss is multiplied wiht flo

  ! NEWTON-RHAPSON TO GET TEMPERATURE AT THE TOP OF THE ICE LAYER

  do iter=1,imax

     B=q1/rhoair*exp(q2/(t+tmelt))		! (saturated) specific humidity over ice
     A2=-d1*ug*t-d3*ug*B &
          -emiss_ice*boltzmann*((t+tmelt)**4)	! sensible and latent heat,and outward radiation
     A3=-d3*ug*B*q2/((t+tmelt)**2)		! gradient coefficient for the latent heat part
     C=con/hice                     		! gradient coefficient for downward heat conductivity
     A3=A3+C+d1*ug & 			! gradient coefficient for sensible heat and radiation 
          +4.0*emiss_ice*boltzmann*((t+tmelt)**3)    
     C=C*(TFrez(S_oc)-t)       		! downward conductivity term

     t=t+(A1+A2+C)/A3 		        ! NEW ICE TEMPERATURE AS THE SUM OF ALL COMPONENTS
  end do

  t=min(0.0,t)

  ! heat fluxes [W/m**2]:

  hfrad= (1.0-alb)*fsh &	        ! absorbed short wave radiation
       +flo &           	        ! long wave radiation coming in  ! in LY2004 emiss is multiplied
       -emiss_ice*boltzmann*((t+tmelt)**4) 	! long wave radiation going out

  hfsen=d1*ug*(ta-t) 			! sensible heat 
  subli=d2*ug*(qa-B) 			! sublimation
  hflat=clhi*subli                     	! latent heat

  hftot=hfrad+hfsen+hflat               ! total heat

  fh= -hftot/cl                         ! growth rate [m ice/sec]
  !                                      	+: ML gains energy, ice melts
  !                                      	-: ML loses energy, ice grows
  subli=subli/rhowat                    ! negative upward

  return
end subroutine budget
!
!======================================================================================

subroutine obudget_rt(t,ug,Ta,Td,qa,p_sl,fsh,flo,A,ch,ce,rhow,evap,hflatow,hfsenow,hflwrdout)
! Ice growth rate for open ocean [m ice/sec]
!
! INPUT:
! t - temperature of open water [C]
! fsh - shortwave radiation
! flo - longwave radiation
! ta - air temperature [C]
! td - dew point temperature [C]
! p_sl - sea level pressure
! ug - wind speed [m/sec]
! ch - transfer coefficient for sensible heat
! ce - transfer coefficient for evaporation
! A  - ice concentration
!
! OUTPUT: rhow - growth rate
!         evap - evaporation
!
! Ralph Timmermann for RG4125, 30.10.2011
!-----------------------------------------------------------------------------------------------
use i_therm_parms
use g_forcing_param

implicit none

real*8 t,Ta,Td,p_sl,fsh,flo,ug,ch,ce,rhow,evap,qa,A
real*8 hfsenow,hfradow,hflatow,hftotow,hflwrdout,b
real*8 q1, q2           ! coefficients for saturated specific humidity
real*8 tqw, qsw         ! coefficients for saturated specific humidity

data q1 /640380./, q2 /-5107.4/     ! for LY2004
data tqw /237.3/, qsw /17.2694/     ! adopted from BRIOS


if(wind_data_source=='GME' .or. wind_data_source=='COSMO' .or. wind_data_source(1:6)=='MPEH5C' .or. wind_data_source(1:6)=='ERAINT') then
 ! prescribe dew point temperature and
 ! calculate specific humidities in the air and over water [kg/kg]
 ! SLP needs to be converted to hPa.
 qa = 0.622/(p_sl/100.)*6.11*exp(qsw*td/(td+tqw))
  b = 0.622/(p_sl/100.)*6.11*exp(qsw*t /(t +tqw))
elseif (wind_data_source=='NCEP'.or.wind_data_source(1:4)=='CORE'.or.wind_data_source(1:6)=='HadCM3'.or.wind_data_source(1:4)=='CFSR') then 
 ! specific humidity qa in 2m already exists; now
 ! calculate (saturated) surface specific humidity
  b = 0.98*q1/rhoair*exp(q2/(t+tmelt))            ! LY2004
else
 write(*,*)'STOP: please adjust treatment of humidity-vs-dew-point issue before calling obudget'
 stop
endif

!heat fluxes [W/m**2]:

hfsenow=rhoair*cpair*ch*ug*(ta-t)             ! sensible heat
hflatow=rhoair*clhw *ce*ug*(qa-b)             ! latent heat W/m2 = latent heat * evap rate

hfradow= (1.0-albw)*fsh &                         ! absorbed short wave radiation
       +flo*emiss_wat                             ! long wave radiation coming in !put emiss/check
hflwrdout=-emiss_wat*boltzmann*((t+tmelt)**4)     ! long wave radiation going out !in LY2004 emiss=1
hfradow=hfradow+hflwrdout

if (precip_data_source=='CORE') then
 evap=rhoair*ce*ug*(qa-b)                     ! evaporation kg/m2/s calculated from latent heat flux
 evap=evap/1000.                              ! evaporation rate [m water/s],negative up !!!
 evap=evap*(1.0-A)                            ! m water/s
elseif (precip_data_source=='NCEP'.or.precip_data_source(1:6)=='HadCM3' &
                                  .or.precip_data_source(1:6)=='MPEH5C'.or.precip_data_source(1:4)=='CFSR'.or.precip_data_source(1:6)=='ERAINT') then
 continue                                     ! use prescribed evaporation from forcing fields
else
 write(*,*)'STOP: Please think about how you want to treat evaporation in subroutine obudget'
 stop
endif

hftotow=hfradow+hfsenow+hflatow                 ! total heat W/m2
rhow   = -hftotow/cl                            ! growth rate [m ice/sec]
  !                                             +: ML gains energy, ice melts
  !                                             -: ML loses energy, ice grows

return
end subroutine obudget_rt

!======================================================================================
!
subroutine obudget (qa,fsh,flo,t,ug,ta,ch,ce,fh,evap,hflatow,hfsenow,hflwrdout)  
  ! Ice growth rate for open ocean [m ice/sec]
  !
  ! INPUT:
  ! t - temperature of open water [C]
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! ta - air temperature [C]
  ! qa  - specific humidity             
  ! ug - wind speed [m/sec]
  ! ch - transfer coefficient for sensible heat
  ! ce - transfer coefficient for evaporation
  !
  ! OUTPUT: fh - growth rate
  !         evap - evaporation

  use i_therm_parms
  implicit none

  real*8 qa,t,Ta,fsh,flo,ug,ch,ce,fh,evap
  real*8 hfsenow,hfradow,hflatow,hftotow,hflwrdout,b
  real*8 q1, q2 		! coefficients for saturated specific humidity
  real*8 c1, c4, c5
  logical :: standard_saturation_shum_formula=.true.


  !data c1, c4, c5 /3.8e-3, 17.67, 243.5/
  data c1, c4, c5 /3.8e-3, 17.27, 237.3/
  data q1 /640380./, q2 /-5107.4/

  ! (saturated) surface specific humidity
  if(standard_saturation_shum_formula) then
     b=c1*exp(c4*t/(t+c5))                      ! a standard one
  else
     b=0.98*q1/rhoair*exp(q2/(t+tmelt)) 	! LY2004 NCAR version 
  end if

  ! heat fluxes [W/m**2]:

  hfradow= (1.0-albw)*fsh &	            	! absorbed short wave radiation
       +flo             	                ! long wave radiation coming in !put emiss/check
  hflwrdout=-emiss_wat*boltzmann*((t+tmelt)**4) 	! long wave radiation going out !in LY2004 emiss=1
  hfradow=hfradow+hflwrdout

  hfsenow=rhoair*cpair*ch*ug*(ta-t)             ! sensible heat 
  evap=rhoair*ce*ug*(qa-b)			! evaporation kg/m2/s
  hflatow=clhw*evap                       	! latent heat W/m2

  hftotow=hfradow+hfsenow+hflatow             	! total heat W/m2

  fh= -hftotow/cl                             	! growth rate [m ice/sec]
  !                                           	+: ML gains energy, ice melts
  !                                           	-: ML loses energy, ice grows
  evap=evap/rhowat 	 			! evaporation rate [m water/s],negative up !!!

  return
end subroutine obudget
!
!======================================================================================
!
subroutine flooding (h,hsn)
  use i_therm_parms

  real*8 h,hsn,hdraft,hflood

  hdraft=(rhosno*hsn+h*rhoice)/rhowat ! Archimedes: displaced water
  hflood=hdraft-min(hdraft,h)         ! Increase in mean ice thickness due to flooding
  h=h+hflood                          ! Add converted snow to ice volume
  hsn=hsn-hflood*rhoice/rhosno        ! Subtract snow from snow layer

  !RT   This is what all AWI sea ice models do, but
  !RT   I wonder whether it really is correct for the heat budget.
  !RT   I suggest we initially keep it to allow for a comparison with BRIOS results
  !RT   and rethink it at a later stage.

  return
end subroutine flooding
!
!======================================================================================
!
function TFrez(S)
  ! Nonlinear correlation for the water freezing temperature.
  ! Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
  implicit none
  real*8 S, TFrez

  TFrez= -0.0575*S+1.7105e-3 *sqrt(S**3)-2.155e-4 *S*S

end function TFrez



