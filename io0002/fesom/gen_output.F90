subroutine init_output
  ! Initialize output files for snapshots
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-----------------------------------------------------------  

  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  character(100)            :: longname
  character(100)            :: filename
  character(1)              :: trind

  !OR needed to initialize every three month for ice sheet coupling
  if(yearnew==yearold) return 
  !if((dayold /= 91 .or. dayold/= 181 .or. dayold/= 271 .or. dayold /=1) .and. &
  !   timeold /= 0.0) return
  !OR-

  if (mype/=0) return

  write(*,*) 'initialize new output files'

  ! first, snapshots

  ! ocean

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'ssh', NF_DOUBLE, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  dimids(1) = dimid_3d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_DOUBLE, 2, dimids, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_DOUBLE, 2, dimids, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w', NF_DOUBLE, 2, dimids, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
  status = nf_def_var(ncid, 'wpot', NF_DOUBLE, 2, dimids, wpot_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  status = nf_def_var(ncid, 'temp', NF_DOUBLE, 2, dimids, tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt', NF_DOUBLE, 2, dimids, tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'ptr'//trind, NF_DOUBLE, 2, dimids, &
             tra_varid(index_passive_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'age'//trind, NF_DOUBLE, 2, dimids, &
             tra_varid(index_age_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='sea surface elevation'
  status = nf_put_att_text(ncid, ssh_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, ssh_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_put_att_text(ncid, u_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, u_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_put_att_text(ncid, v_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, v_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
  longname='vertical velocity potential'
  status = nf_put_att_text(ncid, wpot_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, wpot_varid, 'units', 5, 'm.m/s')
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  longname='potential temperature'
  status = nf_put_att_text(ncid, tra_varid(1), 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(1), 'units', 4, 'degC')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='salinity'
  status = nf_put_att_text(ncid, tra_varid(2), 'description', len_trim(longname), longname) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(2), 'units', 3, 'psu')
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        longname='passive tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        !status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
        !'units', 3, 'NaN')
        !if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        longname='age tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'units', 4, 'Year')
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! ice

#ifdef use_ice
  filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'area', NF_DOUBLE, 2, dimids, area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hice', NF_DOUBLE, 2, dimids, hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hsnow', NF_DOUBLE, 2, dimids, hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'uice', NF_DOUBLE, 2, dimids, uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'vice', NF_DOUBLE, 2, dimids, vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='ice concentration [0 to 1]'
  status = nf_PUT_ATT_TEXT(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective ice thickness'
  status = nf_PUT_ATT_TEXT(ncid, hice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hice_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective snow thickness'
  status = nf_PUT_ATT_TEXT(ncid, hsnow_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hsnow_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_PUT_ATT_TEXT(ncid, uice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, uice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_PUT_ATT_TEXT(ncid, vice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

#endif


  ! second, mean fields
#if defined(allow_calcmeans) || defined(allow_diag)
  call init_output_mean
#endif

  ! third, mesh diagnose
#ifdef allow_diag
  if(diag_mesh) then
     call init_output_mesh
  end if
#endif


  ! initialize the counter for saving results
  save_count=1

end subroutine init_output
!
!----------------------------------------------------------------------------
!
subroutine init_output_mean
  ! initialize output files for diagnose variables
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------

  use o_param
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid
  integer                   :: utemp_varid, vtemp_varid
  integer                   :: usalt_varid, vsalt_varid
  integer                   :: mixlay_varid, Kv_varid
  integer                   :: uu_varid, vv_varid
  integer                   :: rho_varid, urho_varid
  integer                   :: vrho_varid, uv_varid
  integer                   :: sgs_u_varid, sgs_v_varid, sgs_ut_varid
  integer                   :: sgs_vt_varid, sgs_us_varid, sgs_vs_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: thdgr_varid, thdgrsn_varid
  integer                   :: uhice_varid, vhice_varid
  integer                   :: uhsnow_varid, vhsnow_varid
  integer                   :: flice_varid, tair_varid, tdew_varid
  integer                   :: shum_varid, uwind_varid, vwind_varid
  integer                   :: rain_varid, snow_varid, runoff_varid
  integer                   :: evap_varid, lwrd_varid, swrd_varid
  integer                   :: qnet_varid, wnet_varid
  integer                   :: olat_varid, osen_varid, olwout_varid
  integer                   :: virtual_salt_varid, relax_salt_varid
  integer                   :: stress_x_varid, stress_y_varid
  integer                   :: Tsurfmean_varid, Ssurfmean_varid
  character(100)            :: longname
  character(100)            :: filename
  character(1)              :: trind

  if (mype/=0) return

#ifdef allow_calcmeans

  ! ocean

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.mean.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'ssh', NF_FLOAT, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  dimids(1) = dimid_3d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_FLOAT, 2, dimids, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_FLOAT, 2, dimids, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w', NF_FLOAT, 2, dimids, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'temp', NF_FLOAT, 2, dimids, tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt', NF_FLOAT, 2, dimids, tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'ptr'//trind, NF_FLOAT, 2, dimids, &
             tra_varid(index_passive_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'age'//trind, NF_FLOAT, 2, dimids, &
             tra_varid(index_age_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='mean sea surface elevation'
  status = nf_put_att_text(ncid, ssh_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, ssh_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean zonal velocity'
  status = nf_put_att_text(ncid, u_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, u_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean meridional velocity'
  status = nf_put_att_text(ncid, v_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, v_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean potential temperature'
  status = nf_put_att_text(ncid, tra_varid(1), 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(1), 'units', 4, 'degC')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean salinity'
  status = nf_put_att_text(ncid, tra_varid(2), 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(2), 'units', 3, 'psu')
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        longname='passive tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        !status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
        !'units', 3, 'NaN')
        !if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        longname='age tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), & 
             'units', 4, 'Year')
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! ice
#ifdef use_ice

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.mean.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'area', NF_FLOAT, 2, dimids, area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hice', NF_FLOAT, 2, dimids, hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hsnow', NF_FLOAT, 2, dimids, hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'uice', NF_FLOAT, 2, dimids, uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'vice', NF_FLOAT, 2, dimids, vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='ice concentration [0 to 1]'
  status = nf_PUT_ATT_TEXT(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective ice thickness'
  status = nf_PUT_ATT_TEXT(ncid, hice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hice_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective snow thickness'
  status = nf_PUT_ATT_TEXT(ncid, hsnow_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hsnow_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_PUT_ATT_TEXT(ncid, uice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, uice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_PUT_ATT_TEXT(ncid, vice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

#endif
#endif

  ! diagnose
#ifdef allow_diag

  ! ocean
  if(diag_oce) then

     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     if(diag_oce_mix_layer) then
        status = nf_def_var(ncid, 'mixlay', NF_FLOAT, 2, dimids, mixlay_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! Define the netCDF variables for 3D fields
     dimids(1) = dimid_3d
     dimids(2) = dimid_rec


     if(diag_oce_KE) then
        status = nf_def_var(ncid, 'uu', NF_FLOAT, 2, dimids, uu_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'vv', NF_FLOAT, 2, dimids, vv_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_energy_conv) then
        status = nf_def_var(ncid, 'rho', NF_FLOAT, 2, dimids, rho_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'urho', NF_FLOAT, 2, dimids, urho_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'vrho', NF_FLOAT, 2, dimids, vrho_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'uv', NF_FLOAT, 2, dimids, uv_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_transp) then
        status = nf_def_var(ncid, 'utemp', NF_FLOAT, 2, dimids, utemp_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'vtemp', NF_FLOAT, 2, dimids, vtemp_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'usalt', NF_FLOAT, 2, dimids, usalt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'vsalt', NF_FLOAT, 2, dimids, vsalt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(Redi_GM .and. diag_oce_GM_vel) then
        status = nf_def_var(ncid, 'sgs_u', NF_FLOAT, 2, dimids, sgs_u_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_v', NF_FLOAT, 2, dimids, sgs_v_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_SGS_transp) then
        status = nf_def_var(ncid, 'sgs_ut', NF_FLOAT, 2, dimids, sgs_ut_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_vt', NF_FLOAT, 2, dimids, sgs_vt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_us', NF_FLOAT, 2, dimids, sgs_us_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_vs', NF_FLOAT, 2, dimids, sgs_vs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_Kv) then
        status = nf_def_var(ncid, 'Kv', NF_FLOAT, 2, dimids, Kv_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! Assign long_name and units attributes to variables.

     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     if(diag_oce_mix_layer) then
        longname='mixed layer thickness'
        status = nf_put_att_text(ncid, mixlay_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, mixlay_varid, 'units', 1, 'm')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     !KE
     if(diag_oce_KE) then
        longname='u*u'
        status = nf_put_att_text(ncid, uu_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, uu_varid, 'units', 5, 'm2/s2')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='v*v'
        status = nf_put_att_text(ncid, vv_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, vv_varid, 'units', 5, 'm2/s2')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     !energy convert
     if(diag_oce_energy_conv) then
        longname='insitu density'
        status = nf_put_att_text(ncid, rho_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, rho_varid, 'units', 5, 'kg/m3')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='u * insitu density'
        status = nf_put_att_text(ncid, urho_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, urho_varid, 'units', 9, 'm/s*kg/m3')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='v * insitu density'
        status = nf_put_att_text(ncid, vrho_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, vrho_varid, 'units', 9, 'm/s*kg/m3')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='u*v'
        status = nf_put_att_text(ncid, uv_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, uv_varid, 'units', 5, 'm2/s2')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_transp) then
        longname='zonal advective flux of temperature (zonal velocity * temperature)'
        status = nf_put_att_text(ncid, utemp_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, utemp_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='meridional advective flux of temperature (meridional velocity * temperature)'
        status = nf_put_att_text(ncid, vtemp_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, vtemp_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='zonal advective flux of salinity (zonal velocity * salinity)'
        status = nf_put_att_text(ncid, usalt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, usalt_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='meridional advective flux of salinity (meridional velocity * salinity)'
        status = nf_put_att_text(ncid, vsalt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, vsalt_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(Redi_GM .and. diag_oce_GM_vel) then
        longname='SGS (GM) zonal velocity integrated from bottom (k*S_x)'
        status = nf_put_att_text(ncid, sgs_u_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_u_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS (GM) meridional velocity integrated from bottom (k*S_y)'
        status = nf_put_att_text(ncid, sgs_v_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_v_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_SGS_transp) then
        longname='SGS zonal temperature flux'
        status = nf_put_att_text(ncid, sgs_ut_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_ut_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS meridional temperature flux'
        status = nf_put_att_text(ncid, sgs_vt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_vt_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS zonal salinity flux'
        status = nf_put_att_text(ncid, sgs_us_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_us_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS meridional salinity flux'
        status = nf_put_att_text(ncid, sgs_vs_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_vs_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_Kv) then
        longname='Instantaneous vertical diffusivity'
        status = nf_put_att_text(ncid, Kv_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, Kv_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     status = nf_enddef(ncid)  !end def
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)     !close file
     if (status .ne. nf_noerr) call handle_err(status)

  endif  !ocean

  ! ice

#ifdef use_ice
  if(diag_ice) then
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     status = nf_def_var(ncid, 'thdgr', NF_FLOAT, 2, dimids, thdgr_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'thdgrsn', NF_FLOAT, 2, dimids, thdgrsn_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'uhice', NF_FLOAT, 2, dimids, uhice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vhice', NF_FLOAT, 2, dimids, vhice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'uhsnow', NF_FLOAT, 2, dimids, uhsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vhsnow', NF_FLOAT, 2, dimids, vhsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'flice', NF_FLOAT, 2, dimids, flice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Assign long_name and units attributes to variables.
     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     longname='thermodynamic growth rate of eff. ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, thdgr_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, thdgr_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='melting rate of snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, thdgrsn_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, thdgrsn_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='zonal advective flux of eff. ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, uhice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uhice_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional advective flux of eff. ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, vhice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vhice_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='zonal advective flux of eff. snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, uhsnow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uhsnow_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional advective flux of eff. snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, vhsnow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vhsnow_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='rate of flooding snow to ice'
     status = nf_PUT_ATT_TEXT(ncid, flice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, flice_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

  endif
#endif


  ! forcing
#ifdef use_ice
  if(diag_forcing) then
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.forcing.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     status = nf_def_var(ncid, 'tair', NF_FLOAT, 2, dimids, tair_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'tdew', NF_FLOAT, 2, dimids, tdew_varid)   !RT MK44003
     if (status .ne. nf_noerr) call handle_err(status)                    !RT MK44003
     status = nf_def_var(ncid, 'shum', NF_FLOAT, 2, dimids, shum_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'uwind', NF_FLOAT, 2, dimids, uwind_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vwind', NF_FLOAT, 2, dimids, vwind_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'rain', NF_FLOAT, 2, dimids, rain_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'snow', NF_FLOAT, 2, dimids, snow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'runoff', NF_FLOAT, 2, dimids, runoff_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'evap', NF_FLOAT, 2, dimids, evap_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'lwrd', NF_FLOAT, 2, dimids, lwrd_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'swrd', NF_FLOAT, 2, dimids, swrd_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_var(ncid, 'qnet', NF_FLOAT, 2, dimids, qnet_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'olat', NF_FLOAT, 2, dimids, olat_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'osen', NF_FLOAT, 2, dimids, osen_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'olwout', NF_FLOAT, 2, dimids, olwout_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'wnet', NF_FLOAT, 2, dimids, wnet_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'virtual_salt', NF_FLOAT, 2, dimids, virtual_salt_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'relax_salt', NF_FLOAT, 2, dimids, relax_salt_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'stress_x', NF_FLOAT, 2, dimids, stress_x_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'stress_y', NF_FLOAT, 2, dimids, stress_y_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'Tsurf', NF_FLOAT, 2, dimids, Tsurfmean_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'Ssurf', NF_FLOAT, 2, dimids, Ssurfmean_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Assign long_name and units attributes to variables.

     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)

     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     longname='air temperature'
     status = nf_PUT_ATT_TEXT(ncid, tair_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, tair_varid, 'units', 4, 'degC')
     if (status .ne. nf_noerr) call handle_err(status)
!RT MK44003:
     longname='dew point temperature'
     status = nf_PUT_ATT_TEXT(ncid, tdew_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, tdew_varid, 'units', 4, 'degC')
     if (status .ne. nf_noerr) call handle_err(status)
!RT-
     longname='air specific humidity'
     status = nf_PUT_ATT_TEXT(ncid, shum_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, shum_varid, 'units', 5, 'kg/kg')
     if (status .ne. nf_noerr) call handle_err(status)

     longname='zonal wind speed'
     status = nf_PUT_ATT_TEXT(ncid, uwind_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uwind_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional wind speed'
     status = nf_PUT_ATT_TEXT(ncid, vwind_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vwind_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='precipitation rain'
     status = nf_PUT_ATT_TEXT(ncid, rain_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, rain_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='precipitation snow (in m/s water)'
     status = nf_PUT_ATT_TEXT(ncid, snow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, snow_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='runoff'
     status = nf_PUT_ATT_TEXT(ncid, runoff_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, runoff_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='evaporation'
     status = nf_PUT_ATT_TEXT(ncid, evap_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, evap_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='atmosphere longwave radiation'
     status = nf_PUT_ATT_TEXT(ncid, lwrd_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, lwrd_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='atmosphere shortwave radiation'
     status = nf_PUT_ATT_TEXT(ncid, swrd_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, swrd_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)

     longname='net heat flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, qnet_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, qnet_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='latent heat flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, olat_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, olat_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='sensible heat flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, osen_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, osen_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='longwave radiation from ocean, downward positve'
     status = nf_PUT_ATT_TEXT(ncid, olwout_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, olwout_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='net freshwater flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, wnet_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, wnet_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='virtual salt flux to ocean, >0 increase salinity'
     status = nf_PUT_ATT_TEXT(ncid, virtual_salt_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, virtual_salt_varid, 'units', 7, 'psu m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface salinity relaxation, >0 increase salinity'
     status = nf_PUT_ATT_TEXT(ncid, relax_salt_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, relax_salt_varid, 'units', 7, 'psu m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface zonal wind stress'
     status = nf_PUT_ATT_TEXT(ncid, stress_x_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, stress_x_varid, 'units', 5, 'N/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface meridional wind stress'
     status = nf_PUT_ATT_TEXT(ncid, stress_y_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, stress_y_varid, 'units', 5, 'N/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface temperature field'
     status = nf_PUT_ATT_TEXT(ncid, Tsurfmean_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, Tsurfmean_varid, 'units', 5, 'deg C')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface salinity field'
     status = nf_PUT_ATT_TEXT(ncid, Ssurfmean_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, Ssurfmean_varid, 'units', 3, 'psu')
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

  endif
#endif

#endif

end subroutine init_output_mean
!
!--------------------------------------------------------------------------------------------
!
subroutine write_snapshots
  ! write snapshots
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------

  use o_array
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use i_array
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: start(2), count(2), n3
  real(kind=8)              :: sec_in_year
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  allocate(aux2(nod2D), aux3(nod3D)) 
  n3=myDim_nod3D+eDim_nod3D             

  if (mype==0) then 

!     sec_in_year=dt*istep
   sec_in_year=float(daynew)*86400.   !RT  RG4166


     ! ocean

     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'ssh', ssh_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'u', u_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'v', v_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
     status=nf_inq_varid(ncid, 'wpot', wpot_varid)
     if (status .ne. nf_noerr) call handle_err(status)
#endif
     status=nf_inq_varid(ncid, 'temp', tra_varid(1))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'salt', tra_varid(2))
     if (status .ne. nf_noerr) call handle_err(status)

     if(use_passive_tracer) then
        do j=1,num_passive_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'ptr'//trind, tra_varid(index_passive_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     if(use_age_tracer) then
        do j=1,num_age_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'age'//trind, tra_varid(index_age_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)
  end if    !! mype==0     

  ! 2d fields
  call broadcast2D(ssh,aux2)  
  if(mype==0) then            
     start=(/1,save_count/)
     count=(/nod2d, 1/)
     status=nf_put_vara_double(ncid, ssh_varid, start, count, aux2) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if


  ! 3d fields
  call broadcast3D(uf(1:n3), aux3)   
  if (mype==0) then                  
     start=(/1,save_count/)
     count=(/nod3d, 1/)
     status=nf_put_vara_double(ncid, u_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call broadcast3D(uf(1+n3:2*n3),aux3)  
  if(mype==0) then                      
     status=nf_put_vara_double(ncid, v_varid, start, count, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

#ifdef use_non_hydrostatic
  call broadcast3D(uf(1+2*n3:3*n3), aux3)  
  if(mype==0) then                       
     status=nf_put_vara_double(ncid, w_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#else
  call broadcast3D(wrhs,aux3)             
  if(mype==0) then                        
     status=nf_put_vara_double(ncid, w_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast3D(w,aux3)             
  if(mype==0) then                     
     status=nf_put_vara_double(ncid, wpot_varid, start, count, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif
  do j=1,num_tracer
     call broadcast3D(tracer(:,j),aux3)    
     if(mype==0) then                      
        status=nf_put_vara_double(ncid, tra_varid(j), start, count, aux3) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do

  if(mype==0) then
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! ice

#ifdef use_ice
  if(mype==0) then                     
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'area', area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hice', hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'uice', uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'vice', vice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count/)
     count=(/nod2d, 1/)
  end if     !! mype=0                      
  call broadcast2D(a_ice,aux2)                
  if(mype==0) then                           
     status=nf_put_vara_double(ncid, area_varid, start, count, aux2)  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_ice,aux2)              
  if(mype==0) then                          
     status=nf_put_vara_double(ncid, hice_varid, start, count, aux2)   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_snow,aux2)              
  if(mype==0) then                            
     status=nf_put_vara_double(ncid, hsnow_varid, start, count, aux2)  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(u_ice,aux2)              
  if(mype==0) then                           
     status=nf_put_vara_double(ncid, uice_varid, start, count, aux2)    
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(v_ice,aux2)              
  if(mype==0) then                          
     status=nf_put_vara_double(ncid, vice_varid, start, count, aux2)    
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

  deallocate(aux3, aux2)
end subroutine write_snapshots
!
!--------------------------------------------------------------------------------------------
!
subroutine write_means_part1
  ! write mean arrays and diagnose variables
  ! SGS parameterizations are saved by write_means_part2
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------  
  
  use o_mesh
  use o_array
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid
  integer                   :: utemp_varid, vtemp_varid
  integer                   :: usalt_varid, vsalt_varid
  integer                   :: mixlay_varid, Kv_varid
  integer                   :: uu_varid, vv_varid
  integer                   :: rho_varid, urho_varid
  integer                   :: vrho_varid, uv_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: thdgr_varid, thdgrsn_varid
  integer                   :: uhice_varid, vhice_varid
  integer                   :: uhsnow_varid, vhsnow_varid
  integer                   :: flice_varid, tair_varid, tdew_varid
  integer                   :: shum_varid, uwind_varid, vwind_varid
  integer                   :: rain_varid, snow_varid, runoff_varid
  integer                   :: evap_varid, lwrd_varid, swrd_varid
  integer                   :: qnet_varid, wnet_varid
  integer                   :: olat_varid, osen_varid, olwout_varid
  integer                   :: virtual_salt_varid, relax_salt_varid
  integer                   :: stress_x_varid, stress_y_varid
  integer                   :: Tsurfmean_varid, Ssurfmean_varid
  integer                   :: start(2), count(2), n3
  real(kind=8)              :: sec_in_year
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  n3=myDim_nod3D+eDim_nod3D         
  allocate(aux2(nod2D), aux3(nod3D))  

!  sec_in_year=dt*istep
   sec_in_year=float(daynew)*86400.   !RT  RG4166


#ifdef allow_calcmeans 
  if (mype==0) then 

     ! ocean

     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.mean.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'ssh', ssh_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'u', u_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'v', v_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'temp', tra_varid(1))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'salt', tra_varid(2))
     if (status .ne. nf_noerr) call handle_err(status)

     if(use_passive_tracer) then
        do j=1,num_passive_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'ptr'//trind, tra_varid(index_passive_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     if(use_age_tracer) then
        do j=1,num_age_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'age'//trind, tra_varid(index_age_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count/)
     count=(/nod2d, 1/)
  end if    !! mype=0        

  call broadcast2D(sshmean, aux2)       
  if(mype==0) then
     status=nf_put_vara_real(ncid, ssh_varid, start, count, real(aux2,4)) 
     if (status .ne. nf_noerr) call handle_err(status)    
     ! 3d fields
     start=(/1,save_count/)
     count=(/nod3d, 1/)
  end if
  call broadcast3D(ufmean(1:n3),aux3)       
  if(mype==0) then                      
     status=nf_put_vara_real(ncid, u_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast3D(ufmean(1+n3:2*n3),aux3)   
  if(mype==0) then                     
     status=nf_put_vara_real(ncid, v_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#ifdef use_non_hydrostatic
  call broadcast3D(ufmean(1+2*n3:3*n3),aux3)   
  if(mype==0) then                      
     status=nf_put_vara_real(ncid, w_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#else
  call broadcast3D(wrhsmean, aux3)      
  if(mype==0) then
     status=nf_put_vara_real(ncid, w_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

  do j=1,num_tracer
     call broadcast3D(tracermean(:,j),aux3) 
     if(mype==0) then                    
        status=nf_put_vara_real(ncid, tra_varid(j), start, count, real(aux3,4))
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do

  if(mype==0) then                       
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if


  ! ice
#ifdef use_ice
  if(mype==0) then                      
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.mean.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'area', area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hice', hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'uice', uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'vice', vice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count/)
     count=(/nod2d, 1/)
  end if   !! mype=0                        

  call broadcast2D(a_ice_mean, aux2)        
  if(mype==0) then                         
     status=nf_put_vara_real(ncid, area_varid, start, count, real(aux2,4)) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_ice_mean, aux2)         
  if(mype==0) then                           
     status=nf_put_vara_real(ncid, hice_varid, start, count, real(aux2,4))   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_snow_mean, aux2)        
  if(mype==0) then                           
     status=nf_put_vara_real(ncid, hsnow_varid, start, count, real(aux2,4))   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(u_ice_mean, aux2)        
  if(mype==0) then                           
     status=nf_put_vara_real(ncid, uice_varid, start, count, real(aux2,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(v_ice_mean, aux2)       
  if(mype==0) then                          
     status=nf_put_vara_real(ncid, vice_varid, start, count, real(aux2,4))   
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

#endif
  ! mean arrays have been saved


#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(mype==0) then                          
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        if(diag_oce_mix_layer) then
           status=nf_inq_varid(ncid, 'mixlay', mixlay_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        endif

        if(diag_oce_KE) then
           status=nf_inq_varid(ncid, 'uu', uu_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'vv', vv_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        end if

        if(diag_oce_energy_conv) then
           status=nf_inq_varid(ncid, 'rho', rho_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'urho', urho_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'vrho', vrho_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'uv', uv_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        end if

        if(diag_oce_transp) then
           status=nf_inq_varid(ncid, 'utemp', utemp_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'vtemp', vtemp_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'usalt', usalt_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'vsalt', vsalt_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        endif

        if(diag_oce_Kv) then
           status=nf_inq_varid(ncid, 'Kv', Kv_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        end if

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if    !! mype=0                        

     if(diag_oce_mix_layer) then
        call broadcast2D(mixlay_dep_mean,aux2)    
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, mixlay_varid, start, count, real(aux2,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     endif

     ! 3d fields
     start=(/1,save_count/)
     count=(/nod3d, 1/)

     !KE
     if(diag_oce_KE) then
        call broadcast3D(uumean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, uu_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(vvmean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, vv_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     end if
     
     !energy convert
     if(diag_oce_energy_conv) then
        call broadcast3D(rhomean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, rho_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(urhomean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, urho_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(vrhomean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, vrho_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(uvmean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, uv_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     end if

     ! transport
     if(diag_oce_transp) then
        ! the fields which are read
        call broadcast3D(uTFmean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, utemp_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(vTFmean,aux3)          
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, vtemp_varid, start, count,real(aux3,4)) 
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(uSFmean,aux3)           
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, usalt_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(vSFmean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, vsalt_varid, start, count,real(aux3,4)) 
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     endif

     ! Kv
     if(diag_oce_Kv) then
        call broadcast3D(Kv,aux3)                  
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, Kv_varid, start, count, real(aux3,4))    
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     end if

     if(mype==0) then  
        status=nf_close(ncid)  !close file
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     
  endif  ! diag ocean

  ! ice
#ifdef use_ice
  if(diag_ice) then
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'thdgr', thdgr_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'thdgrsn', thdgrsn_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uhice', uhice_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vhice', vhice_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uhsnow', uhsnow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vhsnow', vhsnow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'flice', flice_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if    !! mype=0                        
     call broadcast2D(thdgr_mean,aux2)          
     if(mype==0) then                          
        status=nf_put_vara_real(ncid, thdgr_varid, start, count, real(aux2,4))
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(thdgrsn_mean,aux2)        
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, thdgrsn_varid, start, count, real(aux2,4))
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(uhice_mean,aux2)          
     if(mype==0) then                            
        status=nf_put_vara_real(ncid, uhice_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(vhice_mean,aux2)          
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, vhice_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(uhsnow_mean,aux2)          
     if(mype==0) then                            
        status=nf_put_vara_real(ncid, uhsnow_varid, start, count, real(aux2,4)) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(vhsnow_mean,aux2)         
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, vhsnow_varid, start, count, real(aux2,4)) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(flice_mean,aux2)          
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, flice_varid, start, count, real(aux2,4)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  endif
#endif


  ! forcing 
#ifdef use_ice
  if(diag_forcing) then
     ! open files
     if(mype==0) then 
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.forcing.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'tair', tair_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'tdew', tdew_varid)      !RT MK44003
        if (status .ne. nf_noerr) call handle_err(status)  !RT MK44003
        status=nf_inq_varid(ncid, 'shum', shum_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uwind', uwind_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vwind', vwind_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'rain', rain_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'snow', snow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'runoff', runoff_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'evap', evap_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'lwrd', lwrd_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'swrd', swrd_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'qnet', qnet_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'olat', olat_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'osen', osen_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'olwout', olwout_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'wnet', wnet_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'virtual_salt', virtual_salt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'relax_salt', relax_salt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'stress_x', stress_x_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'stress_y', stress_y_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'Tsurf', Tsurfmean_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'Ssurf', Ssurfmean_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if     ! mype=0                            
     call broadcast2D(tair_mean,aux2)                
     if(mype==0) then                                
        status=nf_put_vara_real(ncid, tair_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                           
!RT MK44003:
     call broadcast2D(tdew_mean,aux2)                
     if(mype==0) then                                
        status=nf_put_vara_real(ncid, tdew_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                           
!RT-
     call broadcast2D(shum_mean,aux2)               
     if(mype==0) then                               
        status=nf_put_vara_real(ncid, shum_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                            
     call broadcast2D(uwind_mean,aux2)               
     if(mype==0) then                               
        status=nf_put_vara_real(ncid, uwind_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(vwind_mean,aux2)                 
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, vwind_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(rain_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, rain_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(snow_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, snow_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(runoff_mean,aux2)                
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, runoff_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(evap_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, evap_varid, start, count, real(aux2,4))     
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(lwrd_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, lwrd_varid, start, count, real(aux2,4))     
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(swrd_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, swrd_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(qnet_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, qnet_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(olat_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, olat_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(osen_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, osen_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(olwout_mean,aux2)                
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, olwout_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(wnet_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, wnet_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(virtual_salt_mean,aux2)          
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, virtual_salt_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(relax_salt_mean,aux2)            
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, relax_salt_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               

     call broadcast2D(stress_x_mean,aux2)              
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, stress_x_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               

     call broadcast2D(stress_y_mean,aux2)              
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, stress_y_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     call broadcast2D(Tsurfmean,aux2)              
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, Tsurfmean_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     call broadcast2D(Ssurfmean,aux2)              
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, Ssurfmean_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(mype==0) then                                    
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

  endif
#endif

#endif

  deallocate(aux3, aux2)
end subroutine write_means_part1
!
!--------------------------------------------------------------------------------------------
!
subroutine write_means_part2
  ! write ocean SGS parameterizations
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------
    
  use o_param
  use o_mesh
  use o_elements
  use o_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: m, elem, elnodes(4)
  integer                   :: status, ncid, sgs_varid
  integer                   :: start(2), count(2)
  real(kind=8)              :: array_3d(nod3d)
  character(100)            :: filename

  ! prepare cluster volume
  ! use wrhs as a temporary array
  wrhs=0.0
  do elem=1,myDim_elem3d                                          
     elnodes=elem3d_nodes(:,elem)
     wrhs(elnodes)=wrhs(elnodes)+voltetra(elem)
  end do

  ! 3d fields
  start=(/1,save_count/)
  count=(/nod3d, 1/)

  if(Redi_GM .and. diag_oce_GM_vel) then
     ! processing
     call process_elem2node(1,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_u', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(2,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_v', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

  end if

  if(diag_oce_SGS_transp) then
     ! processing
     call process_elem2node(3,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_ut', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(4,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_vt', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(5,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_us', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(6,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_vs', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end if

end subroutine write_means_part2
!
!--------------------------------------------------------------------------------------------
!
subroutine output(directionflag)
  ! main output routine
  !
  ! Coded by Ralph Timmermann
  ! Modified by Qiang Wang for more diagnose output
  ! Reviewed by ??
  !--------------------------------------------------------------	

  use g_config
  use g_clock
  use g_diag
  use g_PARFE
  implicit none

  logical :: do_output=.false.
  integer :: directionflag

  !check whether we want to do output
  if (output_length_unit.eq.'y') then
     call annual_output(do_output)
  else if (output_length_unit.eq.'m') then 
     call monthly_output(do_output) 
  else if (output_length_unit.eq.'d') then
     call daily_output(do_output)  
  else if (output_length_unit.eq.'h') then
     call hourly_output(do_output) 
  else if (output_length_unit.eq.'s') then
     call step_output(do_output) 
  else
     write(*,*) 'You did not specify a supported outputflag.'
     write(*,*) 'The program will stop to give you opportunity to do it.'
     call par_ex
     stop
  endif

  if (directionflag.eq.1) do_output=.true.  

  if (.not.do_output) return

  ! write results
  if(mype==0) write(*,*)'Do output (netCDF) ...'
  call write_snapshots


  ! write mean fields
#if defined(allow_calcmeans) || defined(allow_diag)
  call compute_means
  call write_means_part1
#ifdef allow_diag
  if(diag_oce .and. (diag_oce_SGS_transp .or. diag_oce_GM_vel)) call write_means_part2
#endif

  call clean_meanarrays
#endif

  ! write updated mesh (cluster volume by now):
#ifdef allow_diag
  if(diag_mesh) then
     call write_mesh_diag
  end if
#endif

  save_count=save_count+1

  if(mype==0) write(*,*)'output done'
end subroutine output
!
!--------------------------------------------------------------------------------------------
!
subroutine annual_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if ((daynew == ndpyr) .and. (timenew==86400.)) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine annual_output
!
!--------------------------------------------------------------------------------------------
!
subroutine monthly_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (day_in_month==num_day_in_month(fleapyear,month) .and. &
       timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  end if

end subroutine monthly_output
!
!--------------------------------------------------------------------------------------------
!
subroutine daily_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (mod(daynew,output_length)==0 .and. timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine daily_output
!
!--------------------------------------------------------------------------------------------
!
subroutine hourly_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (mod(timenew, 3600.*output_length)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine hourly_output
!
!--------------------------------------------------------------------------------------------
!
subroutine step_output(do_output)
  !decides whether it's time to do output
  use g_config
  implicit none

  logical :: do_output

  if (mod(istep, output_length)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine step_output
!
!--------------------------------------------------------------------------------------------
!
subroutine handle_err(errcode)
  use g_parfe
  implicit none
  
#include "netcdf.inc" 
  
  integer errcode
  
  write(*,*) 'Error: ', nf_strerror(errcode)
  call par_ex
  stop
end subroutine handle_err
!
!--------------------------------------------------------------------------------------------
!
subroutine oce_out
  use o_param
  use o_MESH
  use o_array
  use g_PARFE
  use g_config
  implicit none
  !
  integer                    :: i, j, n3
  real(kind=8), allocatable  :: aux2(:), aux3(:)

  n3=myDim_nod3D+eDim_nod3D              
  allocate(aux2(nod2D), aux3(nod3D))     

  call broadcast2D(ssh,aux2)
  call broadcast3D(uf(1:n3),aux3)
  if (mype==0) then      
     write(*,*) 'writing ocean results (ASCII)'
     open(35,file='ssh.out')
     do i=1,nod3D
        write(35,'(1f9.5)') aux3(i)
        !write(35,'(1e10.3)') aux3(i)
     end do
  end if
  call broadcast3D(uf(1+n3:2*n3), aux3)   
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1f9.5)') aux3(i)
        !write(35,'(1e10.3)') aux3(i)
     end do
     !
     do i=1,nod2d
        write(35,'(1f9.5)') aux2(i)
        !write(35,'(1e10.3)') aux2(i)
     end do
     close(35)
  end if

  if(mype==0) open(36,file='TS.out')
  do j=1,num_tracer
     call broadcast3D(tracer(:,j),aux3)
     if(mype==0) then
        do i=1,nod3D
           write(36,'(1f9.5)') aux3(i)
           !write(36,'(1e12.4)') aux3(i)
        end do
     end if
  end do
  if(mype==0) close(36)

#ifndef use_non_hydrostatic
  call broadcast3D(wrhs,aux3)
  if(mype==0) then
     open(37,file='vvel.out') 
     do i=1,nod3D
        write(37,'(1e10.3)') aux3(i)
     end do
     close(37)
  end if
#else
  call broadcast3D(uf(2*n3+1:3*n3), aux3)
  if(mype==0) then
     open(37,file='vvel.out') 
     do i=1,nod3D
        write(37,'(1e10.3)') aux3(i)
     end do
     close(37)
  end if
#endif

  call broadcast3D(Kv, aux3)
  if(mype==0) then 
     open(38,file='mix_coeff.out')
     do i=1,nod3D
        write(38,'(1e10.3)') aux3(i)
     end do
     close(38)
  end if

  deallocate(aux3, aux2)

end subroutine oce_out
!
!--------------------------------------------------------------------------------------------
!
subroutine ice_out
  use o_MESH
  use i_array
  use g_parfe

  implicit none
  integer :: i
  real(kind=8), allocatable   :: aux2(:) 

  allocate(aux2(nod2D))         

  call broadcast2D(m_ice, aux2) 
  if (mype==0) then 
     write(*,*) 'writing ice results (ASCII)'
     open(35,file='m_ice.out') 
     do i=1,nod2D
        !write(35,'(1f9.5)')  aux2(i)
        write(35,'(1e10.3)') aux2(i)
     end do
     close(35)
  end if

  call broadcast2D(a_ice, aux2)
  if(mype==0) then
     open(36,file='a_ice.out') 
     do i=1,nod2D
        !write(36, '(1f9.5)') aux2(i)
        write(36,'(1e10.3)') aux2(i)
     end do
     close(36)
  end if

  call broadcast2D(m_snow, aux2)
  if(mype==0) then 
     open(37,file='m_snow.out') 
     do i=1,nod2D
        !write(37, '(1f9.5)') aux2(i)
        write(37,'(1e10.3)') aux2(i)
     end do
     close(37)
  end if

  call broadcast2D(u_ice, aux2)
  if(mype==0) then
     open(38,file='u_ice.out')
     do i=1,nod2D
        write(38, '(1f9.5)') aux2(i)
     end do
  end if
  call broadcast2D(v_ice, aux2)
  if(mype==0) then
     do i=1,nod2D
        write(38, '(1f9.5)') aux2(i)
     end do
     close(38)
  end if

  call broadcast2D(net_heat_flux, aux2)
  if(mype==0) then
     open(39,file='heat_water_flux.out')
     do i=1,nod2D
        write(39, *) aux2(i)
     end do
  end if
  call broadcast2D(fresh_wa_flux,aux2)
  if(mype==0) then
     do i=1,nod2D
        write(39, *) aux2(i)
     end do
     close(39)
  end if

  deallocate(aux2)

end subroutine ice_out
!
!--------------------------------------------------------------------------------------------
!
subroutine save_MY_vara
  ! MY output
  ! Better format for MY output/input should be used in the future.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------
  
  use o_MESH
  use o_array
  use o_mixing_my2p5_mod
  use g_PARFE
  use g_config
  implicit none

  integer :: i
  real(kind=8), allocatable  :: aux3(:)

  allocate(aux3(nod3D))  

  if (mype==0) then
     write(*,*) 'writing MY2.5 variables (ascII)'
     Write(*,*) 'Notice: Currently MY variables are only saved once at the end of a run.'
     ! Later we may update to save MY in Netcdf formate when required.

     open(35,file=trim(ResultPath)//runid//'_MY_restart.out')
  end if

  call broadcast3D(Kv,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(Av,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(Kq,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2b,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2l,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2lb,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  if(mype==0) close(35)

  deallocate(aux3)
end subroutine Save_MY_vara
