module g_read_CORE_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid

contains 

  subroutine read_CORE_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    io=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    io=nf_close(ncid)

    return
  end subroutine read_CORE_NetCDF
 ! -------------------------------------------------------------------------------
  
  subroutine read_GME_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=720, ncj=361  ! GME grid dimension
    integer, dimension(4)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    real(kind=4), dimension (nci,ncj)   :: ncdataaux
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,1,itime/)
    icount= (/nci,ncj,1,1/)
    
    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT OPEN FORCING FILE CORRECTLY in read_GME_NetCDF!!!!!'
       print*,'Error in opening netcdf file '//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT GET VARID CORRECTLY in read_GME_NetCDF!!!!!'
       print*, varid
       stop
    endif
    
    io=nf_get_vara_real(ncid,varid,istart,icount,ncdataaux)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY in read_GME_NetCDF!!!!!'
       stop
    endif
    ncdata=real(ncdataaux,8)
    
    io=nf_close(ncid)

    return
  end subroutine read_GME_NetCDF
 ! -------------------------------------------------------------------------------
  
  subroutine read_GMEb_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=720, ncj=361  ! GME grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    real(kind=4), dimension (nci,ncj)   :: ncdataaux
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)
    
    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT OPEN FORCING FILE CORRECTLY in read_GME_NetCDF!!!!!'
       print*,'Error in opening netcdf file '//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT GET VARID CORRECTLY in read_GME_NetCDF!!!!!'
       print*, varid
       stop
    endif
    
    io=nf_get_vara_real(ncid,varid,istart,icount,ncdataaux)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY in read_GME_NetCDF!!!!!'
       stop
    endif
    ncdata=real(ncdataaux,8)
    
    io=nf_close(ncid)

    return
  end subroutine read_GMEb_NetCDF
  
   ! -------------------------------------------------------------------------------
  
  subroutine read_COSMOa_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=600, ncj=180  ! COSMO interpolated grid dimension
    integer, dimension(4)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    real(kind=4), dimension (nci,ncj)   :: ncdataaux
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,1,itime/)
    icount= (/nci,ncj,1,1/)
    
    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT OPEN FORCING FILE CORRECTLY in read_COSMOa_NetCDF!!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT GET VARID CORRECTLY in read_COSMOa_NetCDF!!!!!'
       print*, varid
       stop
    endif
    
    io=nf_get_vara_real(ncid,varid,istart,icount,ncdataaux)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY in read_COSMOa_NetCDF!!!!!'
       stop
    endif
    ncdata=real(ncdataaux,8)
    	
    io=nf_close(ncid)

    return
  end subroutine read_COSMOa_NetCDF
   ! -------------------------------------------------------------------------------
  
  subroutine read_COSMOb_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=600, ncj=180  ! COSMO interpolated grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    real(kind=4), dimension (nci,ncj)   :: ncdataaux
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)
   
    ! open netcdf input file
    io=nf_open(trim(file), nf_nowrite, ncid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT OPEN FORCING FILE CORRECTLY in read_COSMOb_NetCDF!!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT GET VARID CORRECTLY in read_COSMOb_NetCDF!!!!!'
       print*, varid
       stop
    endif
   
    io=nf_get_vara_real(ncid,varid,istart,icount,ncdataaux)
    if (io.ne.nf_noerr)then
       write(*,*) file, itime, vari
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY in read_COSMOb_NetCDF!!!!!'
       stop
    endif
    ncdata=real(ncdataaux,8)
    
    io=nf_close(ncid)

    return
  end subroutine read_COSMOb_NetCDF

end module g_read_CORE_NetCDF
!
!------------------------------------------------------------------------------------
!
module g_read_NCEP_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid

contains 

  subroutine read_NCEP_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    integer, dimension(3)           	:: istart, icount
    integer(kind=2), dimension(nci,ncj) :: iuw 
    real(kind=4)                        :: xscale, xoff, miss
    real(kind=8), dimension(nci,ncj)    :: ncdata
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    ! get att
    io=nf_get_att_real(ncid,varid,'scale_factor',xscale)
    io= nf_get_att_real(ncid,varid,'add_offset',xoff)
    ! io= nf_get_att_int(ncid,varid,'missing_value',miss)

    ! get variable
    io=nf_get_vara_int2(ncid,varid,istart,icount,iuw)

    ! close file 
    io=nf_close(ncid)

    ! ncdata
    do j=1,ncj
       do i=1,nci
          ncdata(i,j)=real(iuw(i,j))*xscale+xoff
       enddo
    end do
    write(*,*) 'variabe '//vari//': ',ncdata(4:7,10)

    return
  end subroutine read_NCEP_NetCDF
  ! --------------------------------------------------------------------

  subroutine read_NCEP2_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    integer, dimension(4)           	:: istart, icount
    integer(kind=2), dimension(nci,ncj) :: iuw 
    real(kind=4)                        :: xscale, xoff, miss
    real(kind=8), dimension(nci,ncj)    :: ncdata
    character(15)                       :: vari
    character(100)                      	:: file

    istart = (/1,1,1,itime/)
    icount= (/nci,ncj,1,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    ! get att
    io=nf_get_att_real(ncid,varid,'scale_factor',xscale)
    io= nf_get_att_real(ncid,varid,'add_offset',xoff)
    ! io= nf_get_att_int(ncid,varid,'missing_value',miss)

    ! get variable
    io=nf_get_vara_int2(ncid,varid,istart,icount,iuw)

    ! close file 
    io=nf_close(ncid)

    ! ncdata
    do j=1,ncj
       do i=1,nci
          ncdata(i,j)=real(iuw(i,j))*xscale+xoff
       enddo
    end do
  !  write(*,*) 'variable '//vari//': ',ncdata(4:7,10)

    return
  end subroutine read_NCEP2_NetCDF
  !--------------------------------------------------------------
  subroutine read_CFSR_NetCDF(file,vari,itime,ncdataout)

    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter			:: nci=1152, ncj=576  ! hadcm3 grid dimension
    integer, dimension(3)           	:: istart, icount
    integer  			        :: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=4), dimension (nci,ncj)   :: ncdata
    real(kind=8), dimension (nci,ncj)   :: ncdataout
    character(25)                       :: vari
    character(80)                       :: file

    istart = (/1,1,itime/)  ! NCEP
    icount = (/nci,ncj,1/)  ! NCEP

!    istart = (/1,1,1,itime/)  ! HadCM3
!    icount = (/nci,ncj,1,1/)  ! HadCM3


    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    io=nf_get_vara_real(ncid,varid,istart,icount,ncdata)

    io=nf_close(ncid)
    
    ncdataout=dble(ncdata)
    write(*,*)file,vari,itime,ncdata(10,10),ncdataout(10,10)

    return

!    use o_param
!    use g_parfe
!    implicit none
!
!#include "netcdf.inc" 
!
!    integer, parameter             	:: nci=1152, ncj=576  ! T62 grid dimension
!    integer				:: ncid, varid, io
!    integer                             :: i, j, itime
!    integer, dimension(4)           	:: istart, icount
!    integer(kind=2), dimension(nci,ncj) :: iuw 
!    real(kind=4)                        :: xscale, xoff, miss
!    real(kind=8), dimension(nci,ncj)    :: ncdata
!    character(15)                       :: vari
!    character(80)                      	:: file
!
!    istart = (/1,1,1,itime/)
!    icount= (/nci,ncj,1,1/)
!
!    write(*,*) 'start read CFSR 1'
!    ! open netcdf input file
!    io=nf_open(file, nf_nowrite, ncid)
!
!    if (io.ne.nf_noerr)then
!       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
!       print*,'Error in opening netcdf file'//file
!       stop 'Fatal error in open_netcdf'
!    endif
!
!    write(*,*) 'start read CFSR 2'
!    ! get handle for data variable    
!    io=nf_inq_varid(ncid, vari, varid)
!
!    ! get att
!    io=nf_get_att_real(ncid,varid,'scale_factor',xscale)
!    io= nf_get_att_real(ncid,varid,'add_offset',xoff)
!    ! io= nf_get_att_int(ncid,varid,'missing_value',miss)
!
!    write(*,*) 'start read CFSR 3'
!    ! get variable
!    io=nf_get_vara_int2(ncid,varid,istart,icount,iuw)
!
!    write(*,*) 'start read CFSR 4'
!    ! close file 
!    io=nf_close(ncid)
!    
!    ! ncdata
!    do j=1,ncj
!       do i=1,nci
!          ncdata(i,j)=real(iuw(i,j))*xscale+xoff
!       enddo
!    end do
!  !  write(*,*) 'variabe '//vari//': ',ncdata(4:7,10)
!
!    return
  end subroutine read_CFSR_NetCDF


  !--------------------------------------------------------------
  !
  subroutine upside_down(array,ni,nj)
    implicit none

    integer                         :: i, j, ji
    integer                         :: ni, nj
    real(kind=8), dimension(ni,nj)  :: array
    real(kind=8), dimension(nj)     :: ywrk

    do i=1,ni  
       do j=1,nj
          ywrk(j)=array(i,j)
       end do
       do j=1,nj
          ji =nj-j+1
          array(i,ji)=ywrk(j)            
       enddo
    enddo

    return
  end subroutine upside_down

end module g_read_NCEP_NetCDF
!
!------------------------------------------------------------------------------------
!
module g_read_other_NetCDF
! Read global data with NetCDF format and interpolate to the model grid.
! Currently used for reading runoff and SSS
! The check_dummy part should be carefully modified in new applications!
! Better interpolation scheme is to be found/updated.
 
contains

  subroutine read_other_NetCDF(file, vari, itime, model_2Darray, check_dummy)
    ! if check_dummy=.true., replace missing value with a meaningful value nearby
    ! if check_dummy=.false., replace missing value with 0.0
    use o_PARAM
    use o_MESH
    use g_rotate_grid
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer			:: i, j, ii, jj, k, n, num, flag, cnt
    integer			:: itime, latlen, lonlen
    integer			:: status, ncid, varid
    integer			:: lonid, latid
    integer			:: istart(3), icount(3)
    real(kind=8)		:: x, y, miss, aux
    real(kind=8), allocatable	:: lon(:), lat(:)
    real(kind=8), allocatable	:: ncdata(:,:), ncdata_temp(:,:)
    real(kind=8), allocatable	:: temp_x(:), temp_y(:), temp_d(:)
    real(kind=8)		:: model_2Darray(myDim_nod2d+eDim_nod2D)   
    character(15)		:: vari
    character(100)              	:: file
    logical                     :: check_dummy

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! lat
    status=nf_inq_dimid(ncid, 'lat', latid)
    status=nf_inq_dimlen(ncid, latid, latlen)
    allocate(lat(latlen))
    status=nf_inq_varid(ncid, 'lat', varid)
    status=nf_get_vara_double(ncid,varid,1,latlen,lat)

    ! lon
    status=nf_inq_dimid(ncid, 'lon', lonid)
    status=nf_inq_dimlen(ncid, lonid, lonlen)
    allocate(lon(lonlen))
    status=nf_inq_varid(ncid, 'lon', varid)
    status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
    ! make sure range 0. - 360.
    do n=1,lonlen
       if(lon(n)<0.0) then
          lon(n)=lon(n)+360.
       end if
    end do

    ! data
    allocate(ncdata(lonlen,latlen), ncdata_temp(lonlen,latlen))
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1,1,itime/)
    icount= (/lonlen,latlen,1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    ! missing value
    status= nf_get_att_double(ncid,varid,'missing_value',miss)
    !write(*,*)'miss', miss
    !write(*,*)'raw',minval(ncdata),maxval(ncdata)
    ncdata_temp=ncdata
    do i=1,lonlen
       do j=1,latlen
          if(ncdata(i,j)==miss .or. ncdata(i,j)==-99.0_8) then  !!
             if(check_dummy) then
                aux=0.0
                cnt=0
                do k=1,30
                   do ii=max(1,i-k),min(lonlen,i+k)
                      do jj=max(1,j-k),min(latlen,j+k)
                         if(ncdata_temp(ii,jj)/=miss .and. ncdata_temp(ii,jj)/=-99.0_8) then  !!
                            aux=aux+ncdata_temp(ii,jj)
                            cnt=cnt+1                         
                         end if
                      end do	!ii
                   end do	!jj
                   if(cnt>0) then
                      ncdata(i,j)=aux/cnt
                      exit
                   end if
                end do  	!k    
             else
                ncdata(i,j)=0.0
             end if
          end if
       end do
    end do
    !write(*,*) 'post',minval(ncdata), maxval(ncdata)

    ! close file
    status=nf_close(ncid)
    ! model grid coordinates
    num=myDim_nod2d+eDim_nod2d
    allocate(temp_x(num), temp_y(num))  
    do n=1, num                        
       if(rotated_grid) then
          call r2g(x, y, coord_nod2d(1,n), coord_nod2d(2,n))
          temp_x(n)=x/rad   ! change unit to degree  
          temp_y(n)=y/rad                             
       else
          temp_x(n)=coord_nod2d(1,n)/rad              
          temp_y(n)=coord_nod2d(2,n)/rad             
       end if
       ! change lon range to [0 360]
       if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
    end do

    ! interpolation
    flag=0
    call interp_2d_field(lonlen, latlen, lon, lat, ncdata, num, temp_x, temp_y, & 
         model_2Darray, flag) 
    deallocate(temp_y, temp_x, ncdata_temp, ncdata, lon, lat)

  end subroutine read_other_NetCDF
end module g_read_other_NetCDF

module g_read_hadcm3_NetCDF

contains 

  subroutine read_hadcm3_vec(file,vari,itime,ncdataout)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=96, ncj=72  ! hadcm3 grid dimension
    integer, dimension(4)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=4), dimension (nci,ncj)   :: ncdata
    real(kind=8), dimension (nci,ncj)   :: ncdataout
    character(15)                       :: vari
    character(100)                      	:: file

!    istart = (/1,1,itime/)  ! NCEP
!    icount = (/nci,ncj,1/)  ! NCEP

    istart = (/1,1,1,itime/)  ! HadCM3
    icount = (/nci,ncj,1,1/)  ! HadCM3


    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    io=nf_get_vara_real(ncid,varid,istart,icount,ncdata)

    io=nf_close(ncid)
    
    ncdataout=dble(ncdata)
    write(*,*)file,vari,itime,ncdata(10,10),ncdataout(10,10)

    return
  end subroutine read_hadcm3_vec

 ! -------------------------------------------------------------------------------
 
  subroutine read_hadcm3_scal(file,vari,itime,ncdataout)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=96, ncj=73  ! hadcm3 grid dimension
    integer, parameter                  :: ncjout=72
    integer, dimension(4)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=4), dimension (nci,ncj)   :: ncdata
    real(kind=8), dimension (nci,ncjout):: ncdataout
    character(15)                       :: vari
    character(100)                      	:: file


    istart = (/1,1,1,itime/)  ! HadCM3
    icount = (/nci,ncj,1,1/)  ! HadCM3


    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in read_hadcm3_scal'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT FIND HANDLE!'
       stop 'Fatal error in read_hadcm3_scal'
    endif

    ! read data
    io=nf_get_vara_real(ncid,varid,istart,icount,ncdata)
    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       stop 'Fatal error in read_hadcm3_scal'
    endif

    io=nf_close(ncid)
    
    do i=1,nci-1
     do j=1,ncj
      ncdataout(i,j)=0.25*(ncdata(i,j)+ncdata(i+1,j)+ncdata(i,j+1)+ncdata(i+1,j+1))
     enddo
    enddo
    i=nci
    do j=1,ncj
     ncdataout(i,j)=0.25*(ncdata(i,j)+ncdata(1,j)+ncdata(i,j+1)+ncdata(1,j+1))
    enddo
    write(*,*)file,vari,itime,ncdata(10,10),ncdataout(10,10)
    

    return
  end subroutine read_hadcm3_scal

!--------------------------------------------------------------------------------------

  subroutine upside_down(array,ni,nj)
    implicit none

    integer                         :: i, j, ji
    integer                         :: ni, nj
    real(kind=8), dimension(ni,nj)  :: array
    real(kind=8), dimension(nj)     :: ywrk

    do i=1,ni  
       do j=1,nj
          ywrk(j)=array(i,j)
       end do
       do j=1,nj
          ji =nj-j+1
          array(i,ji)=ywrk(j)            
       enddo
    enddo

    return
  end subroutine upside_down

 ! -------------------------------------------------------------------------------
end module g_read_hadcm3_NetCDF

!==================================================================================

module g_read_mpeh5c_NetCDF

contains 

  subroutine read_mpeh5c_NetCDF(file,vari,itime,ncdataout)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=96, ncj=48  ! mpeh5c grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=4), dimension (nci,ncj)   :: ncdata
    real(kind=8), dimension (nci,ncj)   :: ncdataout
    character(15)                       :: vari
    character(100)                     	:: file

!    istart = (/1,1,itime/)  ! NCEP
!    icount = (/nci,ncj,1/)  ! NCEP

!    istart = (/1,1,1,itime/)  ! HadCM3
!    icount = (/nci,ncj,1,1/)  ! HadCM3

    istart = (/1,1,itime/)  ! MPEH5C
    icount = (/nci,ncj,1/)  ! MPEH5C


    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    io=nf_get_vara_real(ncid,varid,istart,icount,ncdata)

    io=nf_close(ncid)
    
    ncdataout=dble(ncdata)
    write(*,*)file,vari,itime,ncdata(10,10),ncdataout(10,10)

    return
  end subroutine read_mpeh5c_NetCDF

 ! -------------------------------------------------------------------------------
 
  subroutine upside_down(array,ni,nj)
    implicit none

    integer                         :: i, j, ji
    integer                         :: ni, nj
    real(kind=8), dimension(ni,nj)  :: array
    real(kind=8), dimension(nj)     :: ywrk

    do i=1,ni  
       do j=1,nj
          ywrk(j)=array(i,j)
       end do
       do j=1,nj
          ji =nj-j+1
          array(i,ji)=ywrk(j)            
       enddo
    enddo

    return
  end subroutine upside_down

 ! -------------------------------------------------------------------------------
end module g_read_mpeh5c_NetCDF


!==================================================================================

module g_read_eraint_NetCDF

contains 

  subroutine read_eraint_NetCDF(file,vari,itime,ncdataout)
    use o_param
    use g_parfe
    use g_config
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=480, ncj=241  ! mpeh5c grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io, iu
    integer                             :: i, j, itime
    integer (kind=2), dimension(nci,ncj)::iuw
    real(kind=4), dimension (nci,ncj)   :: ncdata
    real(kind=8), dimension (nci,ncj)   :: ncdataout
    real                                :: xscale, xoff, miss 
    character(15)                       :: vari
    character(100)                     	:: file

    istart = (/1,1,itime/)  ! ERAint
    icount = (/nci,ncj,1/)  ! ERAint


    ! open netcdf input file
    call open_netcdf_new(file,vari,ncid,iu,xscale,xoff,miss)

    io=nf_get_vara_int2(ncid,iu,istart,icount,iuw)
!    io=nf_get_vara_double(ncid, iu, istart, icount, iuw)

    if (io.ne.0) then  !modified SB
       print*,'ERROR in reading wind forcing from 1st netcdf file'
       stop
    endif
    io=nf_close(ncid)
    
   do j=1,ncj
     do i=1,nci
      ncdata(i,j)=0.
      ncdata(i,j)=real(iuw(i,j))*xscale+xoff
     enddo
    enddo

    ncdataout=dble(ncdata)
!    write(*,*)file,vari,itime,ncdata(10,10),ncdataout(10,10)

    return

  end subroutine read_eraint_NetCDF

 ! -------------------------------------------------------------------------------
 
  subroutine upside_down(array,ni,nj)
    implicit none

    integer                         :: i, j, ji
    integer                         :: ni, nj
    real(kind=8), dimension(ni,nj)  :: array
    real(kind=8), dimension(nj)     :: ywrk

    do i=1,ni  
       do j=1,nj
          ywrk(j)=array(i,j)
       end do
       do j=1,nj
          ji =nj-j+1
          array(i,ji)=ywrk(j)            
       enddo
    enddo

    return
  end subroutine upside_down

 ! -------------------------------------------------------------------------------

! ----------------------------------------------------------------------
  subroutine open_netcdf_new(file,vari,ncidu,iu,xscale,xoff,miss)
! ----------------------------------------------------------------------
! opens netcdf data file

IMPLICIT NONE
    
#include "netcdf.inc"
    character(15)                       :: vari
    character(100)                     	:: file
real (kind=4):: xs,xo
real         :: xscale, xoff, miss
integer      :: ncidu,iu
integer      :: io

!  open netcdf input data sets

io=nf_open(file,nf_nowrite,ncidu)

if (io.ne.nf_noerr)then
 print*,'!!!!! ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
 print*,'io <> 0, error in opening netcdf file'//file
 stop 'Fatal error in open_netcdf'
endif

!  get handle for ncep data

io=nf_inq_varid(ncidu,vari,iu)
io=nf_get_att_real(ncidu,iu,'scale_factor',xs)
xscale=real(xs)

io= nf_get_att_real(ncidu,iu,'add_offset',xo)
xoff=real(xo)

io= nf_get_att_int(ncidu,iu,'missing_value',miss)

!    write(*,*)'end subroutine open_netcdf ',file,xscale,xoff
return
end subroutine open_netcdf_new

 ! -------------------------------------------------------------------------------
end module g_read_eraint_NetCDF
