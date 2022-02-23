module o_passive_tracer_mod
  ! Ocean passive tracer module
  !
  ! Coded by Qiang Wang
  ! Adjusted for ice shelf melt water tracer purposes by Ralph Timmermann 08.04.2014
  !-------------------------------------------------------

  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  integer, allocatable, dimension(:)        :: index_passive_tracer
  integer, allocatable, dimension(:,:)      :: passive_tracer_loc_index
  real(kind=8), allocatable, dimension(:,:) :: ptr_sfc_force

contains


  subroutine passive_tracer_init

! Special version originally created for Clara Stolle's master thesis and used to detect water mass pathways in FRIS cavity.
! The original FESOM subroutine passive_tracer_init has been renamed to subroutine passive_tracer_init_original (see below).
! Ralph Timmermann, 24.10.2018

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cptrind
    character(4)         :: tr_name
    character(100)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_passive_tracer(num_passive_tracer))
    do j=1, num_passive_tracer
       write(cptrind,'(i1)') j
       tr_name='ptr'//cptrind
       do i=1, num_tracer
          if(prog_tracer_name(i) == tr_name) index_passive_tracer(j)=i
       end do
    end do

    !--------------------------------------------------------------  
    ! initial values
    do j=1, num_passive_tracer
       tracer(:,index_passive_tracer(j))=ptr_background_value
    end do

    !--------------------------------------------------------------
    ! in case that p.tr is restored in a region
    if(passive_tracer_restore) then

       ! set passive tracer location index: 1 at release, 0 otherwise

       allocate(passive_tracer_loc_index(ToDim_nod3d,num_passive_tracer))
       passive_tracer_loc_index=0


!rt       allocate(temp_arr2d(nod2d))
!rt       temp_arr2d=0
!rt       do n=1, ToDim_nod2D
!rt          temp_arr2d(myList_nod2D(n))=n
!rt       end do

!rt       do j=1, num_passive_tracer
!rt        write(cptrind,'(i1)') j
!rt        tr_name='ptr'//cptrind
!rt        file_name=trim(meshpath)//'passive_tracer_restore_nodes_'//tr_name//'.out'
!rt        fileID=160
!rt        write(*,*)'passive_tracer_init: now open ',file_name
!rt        open(fileID, file=file_name)
!rt        read(fileID,*) num_nod
!rt        allocate(nodes_release(num_nod))
!rt        read(fileID,*) nodes_release
!rt        close(fileID)
!rt          do n=1,num_nod
!rt             n_loc=temp_arr2d(nodes_release(n))
!rt             if(n_loc>0) then
!rt                n_loc=nod3d_below_nod2d(1,n_loc)
!rt                passive_tracer_loc_index(n_loc,j)=1
!rt                tracer(n_loc,index_passive_tracer(j))=ptr_restore_value
!rt             end if
!rt          end do
!rt          deallocate(nodes_release)
!rt       end do

!rt       deallocate(temp_arr2d)

       !--------------------------------------------------------------
       ! in case restore volume
       if(ptr_restore_in_volume) then
        do i=1,ToDim_nod2d
         row=nod3d_below_nod2d(1,i)
         do j=1, num_passive_tracer
          if (j.eq.1) then                ! Filchner Trough RT 23.01.2019
           if (geolon2d(i).gt.-45.0*rad .and. geolon2d(i).lt.-35.0*rad &
         .and. geolat2d(i).gt.-78.0*rad .and. geolat2d(i).lt.-77.5*rad) then
            do k=2,num_layers_below_nod2d(i)+1
             n3=nod3d_below_nod2d(k,i)
             if (abs(coord_nod3d(3,n3)).gt.100.)then
!              write(*,*)'ptr_restore',mype, j, k, n3, abs(coord_nod3d(3,n3)),index_passive_tracer(j)
              passive_tracer_loc_index(n3,j)=1
              tracer(n3,index_passive_tracer(j))=ptr_restore_value
             endif
            enddo
           endif
          else if (j.eq.2) then     ! Ronne Trough  RT 23.01.2019
           if (geolon2d(i).gt.-62.0*rad .and. geolon2d(i).lt.-58.0*rad &
         .and. geolat2d(i).gt.-75.5*rad .and. geolat2d(i).lt.-74.5*rad) then
            do k=2,num_layers_below_nod2d(i)+1
             n3=nod3d_below_nod2d(k,i)
             if (abs(coord_nod3d(3,n3)).gt.400.)then
!              write(*,*)'ptr_restore',mype, j, k, n3, abs(coord_nod3d(3,n3)),index_passive_tracer(j)
              passive_tracer_loc_index(n3,j)=1
              tracer(n3,index_passive_tracer(j))=ptr_restore_value
             endif
            enddo
           endif
          else if (j.eq.3) then     ! Foundation Ice Stream Grounding Line  RT 23.01.2019
           if (geolon2d(i).gt.-63.0*rad .and. geolon2d(i).lt.-57.0*rad &
         .and. geolat2d(i).gt.-84.0*rad .and. geolat2d(i).lt.-82.0*rad) then
!            do k=2,num_layers_below_nod2d(i)+1
!             n3=nod3d_below_nod2d(k,i)
!             if (abs(coord_nod3d(3,n3)).gt.1100.)then
!              write(*,*)'ptr_restore',mype, j, k, n3, abs(coord_nod3d(3,n3)),index_passive_tracer(j)
!              passive_tracer_loc_index(n3,j)=1
!              tracer(n3,index_passive_tracer(j))=ptr_restore_value
!             endif
!            enddo
            k=1                        ! modified RT RG46451 ab Januar 1985
            n3=nod3d_below_nod2d(k,i)
!            write(*,*)'ptr_restore',mype, j, k, n3, abs(coord_nod3d(3,n3)),index_passive_tracer(j)
            passive_tracer_loc_index(n3,j)=1
            tracer(n3,index_passive_tracer(j))=ptr_restore_value
           endif
          else if (j.eq.4) then     ! Support Force Glacier Grounding Line  RT 23.01.2019
           if (geolon2d(i).gt.-48.0*rad .and. geolon2d(i).lt.-42.0*rad &
         .and. geolat2d(i).gt.-83.0*rad .and. geolat2d(i).lt.-81.9*rad) then
!            do k=2,num_layers_below_nod2d(i)+1
!             n3=nod3d_below_nod2d(k,i)
!             if (abs(coord_nod3d(3,n3)).gt.900.)then
!              write(*,*)'ptr_restore',mype, j, k, n3, abs(coord_nod3d(3,n3)),index_passive_tracer(j)
!              passive_tracer_loc_index(n3,j)=1
!              tracer(n3,index_passive_tracer(j))=ptr_restore_value
!             endif
!            enddo
            k=1                        ! modified RT RG46451 ab Januar 1985
            n3=nod3d_below_nod2d(k,i)
!            write(*,*)'ptr_restore',mype, j, k, n3, abs(coord_nod3d(3,n3)),index_passive_tracer(j)
            passive_tracer_loc_index(n3,j)=1
            tracer(n3,index_passive_tracer(j))=ptr_restore_value
           endif
          endif
         enddo
        enddo
       endif
    endif
    !--------------------------------------------------------------
    ! in case that passive tracers enter through surface fluxes
    if(passive_tracer_flux) then
       allocate(ptr_sfc_force(ToDim_nod2d,num_passive_tracer))
    end if

    !backup 
    tracer0(:,index_passive_tracer)=tracer(:,index_passive_tracer)

  end subroutine passive_tracer_init
  !


  subroutine passive_tracer_init_original

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cptrind
    character(4)         :: tr_name
    character(100)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_passive_tracer(num_passive_tracer))
    do j=1, num_passive_tracer
       write(cptrind,'(i1)') j
       tr_name='ptr'//cptrind
       do i=1, num_tracer
          if(prog_tracer_name(i) == tr_name) index_passive_tracer(j)=i
       end do
    end do

    !--------------------------------------------------------------  
    ! initial values
    do j=1, num_age_tracer   ! Ralph thinks that this should be num_passive_tracer
       tracer(:,index_passive_tracer(j))=ptr_background_value
    end do

    !--------------------------------------------------------------
    ! in case that p.tr is restored in a region
    if(passive_tracer_restore) then

       ! set passive tracer location index: 1 at release, 0 otherwise

       allocate(passive_tracer_loc_index(ToDim_nod3d,num_passive_tracer))
       passive_tracer_loc_index=0

       allocate(temp_arr2d(nod2d))
       temp_arr2d=0
       do n=1, ToDim_nod2D
          temp_arr2d(myList_nod2D(n))=n
       end do

       do j=1, num_passive_tracer
          write(cptrind,'(i1)') j
          tr_name='ptr'//cptrind
          file_name=trim(meshpath)//'passive_tracer_restore_nodes_'//tr_name//'.out'
          fileID=160
          write(*,*)'passive_tracer_init: now open ',file_name
          open(fileID, file=file_name)
          read(fileID,*) num_nod
          allocate(nodes_release(num_nod))
          read(fileID,*) nodes_release
          close(fileID)
          do n=1,num_nod
             n_loc=temp_arr2d(nodes_release(n))
             if(n_loc>0) then
                n_loc=nod3d_below_nod2d(1,n_loc)
                passive_tracer_loc_index(n_loc,j)=1
                tracer(n_loc,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
          deallocate(nodes_release)
       end do

       deallocate(temp_arr2d)

       !--------------------------------------------------------------
       ! in case restore volume
       if(ptr_restore_in_volume) then
          do i=1,ToDim_nod2d
             row=nod3d_below_nod2d(1,i)
             do j=1, num_passive_tracer
                if(passive_tracer_loc_index(row,j)==1) then
                   do k=2,num_layers_below_nod2d(i)+1
                      n3=nod3d_below_nod2d(k,i)
                      passive_tracer_loc_index(n3,j)=1
                      tracer(n3,index_passive_tracer(j))=ptr_restore_value
                   end do
                end if
             end do
          end do
       end if

    end if

    !--------------------------------------------------------------
    ! in case that passive tracers enter through surface fluxes
    if(passive_tracer_flux) then
       allocate(ptr_sfc_force(ToDim_nod2d,num_passive_tracer))
    end if

    !backup 
    tracer0(:,index_passive_tracer)=tracer(:,index_passive_tracer)

  end subroutine passive_tracer_init_original
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_sfc_bc

    integer      :: elem, j, elnodes(3), elnodes2(3)
    real(kind=8) :: auxf, entries(3)

    if(passive_tracer_flux) then

       ptr_sfc_force=0.0

!rt       do elem=1,myDim_elem2d             
!rt          elnodes2=elem2D_nodes(:,elem)
!rt          elnodes=nod3D_below_nod2D(1,elnodes2)   
!rt          auxf=voltriangle(elem)/12.0_8    
!rt          do j=1,num_passive_tracer
!rt             entries=-auxf*(tracer(elnodes,2)+tracer(elnodes,index_passive_tracer(j))) &
!rt                  * runoff_landice(elnodes2)*landice_season(month)               !!!!! HIER  !!!!!
!rt             ptr_sfc_force(elnodes2,j)=ptr_sfc_force(elnodes2,j)+sum(entries)+entries
!rt          end do
!rt       end do

    end if

  end subroutine ptr_sfc_bc
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_cutoff_restore
   ! cleaned up and corrected by Ralph Timmermann for RG46421t1

    integer   :: j, row

    if(passive_tracer_restore) then
       do j=1, num_passive_tracer
          do row=1,ToDim_nod3d
             if(passive_tracer_loc_index(row,j)==1) then
                tracer(row,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
       end do
    endif

  end subroutine ptr_cutoff_restore


end module o_passive_tracer_mod
