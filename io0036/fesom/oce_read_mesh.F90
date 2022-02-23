! Reads mesh and communication information in a distributed way.
! Some extra info. (nodal flag for cavity, region type
! of 2d elements, sigma grid slope) is also read in here.

subroutine read_mesh
  ! read mesh
  ! 
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  ! Modified by Qiang Wang, read extra information.
  !------------------------------------------------------------------

  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use o_ARRAY
  use g_config
  use g_PARfe 
  implicit none

  integer        n, m, fileID, ind, nini, nend, n1, n2, n3, n4
  integer        vert_nodes(100)
  real(kind=8)   x, y, z
  character*10   mype_string
  character*90   file_name
  character*80   dist_mesh_dir

  write(mype_string,'(i4.4)') mype  
  dist_mesh_dir=trim(meshpath)//'dist/'

  !=======================
  ! rank partitioning vectors
  !=======================
  file_name=trim(dist_mesh_dir)//'rpart.out' 
  fileID=10+mype
  open(fileID, file=trim(file_name)) 
  allocate(part2D(npes+1), part3D(npes+1))

  read(fileID,*) n
  if (n.ne.npes) then
     write(*,*) 'current NPES does not coincide with that used for mesh pre-partition'
     call par_ex
     stop
  end if
  part2D(1)=1
  read(fileID,*) part2D(2:npes+1)
  do n=2, npes+1
     part2D(n)=part2D(n-1)+part2D(n)
  end do

  part3D(1)=1
  read(fileID,*) part3D(2:npes+1)
  do n=2, npes+1
     part3D(n)=part3D(n-1)+part3D(n)
  end do
  close(fileID)
  !write(*,*) 'rpart is read'

  !===========================
  ! Lists of nodes and elements 
  ! in global indexing. Not everything
  ! is needed
  !===========================

  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
  fileID=10+mype  

  open(fileID, file=trim(file_name))
  read(fileID,*) n

  read(fileID,*) myDim_nod2D
  read(fileID,*) eDim_nod2D
  ToDim_nod2d=myDim_nod2d+eDim_nod2d
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
  read(fileID,*) myList_nod2D

  read(fileID,*) myDim_nod3D
  read(fileID,*) eDim_nod3D 	 
  ToDim_nod3d=myDim_nod3d+eDim_nod3d
  allocate(myList_nod3D(myDim_nod3D+eDim_nod3D)) 	 
  read(fileID,*) myList_nod3D

  read(fileID,*) myDim_elem2D
  allocate(myList_elem2D(myDim_elem2D))
  read(fileID,*) myList_elem2D

  read(fileID,*) myDim_elem3D
  allocate(myList_elem3D(myDim_elem3D))
  read(fileID,*) myList_elem3D ! m

  close(fileID)

  !==============================
  ! Allocate mapping array
  !==============================
  nod3D=part3D(npes+1)-1
  nod2D=part2D(npes+1)-1
  file_name=trim(meshpath)//'elem3d.out'
  open(fileID, file=file_name)
  read(fileID,*) elem3D
  close(fileID)
  allocate(mapping(elem3D))  
  mapping=0 
  !==============================
  ! It will be used for several purposes 
  ! and finally will be filled with a correct mapping 
  !==============================
  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do

  !==============================
  ! read 2d node data
  !==============================

  allocate(coord_nod2D(2,myDim_nod2D+eDim_nod2D))
  allocate(index_nod2D(myDim_nod2D+eDim_nod2D))	 			

  file_name=trim(meshpath)//'nod2d.out'
  open(fileID, file=file_name)
  read(fileID,*) n      ! nod2D, we know it already

  do n=1,nod2D
     read(fileID,*) m, x, y, ind
     if (mapping(n)>0) then
        coord_nod2D(1,mapping(n))=x
        coord_nod2D(2,mapping(n))=y
        index_nod2D(mapping(n))=ind
     end if
  end do
  mapping(1:nod2D)=0
  close(fileID)

  !==============================
  ! read 2d elem data
  !==============================
  file_name=trim(meshpath)//'elem2d.out' 
  open(fileID, file=file_name)

  allocate(elem2D_nodes(3, myDim_elem2D))
  do n=1, myDim_elem2D
     mapping(myList_elem2D(n))=n
  end do
  read(fileID,*) elem2d    
  do n=1,elem2D
     read(fileID,*) n1, n2, n3
     if (mapping(n)>0) then
        elem2D_nodes(1,mapping(n))=n1
        elem2D_nodes(2,mapping(n))=n2
        elem2D_nodes(3,mapping(n))=n3
     end if
  end do
  close(fileID)
  ! nodes in elem2d are in global numbering. convert to local:

  mapping(1:elem2D)=0
  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do
  do n=1, myDim_elem2D
     do m=1,3
        n1=elem2D_nodes(m,n)	 
        elem2D_nodes(m,n)=mapping(n1)	 
     end do
  end do
  mapping(1:nod2D)=0

  !==============================
  ! Ice shelf variables
  !==============================
  allocate(cavity_flag_nod2d(myDim_nod2d+eDim_nod2d))
  cavity_flag_nod2d=0
  ! 1 under cavity, 0 outside cavity
#ifdef use_cavity
  file_name=trim(meshpath)//'cavity_flag_nod2d.out'
  open(fileID, file=file_name)
  do n=1, myDim_nod2D+eDim_nod2d
     mapping(myList_nod2D(n))=n
  end do
  do n=1,nod2D
     read(fileID,*) n1
     if(mapping(n)>0) cavity_flag_nod2d(mapping(n))=n1
  end do
  mapping(1:nod2D)=0
  close(fileID)

  allocate(cavity_flag_extended(toDim_nod2d))
  cavity_flag_extended=0
  do n=1, myDim_elem2D
     if(any(cavity_flag_nod2d(elem2d_nodes(:,n))==1)) then
        cavity_flag_extended(elem2d_nodes(:,n))=1
     end if
  end do
#endif

  !=============================
  ! Region type of 2d elements: z-level or sigma
  !=============================
  allocate(grid_type_elem2d(myDim_elem2d))
  if(grid_type==1) then
     grid_type_elem2d=0
  elseif(grid_type==2) then
     grid_type_elem2d=1
  else
     file_name=trim(meshpath)//'grid_type_elem2d.out'
     open(fileID, file=file_name)
     do n=1, myDim_elem2D
        mapping(myList_elem2D(n))=n
     end do
     do n=1,elem2D
        read(fileID,*) n1
	if(mapping(n)>0) grid_type_elem2d(mapping(n))=n1
     end do
     mapping(1:elem2D)=0
  end if
  close(fileID)

  !==============================
  ! read 3d node data
  !==============================

  allocate(coord_nod3D(3,myDim_nod3D+eDim_nod3D))
  allocate(index_nod3D(myDim_nod3D+eDim_nod3D))	 			
  do n=1, myDim_nod3D+eDim_nod3D
     mapping(myList_nod3D(n))=n
  end do

  file_name=trim(meshpath)//'nod3d.out' 
  open(fileID, file=file_name)
  read(fileID,*) n      ! nod3D, we know it already

  do n=1,nod3D
     read(fileID,*) m, x, y, z, ind
     if (mapping(n)>0) then
        coord_nod3D(1,mapping(n))=x
        coord_nod3D(2,mapping(n))=y
        coord_nod3D(3,mapping(n))=z
        index_nod3D(mapping(n))=ind
     end if
  end do

  mapping(1:nod3D)=0
  close(fileID)

  !==============================
  ! read 3d elem data
  !==============================

  file_name=trim(meshpath)//'elem3d.out' 
  open(fileID, file=file_name)

  allocate(elem3D_nodes(4, myDim_elem3D))
  do n=1, myDim_elem3D
     mapping(myList_elem3D(n))=n
  end do
  read(fileID,*) elem3d    
  do n=1,elem3D
     read(fileID,*) n1, n2, n3, n4
     if (mapping(n)>0) then
        elem3D_nodes(1,mapping(n))=n1
        elem3D_nodes(2,mapping(n))=n2
        elem3D_nodes(3,mapping(n))=n3
        elem3D_nodes(4,mapping(n))=n4
     end if
  end do

  ! nodes in elem3d are in natural numbering. convert to local:

  mapping(1:elem3D)=0
  do n=1, myDim_nod3D+eDim_nod3D
     mapping(myList_nod3D(n))=n
  end do
  do n=1, myDim_elem3D
     do m=1,4
        n1=elem3D_nodes(m,n)	 
        elem3D_nodes(m,n)=mapping(n1)	 
     end do
  end do
  mapping(1:nod3D)=0
  close(fileID)

  !==============================
  ! read aux. arrays 
  !==============================

  file_name=trim(meshpath)//'aux3d.out'  
  open(fileID, file=file_name)

  read(fileID, *) max_num_layers

  !=============================
  ! nod3D_below_nod2D
  !============================= 
  ! ATTENTION: the array is to be stored in slices of max_num_layers
  !            or as a column   
  allocate(nod3D_below_nod2D(max_num_layers,myDim_nod2D+eDim_nod2D))       
  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do

  do n=1, nod2D
     read(fileID, *) vert_nodes(1:max_num_layers)
     if (mapping(n)>0)  then
        nod3D_below_nod2D(:,mapping(n))=vert_nodes(1:max_num_layers)
     end if
  end do
  mapping(1:nod2D)=0

  do n=1, myDim_nod3D+eDim_nod3D
     mapping(myList_nod3D(n))=n
  end do

  do n=1, myDim_nod2D+eDim_nod2D
     do m=1, max_num_layers
        n1=nod3D_below_nod2D(m,n)
        if(n1>0)  then
           nod3D_below_nod2D(m,n)=mapping(n1)
        end if
     end do
  end do

  !=============================
  ! nod2D_corresp_to_nod3D
  !============================= 

  allocate(nod2D_corresp_to_nod3D(myDim_nod3D+eDim_nod3D)) 

  do n=1, nod3D
     read(fileID, *) m
     if (mapping(n)>0)  then
        nod2D_corresp_to_nod3D(mapping(n))=m
     end if
  end do
  mapping(1:nod3D)=0

  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do

  do n=1,myDim_nod3D+eDim_nod3D
     m=nod2D_corresp_to_nod3D(n)
     nod2D_corresp_to_nod3D(n)=mapping(m)
  end do
  mapping(1:nod2D)=0

  !=============================
  ! elem2D_corresp_to_elem3D
  !=============================

  allocate(elem2D_corresp_to_elem3D(myDim_elem3D)) 
  do n=1, myDim_elem3D
     mapping(myList_elem3D(n))=n
  end do
  do n=1,elem3D
     read(fileID,*) m
     if(mapping(n)>0) then
        elem2D_corresp_to_elem3D(mapping(n))=m
     end if
  end do
  mapping(1:elem3D)=0
  do n=1, myDim_elem2D
     mapping(myList_elem2D(n))=n
  end do
  do n=1,myDim_elem3D
     m=elem2D_corresp_to_elem3D(n)
     elem2D_corresp_to_elem3D(n)=mapping(m)
  end do
  mapping(1:elem2D)=0

  close(fileID)

  !correcting 3d nodal indices
#ifndef use_opbnd_restoring
#ifndef use_opbnd_tide
  do n=1,myDim_nod3d+eDim_nod3D
     if (index_nod3D(n)==12) index_nod3D(n)=11
     if (index_nod3D(n)==22) index_nod3D(n)=21
     if (index_nod3D(n)==32) index_nod3D(n)=31
  end do
#endif
#endif

  ! ==============================
  ! read sigma slope and define layer of elem.
  ! ==============================
  if(grid_type/=1) call read_grid_slope	


  if(mype==0) then
     write(*,*) 'mesh (according to pre-partition) is read in'	
     write(*,*) 'configured with nod2D=',nod2D,' nod3D=',nod3D
  end if

  ! ==============================
  ! Communication information
  ! ==============================
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out'  
  fileID=10+mype  
  open(fileID, file=file_name)
  read(fileID,*)  n
  read(fileID,*) com_nod2D%rPEnum
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  read(fileID,*) com_nod2D%rPE
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1))
  read(fileID,*) com_nod2D%rptr
  allocate(com_nod2D%rlist(eDim_nod2D))
  read(fileID,*) com_nod2D%rlist

  read(fileID,*) com_nod2D%sPEnum
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  read(fileID,*) com_nod2D%sPE
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))
  read(fileID,*) com_nod2D%sptr
  n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
  allocate(com_nod2D%slist(n))
  read(fileID,*) com_nod2D%slist

  read(fileID,*) com_nod3D%rPEnum
  allocate(com_nod3D%rPE(com_nod3D%rPEnum))
  read(fileID,*) com_nod3D%rPE
  allocate(com_nod3D%rptr(com_nod3D%rPEnum+1))
  read(fileID,*) com_nod3D%rptr
  allocate(com_nod3D%rlist(eDim_nod3D))
  read(fileID,*) com_nod3D%rlist

  read(fileID,*) com_nod3D%sPEnum
  allocate(com_nod3D%sPE(com_nod3D%sPEnum))
  read(fileID,*) com_nod3D%sPE
  allocate(com_nod3D%sptr(com_nod3D%sPEnum+1))
  read(fileID,*) com_nod3D%sptr
  n=com_nod3D%sptr(com_nod3D%sPEnum+1)-1
  allocate(com_nod3D%slist(n))
  read(fileID,*) com_nod3D%slist
  close(fileID)

  ! mapping
  file_name=trim(dist_mesh_dir)//'mapping'//'.out'  
  fileID=10+mype
  open(fileID, file=trim(file_name))
  read(fileID, *) mapping(1:nod3D+nod2D)
  close(fileID)

  allocate(col_pos(myDim_nod3D+eDim_nod3D))
  col_pos=0

  if(mype==0) write(*,*) 'comm and mapping is read'
end subroutine  read_mesh
!
!=========================================================================
!
subroutine read_grid_slope  
  ! read sigma grid slope
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-----------------------------------------------------------
  
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_config
  use g_parfe

  implicit none

  integer                   :: i, j
  real(kind=8)              :: temp(2)

  ! sigma grid slope
  allocate(grid_slope(2,max_num_layers-1,myDim_elem2d))
  grid_slope=0.0
  do i=1,myDim_elem2D
     mapping(myList_elem2D(i))=i
  end do
  open(13,file=trim(MeshPath)//'sigma_grid_slope_elem.out', status='old')
  do i=1,elem2d
     if(mapping(i)>0) then
        do j=1,max_num_layers-1
           read(13,*) temp(:)
           grid_slope(:,j,mapping(i))=temp(:)
        end do
     else
        do j=1,max_num_layers-1
           read(13,*) temp(:)
        end do
     end if
  end do
  close(13)
  mapping(1:elem2D)=0

end subroutine read_grid_slope
!
!=========================================================================
