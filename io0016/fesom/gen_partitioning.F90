! 1)MPI initialization and finalization;
! 2)Partitioning
!
! Coded by Sergey Danilov
! Reviewed by ??
!---------------------------------------------------------------------------

subroutine par_init      ! initializes MPI
  use g_PARFE
  implicit none

  integer :: i

  call MPI_INIT(i)
  call MPI_Comm_Size(MPI_COMM_WORLD,npes,i)
  call MPI_Comm_Rank(MPI_COMM_WORLD,mype,i)

#ifdef PETSC
  call PETSCInitialize(PETSC_NULL_CHARACTER,i)
#endif

  if(mype==0) write(*,*) 'MPI has been initialized'
end subroutine par_init
!
!============================================================================
!
subroutine par_ex       ! finalizes MPI
  use g_PARFE

  call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call  MPI_Finalize(MPIerr)

end subroutine par_ex
!
!============================================================================
!
subroutine set_par_support_ini
  ! used during preparing distributed memory
  use o_MESH
  use o_elements
  use o_matrices
  use g_PARFE 
  use o_MATRICES
  implicit none
  !
  integer   n, j, k, nini, nend, count, dim_array

  ! Define partitioning vector

  allocate(part2D(nod2D))
  allocate(part3D(nod3D))

  do n=0, npes-1
     nini=(nod2D/npes)*n+1
     nend=(nod2D/npes)*(n+1)
     if (n==npes-1) nend=nod2D
     part2D(nini:nend)=n
  end do
 ! call partit(sshstiff%dim,sshstiff%rowptr,sshstiff%colind,npes, part2D)
 !call partit2(sshstiff%dim,sshstiff%rowptr,sshstiff%colind,num_layers_below_nod2D,npes,part2D)
 call partit3(sshstiff%dim,sshstiff%rowptr,sshstiff%colind,num_layers_below_nod2D,npes,part2D)

  ! propagate partitioning vector down:
  do n=1,nod2D
     do j=1,num_layers_below_nod2D(n)+1
        k = nod3D_below_nod2D(j,n)
        part3D(k)=part2D(n)
     end do
  end do

  call communication_nod
  call mymesh


  if(mype==0) write(*,*) 'Communication arrays are set up'   
end subroutine set_par_support_ini

!=======================================================================

subroutine set_par_support
  ! use during model run
  use o_MESH
  use g_PARFE 
  use o_MATRICES
  implicit none

  integer   n, cnt
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here only communication buffers are set. 

  ! Allocate communication buffers: 

  if (npes>1) then
     allocate(s_buff_2d(com_nod2D%sPEnum),r_buff_2d(com_nod2D%sPEnum))
     do n=1, com_nod2D%sPEnum
        cnt=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)
        allocate(s_buff_2d(n)%array(cnt))
     end do
     do n=1, com_nod2D%rPEnum
        cnt=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)
        allocate(r_buff_2d(n)%array(cnt))
     end do

     allocate(s_buff_3d(com_nod3D%sPEnum),r_buff_3d(com_nod3D%sPEnum))
     do n=1, com_nod3D%sPEnum
        cnt=com_nod3D%sptr(n+1) - com_nod3D%sptr(n)
        allocate(s_buff_3d(n)%array(cnt))
     end do
     do n=1, com_nod3D%rPEnum
        cnt=com_nod3D%rptr(n+1) - com_nod3D%rptr(n)
        allocate(r_buff_3d(n)%array(cnt))
     end do
  end if

  if(mype==0) write(*,*) 'Communication buffer arrays are set up' 
end subroutine set_par_support
!
!============================================================================
!
subroutine communication_nod
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none

  integer n,np, nz, prank, elem, elnodes(3), epe(3), counter, nini, nend
  integer, allocatable :: aux(:,:), pnum(:,:),pmap(:)

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules
  allocate(pnum(npes,npes))

  do elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     epe=part2D(elnodes)+1
     if(epe(1).ne.epe(2)) then
        pnum(epe(1), epe(2))=1
        pnum(epe(2), epe(1))=1
     end if

     if(epe(2).ne.epe(3)) then
        pnum(epe(3), epe(2))=1
        pnum(epe(2), epe(3))=1
     end if
     if(epe(1).ne.epe(3)) then
        pnum(epe(1), epe(3))=1
        pnum(epe(3), epe(1))=1
     end if
  end do
  ! mype interacts with a limited number of PEs
  allocate(pmap(npes))
  pmap=0
  counter=0
  do n=1,npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        pmap(n)=counter
     end if
     if(n==mype+1) then
        counter=counter+1
        pmap(n)=counter
     end if
  end do

  allocate(aux(counter,nod2D)) 
  ! This has a much smaller size if npes is large
  ! but we need pmap to address it
  aux=0
  pnum=0

  do elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     epe=part2D(elnodes)+1
     if(epe(1).ne.epe(2)) then
        if(pmap(epe(1)).ne.0) aux(pmap(epe(1)), elnodes(2))=1
        if(pmap(epe(2)).ne.0) aux(pmap(epe(2)), elnodes(1))=1
     end if

     if(epe(2).ne.epe(3)) then
        if(pmap(epe(3)).ne.0) aux(pmap(epe(3)), elnodes(2))=1
        if(pmap(epe(2)).ne.0) aux(pmap(epe(2)), elnodes(3))=1
     end if
     if(epe(1).ne.epe(3)) then
        if(pmap(epe(1)).ne.0) aux(pmap(epe(1)), elnodes(3))=1
        if(pmap(epe(3)).ne.0) aux(pmap(epe(3)), elnodes(1))=1
     end if
  end do

  do n=1, nod2D
     do np=1, npes
        if(pmap(np).ne.0) then
           if(aux(pmap(np),n).ne.0) then 
              pnum(np,part2D(n)+1)=pnum(np,part2D(n)+1)+1
           end if
	end if
     end do
  end do

  ! We know how many external nodes each PE needs
  ! This is the 'receive' list   
  ! com_nod2D for 2D nodes

  ! The number of external PE I receive information from
  com_nod2D%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_nod2D%rPEnum=com_nod2D%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external nodes ordered by external PE ranks

  counter=0
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
  com_nod2D%rptr(1)=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rptr(counter+1)=com_nod2D%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%rlist(com_nod2D%rptr(com_nod2D%rPEnum+1)-1)) 
  do np=1,com_nod2D%rPEnum
     prank=com_nod2D%rPE(np)
     do n=1, nod2D
        if((aux(pmap(mype+1),n)==1).and.(part2D(n)==prank)) then
           counter=counter+1
           com_nod2D%rlist(counter)=n
        end if
     end do
  end do
  ! Summary of this piece: mype receives
  ! information on external 2D nodes from
  ! comm_nod2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%rptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part, and we need analogous
  ! information for 3D nodes and edges.

  ! SENDING PART
  com_nod2D%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_nod2D%sPEnum=com_nod2D%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1)) 
  com_nod2D%sptr(1)=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sptr(counter+1)=com_nod2D%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
  do np=1,com_nod2D%sPEnum
     prank=com_nod2D%sPE(np)
     do n=1, nod2D
        if(pmap(prank+1).ne.0) then 
           if((aux(pmap(prank+1),n)==1).and.(part2D(n)==mype)) then
              counter=counter+1
              com_nod2D%slist(counter)=n
           end if
	end if
     end do
  end do

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)

  deallocate(pmap, pnum, aux)

  ! 3D part looks similar and simply requres to use larger arrays
  ! to store indices:
  com_nod3D%rPEnum=com_nod2D%rPEnum
  com_nod3D%sPEnum=com_nod2D%sPEnum
  allocate(com_nod3D%sPE(com_nod2D%sPEnum))
  allocate(com_nod3D%rPE(com_nod2D%rPEnum))
  com_nod3D%rPE=com_nod2D%rPE
  com_nod3D%sPE=com_nod2D%sPE
  allocate(com_nod3D%rptr(com_nod2D%rPEnum+1))
  allocate(com_nod3D%sptr(com_nod2D%sPEnum+1))
  com_nod3D%rptr(1)=1
  com_nod3D%sptr(1)=1

  ! Compute pointers and lists of communicating 3D nodes

  do np=1, com_nod3D%rPEnum
     counter=0
     nini=com_nod2D%rptr(np)
     nend=com_nod2D%rptr(np+1)-1
     do n=nini,nend
        counter=counter+num_layers_below_nod2D(com_nod2D%rlist(n))+1
     end do
     com_nod3D%rptr(np+1)=com_nod3D%rptr(np)+counter
  end do
  allocate(com_nod3D%rlist(com_nod3D%rptr(com_nod3D%rPEnum+1)-1)) 
  counter=0
  do np=1, com_nod3D%rPEnum
     nini=com_nod2D%rptr(np)
     nend=com_nod2D%rptr(np+1)-1
     do n=nini,nend
        do nz=1,num_layers_below_nod2D(com_nod2D%rlist(n))+1
           counter=counter+1
           com_nod3D%rlist(counter)=nod3D_below_nod2D(nz,com_nod2D%rlist(n))
        end do
     end do
  end do

  do np=1, com_nod3D%sPEnum
     counter=0
     nini=com_nod2D%sptr(np)
     nend=com_nod2D%sptr(np+1)-1
     do n=nini,nend
        counter=counter+num_layers_below_nod2D(com_nod2D%slist(n))+1
     end do
     com_nod3D%sptr(np+1)=com_nod3D%sptr(np)+counter
  end do
  allocate(com_nod3D%slist(com_nod3D%sptr(com_nod3D%sPEnum+1)-1)) 
  counter=0
  do np=1, com_nod3D%sPEnum
     nini=com_nod2D%sptr(np)
     nend=com_nod2D%sptr(np+1)-1
     do n=nini,nend
        do nz=1,num_layers_below_nod2D(com_nod2D%slist(n))+1
           counter=counter+1
           com_nod3D%slist(counter)=nod3D_below_nod2D(nz,com_nod2D%slist(n))
        end do
     end do
  end do
  ! For the case of distributed memory, numbering of nodes should be 
  ! made consistent before this procedure. Then transition between
  ! global--local is elementary

end subroutine communication_nod
!
!==========================================================================
!
subroutine mymesh
  use o_MESH
  use o_ELEMENTS
  use g_PARFE 
  implicit none
  !
  integer     n, counter, q

  !======= NODES 

  ! 2D nodes
  ! Owned nodes + external nodes which I need:

  counter=0
  do n=1, nod2D
     if (part2D(n)==mype) counter=counter+1
  end do
  myDim_nod2D=counter
  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1   
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
  counter=0   
  do n=1, nod2D
     if (part2D(n)==mype) then
        counter=counter+1
        myList_nod2D(counter)=n
     end if
  end do
  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)=&
       com_nod2D%rlist

  ! Summary:  	     
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs;    
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes

  ! 3D nodes 
  counter=0

  do n=1, nod3D
     if (part3D(n)==mype) counter=counter+1
  end do
  myDim_nod3D=counter
  eDim_nod3D=com_nod3D%rptr(com_nod3D%rPEnum+1)-1   
  allocate(myList_nod3D(myDim_nod3D+eDim_nod3D))
  counter=0   
  do n=1, nod3D
     if (part3D(n)==mype) then
        counter=counter+1
        myList_nod3D(counter)=n
     end if
  end do
  myList_nod3D(myDim_nod3D+1:myDim_nod3D+eDim_nod3D)=&
       com_nod3D%rlist

  ! Summary:  	     
  ! myList_nod3D(myDim_nod3D+1:myDim_nod3D+eDim_nod3D)
  ! contains external nodes which mype needs;    
  ! myList_nod3D(1:myDim_nod3D) contains owned nodes

  !======= ELEMENTS

  ! 2D elements 
  counter=0
  do n=1, elem2D
     do q=1,3
        if(part2D(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_elem2D=counter
  allocate(myList_elem2D(myDim_elem2D))

  counter=0
  do n=1, elem2D
     do q=1,3
        if(part2D(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           myList_elem2D(counter)=n
           exit
        end if
     end do
  end do

  ! 3D elements
  counter=0
  do n=1, elem3D
     do q=1,4
        if(part3D(elem3D_nodes(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_elem3D=counter
  allocate(myList_elem3D(myDim_elem3D))

  counter=0
  do n=1, elem3D
     do q=1,4
        if(part3D(elem3D_nodes(q,n))==mype) then
           counter=counter+1
           myList_elem3D(counter)=n
           exit
        end if
     end do
  end do

  ! Summary: element lists are only needed for parallel assembling.
  ! For this reason, we do not need to distinguish between owned
  ! and shared (but visited) elements --- they should be in the list.
end subroutine mymesh
