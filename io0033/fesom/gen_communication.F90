! 1) Communicate neighbor nodes between partitions
! 2) Exchange full arrays
!-------------------------------------------------------------------------

subroutine com_2D(arr2d)
  ! Coded by Sergey Danilov
  ! Reviewed by ??

  use o_MESH
  use i_ARRAY
  use g_PARFE 
  implicit none

  integer       :: sreq(npes)
  integer       :: rreq(npes)
  integer       :: sstat(MPI_STATUS_SIZE,npes)
  integer       :: rstat(MPI_STATUS_SIZE,npes)
  integer       :: n, sn, rn, dest, nini, nend
  integer       :: offset, count, source
  real(kind=8)	:: arr2d(nod2d)

  ! Put data to be communicated into send buffer 
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum

  do n=1, sn
     nini=com_nod2D%sptr(n)
     nend=com_nod2D%sptr(n+1) - 1
     s_buff_2d(n)%array=arr2d(com_nod2D%slist(nini:nend))
  end do

  do n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     count=com_nod2D%sptr(n+1) - nini

     call MPI_ISEND(s_buff_2d(n)%array, count, MPI_DOUBLE_PRECISION, dest, mype, & 
          MPI_COMM_WORLD, sreq(n), MPIerr)

     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     count=com_nod2D%rptr(n+1) - nini

     call MPI_IRECV(r_buff_2d(n)%array, count, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_WORLD, rreq(n), MPIerr) 
  end do

  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
  call MPI_WAITALL(sn,rreq,rstat, MPIerr)

  ! Put received data to their destination
  do n=1, rn
     nini=com_nod2D%rptr(n)
     nend=com_nod2D%rptr(n+1) - 1
     count=com_nod2D%rptr(n+1) - nini
     arr2d(com_nod2D%rlist(nini:nend))=r_buff_2d(n)%array
  end do
end subroutine com_2D
!
!===================================================================
!
subroutine com_3D(arr3d)
  ! Coded by Sergey Danilov
  ! Reviewed by ??

  use o_MESH
  use o_ARRAY
  use g_PARFE 
  implicit none
  !
  integer  	:: sreq(npes)
  integer  	:: rreq(npes)
  integer  	:: sstat(MPI_STATUS_SIZE,npes)
  integer  	:: rstat(MPI_STATUS_SIZE,npes)
  integer  	:: n, sn, rn, dest, nini, nend
  integer	:: offset, count, source
  real(kind=8)	:: arr3d(nod3d)


  ! Put data to be communicated into send buffer 
  sn=com_nod3D%sPEnum
  rn=com_nod3D%rPEnum

  do n=1, sn
     nini=com_nod3D%sptr(n)
     nend=com_nod3D%sptr(n+1) - 1
     count=com_nod3D%sptr(n+1) - nini
     s_buff_3d(n)%array=arr3d(com_nod3D%slist(nini:nend))
  end do

  do n=1, sn
     dest=com_nod3D%sPE(n)
     nini=com_nod3D%sptr(n)
     count=com_nod3D%sptr(n+1) - nini

     call MPI_ISEND(s_buff_3d(n)%array, count, MPI_DOUBLE_PRECISION, dest, mype, & 
          MPI_COMM_WORLD, sreq(n), MPIerr)

     source=com_nod3D%rPE(n)
     nini=com_nod3D%rptr(n)
     count=com_nod3D%rptr(n+1) - nini

     call MPI_IRECV(r_buff_3d(n)%array, count, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_WORLD, rreq(n), MPIerr) 
  end do

  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
  call MPI_WAITALL(sn,rreq,rstat, MPIerr)

  ! Put received data to their destination
  do n=1, rn
     nini=com_nod3D%rptr(n)
     nend=com_nod3D%rptr(n+1) - 1
     count=com_nod3D%rptr(n+1) - nini
     arr3d(com_nod3D%rlist(nini:nend))=r_buff_3d(n)%array
  end do
end subroutine com_3D
!
!===================================================================
!
subroutine broadcast3D(arr3D, arr3Dglobal)
  ! Makes nodal information available to all PE
  ! arr3d is any array like TF or SF of local size.
  ! arr3Dglobal is an array of nod3D size which 
  ! should be allocated before calling this routine.
  ! It will be filled with information on other PE in 
  ! natural numbering. The routine can be used to organize
  ! output in the same way as in global memory setup  
  !
  ! Coded by Sergey Danilov
  ! Reviewed by ??
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer :: ireals
  integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
  integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

  real(kind=8) ::  arr3D(myDim_nod3D+eDim_nod3D)
  real(kind=8) ::  arr3Dglobal(nod3D)
  real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)

  if ( mype == 0 ) then
     if (npes>1) then
        arr3Dglobal(myList_nod3D(1:myDim_nod3D))=arr3D(1:myDim_nod3D)
     end if
     do  n = 1, npes-1

        call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             0, MPI_COMM_WORLD, status, MPIerr )
        sender = status(MPI_SOURCE)
        allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
        call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
             1, MPI_COMM_WORLD, status, MPIerr )
        call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
             2, MPI_COMM_WORLD, status, MPIerr )

        do i = 1, nTS
           arr3Dglobal(irecvbuf(i)) = recvbuf(i)
        enddo
        deallocate( recvbuf, irecvbuf )

     enddo

  else

     allocate( sendbuf(1:myDim_nod3D), isendbuf(1:myDim_nod3D) )
     do n = 1, myDim_nod3D
        isendbuf(n) = myList_nod3D(n)
        sendbuf(n)  = arr3D(n)
     enddo
     call MPI_SEND( myDim_nod3D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( isendbuf(1), myDim_nod3D, MPI_INTEGER, 0, 1, &
          MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( sendbuf(1), myDim_nod3D, MPI_DOUBLE_PRECISION, &
          0, 2, MPI_COMM_WORLD, MPIerr )
     deallocate( sendbuf, isendbuf )

  endif

  call MPI_BCAST( arr3Dglobal, nod3d, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, MPIerr)

end subroutine broadcast3D
!
!===================================================================
!
subroutine broadcast2D(arr2D, arr2Dglobal)
  ! Makes nodal information available to all PE 
  ! As the preceeding routine, but for 2D arrays
  !
  ! Coded by Sergey Danilov
  ! Reviewed by ??
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer :: ireals
  integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
  integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

  real(kind=8) ::  arr2D(myDim_nod2D+eDim_nod2D)
  real(kind=8) ::  arr2Dglobal(nod2D)
  real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
  if ( mype == 0 ) then
     if (npes>1) then
        arr2Dglobal(myList_nod2D(1:myDim_nod2D))=arr2D(1:myDim_nod2D)
     end if
     do  n = 1, npes-1

        call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             0, MPI_COMM_WORLD, status, MPIerr )
        sender = status(MPI_SOURCE)
        allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
        call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
             1, MPI_COMM_WORLD, status, MPIerr )
        call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
             2, MPI_COMM_WORLD, status, MPIerr )

        do i = 1, nTS
           arr2Dglobal(irecvbuf(i)) = recvbuf(i)
        enddo
        deallocate( recvbuf, irecvbuf )

     enddo

  else

     allocate( sendbuf(1:myDim_nod2D), isendbuf(1:myDim_nod2D) )
     do n = 1, myDim_nod2D
        isendbuf(n) = myList_nod2D(n)
        sendbuf(n)  = arr2D(n)
     enddo
     call MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
          MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
          0, 2, MPI_COMM_WORLD, MPIerr )
     deallocate( sendbuf, isendbuf )

  endif

  call MPI_BCAST( arr2Dglobal, nod2d, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, MPIerr)

end subroutine broadcast2D
!===================================================================
