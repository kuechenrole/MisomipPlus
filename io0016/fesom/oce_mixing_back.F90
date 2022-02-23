subroutine oce_mixing_back
  ! set up background mixing for Av and Kv
  ! only needed for the "no" mixing option; all others (e.g. kpp) do this
  ! already
  !
  ! Coded by ole Richter
  ! Reviewed by ??
  !---------------------------------------------------
  
  use o_MESH
  use o_ELEMENTS
  use g_config
  use o_PARAM
  use o_array
  use g_forcing_arrays
  use g_PARFE


  implicit none

  integer :: i

  do i=1,ToDim_nod3d
    Kv(i,1)=Kv0
    Av(i)=Av0
  end do

end subroutine oce_mixing_back

