!
!--------------------------------------------------------------------
subroutine init_pot ( zeta, r, mesh, vpot )
!--------------------------------------------------------------------
  !
  ! initialize potential
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: mesh
  real (dp), intent(in) :: zeta, r(0:mesh)
  real (dp), intent(out):: vpot(0:mesh)
  integer :: i

  open (7,file='pot.out',status='unknown',form='formatted')
  write(7,'("#       r             V(r)")')
  do i =0,mesh
     vpot (i) = - 2.0_dp*zeta/r(i)
     write (7,*) r(i),vpot(i)
  end do
  close(7)
  
  return
end subroutine init_pot
