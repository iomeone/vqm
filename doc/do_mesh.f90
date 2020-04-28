!
!--------------------------------------------------------------------
subroutine do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer, intent (in) :: mesh
  real (dp), intent (in) :: zmesh, xmin, dx
  real (dp), intent (out) :: r(0:mesh), sqr(0:mesh), r2(0:mesh)
  !
  integer :: i
  real(dp) :: x
  !
  do i=0,mesh
     x = xmin + dx * i
     r(i)  = exp(x)/zmesh
     sqr(i)= sqrt(r(i))
     r2(i) = r(i) * r(i)
  end do
  write(*,'(/" radial grid information:")')
  write(*,'(" dx =",f10.6,", xmin =",f10.6,", zmesh =",f10.6)') &
       dx,xmin,zmesh
  write(*,'(" mesh =",i6,", r(0) =",f10.6,", r(mesh) =",f10.6)') &
       mesh, r(0), r(mesh)
  write(*,*) 
  !
  return
end subroutine do_mesh
