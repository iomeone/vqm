!---------------------------------------------------------------
program helium_hf
!---------------------------------------------------------------
  !
  ! Solve the Hartree(-Fock) equations for He in the ground state
  ! using self-consistency on the potential and the Numerov algorithm 
  ! to integrate the radial Schroedinger equation under an effective
  ! potential (assumed to have spherical symmetry)
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi
  integer :: mesh
  integer :: n, l, i, iter
  real (dp) :: zeta=2.0_dp, zmesh, rmax, xmin, dx, de, beta, tol
  real (dp) :: e, ehx, evion, ekin, etot, etot1
  real (dp), allocatable :: r(:), sqr(:), r2(:), y(:), vpot(:), &
       rho(:), vscf(:), vhx(:), deltav(:)
  !
  ! read atomic charge and other parameters
  !
  write(*,*) 
  write(*,'(" Hartree-Fock calculation for Helium-like atom")')
  write(*,*) 
  ! write (*,'(a,f10.6)') 'Atomic Charge = ', zeta
  write (*,'(" Mixing parameter beta [0.0-1.0] > ",$)')
  read(*,*) beta
  if ( beta < 0.0_dp .or. beta > 1.0_dp ) stop 'beta not in [0.0-1.0]'
  write (*,'(" SCF accuracy > ",$)')
  read (*,*) tol
  if ( tol < 0.0_dp ) stop 'tol should be positive'
  !
  ! initialize logarithmic mesh
  !
  zmesh=  zeta 
  rmax = 100.0_dp
  xmin = -8.0_dp
  dx   =  0.01_dp
  !
  mesh = (log(zmesh*rmax) - xmin)/dx
  !
  allocate ( r(0:mesh), sqr(0:mesh), r2(0:mesh), vpot(0:mesh), y(0:mesh), &
       rho(0:mesh), vscf(0:mesh), vhx(0:mesh), deltav(0:mesh) )
  !
  call do_mesh ( zmesh, xmin, dx, mesh, r, sqr, r2 )
  !
  ! initialize the potential
  !
  call init_pot ( zeta, r, mesh, vpot )
  vscf = vpot
  !
  ! The GS configuration for helium-like atoms is (1s)**2
  !
  n = 1
  l = 0
  !
  iter = 0
1 continue ! this is the entry point of the SCF cycle
  iter = iter + 1
  write(*,'("####################################################")')
  write(*,'(" SCF iteration # ",i3)') iter
  !
  ! solve the schroedinger equation in radial coordinates by Numerov method
  !
  call solve_sheq ( n, l, e, mesh, dx, r, sqr, r2, vscf, zeta, y )
  !
  ! calculate the charge density from the wfc
  !
  call rho_of_r ( mesh, r, r2, y, rho )
  !
  ! calculate the Hartree + Exchange potential and energy
  !
  call v_of_rho ( mesh, dx, r, r2, rho, vhx, ehx )
  !
  ! calculate the kinetic energy and the energy in the external potential
  ! 
  evion = 0.0_dp
  ekin = 2.0_dp * e
  do i=0,mesh
     evion = evion + vpot(i) * rho(i) * fpi * r2(i) * r(i) * dx
     ekin  = ekin - vscf(i) * rho(i) * fpi * r2(i) * r(i) * dx
  end do
  !
  deltav = vpot + vhx -vscf
  vscf = vscf + beta * deltav
  !
  de = 0.0_dp
  do i=0,mesh
     de = de + deltav(i) * rho(i) * fpi * r2(i) * r(i) * dx
  end do
  !
  ! write out the eigenvalue energy to be compared with the external potential
  !
  etot = 2 * e - ehx + de
  etot1= ekin + evion + ehx
  write (*,'(" eigenvalue   = ",f15.6)') e
  write (*,'(" Eigenvalue energy       ",f15.6)' )  2.0_dp * e 
  write (*,'(" Kinetic energy          ",f15.6)' )  ekin
  write (*,'(" External pot. energy    ",f15.6)' )  evion
  write (*,'(" Hartree+Exchange energy ",f15.6)' )  ehx
  write (*,'(" Variational correction  ",f15.6)' )  de
  write (*,'(" Total energy            ",2f15.6)')  etot, etot1
  write (*,'(" Virial check            ",f15.6)' ) -(evion+ehx)/ekin
  !
  ! The variational correction "de" is used to check for convergence
  !
  if ( abs(de) > tol ) go to 1
  
  write (*,*) 
  write (*,'(" SCF Convergence has been achieved")')
  write (*,*) 
  write (*,'(" compute additional single-particle states")')
  write (*,*) 
  !
  ! write potentials to file pot.out
  !
  open (7,file='pot.out',status='unknown',form='formatted')
  do i =0,mesh
     write (7,*) r(i),vpot(i),vhx(i),vscf(i)
  end do
  close(7)
  !
  ! open wfc file
  !
  open (7,file='wfc.out',status='unknown',form='formatted')
  !
2 continue     
  !
  ! read principal and angular quantum numbers
  !
  write (*,'(" n, l >" ,$)') 
  read (*,*) n, l
  if ( n < 1 ) stop
  if ( n < l+1 ) then
     write(*,*) 'error in main: n < l+1 -> wrong number of nodes '
     go to 2
  end if
  call solve_sheq ( n, l, e, mesh, dx, r, sqr, r2, vscf, zeta, y )
  write (*,'(" eigenvalue ",2i2," = ",f15.6)') n,l,e
  do i=0,mesh
     write (7,*) r(i),y(i)/sqr(i), y(i)*sqr(i), e
  end do
  write (7,'(/)')
  go to 2
  
end program helium_hf
!
!--------------------------------------------------------------------
subroutine rho_of_r ( mesh, r, r2, y, rho )
  !--------------------------------------------------------------------
  !
  ! compute the charge density of an He-like atom
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi, &
       nelec=2.0_dp
  
  integer, intent(in) :: mesh
  real (dp), intent(in) :: r(0:mesh), r2(0:mesh), y(0:mesh)
  real (dp), intent(out):: rho(0:mesh)
  
  integer :: i
  
  rho = nelec * y**2 * r / (fpi*r2)
  
  return
end subroutine rho_of_r
!
!--------------------------------------------------------------------
subroutine v_of_rho ( mesh, dx, r, r2, rho, vhx, ehx )
  !--------------------------------------------------------------------
  !
  ! compute the Hartree + Exchange potential for an He-like atom
  ! Exchange cancels exactly half of the Hartree result for both 
  ! potential and energy
  !
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, fpi=4.0_dp*pi, &
       e2=2.0_dp
  integer, intent(in) ::mesh
  real (dp), intent (in) :: dx, r(0:mesh), r2(0:mesh), rho(0:mesh)
  real (dp), intent (out):: vhx(0:mesh), ehx
  
  integer :: i
  real (dp) :: charge
  !
  ! calculate the Hartree potential and energy by integrating the 
  ! electric field generated by the electronic charge. 
  ! This is done in 2 steps
  !
  ! 1) calculate the charge inside a sphere and fill vhx with the 
  !    electric field generated by this charge
  !
  charge = 0.0_dp
  do i=0,mesh
     charge = charge + rho(i) * fpi * r2(i) * r(i) * dx
     vhx(i) = e2*charge/r2(i)
  end do
  !
  !   (the total charge is writen in output as a check)
  ! 
  write (*,'(" Check: total charge = ",f15.6)') charge
  !
  ! 2) integrate the electric field from +\infty to r to get the potential
  !    and integrate V_{Hartree}*rho to get the energy
  !
  ehx = 0.d0
  vhx(mesh) = e2 * charge / r(mesh)
  do i = mesh-1, 0, -1
     vhx(i) = vhx(i+1) + vhx(i) * r(i) * dx
     ehx = ehx + vhx(i) * rho(i) * fpi * r2(i) * r(i) * dx
  end do
  ehx = ehx/2.0_dp
  !
  ! Exchange cancels exactly half of the Hartree result for both 
  ! potential and energy
  !
  ehx = 0.5_dp * ehx 
  vhx = 0.5_dp * vhx
  !
  return
end subroutine v_of_rho
