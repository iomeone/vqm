!---------------------------------------------------------------
program helium_hf_gauss
  !---------------------------------------------------------------
  !
  !    Hartree-Fock ground-state solution for the Helium atom
  !    S=0 state - equivalent to Hartree approximation
  !    expansion on a gaussian basis set and self-consistency
  !    with diagonalization of the Fock matrix
  !
  !    Requires lapack and blas
  !
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  real (dp), parameter :: pi=3.14159265358979_dp, tpi=2.0_dp*pi, &
       e2=2.0_dp, zeta=2.0_dp
  integer :: ia,ib,ic,id, i, n_alpha, iter, maxter=100
  real (dp) :: aa, eold, enew
  real (dp), allocatable ::  alpha(:), e(:), c(:)
  real (dp), allocatable ::  h(:,:), f(:,:), s(:,:), v(:,:)
  real (dp), allocatable ::  q(:,:,:,:)
  character (len=80) :: filename
  integer :: iunit, ios
  !
  !       Input data
  !
  ! write (*,'(" Carica atomica = ", f10.6)')
  write (*,'(" Parameters of the Gaussians from file >> ",$)')
  read (*,'(a)') filename
  if ( trim(filename) == '') then
     !  read from terminal
     iunit = 5
  else
     write (*,*)
     !  read from specified file
     iunit = 10
     open (unit=iunit, file=trim(filename), status='old', &
           form='formatted', iostat=ios )
     if ( ios /= 0 ) stop 'file not found'
  end if 
  if ( iunit == 5 ) write (*,'(" Number of Gaussians >> ",$)')
  read (iunit,*) n_alpha
  allocate (alpha (n_alpha))
  do i=1,n_alpha
     if ( iunit == 5 ) then
        write (*,'(" Gaussian # ",i3,",  coefficient >> ",$)') i
        read (iunit,*) alpha(i)
     else
        read (iunit,*) alpha(i)
        write (*,'(" Gaussian # ",i3,",   coefficient >> ",f10.6)') i, alpha(i)
     end if
     if ( alpha(i) <= 0.0_dp ) stop 'wrong coefficient'
  end do
  if ( iunit /= 5 ) close (unit = iunit, status='keep')
  !
  !     Fill the Q matrix with the matrix elements of the e-e interaction:
  !       q(i,j,k,l) = \int b_i(r) b_j(r') 1/|r-r'| b_k(r) b_l(r') dr dr'
  !     where b_i(r) = exp (-alpha_i r^2) are the gaussians  
  !
  allocate ( q(n_alpha, n_alpha, n_alpha, n_alpha) )
  !
  do id=1,n_alpha
     do ic=1,n_alpha
        do ib=1,n_alpha
           do ia=1,n_alpha
              q(ia,ib,ic,id) = 4.0_dp*pi**(2.5_dp) / &
                   (alpha(ia)+alpha(ic))/(alpha(ib)+alpha(id)) / &
                    sqrt(alpha(ia)+alpha(ib)+alpha(ic)+alpha(id))
           end do
        end do
     end do
  end do
  !
  !       Assign values of overlap integrals S and of matrix elements 
  !       of the one-electron hamiltonian H on the gaussian basis:
  !
  allocate ( s(n_alpha, n_alpha), h(n_alpha, n_alpha), f (n_alpha, n_alpha), &
       v(n_alpha, n_alpha), c(n_alpha), e (n_alpha)  )
  !
  do ib=1,n_alpha
     do ia=1,n_alpha
        aa = alpha(ia)+alpha(ib)
        s(ia,ib) = (pi/aa)**1.5d0
        h(ia,ib) = s(ia,ib)*6.d0*alpha(ia)*alpha(ib)/aa - e2*zeta*tpi/aa 
     end do
  end do
  !
  !       Starting solution (very lousy guess)
  !
  do ia=1,n_alpha
     c(ia) = 0.d0
  end do
  c(1) = 1.d0
  !
  !       Self-consistency iteration
  !
  enew = 0.d0
  write (*,*)
  !
  do iter=1,maxter
     !
     !       Fill the Fock matrix
     !
     do ia=1,n_alpha
        do ib=1,n_alpha
           f(ia,ib) = h(ia,ib)
           do id=1,n_alpha
              do ic=1,n_alpha
                 ! Hartree term only:
                 ! f(ia,ib) = f(ia,ib) + q(ia,ic,ib,id)*c(ic)*c(id)
                 ! Hartree-Fock (yields same results in this case:
                 f(ia,ib) = f(ia,ib) + ( q(ia,ic,ib,id)*2.0_dp - &
                                         q(ia,ib,ic,id)) * c(ic)*c(id)
              end do
           end do
        end do
     end do
     !
     !       Solution [expansion coefficients are stored into v(j,i)
     !                 j=basis function index, i= eigenvalue index]
     !
     call diag ( n_alpha, n_alpha, f, s,  e, v )
     !
     c(:) =  v(:,1)
     !
     eold = enew
     enew = 0.d0
     do ia=1,n_alpha
        do ib=1,n_alpha
           enew = enew + 2.0_dp*h(ia,ib)*c(ia)*c(ib)
           do ic=1,n_alpha
              do id=1,n_alpha
                 enew = enew + q(ia,ic,ib,id)*c(ia)*c(ib)*c(ic)*c(id)
              end do
           end do
        end do
     end do
     write (*,'(" Iteration # ",i3,":  HF eigenvalue, energy = ",2f12.6)') &
        iter, e(1), enew
     !
     if ( abs (enew-eold) < 1.d-8 ) then
        write (*,'(/" Convergence achieved, stopping")')
        deallocate ( e, c, v, f, h, s, alpha )
        stop 
     end if
     !
  end do
  !
  write (*,'(/" Convergence not reached, stopping")')
  deallocate ( e, c, v, f, h, s, alpha )
  stop 
  !
end program helium_hf_gauss
