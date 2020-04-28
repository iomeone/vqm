!---------------------------------------------------------------
program periodicwell
  !---------------------------------------------------------------
  !
  !     Band structure of a 1d model periodic system (Kronig-Penney)
  !     Expansion on a plane-wave basis set and diagonalization
  !     Units: hbar^2/2m = 1
  !     Requires lapack dsyev, cern library cft 

  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  real(dp), parameter :: pi=3.14159265358979_dp
  integer :: n, npw, nfft, ifft
  real(dp) :: v0, a, b, ecut, k, x, gx
  real(dp), allocatable :: g(:), e(:), h(:,:), work (:), vr(:), vi(:)
  complex(dp) :: f
  integer :: i,j,m, lwork, info
  integer, external :: good_fft_dimension
  !
  !       Input data
  !
  !     Potential well: V(x)=-V_0 for |x|<b/2, V(x)=0 for |x|>b/2
  !     Periodicity:    V(x+a)=V(x)
  !
  write (*,"('Parameters for potential well: V_0, a, b > ')",advance='no')
  read (*,*) v0, a, b
  if ( v0 <= 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp .or. a <= b ) &
     stop ' wrong input parameters '
  write (*,"('   V_0, a, b =',3f10.4)") v0, a, b
  !
  !     Plane waves basis set: G_n=n*2pi/a, \hbar^2/2m*G^2 < Ecut
  !
  write (*,"('Cutoff for plane waves: ecut > ',$)")
  read (*,*) ecut
  if ( ecut <= 0.0_dp) stop ' wrong input parameter '
  !
  !     Number of plane waves
  !
  npw = nint ( sqrt ( ecut/(2.0_dp*pi/a)**2 ) + 0.5_dp )
  npw = 2*npw+1
  write (*,"('Ecut=',f8.2,'  # of PWs=',i4)") ecut,npw
  allocate (g(npw), e(npw), work(3*npw), h(npw,npw) )
  !
  !       Assign values of G_n: n=0,+1,-1,+2,-2, etc
  !
  g(1) = 0.0_dp
  do i=2,npw-1,2
     g(i  ) = (i/2)*2.0_dp*pi/a
     g(i+1) =-(i/2)*2.0_dp*pi/a
  end do
  !
  !     Compute V(G) with FFT
  !     nfft = number of points in real and reciprocal space
  !            must be at least four times the number of plane waves
  !            because we need V(G) for G=G_i-G_j
  !
  nfft = 4*npw
  write (*,"('FFT dimension: min=',i4,'  Your value >')",advance='no') nfft
  read (*,*) nfft
  allocate (vr(nfft), vi(nfft))
  ! compute v(x_i), x_i=(i-1)*dx, dx=a/nfft
  do i=1,nfft
     x=(i-1)*(a/nfft)
     if ( x <= b/2.0_dp .or. x > (a-b/2.0_dp) ) then
        vr(i) = -v0
     else
        vr(i) = 0.0_dp
     end if
     ! imaginary part is zero in this case, for both real and reciprocal space
     vi(i) = 0.0_dp
  end do
  call cft(vr,vi,nfft,nfft,nfft,-1)
  vr = vr/nfft
  ! test: compare v(g) with analytical results
  !write(33,*) 1,0.0_dp,vr(1),-v0/a*b
  !do i=2,nfft
  !   x = (i-1)*(2.0_dp*pi/a)
  !   write(33,*) i,x,vr(i),-v0/a*sin(x*b/2.0_dp)/x*2.0_dp
  !end do
  !
  !     Loop on k-vectors: k runs from -pi/a to pi/a
  !
  open (7,file='bands.out',status='unknown',form='formatted')
  n = 20
  do m=-n,n
     !
     k = m*pi/n/a
     !       cleanup
     h(:,:) = 0.0_dp
     !
     !       Assign values of the matrix elements of the hamiltonian 
     !       on the plane wave basis
     !
     do i=1,npw
        do j=1,npw
           if ( i ==j ) then
              h(i,j) = (k+g(i))**2 + vr(1)
              ! vr(1) = v(G=0) = -v0/a*b
           else
              !
              ! ifft points to the component G=g(j)-g(i) 
              !
              ifft = nint( (g(j)-g(i))*a/pi ) / 2 + 1
              if ( ifft < 1 ) ifft = ifft + nfft
              h(i,j) = vr(ifft)
              ! h(i,j) = -v0/a * sin((g(j)-g(i))*b/2.0_dp) / (g(j)-g(i))*2.0_dp
           end if
           ! print  '(2i4,f12.6)', i,j, h(i,j)
        end do
     end do
     !
     !       Solution [expansion coefficients are stored into h(j,i)
     !                 j=basis function index, i= eigenvalue index]
     !
     lwork = 3*npw
     call dsyev ( 'V', 'U', npw, h, npw, e, work, lwork, info )
     if (info /= 0) stop 'H-matrix diagonalization failed '
     !
     write (*,"('k=',f10.4,'    Lowest eigenvalues:',3f12.6)") k,e(1),e(2),e(3)
     !
     !       Write to output file the band dispersion e(k)
     !
     write (7,'(4f12.6)') k, e(1), e(2), e(3)
  end do
  close(7)
  deallocate ( h, work, e, g)
  !
end program periodicwell
