!-----------------------------------------------------------------------
subroutine diag ( n, ldh, h, s, e, v )
!     =================
!
!       Finds eigenvalues and eigenvectors of the generalized problem
!       Hv=eSv, where H=hermitian matrix, S=overlap matrix
!       On return S and H are unchanged.
!       Uses level-2 blas dgemm for matrix-matrix multiplication
!       Uses lapack dsyev for matrix diagonalization, checks for
!       too small eigenvalues of S, removes corresponding eigenvectors
!
!-----------------------------------------------------------------------
!
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)

  integer, intent(in) :: n, &! dimension of the matrix to be diagonalized
                         ldh ! leading dimension of h,s and v, as declared
                             ! in the calling pgm unit
  real (dp), intent(in) :: h(ldh,n), & ! matrix to be diagonalized
                           s(ldh,n)    ! overlap matrix
  real (dp), intent(out):: e(n), &     ! eigenvalues
                           v(ldh,n)    ! eigenvectors (column-wise)
  !
  ! local variables
  integer lwork, info, i, j, nn
  !
  real (dp), parameter :: small=1.d-10
  real (dp), allocatable :: work(:), aux(:,:), h1(:,:)
  !
  info = 0
  lwork=3*n
  allocate (work(lwork), aux(ldh,n))
  !
  !       Copy S into an auxiliary matrix because dsyev destroys the matrix
  !
  aux = s
  !
  !       Diagonalize S
  !
  call dsyev ( 'V', 'U', n, aux, ldh, e, work, lwork, info )
  if (info /= 0) stop 'S-matrix diagonalization failed '
  !
  !   Keep only linearly independent combinations (within a given threshold)
  !   store into matrix "aux" the eigenvectors of S divided by the squares
  !   of the eigenvectors
  !
  nn = 0
  do i=1,n
     if (e(i) > small) then
        nn = nn + 1
        aux(:,nn) = aux(:,i) / sqrt(e(i))
     end if
  end do
  if ( nn < n) write(*,*) " # of linearly independent vectors =", nn, n
  !
  !   Transform H using the "aux" matrix
  !
  !   V(i,j) = \sum_{k=1}^{n} H(i,k) aux(k,j),  i=1,n, j=1,nn
  !
  call dgemm ( 'N', 'N', n, nn, n, 1.0_dp, h, ldh, aux, ldh, 0.0_dp, v, ldh )
  !
  !    h1(i,j) = \sum_{k=1}^{n} aux(k,i) v(k,j),  i=1,nn, j=1,nn
  !    H' = transpose(aux)*H*aux
  !
  allocate (h1(nn,nn) )
  call dgemm ( 'T', 'N', nn, nn, n, 1.0_dp, aux, ldh, v, ldh, 0.0_dp, h1, nn )
  !
  !    Diagonalize transformed H
  !
  info = 0
  call dsyev ('V', 'U', nn, h1, nn, e, work, lwork, info)
  if (info /= 0) stop 'H-matrix diagonalization failed '
  !
  !    Back-transform eigenvectors
  !
  call dgemm ( 'N', 'N', n, nn, nn, 1.0_dp, aux, ldh, h1, nn, 0.0_dp, v, ldh )
  !
  deallocate (h1, aux, work)
  !
end subroutine diag

