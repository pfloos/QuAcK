subroutine rotate_orbitals(N, Kap, c)
  
  implicit none

  ! Input
  integer, intent(in) :: N
  double precision, intent(in) :: Kap(N,N)
  double precision, intent(inout) :: c(N,N)
  
  ! Local
  double precision, allocatable :: ExpKap(:,:)

  allocate(ExpKap(N,N))

  ! Compute matrix exponential of Kap
  call matrix_exp(N, Kap, ExpKap)

  ! Update c <- c * Exp(Kap)
  c = matmul(c, ExpKap)

  ! Free temporary array
  deallocate(ExpKap)

end subroutine
