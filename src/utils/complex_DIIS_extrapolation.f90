subroutine complex_DIIS_extrapolation(rcond,n_err,n_e,n_diis,error,e,error_in,e_inout)

! Perform DIIS extrapolation

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: n_err,n_e
  complex*16,intent(in)         :: error_in(n_err)
  complex*16,intent(in)         :: error(n_err,n_diis)
  complex*16,intent(in)         :: e(n_e,n_diis)

! Local variables

  complex*16,allocatable        :: A(:,:)
  complex*16,allocatable        :: b(:)
  complex*16,allocatable        :: w(:)

! Output variables

  double precision,intent(out)  :: rcond
  integer,intent(inout)         :: n_diis
  complex*16,intent(inout)      :: e_inout(n_e)

! Memory allocation
  allocate(A(n_diis+1,n_diis+1),b(n_diis+1),w(n_diis+1))

! Update DIIS "history"
  call complex_prepend(n_err,n_diis,error,error_in)
  call complex_prepend(n_e,n_diis,e,e_inout)
  
!  Build A matrix

  A(1:n_diis,1:n_diis) = matmul(transpose(error),error)

  A(1:n_diis,n_diis+1) = cmplx(-1d0,0d0,kind=8)
  A(n_diis+1,1:n_diis) = cmplx(-1d0,0d0,kind=8)
  A(n_diis+1,n_diis+1) = cmplx(0d0,0d0,kind=8)

! Build x matrix

  b(1:n_diis) = cmplx(+0d0,0d0,kind=8)
  b(n_diis+1) = cmplx(-1d0,0d0,kind=8)

  ! Solve linear system
  call complex_linear_solve(n_diis+1,A,b,w,rcond)
  ! Perform extrapolation only if the system is well-conditioned
  e_inout(:) = matmul(w(1:n_diis),transpose(e(:,1:n_diis)))

! Extrapolate

  e_inout(:) = matmul(w(1:n_diis),transpose(e(:,1:n_diis)))
  deallocate(A,b,w)
end subroutine 
