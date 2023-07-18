subroutine DIIS_extrapolation(rcond,n_err,n_e,n_diis,error,e,error_in,e_inout)

! Perform DIIS extrapolation

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: n_err,n_e
  double precision,intent(in)   :: error_in(n_err),error(n_err,n_diis),e(n_e,n_diis)

! Local variables

  double precision,allocatable  :: A(:,:),b(:),w(:)

! Output variables

  double precision,intent(out)  :: rcond
  integer,intent(inout)         :: n_diis
  double precision,intent(inout):: e_inout(n_e)

! Memory allocaiton

  allocate(A(n_diis+1,n_diis+1),b(n_diis+1),w(n_diis+1))

! Update DIIS "history"

  call prepend(n_err,n_diis,error,error_in)
  call prepend(n_e,n_diis,e,e_inout)

!  Build A matrix

  A(1:n_diis,1:n_diis) = matmul(transpose(error),error)

  A(1:n_diis,n_diis+1) = -1d0
  A(n_diis+1,1:n_diis) = -1d0
  A(n_diis+1,n_diis+1) = +0d0

! Build x matrix

  b(1:n_diis) = +0d0
  b(n_diis+1) = -1d0

! Solve linear system

! call easy_linear_solve(n_diis+1,A,b,w)
  call linear_solve(n_diis+1,A,b,w,rcond)

! Extrapolate

  e_inout(:) = matmul(w(1:n_diis),transpose(e(:,1:n_diis)))

end subroutine 
