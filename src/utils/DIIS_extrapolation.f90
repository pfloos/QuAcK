subroutine DIIS_extrapolation(rcond,n_err,n_e,n_diis,error,e,error_in,e_inout)

! Perform DIIS extrapolation

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: n_err,n_e
  double precision,intent(in)   :: error_in(n_err)
  double precision,intent(in)   :: error(n_err,n_diis)
  double precision,intent(in)   :: e(n_e,n_diis)

! Local variables

  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: b(:)
  double precision,allocatable  :: w(:)

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

!  A(1:n_diis,1:n_diis) = matmul(transpose(error),error)
  call dgemm('T','N',n_diis,n_diis,n_err,1d0,error(1,1),n_err,error(1,1),n_err,0d0,A(1:n_diis,1:n_diis),n_diis)

  A(1:n_diis,n_diis+1) = -1d0
  A(n_diis+1,1:n_diis) = -1d0
  A(n_diis+1,n_diis+1) = +0d0

! Build x matrix

  b(1:n_diis) = +0d0
  b(n_diis+1) = -1d0

! Solve linear system
  call linear_solve(n_diis+1,A,b,w,rcond)

! Extrapolate

!  e_inout(:) = matmul(w(1:n_diis),transpose(e(:,1:n_diis)))
  call dgemm('N','T',1,n_e,n_diis,1d0,w(1),1,e(1,1),n_e,0d0,e_inout(1),1)

  
  deallocate(A,b,w)
end subroutine 
