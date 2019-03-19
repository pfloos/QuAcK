subroutine Newton(nWSq,gradient,hessian,cWeight)

! Calculate the Green functions

  implicit none

! Input variables

  integer,intent(in)            :: nWSq
  double precision,intent(in)   :: gradient(nWSq),hessian(nWSq,nWSq)

! Local variables

  integer                       :: info
  integer,allocatable           :: ipiv(:)
  double precision,allocatable  :: scr(:),eigval(:),eigvec(:,:)

! Output variables
  
  double precision,intent(inout):: cWeight(nWSq)

! Memory allocation

  allocate(ipiv(nWSq),scr(3*nWsq),eigval(nWSq),eigvec(nWSq,nWSq))

! Compute eigenvectors and eigenvalues

  eigvec = hessian
  call dsyev('V','U',nWSq,eigvec,nWSq,eigval,scr,3*nWSq,info)

  if(info /= 0)then
    write(*,*) ' Problem with dsyev!'
    stop
  endif

  write(*,*)
  write(*,*) 'Eigenvalues of hessian'
  call matout(nWSq,1,eigval)
  write(*,*)
!  write(*,*) 'Eigenvectors of hessian'
!  call matout(nWSq,1,eigval)
!  write(*,*)

! Compute inverse of the hessian

  call dgetrf(nWSq,nWSq,hessian,nWSq,ipiv,info)

  if(info /= 0) then 
    write(*,*) ' Problem in dgetrf!'
    stop
  endif

  call dgetri(nWSq,hessian,nWSq,ipiv,scr,nWSq,info)

  if(info /= 0) then
    write(*,*) ' Problem in dgetri!'
    stop
  endif

  print*,'inverse hessian'
  call matout(nWSq,nWSq,hessian)

! Compute new coefficients

  cWeight = cWeight - matmul(hessian,gradient)

end subroutine Newton
