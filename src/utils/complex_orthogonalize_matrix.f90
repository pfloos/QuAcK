subroutine complex_orthogonalize_matrix(nBas,S,X)

! Compute the orthogonalization matrix X

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  complex*16,intent(in)         :: S(nBas,nBas)

! Local variables

  logical                       :: debug
  complex*16,allocatable        :: UVec(:,:),Uval(:)
  double precision,parameter    :: thresh = 1d-6
  complex*16,allocatable        :: D(:,:)

  integer                       :: i

! Output variables

  complex*16,intent(out)        :: X(nBas,nBas)

  debug = .false.

  allocate(Uvec(nBas,nBas),Uval(nBas),D(nBas,nBas))


!   write(*,*)
!   write(*,*) ' Lowdin orthogonalization'
!   write(*,*)

    Uvec = S
    call complex_diagonalize_matrix(nBas,Uvec,Uval)

    do i=1,nBas

      if(abs(Uval(i)) < thresh) then 

        write(*,*) 'Eigenvalue',i,' is very small in Lowdin orthogonalization = ',Uval(i)

      endif
      Uval(i) = 1d0/sqrt(Uval(i))

    enddo
    call diag(nBas,Uval, D)
    X = matmul(Uvec,matmul(D,conjg(transpose(Uvec))))
    deallocate(Uvec,Uval,D)
end subroutine 

subroutine diag(n,vec,M)
! Create diag matrix from vector

implicit none

! input variables

integer,intent(in) :: n
complex*16,intent(in) :: vec(n)

! local variables
integer :: i

! Output variables
complex*16   :: M(n,n)

M(:,:) = 0d0
do i=1,n
  M(i,i) = vec(i)
enddo
end subroutine
