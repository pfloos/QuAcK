subroutine orthogonalization_matrix(nBas,S,X)

! Compute the orthogonalization matrix X = S^(-1/2)

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)

! Local variables

  logical                       :: debug
  double precision,allocatable  :: UVec(:,:),Uval(:)
  double precision,parameter    :: thresh = 1d-6

  integer                       :: i

! Output variables

  double precision,intent(out)  :: X(nBas,nBas)

  debug = .false.

  allocate(Uvec(nBas,nBas),Uval(nBas))

  write(*,*)
  write(*,*) ' *** Lowdin orthogonalization X = S^(-1/2)  *** '
  write(*,*)

  Uvec = S
  call diagonalize_matrix(nBas,Uvec,Uval)

  do i=1,nBas

    if(Uval(i) > thresh) then 

      Uval(i) = 1d0/sqrt(Uval(i))

    else

      write(*,*) 'Eigenvalue',i,'too small for Lowdin orthogonalization'

    endif

  enddo
  
  call ADAt(nBas,Uvec,Uval,X)

! Print results

  if(debug) then

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Orthogonalization matrix'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,X)
    write(*,*)

  endif

end subroutine orthogonalization_matrix
