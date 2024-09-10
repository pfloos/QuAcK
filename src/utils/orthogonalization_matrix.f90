subroutine orthogonalization_matrix(nBas,nOrb,S,X)

! Compute the orthogonalization matrix X

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)

! Local variables

  logical                       :: debug
  double precision,allocatable  :: UVec(:,:),Uval(:)
  double precision,parameter    :: thresh = 1d-6

  integer                       :: i,j,j0

! Output variables

  integer                       :: nOrb
  double precision,intent(out)  :: X(nBas,nBas)

  debug = .false.

  allocate(Uvec(nBas,nBas),Uval(nBas))

  Uvec(1:nBas,1:nBas) = S(1:nBas,1:nBas)
  call diagonalize_matrix(nBas,Uvec,Uval)

  nOrb = 0
  do i = 1,nBas
    if(Uval(i) > thresh) then
        Uval(i) = 1d0 / dsqrt(Uval(i))
        nOrb = nOrb + 1
    else
      write(*,*) ' Eigenvalue',i,'too small for canonical orthogonalization'
    end if
  end do

  write(*,'(A50)') '------------------------------------------------'
  write(*,'(A40,1X,I5)') 'Number of basis functions     = ',nBas
  write(*,'(A40,1X,I5)') 'Number of spatial orbitals    = ',nOrb
  write(*,'(A40,1X,I5)') 'Number of discarded functions = ',nBas - nOrb
  write(*,'(A50)') '------------------------------------------------'
  write(*,*)

  j0 = nBas - nOrb

  do j = j0+1,nBas
    do i = 1,nBas
      X(i,j-j0) = Uvec(i,j) * Uval(j)
    enddo
  enddo

  deallocate(Uvec,Uval)

! Print results

  if(debug) then

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Orthogonalization matrix'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,X)
    write(*,*)

  end if

end subroutine 

! ---

