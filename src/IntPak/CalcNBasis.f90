subroutine CalcNBasis(nShell,atot,NBasis)

  implicit none

! Input variables

  integer,intent(in)   :: nShell
  integer,intent(in)   :: atot(nShell)

! Local variables

  integer              :: iShell

! Output variables

  integer,intent(out)  :: NBasis

  NBasis = 0
  do iShell=1,nShell
    NBasis = NBasis + (atot(iShell)*atot(iShell) + 3*atot(iShell) + 2)/2
  enddo

  write(*,'(A28)') '------------------'
  write(*,'(A28,1X,I16)') 'Number of basis functions',NBasis
  write(*,'(A28)') '------------------'
  write(*,*)

end subroutine CalcNBasis
