subroutine read_basis_pyscf(nBas,nO,nV)

! Read basis set information from PySCF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)           :: nO(nspin)

! Local variables

! Output variables

  integer,intent(out)           :: nV(nspin)
  integer,intent(out)           :: nBas

!------------------------------------------------------------------------
! Primary basis set information
!------------------------------------------------------------------------

  open(unit=3,file='int/nBas.dat')
    read(3, *) nBas
  close(unit=3)

!  write(*,'(A38)') '--------------------------------------'
!  write(*,'(A38,1X,I16)') 'Number of basis functions (AOs)', nBas
!  write(*,'(A38)') '--------------------------------------'
!  write(*,*)

! Number of virtual orbitals

  nV(:) = nBas - nO(:)

end subroutine 
