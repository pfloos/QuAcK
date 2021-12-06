subroutine unrestricted_auxiliary_energy(nBas,nEns,eps,occnum,Eaux)

! Compute the auxiliary KS energies 

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: eps(nBas,nspin)
  double precision,intent(in)   :: occnum(nBas,nspin,nEns)

! Local variables

  integer                       :: iEns,iBas
  integer                       :: ispin
  integer                       :: p
  double precision,allocatable  :: nEl(:) 


! Output variables

  double precision,intent(out)  :: Eaux(nspin,nEns)

! Memory allocation

  allocate(nEl(nEns))

! Compute the number of electrons

  nEl(:) = 0d0
  do iEns=1,nEns
    do iBas=1,nBas
      nEl(iEns) = nEl(iEns) + occnum(iBas,1,iEns) + occnum(iBas,2,iEns)
    end do
  end do

! Compute auxiliary energies for each state of the ensemble based on occupation numbers

  Eaux(:,:) = 0d0
  do iEns=1,nEns
    do ispin=1,nspin
      do p=1,nBas

        Eaux(ispin,iEns) = Eaux(ispin,iEns) + occnum(p,ispin,iEns)*eps(p,ispin) 

      end do
    end do

  end do

end subroutine unrestricted_auxiliary_energy
