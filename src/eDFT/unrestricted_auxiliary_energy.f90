subroutine unrestricted_auxiliary_energy(nBas,nEns,eps,Eaux,occnum)

! Compute the auxiliary KS energies 

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: eps(nBas,nspin)
  double precision,intent(in)   :: occnum(nBas,nspin,nEns)


! Local variables

  integer                       :: iEns
  integer                       :: ispin
  integer                       :: p

! Output variables

  double precision,intent(out)  :: Eaux(nspin,nEns)

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
