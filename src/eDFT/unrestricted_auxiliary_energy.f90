subroutine unrestricted_auxiliary_energy(nBas,nEns,nO,eps,Eaux)

! Compute the auxiliary KS energies 

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: eps(nBas,nspin)

! Local variables

  integer                       :: iEns,ispin

! Output variables

  double precision,intent(out)  :: Eaux(nspin,nEns)

! N-electron ground state

  iEns = 1
  do ispin=1,nspin
    Eaux(ispin,iEns) = sum(eps(1:nO(ispin),ispin))
  end do

! (N-1)-electron ground state

  iEns = 2 
 
  Eaux(1,iEns) = sum(eps(1:nO(1),1))

  if(nO(2) > 1) then 
    Eaux(2,iEns) = sum(eps(1:nO(2)-1,2))
  else
    Eaux(2,iEns) = 0d0
  end if

! (N+1)-electron ground state

  iEns = 3 
 
  Eaux(1,iEns) = sum(eps(1:nO(1)+1,1))
  Eaux(2,iEns) = sum(eps(1:nO(2),2)) 

end subroutine unrestricted_auxiliary_energy
