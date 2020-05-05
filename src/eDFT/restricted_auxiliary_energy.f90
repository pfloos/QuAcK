subroutine restricted_auxiliary_energy(nBas,nEns,nO,eps,Eaux)

! Compute the auxiliary KS energies 

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  integer,intent(in)            :: nO
  double precision,intent(in)   :: eps(nBas)

! Local variables

  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: Eaux(nEns)

! Ground state density matrix

  iEns = 1

  Eaux(iEns) = 2d0*sum(eps(1:nO))

! Singly-excited state density matrix

  iEns = 2 
  
  if(nO > 1) then 
    Eaux(iEns) = 2d0*sum(eps(1:nO-1))
  else
    Eaux(iEns) = 0d0
  end if

  Eaux(iEns) = Eaux(iEns) + eps(nO) + eps(nO+2)

! Doubly-excited state density matrix

  iEns = 3 
  
  if(nO > 1) then 
    Eaux(iEns) = 2d0*sum(eps(1:nO-1))
  else
    Eaux(iEns) = 0d0
  end if

  Eaux(iEns) = Eaux(iEns) + 2d0*eps(nO+1)

end subroutine restricted_auxiliary_energy
