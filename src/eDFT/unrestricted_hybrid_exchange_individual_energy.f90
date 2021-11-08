subroutine unrestricted_hybrid_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,nBas,ERI,Pw,P,rhow,drhow,rho,drho,Ex)

! Compute the hybrid exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)

! Output variables

  double precision              :: Ex

! Select correlation functional

  select case (DFA)

    case (1)

      call unrestricted_fock_exchange_individual_energy(nBas,Pw,P,ERI,Ex)

    case default

      call print_warning('!!! Hybrid exchange individual energy not available !!!')
      stop

  end select

end subroutine unrestricted_hybrid_exchange_individual_energy
