subroutine hybrid_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,nBas,ERI,Pw,rhow,drhow, &
                                                          P,rho,drho,LZx,Ex)

! Compute the hybrid exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin,nEns)

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas,nEns)

! Output variables

  double precision              :: LZx(nspin)
  double precision              :: Ex(nspin,nEns)

! Select correlation functional

  select case (DFA)

    case (1)

      call fock_exchange_individual_energy(nBas,nEns,Pw,P,ERI,LZx,Ex)

    case default

      call print_warning('!!! Hybrid exchange individual energy not available !!!')
      stop

  end select

end subroutine hybrid_exchange_individual_energy
