subroutine mgga_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,drhow,rho,drho,LZx,Ex)

! Compute MGGA exchange energy for individual states

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

! Output variables

  double precision              :: LZx(nspin)
  double precision              :: Ex(nspin,nEns)

! Select correlation functional

  select case (DFA)

    case default

      call print_warning('!!! MGGA exchange individual energy not available !!!')
      stop

  end select

end subroutine mgga_exchange_individual_energy
