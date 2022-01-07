subroutine mgga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

! Select MGGA exchange functional for energy calculation

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Output variables

  double precision              :: Ex

  select case (DFA)

    case default

      call print_warning('!!! MGGA exchange energy not available !!!')
      stop

  end select

end subroutine mgga_exchange_energy
