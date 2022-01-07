subroutine gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

! Select GGA exchange functional for energy calculation

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

    case (1)

      call G96_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

    case (2)

      call B88_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

    case (3)

      call PBE_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

    case default

      call print_warning('!!! GGA exchange energy not available !!!')
      stop

  end select

end subroutine gga_exchange_energy
