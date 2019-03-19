subroutine gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

! Select GGA exchange functional for energy calculation

  implicit none

  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Output variables

  double precision              :: Ex

  select case (DFA)

!   Gill's 96 exchange functional 

    case ('G96')

      call G96_gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

!   Becke's 88 exchange functional 

    case ('B88')

      call B88_gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

    case default

      call print_warning('!!! GGA exchange functional not available !!!')
      stop

  end select

end subroutine gga_exchange_energy
