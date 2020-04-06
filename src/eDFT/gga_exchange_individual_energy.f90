subroutine gga_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,drhow,rho,drho,Ex)

! Compute GGA exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Output variables

  double precision              :: Ex

! Select correlation functional

  select case (DFA)

    case ('RB88')

      call RB88_gga_exchange_individual_energy(nGrid,weight(:),rhow(:),drhow(:,:),rho(:),drho(:,:),Ex)

    case default

      call print_warning('!!! GGA exchange individual energy not available !!!')
      stop

  end select

end subroutine gga_exchange_individual_energy
