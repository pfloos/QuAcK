subroutine unrestricted_mgga_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,drhow,rho,drho,Ex)

! Compute MGGA exchange energy for individual states

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

    case default

      call print_warning('!!! MGGA exchange individual energy not available !!!')
      stop

  end select

end subroutine unrestricted_mgga_exchange_individual_energy
