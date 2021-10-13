subroutine unrestricted_mgga_exchange_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

! Select MGGA exchange functional for potential calculation

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Select MGGA exchange functional

  select case (DFA)

    case default

      call print_warning('!!! MGGA exchange potential not available !!!')
      stop

  end select

end subroutine unrestricted_mgga_exchange_potential
