subroutine unrestricted_gga_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute unrestricted GGA correlation potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(3,nGrid,nspin)

! Local variables

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Select GGA exchange functional

  select case (DFA)

    case (1)

      call ULYP_gga_correlation_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

    case (2)

      call UPBE_gga_correlation_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

    case default

      call print_warning('!!! GGA correlation potential not available !!!')
      stop

  end select

end subroutine unrestricted_gga_correlation_potential
