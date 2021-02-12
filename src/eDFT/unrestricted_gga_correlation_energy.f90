subroutine unrestricted_gga_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute unrestricted GGA correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,ga,gb

! Output variables

  double precision              :: Ec(nsp)

  select case (DFA)

    case ('LYP')

      call ULYP_gga_correlation_energy(nGrid,weight,rho,drho,Ec)

    case default

      call print_warning('!!! GGA correlation energy not available !!!')
      stop

  end select

end subroutine unrestricted_gga_correlation_energy
