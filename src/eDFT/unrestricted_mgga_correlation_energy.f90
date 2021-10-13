subroutine unrestricted_mgga_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute unrestricted MGGA correlation energy

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

    case default

      call print_warning('!!! MGGA correlation energy not available !!!')
      stop

  end select

end subroutine unrestricted_mgga_correlation_energy
