subroutine gga_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute GGA correlation energy

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

! Coefficients for ??? GGA exchange functional

! Compute GGA exchange energy

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = rho(iG,1)
    rb = rho(iG,2)
    ga = drho(1,iG,1)**2 + drho(2,iG,1)**2 + drho(3,iG,1)**2
    gb = drho(1,iG,2)**2 + drho(2,iG,2)**2 + drho(3,iG,2)**2

  enddo

end subroutine gga_correlation_energy
