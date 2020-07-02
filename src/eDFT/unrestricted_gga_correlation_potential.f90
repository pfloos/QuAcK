subroutine unrestricted_gga_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute unrestricted GGA correlation potential

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
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(3,nGrid,nspin)

! Local variables

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Coefficients for GGA correlation functional

! Compute GGA correlation matrix in the AO basis

end subroutine unrestricted_gga_correlation_potential
