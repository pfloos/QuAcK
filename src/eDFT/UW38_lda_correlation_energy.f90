subroutine UW38_lda_correlation_energy(nGrid,weight,rho,Ec)

! Compute the unrestricted version of the Wigner's LDA correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r
  double precision              :: a,d,epsc

! Output variables

  double precision              :: Ec(nsp)

! Coefficients for Wigner's LDA correlation

  a = 0.04918d0
  d = 0.349d0

! Compute LDA correlation energy

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))

    if(ra > threshold .and. rb > threshold) then

      r = ra + rb

      epsc = ra*rb/(r + d*r**(2d0/3d0))

      Ec(2) = Ec(2) + weight(iG)*epsc

    endif

  enddo

  Ec(2) = -4d0*a*Ec(2)

end subroutine UW38_lda_correlation_energy
