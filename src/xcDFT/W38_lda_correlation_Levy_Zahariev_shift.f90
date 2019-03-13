subroutine W38_lda_correlation_Levy_Zahariev_shift(nGrid,weight,rho,EcLZ)

! Compute Wigner's LDA contribution to Levy-Zahariev shift

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
  double precision              :: dFcdra,dFcdrb

! Output variables

  double precision,intent(out)  :: EcLZ(nsp)

! Coefficients for Wigner's LDA correlation

  a = 0.04918d0
  d = 0.349d0

! Compute LDA correlation matrix in the AO basis

  EcLZ(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))

    r  = ra + rb

    if(r > threshold) then

      epsc  = ra*rb/(r + d*r**(2d0/3d0))
      dFcdra = epsc*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 2d0/r + 1d0/ra)
      dFcdrb = epsc*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 2d0/r + 1d0/rb)
  
      EcLZ(2) = EcLZ(2) - weight(iG)*r*r*(dFcdra + dFcdrb)

    endif

  enddo

  EcLZ(:) = -4d0*a*EcLZ(:)

end subroutine W38_lda_correlation_Levy_Zahariev_shift
