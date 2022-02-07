subroutine W38_lda_correlation_individual_energy(nGrid,weight,rhow,rho,Ec)

! Compute the unrestricted version of the Wigner's LDA individual energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r
  double precision              :: raI,rbI,rI
  double precision              :: a,d,epsc
  double precision              :: dFcdra,dFcdrb

! Output variables

  double precision,intent(out)  :: Ec(nsp)

! Coefficients for Wigner's LDA correlation

  a = 0.04918d0
  d = 0.349d0

! Compute LDA correlation individual energy

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rhow(iG,1))
    rb = max(0d0,rhow(iG,2))

    raI = max(0d0,rho(iG,1))
    rbI = max(0d0,rho(iG,2))
    
    r  = ra + rb
    rI = raI + rbI

    if(r > threshold .or. rI > threshold) then

      epsc  = ra*rb/(r + d*r**(2d0/3d0))
      dFcdra = epsc*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 1d0/r + 1d0/ra)
      dFcdrb = epsc*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 1d0/r + 1d0/rb)
  
      Ec(2) = Ec(2) + weight(iG)*rI*0.5d0*(dFcdra + dFcdrb)

    endif

  enddo

  Ec(2) = -4d0*a*Ec(2)

end subroutine W38_lda_correlation_individual_energy
