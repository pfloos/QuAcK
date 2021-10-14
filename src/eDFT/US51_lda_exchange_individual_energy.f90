subroutine US51_lda_exchange_individual_energy(nGrid,weight,rhow,rho,doNcentered,kappa,Ex)

! Compute the restricted version of Slater's LDA exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  logical,intent(in)            :: doNcentered
  double precision,intent(in)   :: kappa

! Local variables

  integer                       :: iG
  double precision              :: r,rI
  double precision              :: e,dedr
  double precision              :: Exrr,ExrI,ExrrI

! Output variables

  double precision,intent(out)  :: Ex

! Compute LDA exchange matrix in the AO basis

  Exrr  = 0d0
  ExrI  = 0d0
  ExrrI = 0d0

  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold) then

      e    =         CxLSDA*r**(1d0/3d0)
      dedr = 1d0/3d0*CxLSDA*r**(-2d0/3d0)

      Exrr = Exrr - weight(iG)*dedr*r*r      

      if(rI > threshold) then

        ExrI  = ExrI  + weight(iG)*e*rI
        ExrrI = ExrrI + weight(iG)*dedr*r*rI

      endif

    endif

  enddo

  if(doNcentered) then
  
    Exrr  = kappa*Exrr
    ExrI  = kappa*ExrI

  endif

  Ex = Exrr + ExrI + ExrrI

end subroutine US51_lda_exchange_individual_energy
