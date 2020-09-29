subroutine US51_lda_exchange_individual_energy(nGrid,weight,rhow,rho,doNcentered,Ex)

! Compute the restricted version of Slater's LDA exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: doNcentered


! Local variables

  integer                       :: iG
  double precision              :: r,rI,alpha
  double precision              :: e,dedr
  double precision              :: nEli,nElw

! Output variables

  double precision,intent(out)  :: Ex

! External variable

  double precision,external     :: electron_number

  nEli = electron_number(nGrid,weight,rho)

  nElw = electron_number(nGrid,weight,rhow)


! Compute LDA exchange matrix in the AO basis

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold) then

      e    =         alpha*r**(1d0/3d0)
      dedr = 1d0/3d0*alpha*r**(-2d0/3d0)

      if (doNcentered == 0) then
        Ex = Ex - weight(iG)*dedr*r*r      
      else
        Ex = Ex - weight(iG)*dedr*r*r*(nEli/nElw)
      end if

      if(rI > threshold) then

        if (doNcentered == 0) then
          Ex = Ex + weight(iG)*(e*rI + dedr*r*rI)
        else
          Ex = Ex + weight(iG)*((nEli/nElw)*e*rI + dedr*r*rI)
        end if

      endif

    endif

  enddo

end subroutine US51_lda_exchange_individual_energy
