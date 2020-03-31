subroutine restricted_elda_correlation_energy(aMFL,nGrid,weight,rho,Ec)

! Compute the restricted LDA correlation energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: aMFL(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r,e

! Output variables

  double precision,intent(out)  :: Ec


! Compute eLDA correlation energy

  Ec = 0d0
  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then
      e = aMFL(1)/(1d0 + aMFL(2)*r**(-1d0/6d0) + aMFL(3)*r**(-1d0/3d0))
      Ec = Ec + weight(iG)*e*r
    end if

  end do

end subroutine restricted_elda_correlation_energy
