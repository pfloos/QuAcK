subroutine restricted_elda_correlation_individual_energy(aMFL,nGrid,weight,rhow,rho,Ec)

! Compute LDA correlation individual energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: aMFL(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r
  double precision              :: rI
  double precision              :: ec_p,dFcdr

! Output variables

  double precision,intent(out)  :: Ec

! Compute eLDA correlation potential

  Ec = 0d0

  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .and. rI > threshold) then

      ec_p  = aMFL(1)/(1d0 + aMFL(2)*r**(-1d0/6d0) + aMFL(3)*r**(-1d0/3d0))

      dFcdr = aMFL(2)*r**(-1d0/6d0) + 2d0*aMFL(3)*r**(-1d0/3d0)
      dFcdr = dFcdr/(1d0 + aMFL(2)*r**(-1d0/6d0) + aMFL(3)*r**(-1d0/3d0))
      dFcdr = ec_p*dFcdr/(6d0*r)

      Ec = Ec + weight(iG)*(ec_p*rI + dFcdr*r*rI - dFcdr*r*r)

     end if

  end do

end subroutine restricted_elda_correlation_individual_energy
