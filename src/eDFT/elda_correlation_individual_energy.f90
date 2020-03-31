subroutine elda_correlation_individual_energy(nEns,aLF,nGrid,weight,rhow,rho,Ec)

! Compute LDA correlation individual energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: aLF(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra, rb, r
  double precision              :: raI,rbI,rI
  double precision              :: ec_p,dFcdr

! Output variables

  double precision,intent(out)  :: Ec(nsp)

! Compute eLDA correlation potential

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rhow(iG,1))
    rb = max(0d0,rhow(iG,2))

    raI = max(0d0,rho(iG,1))
    rbI = max(0d0,rho(iG,2))

    r  = ra  + rb
    rI = raI + rbI

    if(r > threshold .and. rI > threshold) then

      ec_p  = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))

      dFcdr = aLF(2)*r**(-1d0/6d0) + 2d0*aLF(3)*r**(-1d0/3d0)
      dFcdr = dFcdr/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdr = ec_p*dFcdr/(6d0*r)
      dFcdr = ec_p + dFcdr*r

      Ec(2) = Ec(2) + weight(iG)*rI*dFcdr

     end if

  end do

end subroutine elda_correlation_individual_energy
