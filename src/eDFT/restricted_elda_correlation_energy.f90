subroutine restricted_elda_correlation_energy(nEns,aLF,nGrid,weight,rho,Ec)

! Compute the restricted LDA correlation energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: aLF(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r,ec_p

! Output variables

  double precision,intent(out)  :: Ec


! Compute eLDA correlation energy

  Ec = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      ec_p = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))

      Ec = Ec + weight(iG)*ec_p*r

    end if

  end do

end subroutine restricted_elda_correlation_energy