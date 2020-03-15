subroutine restricted_elda_correlation_Levy_Zahariev_shift(nEns,aLF,nGrid,weight,rho,EcLZ)

! Compute the restricted Levy-Zahariev LDA correlation shift of 2-glomium for various states

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
  double precision              :: dFcdr

! Output variables

  double precision,intent(out)  :: EcLZ

! Compute Levy-Zahariev eLDA correlation shift

  EcLZ = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      ec_p   = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdr = aLF(2)*r**(-1d0/6d0) + 2d0*aLF(3)*r**(-1d0/3d0)
      dFcdr = dFcdr/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdr = ec_p*dFcdr/6d0
      dFcdr = ec_p + dFcdr

      EcLZ = EcLZ - weight(iG)*r*r*dFcdr

     end if

  end do

end subroutine restricted_elda_correlation_Levy_Zahariev_shift
