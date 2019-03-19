subroutine elda_correlation_Levy_Zahariev_shift(nEns,aLF,nGrid,weight,rho,EcLZ)

! Compute Levy-Zahariev LDA correlation shift of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: aLF(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r,ec_p
  double precision              :: dFcdra,dFcdrb

! Output variables

  double precision,intent(out)  :: EcLZ

! Compute Levy-Zahariev eLDA correlation shift

  EcLZ = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))
    r = ra + rb

    if(ra > threshold) then

      ec_p   = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdra = aLF(2)*r**(-1d0/6d0) + 2d0*aLF(3)*r**(-1d0/3d0)
      dFcdra = dFcdra/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdra = ec_p*dFcdra/6d0
      dFcdra = ec_p + dFcdra

      EcLZ = EcLZ - weight(iG)*r*r*dFcdra

     end if

    if(rb > threshold) then

      ec_p   = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdrb = aLF(2)*r**(-1d0/6d0) + 2d0*aLF(3)*r**(-1d0/3d0)
      dFcdrb = dFcdrb/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
      dFcdrb = ec_p*dFcdrb/6d0
      dFcdrb = ec_p + dFcdrb

      EcLZ = EcLZ - weight(iG)*r*r*dFcdrb

     end if

  end do

end subroutine elda_correlation_Levy_Zahariev_shift
