subroutine restricted_elda_exchange_energy(nEns,Cx,nGrid,weight,rho,Ex)

! Compute the restricted LDA exchange energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: Cx
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r

! Output variables

  double precision,intent(out)  :: Ex


! Compute eLDA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex + weight(iG)*Cx*r**(4d0/3d0)

    end if

  end do

end subroutine restricted_elda_exchange_energy
