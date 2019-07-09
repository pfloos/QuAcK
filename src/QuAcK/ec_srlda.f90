subroutine ec_srlda(nGrid,weight,rho,mu)

! Compute sr-lda ec

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: mu(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r
  double precision              :: rs
  double precision              :: ecsr
  double precision              :: ec
  double precision,parameter    :: thres = 1d-15

! Initialization

  ec = 0d0

  do iG=1,ngrid

    r = max(0d0,rho(iG))

    if(r > thres) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

      call srlda(rs,mu(iG),ecsr)

      ec = ec + weight(iG)*ecsr*rho(iG)

    end if

  end do

  print*, 'ec = ',ec

end subroutine ec_srlda
