subroutine fc_srlda(nEl,nBas,nGrid,weight,MO,rho,mu)

! Compute sr-lda ec

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEl
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: MO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: mu(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r
  double precision              :: rs
  double precision              :: ecsr,vcup,vcdw
  double precision              :: IP
  double precision,parameter    :: thres = 1d-15

! Initialization

  IP = 0d0

  do iG=1,ngrid

    r = max(0d0,rho(iG))

    if(r > thres) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

      call lsdsr(rs,0d0,mu(iG),ecsr,vcup,vcdw)

      IP = IP + weight(iG)*vcup*MO(nEl,iG)**2

    end if

  end do

  print*, 'IP = ',IP*HaToeV

end subroutine fc_srlda
