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
  double precision              :: ec,vcup,vcdw

! Initialization

  ecsr = 0d0

  do iG=1,ngrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

      call lsdsr(rs,0d0,mu(iG),ec,vcup,vcdw)

      ecsr = ecsr + weight(iG)*ec*r

    end if

  end do

  write(*,'(A32,1X,F16.10)') 'ecsr = ',ecsr

end subroutine ec_srlda
