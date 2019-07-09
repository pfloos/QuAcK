subroutine mu_grid(nGrid,rho,f,mu)

! Compute mu

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: f(nGrid)

! Local variables

  integer                       :: iG
  double precision,parameter    :: thres = 1d-10
  double precision              :: n2

! Output variables

  double precision,intent(out)  :: mu(nGrid)

! Initialization

  mu(:) = 0d0

  do iG=1,ngrid

    n2 = rho(iG)**2

    if(abs(n2) > thres) then

      mu(iG) =  f(iG)/n2
   
    else
    
      mu(iG) = 1d0/thres

    end if

  end do

  mu(:) = 0.5d0*sqrt(pi)*mu(:)

end subroutine mu_grid
