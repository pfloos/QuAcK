subroutine allocate_grid(nNuc,ZNuc,rNuc,nShell,TotAngMomShell,ExpShell,max_ang_mom,min_exponent,max_exponent, & 
                         radial_precision,nRad,nAng,nGrid)

! Allocate quadrature grid with numgrid (Radovan Bast)

  use numgrid
  use, intrinsic :: iso_c_binding, only: c_ptr

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)

  integer,intent(in)            :: nShell
  integer,intent(in)            :: TotAngMomShell(maxShell)
  double precision,intent(in)   :: ExpShell(maxShell,maxK)

  integer,intent(in)            :: max_ang_mom(nNuc)
  double precision,intent(in)   :: min_exponent(nNuc,maxL+1)
  double precision,intent(in)   :: max_exponent(nNuc)

  double precision              :: radial_precision
  integer,intent(in)            :: nRad
  integer,intent(in)            :: nAng

! Local variables
 
  integer                       :: iNuc

  integer                       :: min_num_angular_points
  integer                       :: max_num_angular_points
  integer                       :: num_points

  integer                       :: center_index
  type(c_ptr)                   :: context

! Output variables

  integer,intent(out)           :: nGrid

! Set useful variables

  min_num_angular_points = nAng
  max_num_angular_points = nAng

! Get total number of grid points

  nGrid = 0

  do iNuc=1,nNuc

    context = numgrid_new_atom_grid(radial_precision,min_num_angular_points,max_num_angular_points, &
                                    int(ZNuc(iNuc)),max_exponent(iNuc),max_ang_mom(iNuc),           &
                                    min_exponent(iNuc,1:max_ang_mom(iNuc)+1))

    nGrid = nGrid + numgrid_get_num_grid_points(context)

    call numgrid_free_atom_grid(context)

  end do

end subroutine allocate_grid
