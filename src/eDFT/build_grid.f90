subroutine build_grid(nNuc,ZNuc,rNuc,nShell,TotAngMomShell,ExpShell,max_ang_mom,min_exponent,max_exponent, & 
                      radial_precision,nRad,nAng,nGrid,weight,root)

! Compute quadrature grid with numgrid (Radovan Bast)

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

  double precision,intent(in)   :: radial_precision
  integer,intent(in)            :: nRad
  integer,intent(in)            :: nAng
  integer,intent(in)            :: nGrid

! Local variables
 
  logical                       :: dump_grid = .false.
  integer                       :: iNuc
  integer                       :: iG

  integer                       :: min_num_angular_points
  integer                       :: max_num_angular_points
  integer                       :: num_points

  integer                       :: center_index
  type(c_ptr)                   :: context

! Output variables

  double precision,intent(out)  :: root(ncart,nGrid)
  double precision,intent(out)  :: weight(nGrid)

! Set useful variables

  min_num_angular_points = nAng
  max_num_angular_points = nAng

!------------------------------------------------------------------------
! Main loop over atoms
!------------------------------------------------------------------------

  iG = 0

  do iNuc=1,nNuc

    context = numgrid_new_atom_grid(radial_precision,min_num_angular_points,max_num_angular_points, &
                                    int(ZNuc(iNuc)),max_exponent(iNuc),max_ang_mom(iNuc),           &
                                    min_exponent(iNuc,1:max_ang_mom(iNuc)+1))

    center_index = iNuc - 1
    num_points = numgrid_get_num_grid_points(context)

    call numgrid_get_grid(context,nNuc,center_index,rNuc(:,1),rNuc(:,2),rNuc(:,3),int(ZNuc(:)),    &
                          root(1,iG+1:iG+num_points),root(2,iG+1:iG+num_points),root(3,iG+1:iG+num_points), &
                          weight(iG+1:iG+num_points))

    iG = iG + num_points

    call numgrid_free_atom_grid(context)

  end do

!------------------------------------------------------------------------
! End main loop over atoms
!------------------------------------------------------------------------

! Print grid

  if(dump_grid) then 

    write(*,*) ' ***********************'
    write(*,*) ' *** QUADRATURE GRID ***'
    write(*,*) ' ***********************'
    write(*,'(A10,3X,3A15)') 'Grid point','X','Y','Z'
    do iG=1,nGrid
       write(*,'(I10,3X,4F15.10)') iG,weight(iG),root(:,iG)
    end do
    write(*,*)

  end if

end subroutine build_grid
