subroutine rij(nWalk,r,r12)

! Compute the interelectronic distances

  implicit none

! Input variables

  integer,intent(in)            :: nWalk
  double precision,intent(in)   :: r(nWalk,1:2,1:3)

! Output variables

  double precision,intent(out)  :: r12(nWalk)

! Compute

  r12(1:nWalk) = (r(1:nWalk,1,1)-r(1:nWalk,2,1))**2 &
               + (r(1:nWalk,1,2)-r(1:nWalk,2,2))**2 &
               + (r(1:nWalk,1,3)-r(1:nWalk,2,3))**2

  r12 = sqrt(r12)

end subroutine rij
