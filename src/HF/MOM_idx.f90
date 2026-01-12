subroutine MOM_idx(nO, nOrb, projO, idx)
  ! Returns the indices of the nO orbitals with maximum overlap

  implicit none
  include 'parameters.h'

  ! Input variables
  integer, intent(in)          :: nOrb, nO
  double precision, intent(in) :: projO(nOrb)

  ! Output variables
  integer, intent(out)         :: idx(nO)

  ! Local variables
  integer                      :: i, j, imax
  double precision             :: pmax
  logical                      :: used(nOrb)

  ! Initialize
  used(:) = .false.

  ! Select nO largest overlaps
  do i = 1, nO
     pmax = -1.0d0
     imax = -1

     do j = 1, nOrb
        if (.not. used(j)) then
           if (projO(j) > pmax) then
              pmax = projO(j)
              imax = j
           end if
        end if
     end do

     idx(i)   = imax
     used(imax) = .true.
  end do

end subroutine
