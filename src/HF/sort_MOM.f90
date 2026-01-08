
subroutine sort_MOM(nBas, nOrb, projO, c,eHF)

! Sort indices according to projO (descending)
! Simple selection sort (nOrb is small in practice)

  implicit none
  include 'parameters.h'

! Input variables
  integer, intent(in)             :: nBas, nOrb
  double precision, intent(in)    :: projO(nOrb)

! Input / output variables
  double precision, intent(inout) :: c(nBas, nOrb)
  double precision, intent(inout) :: eHF(nOrb) 

! Local variables
  integer                         :: i, j, imax
  integer                         :: idx(nOrb)

! Initialize index array
  do i = 1, nOrb
     idx(i) = i
  end do

! Sort indices
  do i = 1, nOrb - 1
     imax = i
     do j = i + 1, nOrb
        if (projO(idx(j)) > projO(idx(imax))) then
           imax = j
        end if
     end do

     if (imax /= i) then
        ! swap indices
        j        = idx(i)
        idx(i)   = idx(imax)
        idx(imax)= j
     end if
  end do

! Reorder columns of c and eHF according to sorted indices
  c(:,1:nOrb) = c(:,idx(:))
  eHF(:)      = eHF(idx(:))
end subroutine
