subroutine f_grid(nBas,nO,nGrid,MO,ERI,f)

! Compute f

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: MO(nBas,nGrid)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p,q
  integer                       :: i,j
  integer                       :: iG

! Output variables

  double precision,intent(out)  :: f(nGrid)

! Initialization

  f(:) = 0d0

  do p=1,nBas
    do q=1,nBas
      do i=1,nO
        do j=1,nO
          do iG=1,ngrid

            f(iG) = f(iG) + MO(p,iG)*MO(q,iG)*ERI(p,q,i,j)*MO(i,iG)*MO(j,iG)

          end do
        end do
      end do
    end do
  end do

end subroutine f_grid
