subroutine form_delta_OV(nC,nO,nV,nR,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eO(nO-nC)
  double precision,intent(in)   :: eV(nV-nR)

! Local variables

  integer                       :: i,a

! Output variables

  double precision,intent(out)  :: delta(nO-nC,nV-nR)

    do i=1,nO-nC
      do a=1,nV-nR
        delta(i,a) = eV(a) - eO(i)
      end do
    end do

end subroutine 
