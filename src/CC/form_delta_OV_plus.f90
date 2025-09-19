subroutine form_delta_OV_plus(nC,nO,nV,nR,eO,eV,delta,OVOV)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eO(nO-nC)
  double precision,intent(in)   :: eV(nV-nR)
  double precision,intent(in)   :: OVOV(nO,nV,nO,nV)

! Local variables

  integer                       :: i,a

! Output variables

  double precision,intent(inout):: delta(nO-nC,nV-nR)

    do i=1,nO-nC
      do a=1,nV-nR
        delta(i,a) = delta(i,a) - OVOV(i,a,i,a)
      end do
    end do

end subroutine 
