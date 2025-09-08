subroutine form_delta_OOVV_plus(nC,nO,nV,nR,eO,eV,delta,OVOV,OOOO,VVVV)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: eO(nO-nC)
  double precision,intent(in)   :: eV(nV-nR)
  double precision,intent(in)   :: OVOV(nO,nV,nO,nV)
  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)


! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(inout):: delta(nO-nC,nO-nC,nV-nR,nV-nR)

    do i=1,nO-nC
      do j=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR

            delta(i,j,a,b) = delta(i,j,a,b) &
                           - OVOV(i,a,i,a) - OVOV(j,a,j,a) - OVOV(i,b,i,b) - OVOV(j,b,j,b) &
                           + OOOO(i,j,i,j) + VVVV(a,b,a,b)

          end do
        end do
      end do
    end do

end subroutine 
