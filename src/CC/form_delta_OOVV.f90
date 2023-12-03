subroutine form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: eO(nO-nC)
  double precision,intent(in)   :: eV(nV-nR)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: delta(nO-nC,nO-nC,nV-nR,nV-nR)

    do i=1,nO-nC
      do j=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR

            delta(i,j,a,b) = eV(a) + eV(b) - eO(i) - eO(j)

          end do
        end do
      end do
    end do

end subroutine 
