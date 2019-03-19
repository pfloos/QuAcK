subroutine form_delta_OOVV(nO,nV,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: delta(nO,nO,nV,nV)

    do i=1,nO
      do j=1,nO
        do a=1,nV
          do b=1,nV

            delta(i,j,a,b) = eV(a) + eV(b) - eO(i) - eO(j)

          enddo
        enddo
      enddo
    enddo

end subroutine form_delta_OOVV
