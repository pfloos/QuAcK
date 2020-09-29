subroutine form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: delta(nO,nO,nV,nV)

    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR

            delta(i,j,a,b) = eV(a) + eV(b) - eO(i) - eO(j)

          enddo
        enddo
      enddo
    enddo

end subroutine form_delta_OOVV
