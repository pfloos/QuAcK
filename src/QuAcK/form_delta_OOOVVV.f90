subroutine form_delta_OOOVVV(nO,nV,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)

! Local variables

  integer                       :: i,j,k,a,b,c

! Output variables

  double precision,intent(out)  :: delta(nO,nO,nO,nV,nV,nV)

    do i=1,nO
      do j=1,nO
        do k=1,nO
          do a=1,nV
            do b=1,nV
              do c=1,nV

                delta(i,j,k,a,b,c) = eV(a) + eV(b) + eV(c) - eO(i) - eO(j) - eO(k)

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

end subroutine form_delta_OOOVVV
