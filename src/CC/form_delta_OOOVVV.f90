subroutine form_delta_OOOVVV(nC,nO,nV,nR,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)

! Local variables

  integer                       :: i,j,k,a,b,c

! Output variables

  double precision,intent(out)  :: delta(nO,nO,nO,nV,nV,nV)

    do i=nC+1,nO
      do j=nC+1,nO
        do k=nC+1,nO
          do a=1,nV-nR
            do b=1,nV-nR
              do c=1,nV-nR

                delta(i,j,k,a,b,c) = eV(a) + eV(b) + eV(c) - eO(i) - eO(j) - eO(k)

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

end subroutine form_delta_OOOVVV
