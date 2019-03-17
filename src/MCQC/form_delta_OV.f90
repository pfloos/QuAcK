subroutine form_delta_OV(nO,nV,eO,eV,delta)

! Form energy denominator for CC

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)

! Local variables

  integer                       :: i,a

! Output variables

  double precision,intent(out)  :: delta(nO,nV)

    do i=1,nO
      do a=1,nV
        delta(i,a) = eV(a) - eO(i)
      enddo
    enddo

end subroutine form_delta_OV
