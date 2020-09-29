subroutine antisymmetrize_ERI(ispin,nBas,ERI,db_ERI)

! Antisymmetrize ERIs

  implicit none

! Input variables

  integer,intent(in)            :: ispin,nBas
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,k,l

! Output variables

  double precision,intent(out)  :: db_ERI(nBas,nBas,nBas,nBas)

  if(ispin == 1) then

    do i=1,nBas
      do j=1,nBas
        do k=1,nBas
          do l=1,nBas
            db_ERI(i,j,k,l) = 2d0*ERI(i,j,k,l) - ERI(i,j,l,k)
          enddo
        enddo
      enddo
    enddo

  elseif(ispin == 2) then

    do i=1,nBas
      do j=1,nBas
        do k=1,nBas
          do l=1,nBas
            db_ERI(i,j,k,l) = ERI(i,j,k,l) - ERI(i,j,l,k)
          enddo
        enddo
      enddo
    enddo

  endif

end subroutine antisymmetrize_ERI
