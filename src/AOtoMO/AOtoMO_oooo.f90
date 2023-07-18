subroutine AOtoMO_oooo(nBas,nO,cO,O,ooOoo)

! AO to MO transformation of two-electron integrals for the block oooo
! Semi-direct O(N^5) algorithm

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nO
  double precision,intent(in)   :: cO(nBas,nO),O(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: scr1(:,:,:,:),scr2(:,:,:,:)
  integer                       :: mu,nu,la,si,i,j,k,l

! Output variables

  double precision,intent(out)  :: ooOoo(nO,nO,nO,nO)

! Memory allocation
  allocate(scr1(nBas,nBas,nBas,nBas),scr2(nBas,nBas,nBas,nBas))

  write(*,*)
  write(*,'(A42)')           '----------------------------------------'
  write(*,'(A42)')           ' AO to MO transformation for oooo block '
  write(*,'(A42)')           '----------------------------------------'
  write(*,*)

  scr1 = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          do l=1,nO
            scr1(mu,nu,la,l) = scr1(mu,nu,la,l) + O(mu,nu,la,si)*cO(si,l)
          enddo
        enddo
      enddo
    enddo
  enddo
 
  scr2 = 0d0   
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do i=1,nO
          do l=1,nO
            scr2(i,nu,la,l) = scr2(i,nu,la,l) + cO(mu,i)*scr1(mu,nu,la,l)
          enddo
        enddo
      enddo
    enddo
  enddo

  scr1 = 0d0
  do nu=1,nBas
    do la=1,nBas
      do i=1,nO
        do k=1,nO
          do l=1,nO
            scr1(i,nu,k,l) = scr1(i,nu,k,l) + scr2(i,nu,la,l)*cO(la,k)
          enddo
        enddo
      enddo
    enddo
  enddo

  ooOoo = 0d0
  do nu=1,nBas
    do i=1,nO
      do j=1,nO
        do k=1,nO
          do l=1,nO
            ooOoo(i,j,k,l) = ooOoo(i,j,k,l) + cO(nu,j)*scr1(i,nu,k,l)
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(scr1,scr2)

end subroutine 
