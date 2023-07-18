subroutine AOtoMO_oooa(nBas,nO,nA,cO,cA,O,ooOoa)

! AO to MO transformation of two-electron integrals for the block oooa
! Semi-direct O(N^5) algorithm

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nO,nA
  double precision,intent(in)   :: cO(nBas,nO),cA(nBas,nA),O(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: scr1(:,:,:,:),scr2(:,:,:,:)
  integer                       :: mu,nu,la,si,i,j,k,x

! Output variables

  double precision,intent(out)  :: ooOoa(nO,nO,nO,nA)

! Memory allocation
  allocate(scr1(nBas,nBas,nBas,nBas),scr2(nBas,nBas,nBas,nBas))

  write(*,*)
  write(*,'(A42)')           '----------------------------------------'
  write(*,'(A42)')           ' AO to MO transformation for oooa block '
  write(*,'(A42)')           '----------------------------------------'
  write(*,*)

  scr1 = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          do x=1,nA
            scr1(mu,nu,la,x) = scr1(mu,nu,la,x) + O(mu,nu,la,si)*cA(si,x)
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
          do x=1,nA
            scr2(i,nu,la,x) = scr2(i,nu,la,x) + cO(mu,i)*scr1(mu,nu,la,x)
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
          do x=1,nA
            scr1(i,nu,k,x) = scr1(i,nu,k,x) + scr2(i,nu,la,x)*cO(la,k)
          enddo
        enddo
      enddo
    enddo
  enddo

  ooOoa = 0d0
  do nu=1,nBas
    do i=1,nO
      do j=1,nO
        do k=1,nO
          do x=1,nA
            ooOoa(i,j,k,x) = ooOoa(i,j,k,x) + cO(nu,j)*scr1(i,nu,k,x)
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(scr1,scr2)

end subroutine 
