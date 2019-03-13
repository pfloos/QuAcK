subroutine AOtoMO_oooooo(nBas,nO,cO,O,oooOooo)

! AO to MO transformation of three-electron integrals for the block oooooo
! Semi-direct O(N^7) algorithm

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nO
  double precision,intent(in)   :: cO(nBas,nO),O(nBas,nBas,nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: scr1(:,:,:,:,:,:),scr2(:,:,:,:,:,:)
  integer                       :: mu,nu,la,si,ka,ta,i,j,k,l,m,n

! Output variables

  double precision,intent(out)  :: oooOooo(nO,nO,nO,nO,nO,nO)

! Memory allocation
  allocate(scr1(nBas,nBas,nBas,nBas,nBas,nBas),scr2(nBas,nBas,nBas,nBas,nBas,nBas))

  write(*,*)
  write(*,'(A42)')           '------------------------------------------'
  write(*,'(A42)')           ' AO to MO transformation for oooooo block '
  write(*,'(A42)')           '------------------------------------------'
  write(*,*)

  scr1 = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          do ka=1,nBas
            do ta=1,nBas
              do n=1,nO
                scr1(mu,nu,la,si,ka,n) = scr1(mu,nu,la,si,ka,n) + O(mu,nu,la,si,ka,ta)*cO(ta,n)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
 
  scr2 = 0d0   
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          do ka=1,nBas
            do i=1,nO
              do n=1,nO
                scr2(i,nu,la,si,ka,n) = scr2(i,nu,la,si,ka,n) + cO(mu,i)*scr1(mu,nu,la,si,ka,n)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  scr1 = 0d0
  do nu=1,nBas
    do la=1,nBas
      do si=1,nBas
        do ka=1,nBas
          do i=1,nO
            do m=1,nO
              do n=1,nO
                scr1(i,nu,la,si,m,n) = scr1(i,nu,la,si,m,n) + scr2(i,nu,la,si,m,n)*cO(ka,m)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  scr2 = 0d0   
  do nu=1,nBas
    do la=1,nBas
      do si=1,nBas
        do i=1,nO
          do j=1,nO
            do m=1,nO
              do n=1,nO
                scr2(i,j,la,si,m,n) = scr2(i,j,la,si,m,n) + cO(nu,j)*scr1(i,nu,la,si,m,n)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  scr1 = 0d0   
  do la=1,nBas
    do si=1,nBas
      do i=1,nO
        do j=1,nO
          do l=1,nO
            do m=1,nO
              do n=1,nO
                scr1(i,j,la,l,m,n) = scr1(i,j,la,l,m,n) + cO(si,l)*scr2(i,j,la,si,m,n)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  oooOooo = 0d0
  do si=1,nBas
    do i=1,nO
      do j=1,nO
        do k=1,nO
          do l=1,nO
            do m=1,nO
              do n=1,nO
                oooOooo(i,j,k,l,m,n) = oooOooo(i,j,k,l,m,n) + cO(la,k)*scr1(i,j,la,k,m,n)
               enddo
             enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  deallocate(scr1,scr2)

end subroutine AOtoMO_oooooo
