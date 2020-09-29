subroutine AOtoMO_oovv(nBas,nO,nV,cO,cV,O,ooOvv)

! AO to MO transformation of two-electron integrals for the block oovv
! Semi-direct O(N^5) algorithm

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nO,nV
  double precision,intent(in)   :: cO(nBas,nO),cV(nBas,nV),O(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: scr1(:,:,:,:),scr2(:,:,:,:)
  integer                       :: mu,nu,la,si,i,j,a,b

! Output variables

  double precision,intent(out)  :: ooOvv(nO,nO,nV,nV)

! Memory allocation
  allocate(scr1(nBas,nBas,nBas,nBas),scr2(nBas,nBas,nBas,nBas))

  scr1 = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          do b=1,nV
            scr1(mu,nu,la,b) = scr1(mu,nu,la,b) + O(mu,nu,la,si)*cV(si,b)
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
          do b=1,nV
            scr2(i,nu,la,b) = scr2(i,nu,la,b) + cO(mu,i)*scr1(mu,nu,la,b)
          enddo
        enddo
      enddo
    enddo
  enddo

  scr1 = 0d0
  do nu=1,nBas
    do la=1,nBas
      do i=1,nO
        do a=1,nV
          do b=1,nV
            scr1(i,nu,a,b) = scr1(i,nu,a,b) + scr2(i,nu,la,b)*cV(la,a)
          enddo
        enddo
      enddo
    enddo
  enddo

  ooOvv = 0d0
  do nu=1,nBas
    do i=1,nO
      do j=1,nO
        do a=1,nV
          do b=1,nV
            ooOvv(i,j,a,b) = ooOvv(i,j,a,b) + cO(nu,j)*scr1(i,nu,a,b)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine AOtoMO_oovv
