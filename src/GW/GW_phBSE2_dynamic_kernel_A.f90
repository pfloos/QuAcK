subroutine GW_phBSE2_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,eGW,W,OmBSE,KA_dyn,ZA_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: W(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision              :: num,dem
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc,k,c,l,d

! Output variables

  double precision,intent(out)  :: KA_dyn(nS,nS)
  double precision,intent(out)  :: ZA_dyn(nS,nS)

! Build dynamic A matrix

  jb = 0

!$omp parallel do default(private) shared(KA_dyn,ZA_dyn,OmBSE,W,num,dem,eGW,nO,nBas,eta,nC,nR)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO)

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO)

          do k=nC+1,nO
            do c=nO+1,nBas-nR

              dem =  OmBSE - eGW(a) - eGW(c) + eGW(k) + eGW(j)
              num = 2d0*W(j,k,i,c)*W(a,c,b,k)

              KA_dyn(ia,jb) = KA_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
              ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE + eGW(i) - eGW(c) + eGW(k) - eGW(b)
              num = 2d0*W(j,c,i,k)*W(a,k,b,c)

              KA_dyn(ia,jb) = KA_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
              ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

            enddo
          enddo

          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR 

              dem = OmBSE + eGW(i) + eGW(j) - eGW(c) - eGW(d)
              num = 2d0*W(a,j,c,d)*W(c,d,i,b)

              KA_dyn(ia,jb) = KA_dyn(ia,jb) + num*dem/(dem**2 + eta**2)
              ZA_dyn(ia,jb) = ZA_dyn(ia,jb) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

            enddo
          enddo

          do k=nC+1,nO
            do l=nC+1,nO

              dem = OmBSE - eGW(a) - eGW(b) + eGW(k) + eGW(l)
              num = 2d0*W(a,j,k,l)*W(k,l,i,b)

              KA_dyn(ia,jb) = KA_dyn(ia,jb) + num*dem/(dem**2 + eta**2)
              ZA_dyn(ia,jb) = ZA_dyn(ia,jb) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

            end do
          end do 

        enddo
      enddo
    enddo
  enddo
!$omp end parallel do

end subroutine 
