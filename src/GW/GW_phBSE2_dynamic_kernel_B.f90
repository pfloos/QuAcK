subroutine GW_phBSE2_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,eGW,W,KB_dyn)

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
  
! Local variables

  double precision              :: num,dem
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc,k,c,l,d

! Output variables

  double precision,intent(out)  :: KB_dyn(nS,nS)

! Build dynamic B matrix

  jb = 0

!$omp parallel do default(private) shared(KB_dyn,W,num,dem,eGW,nO,nBas,eta,nC,nR)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO)

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO)

          do k=nC+1,nO
            do c=nO+1,nBas-nR

              dem =        - eGW(a) - eGW(c) + eGW(k) + eGW(j)
              num = 2d0*W(b,k,i,c)*W(a,c,j,k)

              KB_dyn(ia,jb) = KB_dyn(ia,jb) - num*dem/(dem**2 + eta**2)

              dem =       + eGW(i) - eGW(c) + eGW(k) - eGW(b)
              num = 2d0*W(b,c,i,k)*W(a,k,j,c)

              KB_dyn(ia,jb) = KB_dyn(ia,jb) - num*dem/(dem**2 + eta**2)

            end do
          end do

          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR 

              dem =       + eGW(i) + eGW(j) - eGW(c) - eGW(d)
              num = 2d0*W(a,b,c,d)*W(c,d,i,j)

              KB_dyn(ia,jb) =  KB_dyn(ia,jb) + num*dem/(dem**2 + eta**2)

            end do
          end do

          do k=nC+1,nO
            do l=nC+1,nO

              dem =       - eGW(a) - eGW(b) + eGW(k) + eGW(l)
              num = 2d0*W(a,b,k,l)*W(k,l,i,j)

              KB_dyn(ia,jb) =  KB_dyn(ia,jb) + num*dem/(dem**2 + eta**2)

            end do
          end do 

        end do
      end do
    end do
  end do
!$omp end parallel do

end subroutine 
