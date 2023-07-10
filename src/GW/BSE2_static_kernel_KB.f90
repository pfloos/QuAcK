subroutine BSE2_static_kernel_KB(eta,nBas,nC,nO,nV,nR,nS,lambda,eW,W,KB2_sta)

! Compute the second-order static BSE kernel for the coupling block (only for singlets!)

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda

  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: W(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: dem,num
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb

! Output variables

  double precision,intent(out)   :: KB2_sta(nS,nS)

!------------------------------------------------
! Compute BSE2 kernel
!------------------------------------------------

  jb = 0
!$omp parallel do default(private) shared(KB2_sta,W,num,dem,eW,nO,nBas,eta,nC,nR)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO) 

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO) 

          do k=nC+1,nO
            do c=nO+1,nBas-nR

              dem = - (eW(c) - eW(k))
              num = 2d0*W(b,k,i,c)*W(a,c,j,k)

              KB2_sta(ia,jb) =  KB2_sta(ia,jb) - num*dem/(dem**2 + eta**2)

              dem = + (eW(c) - eW(k))
              num = 2d0*W(b,c,i,k)*W(a,k,j,c) 

              KB2_sta(ia,jb) =  KB2_sta(ia,jb) + num*dem/(dem**2 + eta**2)

            end do
          end do

          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR

              dem = - (eW(c) + eW(d))
              num = 2d0*W(a,b,c,d)*W(c,d,i,j)

              KB2_sta(ia,jb) =  KB2_sta(ia,jb) + num*dem/(dem**2 + eta**2)

            end do
          end do

          do k=nC+1,nO
            do l=nC+1,nO

              dem = - (eW(k) + eW(l))
              num = 2d0*W(a,b,k,l)*W(k,l,i,j) 

              KB2_sta(ia,jb) =  KB2_sta(ia,jb) - num*dem/(dem**2 + eta**2)

            end do
          end do

        end do
      end do

    end do
  end do

!$omp end parallel do

end subroutine BSE2_static_kernel_KB
