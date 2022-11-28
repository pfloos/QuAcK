subroutine BSE2_static_kernel(eta,nBas,nC,nO,nV,nR,nS,lambda,eW,ERI,Om,rho,A_sta)

! Compute the second-order static BSE kernel (only for singlets!)

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

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: Om(nS)

  double precision,intent(in)   :: rho(nBas,nBas,nS)


! Local variables

  double precision              :: chi
  integer                       :: p,q,r,s
  integer                       :: m

  double precision              :: dem,num
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb

  double precision,allocatable   :: W(:,:,:,:)

! Output variables

  double precision,intent(inout)   :: A_sta(nS,nS)

! memory allocation

  allocate(W(nBas,nBas,nBas,nBas))

!------------------------------------------------
! Compute static screening (physicist's notation)
!------------------------------------------------

  do p=1,nBas
    do q=1,nBas
      do r=1,nBas
        do s=1,nBas

          chi = 0d0
          do m=1,nS
            dem = Om(m)**2 + eta**2
            chi = chi + rho(p,q,m)*rho(r,s,m)*Om(m)/dem
          enddo

          W(p,s,q,r) = - ERI(p,s,q,r) + 4d0*chi

        enddo
      enddo
    enddo
  enddo

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1

      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          do k=nC+1,nO
            do c=nO+1,nBas-nR

              dem = - (eW(c) - eW(k))
              num = 2d0*W(j,k,i,c)*W(a,c,b,k)

              A_sta(ia,jb) =  A_sta(ia,jb) - num*dem/(dem**2 + eta**2)

              dem = + (eW(c) - eW(k))
              num = 2d0*W(j,c,i,k)*W(a,k,b,c) 

              A_sta(ia,jb) =  A_sta(ia,jb) + num*dem/(dem**2 + eta**2)

            end do
          end do

          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR

              dem = - (eW(c) + eW(d))
              num = 2d0*W(a,j,c,d)*W(c,d,i,b)

              A_sta(ia,jb) =  A_sta(ia,jb) + num*dem/(dem**2 + eta**2)

            end do
          end do

          do k=nC+1,nO
            do l=nC+1,nO

              dem = - (eW(k) + eW(l))
              num = 2d0*W(a,j,k,l)*W(k,l,i,b) 

              A_sta(ia,jb) =  A_sta(ia,jb) - num*dem/(dem**2 + eta**2)

            end do
          end do

        end do
      end do

    end do
  end do

end subroutine BSE2_static_kernel
