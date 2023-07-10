subroutine GF2_phBSE2_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,eGF,A_sta)

! Compute the resonant part of the static BSE2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision              :: dem,num
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb

! Output variables

  double precision,intent(out)  :: A_sta(nS,nS)

! Initialization

   A_sta(:,:) = 0d0

! Second-order correlation kernel for the block A of the singlet manifold

  if(ispin == 1) then

    jb = 0
!$omp parallel do default(private) shared(A_sta,ERI,num,dem,eGF,nO,nBas,eta,nC,nR)
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = (b-nO) + (j-1)*(nBas-nO)
 
        ia = 0
        do i=nC+1,nO
          do a=nO+1,nBas-nR
            ia = (a-nO) + (i-1)*(nBas-nO)
            
  
            do k=nC+1,nO
              do c=nO+1,nBas-nR
     
                dem = - (eGF(c) - eGF(k))
                num = 2d0*ERI(j,k,i,c)*ERI(a,c,b,k) -     ERI(j,k,i,c)*ERI(a,c,k,b) & 
                    -     ERI(j,k,c,i)*ERI(a,c,b,k) + 2d0*ERI(j,k,c,i)*ERI(a,c,k,b)

                A_sta(ia,jb) =  A_sta(ia,jb) - num*dem/(dem**2 + eta**2)
            
                dem = + (eGF(c) - eGF(k))
                num = 2d0*ERI(j,c,i,k)*ERI(a,k,b,c) -     ERI(j,c,i,k)*ERI(a,k,c,b) & 
                    -     ERI(j,c,k,i)*ERI(a,k,b,c) + 2d0*ERI(j,c,k,i)*ERI(a,k,c,b)

                A_sta(ia,jb) =  A_sta(ia,jb) + num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
     
                dem = - (eGF(c) + eGF(d))
                num = 2d0*ERI(a,j,c,d)*ERI(c,d,i,b) -     ERI(a,j,c,d)*ERI(c,d,b,i) & 
                    -     ERI(a,j,d,c)*ERI(c,d,i,b) + 2d0*ERI(a,j,d,c)*ERI(c,d,b,i)

                A_sta(ia,jb) =  A_sta(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do k=nC+1,nO
              do l=nC+1,nO

                dem = - (eGF(k) + eGF(l))
                num = 2d0*ERI(a,j,k,l)*ERI(k,l,i,b) -     ERI(a,j,k,l)*ERI(k,l,b,i) & 
                    -     ERI(a,j,l,k)*ERI(k,l,i,b) + 2d0*ERI(a,j,l,k)*ERI(k,l,b,i)

                A_sta(ia,jb) =  A_sta(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do
 
          end do
        end do

      end do
    end do
!$omp end parallel do

  end if

! Second-order correlation kernel for the block A of the triplet manifold

  if(ispin == 2) then

    jb = 0
!$omp parallel do default(private) shared(A_sta,ERI,num,dem,eGF,nO,nBas,eta,nC,nR)
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = (b-nO) + (j-1)*(nBas-nO)
 
        ia = 0
        do i=nC+1,nO
          do a=nO+1,nBas-nR
            ia = (a-nO) + (i-1)*(nBas-nO)
  
            do k=nC+1,nO
              do c=nO+1,nBas-nR
     
                dem = - (eGF(c) - eGF(k))
                num = 2d0*ERI(j,k,i,c)*ERI(a,c,b,k) - ERI(j,k,i,c)*ERI(a,c,k,b) - ERI(j,k,c,i)*ERI(a,c,b,k) 

                A_sta(ia,jb) =  A_sta(ia,jb) - num*dem/(dem**2 + eta**2)
            
                dem = + (eGF(c) - eGF(k))
                num = 2d0*ERI(j,c,i,k)*ERI(a,k,b,c) - ERI(j,c,i,k)*ERI(a,k,c,b) - ERI(j,c,k,i)*ERI(a,k,b,c)

                A_sta(ia,jb) =  A_sta(ia,jb) + num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
     
                dem = - (eGF(c) + eGF(d))
                num = ERI(a,j,c,d)*ERI(c,d,b,i) + ERI(a,j,d,c)*ERI(c,d,i,b)

                A_sta(ia,jb) =  A_sta(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do k=nC+1,nO
              do l=nC+1,nO

                dem = - (eGF(k) + eGF(l))
                num = ERI(a,j,k,l)*ERI(k,l,b,i) + ERI(a,j,l,k)*ERI(k,l,i,b)

                A_sta(ia,jb) =  A_sta(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do
 
          end do
        end do

      end do
    end do
!$omp end parallel do

  end if


end subroutine 
