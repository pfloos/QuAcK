subroutine GF2_phBSE2_dynamic_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,eGF,OmBSE,KA_dyn,ZA_dyn)

! Compute the resonant part of the dynamic BSE2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision              :: dem,num
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb

! Output variables

  double precision,intent(out)  :: KA_dyn(nS,nS)
  double precision,intent(out)  :: ZA_dyn(nS,nS)

! Initialization

  KA_dyn(:,:) = 0d0
  ZA_dyn(:,:) = 0d0

! Second-order correlation kernel for the block A of the singlet manifold

  if(ispin == 1) then

    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = (b-nO) + (j-1)*(nBas-nO)
 
        ia = 0
        do i=nC+1,nO
          do a=nO+1,nBas-nR
            ia = (a-nO) + (i-1)*(nBas-nO) 
  
            do k=nC+1,nO
              do c=nO+1,nBas-nR
     
                dem = OmBSE - eGF(a) + eGF(k) - eGF(c) + eGF(j)
                num = 2d0*ERI(j,k,i,c)*ERI(a,c,b,k) -     ERI(j,k,i,c)*ERI(a,c,k,b) & 
                    -     ERI(j,k,c,i)*ERI(a,c,b,k) + 2d0*ERI(j,k,c,i)*ERI(a,c,k,b)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE + eGF(i) - eGF(c) + eGF(k) - eGF(b)
                num = 2d0*ERI(j,c,i,k)*ERI(a,k,b,c) -     ERI(j,c,i,k)*ERI(a,k,c,b) & 
                    -     ERI(j,c,k,i)*ERI(a,k,b,c) + 2d0*ERI(j,c,k,i)*ERI(a,k,c,b)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

            do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
     
                dem = OmBSE + eGF(i) + eGF(j) - eGF(c) - eGF(d)
                num = 2d0*ERI(a,j,c,d)*ERI(c,d,i,b) -     ERI(a,j,c,d)*ERI(c,d,b,i) & 
                    -     ERI(a,j,d,c)*ERI(c,d,i,b) + 2d0*ERI(a,j,d,c)*ERI(c,d,b,i)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

            do k=nC+1,nO
              do l=nC+1,nO

                dem = OmBSE - eGF(a) - eGF(b) + eGF(k) + eGF(l)
                num = 2d0*ERI(a,j,k,l)*ERI(k,l,i,b) -     ERI(a,j,k,l)*ERI(k,l,b,i) & 
                    -     ERI(a,j,l,k)*ERI(k,l,i,b) + 2d0*ERI(a,j,l,k)*ERI(k,l,b,i)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
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
!$omp parallel do default(private) shared(A_dyn,ZA_dyn,ERI,OmBSE,num,dem,eGF,nO,nBas,eta,nC,nR)
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = (b-nO) + (j-1)*(nBas-nO)

        ia = 0
        do i=nC+1,nO
          do a=nO+1,nBas-nR
            ia = (a-nO) + (i-1)*(nBas-nO)
  
            do k=nC+1,nO
              do c=nO+1,nBas-nR
     
                dem = OmBSE - eGF(a) + eGF(k) - eGF(c) + eGF(j)
                num = 2d0*ERI(j,k,i,c)*ERI(a,c,b,k) - ERI(j,k,i,c)*ERI(a,c,k,b) - ERI(j,k,c,i)*ERI(a,c,b,k) 

                KA_dyn(ia,jb) = KA_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE + eGF(i) - eGF(c) + eGF(k) - eGF(b)
                num = 2d0*ERI(j,c,i,k)*ERI(a,k,b,c) - ERI(j,c,i,k)*ERI(a,k,c,b) - ERI(j,c,k,i)*ERI(a,k,b,c)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

            do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
     
                dem = OmBSE + eGF(i) + eGF(j) - eGF(c) - eGF(d)
                num = ERI(a,j,c,d)*ERI(c,d,b,i) + ERI(a,j,d,c)*ERI(c,d,i,b)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

            do k=nC+1,nO
              do l=nC+1,nO

                dem = OmBSE - eGF(a) - eGF(b) + eGF(k) + eGF(l)
                num = ERI(a,j,k,l)*ERI(k,l,b,i) + ERI(a,j,l,k)*ERI(k,l,i,b)

                KA_dyn(ia,jb) = KA_dyn(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do
 
          end do
        end do

      end do
    end do
!$omp end parallel do

  end if


end subroutine 
