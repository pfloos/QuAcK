subroutine GF2_phBSE2_dynamic_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,eGF,B_dyn)

! Compute the anti-resonant part of the dynamic BSE2 matrix

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

  double precision,intent(out)  :: B_dyn(nS,nS)

! Initialization

   B_dyn(:,:) = 0d0

! Second-order correlation kernel for the block A of the singlet manifold

  if(ispin == 1) then

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
     
                dem = - eGF(a) + eGF(k) - eGF(c) + eGF(j)
                num = 2d0*ERI(b,k,i,c)*ERI(a,c,j,k) -     ERI(b,k,i,c)*ERI(a,c,k,j) & 
                    -     ERI(b,k,c,i)*ERI(a,c,j,k) + 2d0*ERI(b,k,c,i)*ERI(a,c,k,j)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
            
                dem = + eGF(i) - eGF(c) + eGF(k) - eGF(b)
                num = 2d0*ERI(b,c,i,k)*ERI(a,k,j,c) -     ERI(b,c,i,k)*ERI(a,k,c,j) & 
                    -     ERI(b,c,k,i)*ERI(a,k,j,c) + 2d0*ERI(b,c,k,i)*ERI(a,k,c,j)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
     
                dem = + eGF(i) + eGF(j) - eGF(c) - eGF(d)
                num = 2d0*ERI(a,b,c,d)*ERI(c,d,i,j) -     ERI(a,b,c,d)*ERI(c,d,j,i) & 
                    -     ERI(a,b,d,c)*ERI(c,d,i,j) + 2d0*ERI(a,b,d,c)*ERI(c,d,j,i)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do k=nC+1,nO
              do l=nC+1,nO

                dem = - eGF(a) - eGF(b) + eGF(k) + eGF(l)
                num = 2d0*ERI(a,b,k,l)*ERI(k,l,i,j) -     ERI(a,b,k,l)*ERI(k,l,j,i) & 
                    -     ERI(a,b,l,k)*ERI(k,l,i,j) + 2d0*ERI(a,b,l,k)*ERI(k,l,j,i)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do
 
          end do
        end do

      end do
    end do

  end if

! Second-order correlation kernel for the block A of the triplet manifold

  if(ispin == 2) then

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
     
                dem = - eGF(a) + eGF(k) - eGF(c) + eGF(j)
                num = 2d0*ERI(b,k,i,c)*ERI(a,c,j,k) - ERI(b,k,i,c)*ERI(a,c,k,j) - ERI(b,k,c,i)*ERI(a,c,j,k) 

                 B_dyn(ia,jb) =  B_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
            
                dem = + eGF(i) - eGF(c) + eGF(k) - eGF(b)
                num = 2d0*ERI(b,c,i,k)*ERI(a,k,j,c) - ERI(b,c,i,k)*ERI(a,k,c,j) - ERI(b,c,k,i)*ERI(a,k,j,c)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) - num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
     
                dem = + eGF(i) + eGF(j) - eGF(c) - eGF(d)
                num = ERI(a,b,c,d)*ERI(c,d,j,i) + ERI(a,b,d,c)*ERI(c,d,i,j)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

            do k=nC+1,nO
              do l=nC+1,nO

                dem = - eGF(a) - eGF(b) + eGF(k) + eGF(l)
                num = ERI(a,b,k,l)*ERI(k,l,j,i) + ERI(a,b,l,k)*ERI(k,l,i,j)

                 B_dyn(ia,jb) =  B_dyn(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do
 
          end do
        end do

      end do
    end do

  end if


end subroutine 
