subroutine RGF2_ppBSE2_dynamic_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,eGF,KB_dyn)

! Compute the resonant part of the dynamic BSE2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision              :: dem,num
  integer                       :: i,j,k,a,b,c
  integer                       :: ab,ij

! Output variables

  double precision,intent(out)  :: KB_dyn(nVV,nOO)

! Initialization

  KB_dyn(:,:) = 0d0

! Second-order correlation kernel for the block B of the singlet manifold

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=i,nO
            ij = ij + 1
  
            do k=nC+1,nO
              do c=nO+1,nBas-nR
     
                dem = eGF(j) + eGF(k) - eGF(c) - eGF(b)
                num = 2d0*ERI(a,k,i,c)*ERI(b,c,j,k) -     ERI(a,k,i,c)*ERI(b,c,k,j) & 
                    -     ERI(a,k,c,i)*ERI(b,c,j,k) + 2d0*ERI(a,k,c,i)*ERI(b,c,k,j)

                KB_dyn(ab,ij) = KB_dyn(ab,ij) + 0.5d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(j) + eGF(k) - eGF(c) - eGF(a)
                num = 2d0*ERI(b,k,i,c)*ERI(a,c,j,k) -     ERI(b,k,i,c)*ERI(a,c,k,j) & 
                    -     ERI(b,k,c,i)*ERI(a,c,j,k) + 2d0*ERI(b,k,c,i)*ERI(a,c,k,j)

                KB_dyn(ab,ij) = KB_dyn(ab,ij) - 0.5d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(i) + eGF(k) - eGF(c) - eGF(a)
                num = 2d0*ERI(a,c,i,k)*ERI(b,k,j,c) -     ERI(a,c,i,k)*ERI(b,k,c,j) & 
                    -     ERI(a,c,k,i)*ERI(b,k,j,c) + 2d0*ERI(a,c,k,i)*ERI(b,k,c,j)

                KB_dyn(ab,ij) = KB_dyn(ab,ij) + 0.5d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(i) + eGF(k) - eGF(c) - eGF(b)
                num = 2d0*ERI(b,c,i,k)*ERI(a,k,j,c) -     ERI(b,c,i,k)*ERI(a,k,c,j) & 
                    -     ERI(b,c,k,i)*ERI(a,k,j,c) + 2d0*ERI(b,c,k,i)*ERI(a,k,c,j)

                KB_dyn(ab,ij) = KB_dyn(ab,ij) - 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

! Second-order correlation kernel for the block B of the triplet manifold

  if(ispin == 2) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a+1,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=i+1,nO
            ij = ij + 1
  
            do k=nC+1,nO
              do c=nO+1,nBas-nR
     
                dem = eGF(j) + eGF(k) - eGF(c) - eGF(b)
                num = 2d0*ERI(a,k,i,c)*ERI(b,c,j,k) - ERI(a,k,i,c)*ERI(b,c,k,j) - ERI(a,k,c,i)*ERI(b,c,j,k) 

                KB_dyn(ab,ij) = KB_dyn(ab,ij) + 0.5d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(j) + eGF(k) - eGF(c) - eGF(a)
                num = 2d0*ERI(b,k,i,c)*ERI(a,c,j,k) - ERI(b,k,i,c)*ERI(a,c,k,j) - ERI(b,k,c,i)*ERI(a,c,j,k) 

                KB_dyn(ab,ij) = KB_dyn(ab,ij) - 0.5d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(i) + eGF(k) - eGF(c) - eGF(a)
                num = 2d0*ERI(a,c,i,k)*ERI(b,k,j,c) - ERI(a,c,i,k)*ERI(b,k,c,j) - ERI(a,c,k,i)*ERI(b,k,j,c) 

                KB_dyn(ab,ij) = KB_dyn(ab,ij) + 0.5d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(i) + eGF(k) - eGF(c) - eGF(b)
                num = 2d0*ERI(b,c,i,k)*ERI(a,k,j,c) - ERI(b,c,i,k)*ERI(a,k,c,j) - ERI(b,c,k,i)*ERI(a,k,j,c) 

                KB_dyn(ab,ij) = KB_dyn(ab,ij) - 0.5d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
