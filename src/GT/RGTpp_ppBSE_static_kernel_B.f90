subroutine RGTpp_ppBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,lambda,eGF,Taaaa,Tabab,Tbaab,KB_sta)

! Compute the VVOO block of the static T-matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: Taaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Tabab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Tbaab(nBas,nBas,nBas,nBas)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: i,j,a,b,ij,ab,cd,kl,m,e

! Output variables

  double precision,intent(out)  :: KB_sta(nVV,nOO)

! Initialization

  KB_sta(:,:) = 0d0
  
!===============!
! singlet block !
!===============!

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=i,nO
            ij = ij + 1
 
            do m=nC+1,nO
              do e=nO+1,nBas-nR
                 dem = eGF(m) - eGF(e)
                 ! Wabab_{ijkl}
                 num = Taaaa(a,m,i,e)*Tabab(e,b,m,j) + Tabab(a,m,i,e)*Taaaa(e,b,m,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(a,e,i,m)*Tabab(m,b,e,j) + Tabab(a,e,i,m)*Taaaa(m,b,e,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)
                 
                 num = Taaaa(b,m,i,e)*Tabab(e,a,m,j) + Tabab(b,m,i,e)*Taaaa(e,a,m,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(b,e,i,m)*Tabab(m,a,e,j) + Tabab(b,e,i,m)*Taaaa(m,a,e,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)

                 num = Tbaab(a,m,i,e)*Tbaab(e,b,m,j) + Tbaab(a,e,i,m)*Tbaab(m,b,e,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) - num*dem/(dem**2 + eta**2)
        
                 num = Tbaab(b,m,i,e)*Tbaab(e,a,m,j) + Tbaab(b,e,i,m)*Tbaab(m,a,e,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) - num*dem/(dem**2 + eta**2)
              end do
            end do
 
            KB_sta(ab,ij) = KB_sta(ab,ij) / sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))
 
          end do
        end do

      end do
    end do

  end if

!===============!
! triplet block !
!===============!

  if(ispin == 2) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a+1,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=i+1,nO
            ij = ij + 1
 
            do m=nC+1,nO
              do e=nO+1,nBas-nR
                 dem = eGF(m) - eGF(e)

                 num = Taaaa(a,m,i,e)*Taaaa(e,b,m,j) + Tabab(a,m,i,e)*Tabab(e,b,m,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(a,e,i,m)*Taaaa(m,b,e,j) + Tabab(a,e,i,m)*Tabab(m,b,e,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)
                 
                 num = Taaaa(b,m,i,e)*Taaaa(e,a,m,j) + Tabab(b,m,i,e)*Tabab(e,a,m,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) - num*dem/(dem**2 + eta**2)

                 num = Taaaa(b,e,i,m)*Taaaa(m,a,e,j) + Tabab(b,e,i,m)*Tabab(m,a,e,j)
                 KB_sta(ab,ij) = KB_sta(ab,ij) - num*dem/(dem**2 + eta**2)

              end do
            end do
 
          end do
        end do

      end do
    end do

  end if

end subroutine 
