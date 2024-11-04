subroutine RGTpp_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,lambda,eGF,Taaaa,Tabab,Tbaab,KC_sta)

! Compute the VVVV block of the static T-matrix

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
  integer                       :: a,b,c,d,ab,cd,ef,mn,m,e

! Output variables

  double precision,intent(out)  :: KC_sta(nVV,nVV)

! Initialization

  KC_sta(:,:) = 0d0
  
!===============!
! singlet block !
!===============!

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a,nBas-nR
        ab = ab + 1

        cd = 0
        do c=nO+1,nBas-nR
          do d=c,nBas-nR
            cd = cd + 1
 
            do m=nC+1,nO
              do e=nO+1,nBas-nR
                 dem = eGF(m) - eGF(e)
                 ! Wabab_{ijkl}
                 num = Taaaa(a,m,c,e)*Tabab(e,b,m,d) + Tabab(a,m,c,e)*Taaaa(e,b,m,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(a,e,c,m)*Tabab(m,b,e,d) + Tabab(a,e,c,m)*Taaaa(m,b,e,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
                 
                 num = Taaaa(b,m,c,e)*Tabab(e,a,m,d) + Tabab(b,m,c,e)*Taaaa(e,a,m,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(b,e,c,m)*Tabab(m,a,e,d) + Tabab(b,e,c,m)*Taaaa(m,a,e,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                 num = Tbaab(a,m,c,e)*Tbaab(e,b,m,d) + Tbaab(a,e,c,m)*Tbaab(m,b,e,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)
        
                 num = Tbaab(b,m,c,e)*Tbaab(e,a,m,d) + Tbaab(b,e,c,m)*Tbaab(m,a,e,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)
              end do
            end do
 
            KC_sta(ab,cd) = KC_sta(ab,cd) / sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
 
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

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
 
            do m=nC+1,nO
              do e=nO+1,nBas-nR
                 dem = eGF(m) - eGF(e)

                 num = Taaaa(a,m,c,e)*Taaaa(e,b,m,d) + Tabab(a,m,c,e)*Tabab(e,b,m,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(a,e,c,m)*Taaaa(m,b,e,d) + Tabab(a,e,c,m)*Tabab(m,b,e,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
                 
                 num = Taaaa(b,m,c,e)*Taaaa(e,a,m,d) + Tabab(b,m,c,e)*Tabab(e,a,m,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)

                 num = Taaaa(b,e,c,m)*Taaaa(m,a,e,d) + Tabab(b,e,c,m)*Tabab(m,a,e,d)
                 KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)

              end do
            end do
 
          end do
        end do

      end do
    end do

  end if

end subroutine 
